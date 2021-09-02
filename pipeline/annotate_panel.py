import json
import pdb 
import os
import csv
import sys
import argparse
import math
import glob
from collections import defaultdict
from itertools import chain
from scipy.stats import fisher_exact

from pipeline import run, pipe
from covermi import Panel, appris

BIOTYPE = defaultdict(int, (("protein_coding", 1), ("pseudogene", -1)))
REFSEQ = defaultdict(int, (("NM", 4), ("NR", 3), ("XM", 2), ("XR", 1)))
DELETE_NON_DIGIT = str.maketrans("", "", "ABCDEFGHIJKLMNOPQRSTUVWXYZ_")

CHROM = 0
QUAL = 5
FILTERS = 6



def vardict_read_data(row):
    fmt = dict(zip(row[8].split(":"), row[9].split(":")))
    return {"vaf": fmt["AF"], "depth": fmt["DP"], "alt_depth": fmt["VD"], "ref_fr": fmt["RD"], "alt_fr": fmt["ALD"]}



def varscan2_read_data(row):
    fmt = dict(zip(row[8].split(":"), row[9].split(":")))
    vaf = float(fmt["FREQ"].rstrip("%"))/100
    ref_fr = "{}:{}".format(fmt["RDF"], fmt["RDR"])
    alt_fr = "{}:{}".format(fmt["ADF"], fmt["ADR"])
    return {"vaf": f"{vaf:.4f}", "depth": fmt["DP"], "alt_depth": fmt["AD"], "ref_fr": ref_fr, "alt_fr": alt_fr}



def mutect2_read_data(row):
    fmt = dict(zip(row[8].split(":"), row[9].split(":")))
    vaf = float(fmt["AF"])
    alt_depth = fmt["AD"].split(",")[1]
    d4 = fmt["SB"].split(",")
    ref_fr = "{}:{}".format(*d4[:2])
    alt_fr = "{}:{}".format(*d4[2:])
    return {"vaf": f"{vaf:.4f}", "depth": fmt["DP"], "alt_depth": alt_depth, "ref_fr": ref_fr, "alt_fr": alt_fr}



def chrom2int(chrom):
    try:
        return int(chrom[3:])
    except ValueError:
        return math.inf



# reference should be a required argument as fails in vep 104 without but is made optional here for compatibility with legacy code in cfpipeline.py for vep 101
def annotate_panel(vcf, vep, reference=None, threads=None, output="output.tsv", panel="", buffer_size=None):
    if threads is None:
        threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()

    vepjson = "{}.vep.json".format(output[:-4])
    vep_options = ["--no_stats",
                   "--dir", vep,
                   "--format", "vcf",
                   "--json",
                   "--offline",
                   "--everything",
                   "--warning_file", "STDERR",
                   "--force_overwrite"]
    if reference is not None:
        vep_options += ["--fasta", reference]
    if int(threads) > 1:
        vep_options += ["--fork", threads]
    if "refseq" in vep:
        vep_options += ["--refseq"]
    if buffer_size is not None:
        vep_options += ["--buffer_size", buffer_size]
    
    pipe(["vep", "-i", vcf, "-o", vepjson] + vep_options)
    
    
    read_data = None
    with open(vcf, "rt") as f:
        for row in f:
            if not row.startswith("#"):
                break
            if row.startswith("##source="):
                source = row[9:].strip()
                #if source == "strelka":
                if source.startswith("VarDict")
                    get_read_data = vardict_read_data                    
                elif source == "VarScan2":
                    get_read_data = varscan2_read_data
                elif source == "Mutect2"
                    get_read_data = mutect2_read_data                    
            headings = row
    
    if read_data is None:
        sys.exit(f"Unsupported variant caller {source}", file=sys.stderr)
    if len(headings) > 10:
        sys.exit("Multi-sample vcfs not suppored", file=sys.stderr)
    
    
    targets = None
    principal = {}
    needed_genes = set()
    needed_transcripts = set()
    if panel:
        panel = Panel(panel)
        if "targets" in panel:
            targets = panel.targets
        
        if "names" in panel:
            for name in panel.names:
                name = name.split()
                needed_genes.add(name[0])
                if len(name) > 1:
                    needed_transcripts.add(name[1])
        if "principal" in panel.paths:
            principal = appris(panel.paths["principal"])

    if "refseq" in vep:
        def consequence_sort(cons):
            transcript, minor = cons["transcript_id"].split(".")
            prefix = transcript[:2]
            major = transcript[3:]
            return [BIOTYPE[cons["biotype"]],
                    REFSEQ[prefix],
                    -int(cons["gene_id"]),
                    principal.get(transcript, 0),
                    "canonical" in cons,
                    -int(major),
                    int(minor)]
    
    else: # ensembl transcripts
        def consequence_sort(cons):
            # Version numbers not in vep as of version 101, but who knows the future ...
            return [BIOTYPE[cons["biotype"]],
                    -int(cons["gene_id"].translate(DELETE_NON_DIGIT)),
                    principal.get(cons["transcript_id"], 0),
                    "canonical" in cons,
                    -int(cons["transcript_id"].translate(DELETE_NON_DIGIT))]
    
    annotations = []
    with open(vepjson) as f:
        for line in f:
            vep_output = json.loads(line)

            consequences = vep_output.get("transcript_consequences", ())
            consequences = [cons for cons in consequences if "gene_symbol" in cons and "transcript_id" in cons]
            if not consequences:
                continue
            
            # create variant - assume one variant per row ?always true
            row = vep_output["input"].rstrip().split("\t")
            read_data = get_read_data(row)
            
            if read_data["alt_depth"] == "0":
                continue
            
            # https://gatk.broadinstitute.org/hc/en-us/articles/360035532152-Fisher-s-Exact-Test
            ref_fr = [int(n) for n in data["ref_fr"].split(":")]
            alt_fr = [int(n) for n in data["alt_fr"].split(":")]
            fisher_strand = -10 * math.log10(fisher_exact([ref_fr, alt_fr])[1])
            
            demographics = defaultdict(list)
            for colocated in vep_output.get("colocated_variants", ()):
                allele_string = colocated["allele_string"]
                identifier = colocated["id"]
                if allele_string == "COSMIC_MUTATION":
                    demographics["cosmic"] += [identifier]

                elif allele_string == "HGMD_MUTATION":
                    demographics["hgmd"] += [identifier]

                else:
                    if identifier.startswith("rs"):
                        demographics["dbsnp"] += [identifier]
                    
                    if "pubmed" in colocated:
                        pubmed = colocated["pubmed"] if isinstance(colocated["pubmed"], list) else [colocated["pubmed"]]
                        demographics["pubmed"] += [str(ref) for ref in pubmed]
                        if "clin_sig" in colocated:
                            demographics["clin_sig"] += colocated["clin_sig"]
                        
                    frequencies = colocated.get("frequencies"][vep_output["allele_string"].split("/")[1])
                    if frequencies:
                        statistics = defaultdict(list)
                        lookup = {"gnomad_afr": "gnomad", "gnomad_amr": "gnomad", "gnomad_asj": "gnomad", "gnomad_eas": "gnomad", "gnomad_fin": "gnomad", \
                                "gnomad_nfe": "gnomad", "gnomad_oth": "gnomad", "gnomad_sas": "gnomad", "ea": "nhlbi", "aa": "nhlbi", \
                                "afr": "1k", "amr": "1k", "asn": "1k", "eur": "1k", "eas": "1k", "sas": "1k", "gnomad": "gnomad"}
                        for population, frequency in frequencies.items():
                            statistics[lookup[population]] += [frequency]
                        if statistics:
                            demographics["maf"] = max(chain(*statistics.values()))
            
            if "clin_sig" in demographics:
                demographics["clin_sig"] = sorted(set(demographics["clin_sig"]))
            

            all_genes = set()
            potentials = defaultdict(list)
            genes_with_transcripts = set()
            for cons in consequences:
                gene_symbol = cons["gene_symbol"]
                all_genes.add(gene_symbol)
                transcript_id = cons["transcript_id"].split(".")[0]
                if transcript_id in needed_transcripts:
                    potentials[transcript_id].append(cons)
                    genes_with_transcripts.add(gene_symbol)
                elif gene_symbol in needed_genes:
                    potentials[gene_symbol].append(cons)
            
            if potentials:
                selected = [sorted(val, key=consequence_sort)[-1] for key, val in potentials.items() if key not in genes_with_transcripts]
            else:
                selected = [sorted(consequences, key=consequence_sort)[-1]]
            

            for cons in selected:
                gene_symbol = cons.get("gene_symbol", "")
                annotations.append([gene_symbol,
                                    cons.get("transcript_id", ""),
                                    row[CHROM],
                                    vep_output["start"],
                                    vep_output["allele_string"],
                                    row[QUAL],
                                    row[FILTERS],
                                    read_data["vaf"],
                                    read_data["depth"],
                                    read_data["alt_depth"],
                                    read_data["alt_fr"],
                                    read_data["ref_fr"],
                                    fisher_strand,
                                    cons.get("hgvsc", ""),
                                    cons.get("hgvsp", ""),
                                    cons.get("biotype", ""),
                                    cons.get("impact", ""),
                                    ", ".join(demographics.get("clin_sig", ())),
                                    ", ".join(cons.get("consequence_terms", ())),
                                    cons.get("sift_prediction", ""),
                                    cons.get("polyphen_prediction", ""),
                                    "{:.10f}".format(demographics["maf"]) if "maf" in demographics else "",
                                    ", ".join(sorted(all_genes - set([gene_symbol]))),
                                    ", ".join(demographics.get("dbsnp", ())),
                                    ", ".join(demographics.get("hgmd", ())),
                                    ", ".join(demographics.get("cosmic", ())),
                                    ", ".join(demographics.get("pubmed", ()))])
    
    
    os.unlink(vepjson)
    annotations.sort(key=lambda r:(chrom2int(r[2]), r[2], int(r[3]), r[4]))
    with open(output, "wt") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Gene", "Transcript", "Chrom", "Pos", "Change", "Quality", "Filters", "VAF", "Depth", "Alt Depth", "Alt Depth F:R", "Ref Depth F:R",
                         "FisherStrand", "HGVSc", "HGVSp", "Biotype", "Impact", "Clinical Significance (Pubmed)", "Consequences", "Sift", "Polyphen", "MAF",
                         "Other Genes", "dbSNP",  "HGMD", "COSMIC", "Pubmed"])
        writer.writerows(annotations)
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', help="Input vcf file.")
    parser.add_argument("-v", "--vep", help="Directory containing vep data.", required=True)
    parser.add_argument("-o", "--output", help="Output annotated tsv.", default=argparse.SUPPRESS)
    parser.add_argument("-r", "--reference", help="Fasta that the sample was aligned against.", default=argparse.SUPPRESS) # Make required once vep 101 code  no longer required
    parser.add_argument("-p", "--panel", help="Directory containing panel data.", default=argparse.SUPPRESS)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", type=int, default=argparse.SUPPRESS)
    parser.add_argument("-b", "--buffer-size", help="Number of variants to read into memory simultaneously. " + \
                                                    "Only needed if vep is being killed for running out of memory!", type=int, default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        annotate_panel(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

