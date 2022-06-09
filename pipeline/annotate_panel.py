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
IMPACT = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
REFSEQ = defaultdict(int, (("NM", 4), ("NR", 3), ("XM", 2), ("XR", 1)))
DELETE_NON_DIGIT = str.maketrans("", "", "ABCDEFGHIJKLMNOPQRSTUVWXYZ_")

CHROM = 0
QUAL = 5
FILTERS = 6
FMT_KEYS = 8
FMT_VALS = 9


def vardict_read_data(row):
    fmt = dict(zip(row[FMT_KEYS].split(":"), row[FMT_VALS].split(":")))
    ref_fr = [int(n) for n in fmt["RD"].split(",")]
    alt_fr = [int(n) for n in fmt["ALD"].split(",")]
    return {"vaf": fmt["AF"], "depth": fmt["DP"], "alt_depth": fmt["VD"], "ref_fr": ref_fr, "alt_fr": alt_fr}



def varscan2_read_data(row):
    fmt = dict(zip(row[FMT_KEYS].split(":"), row[FMT_VALS].split(":")))
    vaf = float(fmt["FREQ"].rstrip("%"))/100
    ref_fr = [int(fmt["RDF"]), int(fmt["RDR"])]
    alt_fr = [int(fmt["ADF"]), int(fmt["ADR"])]
    return {"vaf": f"{vaf:.4f}", "depth": fmt["DP"], "alt_depth": fmt["AD"], "ref_fr": ref_fr, "alt_fr": alt_fr}



def mutect2_read_data(row):
    fmt = dict(zip(row[FMT_KEYS].split(":"), row[FMT_VALS].split(":")))
    vaf = float(fmt["AF"])
    alt_depth = fmt["AD"].split(",")[1]
    d4 = fmt["SB"].split(",")
    ref_fr = [int(n) for n in d4[:2]]
    alt_fr = [int(n) for n in d4[2:]]
    return {"vaf": f"{vaf:.4f}", "depth": fmt["DP"], "alt_depth": alt_depth, "ref_fr": ref_fr, "alt_fr": alt_fr}



def chrom2int(chrom):
    try:
        return int(chrom[3:])
    except ValueError:
        return math.inf



def parse_colocated(vep_output):
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
                
            frequencies = colocated.get("frequencies")
            if frequencies:
                frequencies = frequencies.get(vep_output["allele_string"].split("/")[1], {})
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
    
    return demographics



# reference should be a required argument as fails in vep 104 without but is made optional here for compatibility with legacy code in cfpipeline.py for vep 101
def annotate_panel(vcf, vep, reference=None, threads=None, output="", panel="", buffer_size=None):
    if threads is None:
        threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    
    if not output:
        output = "."
    if os.path.isdir(output):
        output = os.path.join(output, "{}.annotation.tsv".format(vcf[:-4] if vcf.endswith(".vcf") else vcf))
    
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
        reference = (glob.glob(f"{reference}/*.fna") + glob.glob(f"{reference}/*.fa") + glob.glob(f"{reference}/*.fasta") + [reference])[0]
        vep_options += ["--fasta", reference]
    if int(threads) > 1:
        vep_options += ["--fork", threads]
    if "refseq" in vep:
        vep_options += ["--refseq"]
    if buffer_size is not None:
        vep_options += ["--buffer_size", buffer_size]
    
    pipe(["vep", "-i", vcf, "-o", vepjson] + vep_options)
    
    
    get_read_data = None
    with open(vcf, "rt") as f:
        for row in f:
            if not row.startswith("#"):
                break
            if row.startswith("##source="):
                source = row[9:].strip()
                #if source == "strelka":
                if source.startswith("VarDict"):
                    get_read_data = vardict_read_data                    
                elif source == "VarScan2":
                    get_read_data = varscan2_read_data
                elif source == "Mutect2":
                    get_read_data = mutect2_read_data                    
            headings = row
    
    if get_read_data is None:
        sys.exit(f"Unsupported variant caller {source}")
    if len(headings.split("\t")) > 10:
        sys.exit("Multi-sample vcfs not suppored")
    
    
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
            return [transcript in needed_transcripts,
                    cons["gene_symbol"] in needed_genes,
                    BIOTYPE[cons["biotype"]],
                    REFSEQ[prefix],
                    -int(cons["gene_id"]),
                    principal.get(transcript, 0),
                    "canonical" in cons,
                    -int(major),
                    int(minor)]
    
    else: # ensembl transcripts
        def consequence_sort(cons):
            # Version numbers not in vep as of version 101, but who knows the future ...
            transcript = cons["transcript_id"]
            return [transcript in needed_transcripts,
                    cons["gene_symbol"] in needed_genes,
                    BIOTYPE[cons["biotype"]],
                    -int(cons["gene_id"].translate(DELETE_NON_DIGIT)),
                    principal.get(transcript, 0),
                    "canonical" in cons,
                    -int(transcript.translate(DELETE_NON_DIGIT))]
    
    annotations = []
    with open(vepjson) as f:
        for line in f:
            vep_output = json.loads(line)
            
            consequences = vep_output.get("transcript_consequences") 
            if consequences:
                cons = sorted(consequences, key=consequence_sort)[-1]
                other_genes = set(c["gene_symbol"] for c in consequences) - set([cons["gene_symbol"]])
            
            else:
                most_severe_consequence = vep_output["most_severe_consequence"]
                for cons in sorted(chain(*[v for k, v in vep_output.items() if k.endswith("_consequences")]), key=lambda x:x.get("biotype", ""), reverse=True):
                    # We are only going to use biotype and impact so probably does not matter which one we choose so long as we are consistent
                    if most_severe_consequence in cons["consequence_terms"]:
                        break
                other_genes = ()
            
            row = vep_output["input"].rstrip().split("\t")
            read_data = get_read_data(row)
            
            if read_data["alt_depth"] == "0":
                continue
            
            # https://gatk.broadinstitute.org/hc/en-us/articles/360035532152-Fisher-s-Exact-Test
            exact = fisher_exact([read_data["ref_fr"], read_data["alt_fr"]])[1]
            if exact:
                fisher_strand = "{:.1f}".format(-10 * math.log10(exact))
            else:
                fisher_strand = ""
            
            demographics = parse_colocated(vep_output)

            annotations.append([cons.get("gene_symbol", ""),
                                cons.get("transcript_id", ""),
                                row[CHROM],
                                vep_output["start"],
                                vep_output["allele_string"],
                                row[QUAL],
                                row[FILTERS],
                                read_data["vaf"],
                                read_data["depth"],
                                read_data["alt_depth"],
                                "{}:{}".format(*read_data["alt_fr"]),
                                "{}:{}".format(*read_data["ref_fr"]),
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
                                ", ".join(sorted(other_genes)),
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

