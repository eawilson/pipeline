import json, pdb, os, csv
from collections import defaultdict
from itertools import chain

from covermi.gr import load_principal, load_targets, SequencedVariant, Chrom
from covermi import Panel, Gr

__all__ = ["create_report"]

BIOTYPE = defaultdict(lambda:1, (("protein_coding", 0), ("pseudogene", 2)))


def transcript_version_sort(cons):
    transcript_id = cons["transcript_id"]
    version = transcript_id.split(".")[-1]
    return (-int(version) if version.isnumeric() else 0, transcript_id)



def create_report(vep_json_file, panel):

    transcript_ids, gene_symbols, _ = load_targets(panel.targets if "targets" in panel else None)
    transcript_ids = set(transcript_id.split(".")[0] for transcript_id in transcript_ids)
    principal = load_principal(panel.principal if "principal" in panel else None)

    def refseq_sort(cons):
        return [BIOTYPE[cons["biotype"]], int(cons["gene_id"]), \
                not(cons["transcript_id"].startswith("N")), not(principal[cons["transcript_id"]]),"canonical" not in cons, cons["transcript_id"]]

    def ensembl_sort(cons):
        return [BIOTYPE[cons["biotype"]], int(cons["gene_id"]), not(principal[cons["transcript_id"]]),"canonical" not in cons, cons["transcript_id"]]

    if not isinstance(panel, Panel):
        panel = Panel(panel)
    consequence_sort = refseq_sort if panel.properties.get("transcript_source", "refseq") == "refseq" else ensembl_sort
        
    annotations = []
    with open(vep_json_file) as f:
        for line in f:
            vep_output = json.loads(line)

            # create variant
            row = vep_output["input"].strip().split("\t")
            try:
                chrom = Chrom(row[0])
            except KeyError:
                #print("Unknown contig {} skipping.".format(row[0]))
                continue
            
            formatdict = dict(zip(row[8].split(":"), row[9].split(":")))
            filters = ";".join(code for code in row[6].split(";") if code != "LowVariantFreq")
            if filters == "":
                filters = "PASS"
            variant = SequencedVariant(chrom, int(row[1]), row[3], row[4], name=row[2])
            try:
                variant.qual = float(row[5])
            except ValueError:
                variant.qual = "."
            variant.filters = filters
            variant.depth = int(formatdict["AD"]) + int(formatdict["RD"])
            variant.alt_depth = int(formatdict["AD"])
            if panel.properties["homo_ll"] <= variant.vaf:
                variant.zygosity = "hemizygous" if chrom in (Chrom(23), Chrom(24), Chrom(25)) else "homozygous"
            elif panel.properties["hetero_ll"] <= variant.vaf <= panel.properties["hetero_ul"]:
                variant.zygosity = "heterozygous"
            else:
                variant.zygosity = "unknown"

            if not Gr(variant).touched_by(panel.amplicons):
                #print("Offtarget skipping.")
                continue

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
                        
                    try:
                        frequencies = colocated["frequencies"][variant.alt or "-"]
                    except KeyError:
                        continue
                        
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
                
            consequences = vep_output.get("transcript_consequences", ())
            consequences = [cons for cons in consequences if "gene_symbol" in cons and "transcript_id" in cons]
            
            if not consequences:
                #print("No transcript consequences skipping.")
                continue
            
            gene_cons = [cons for cons in consequences if cons["gene_symbol"] in gene_symbols]
            transcript_cons = [cons for cons in consequences if cons["transcript_id"].split(".")[0] in transcript_ids]
            consequences = set(gene_cons) & set(transcript_cons) if (gene_cons and transcript_cons) else gene_cons or transcript_cons or consequences
            
            deduplicate = defaultdict(list)
            for cons in consequences:
                deduplicate[cons["transcript_id"].split(".")[0]] += [cons]
            
            consequences = [sorted(values, key=transcript_version_sort)[0] for values in deduplicate.values()]
            consequence = sorted(consequences, key=consequence_sort)[0]
            #print(consequence)
            annotation = {"chrom": variant.chrom, 
                          "pos": int(vep_output["start"]), 
                          "allele_string": vep_output["allele_string"],
                          "quality": variant.qual,
                          "filters": variant.filters,
                          "alt_depth": variant.alt_depth,
                          "depth": variant.depth,
                          "vaf": variant.vaf,
                          "zygosity": variant.zygosity,
                          "depth": variant.depth,
                          "alt_depth": variant.alt_depth,
                          "gene_symbol": consequence.get("gene_symbol", ""),
                          "transcript_id": consequence.get("transcript_id", ""),
                          "hgvsc": consequence.get("hgvsc", ""),
                          "hgvsp": consequence.get("hgvsp", ""),
                          "biotype": consequence["biotype"],
                          "impact": consequence["impact"],
                          "consequence_terms": consequence["consequence_terms"],
                          "sift_score": consequence.get("sift_score", ""),
                          "sift_prediction": consequence.get("sift_prediction", ""),
                          "polyphen_score": consequence.get("polyphen_score", ""),
                          "polyphen_prediction": consequence.get("polyphen_prediction", ""),
                          "other_genes": sorted(set(cons.get("gene_symbol", "") for cons in vep_output["transcript_consequences"]) - set([consequence["gene_symbol"], ""])),
                          **demographics}
            annotations += [annotation]


    annotation_file = "{}.annotation.tsv".format(os.path.splitext(vep_json_file)[0])
    with open(annotation_file, "wt") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Gene", "Transcript", "Chrom", "Pos", "Change", "Quality", "Filters", "Depth", "Alt Depth", "VAF", "Zygosity", "HGVSc", "HGVSp", "Impact", \
                            "Clinical Significance (Pubmed)", "Consequences", "Sift", "Polyphen", "MAF", "COSMIC", "dbSNP", "Pubmed",  "HGMD", "Other_Genes"])
        for a in sorted(annotations, key=lambda a:(a["chrom"], a["pos"], a["allele_string"])):
           writer.writerow([a["gene_symbol"],
                             a["transcript_id"],
                             a["chrom"],
                             a["pos"],
                             a["allele_string"],
                             a["quality"],
                             a["filters"],
                             a["depth"],
                             a["alt_depth"],
                             a["vaf"],
                             a["zygosity"],
                             a["hgvsc"],
                             a["hgvsp"],
                             a["impact"],
                             ", ".join(a.get("clin_sig", ())),
                             ", ".join(a["consequence_terms"]),
                             a["sift_prediction"],
                             a["polyphen_prediction"],
                             a.get("maf", ""),
                             ", ".join(a.get("cosmic", ())),
                             ", ".join(a.get("dbsnp", ())),
                             ", ".join(a.get("pubmed", ())),
                             ", ".join(a.get("hgmd", ())),
                             ", ".join(a.get("other_genes", ())),
                             ])
    return annotation_file
    

