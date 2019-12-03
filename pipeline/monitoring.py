import csv, pdb
from collections import defaultdict, Counter



missing_patients = set()
patients = defaultdict(lambda:defaultdict(lambda:["", "", "", ""]))
with open("accept.annotation.combined.reviewed.tsv", "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        missing_patients.add(row["Patient"])
        if row["Relevant"] == "Y":
            patients[row["Patient"]][(row["Gene"], row["HGVSc"])][int(row["Timepoint"])] = row["VAF"]

with open("accept.monitoring.tsv", "wt") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Patient", "Gene", "Variant", "Timepoint 0", "Timepoint 1", "Timepoint 2", "Timepoint 3"])
    for patient, variants in patients.items():
        for variant, vafs in sorted(variants.items()):
            vafs = ["{:.3f}".format(float(vaf)) if vaf else "" for vaf in vafs]
            writer.writerow([patient] + list(variant) + vafs)
        writer.writerow([])

with open("accept.monitoring.long.tsv", "wt") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Patient", "Gene", "Variant", "Timepoint", "VAF"])
    for patient, variants in patients.items():
        for variant, vafs in sorted(variants.items()):
            for timepoint, vaf in enumerate(["{:.3f}".format(float(vaf)) if vaf else "0" for vaf in vafs]):
                writer.writerow([patient] + list(variant) + [timepoint, vaf])

    for patient in missing_patients:
        if patient not in patients:
            writer.writerow([patient, "", "",  "0", "0"])



pathways = {}
with open("PMAL.tsv", "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        pathways[row["GENE"]] = row["PATHWAY"]

counts = defaultdict(Counter)
for patient, variants in patients.items():
    for variant, vafs in sorted(variants.items()):
        counts[patient][variant[0]] += 1

rows = []
for patient, genes in counts.items():
    for gene, count in sorted(genes.items()):
        rows += [[str(patient), gene, pathways[gene], count]]
            
            
with open("accept.pathways.tsv", "wt") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Patient", "Gene", "Pathway", "Count"])
    for row in sorted(rows, key=lambda r:(r[2], r[1])):
        writer.writerow(row)


genelist = set()
for patient, genes in counts.items():
    for gene, count in sorted(genes.items()):
        genelist.add(gene)
genelist = sorted(genelist)
patientlist = sorted(counts)

with open("accept.heatmap.tsv", "wt") as f:
    writer = csv.writer(f, delimiter="\t")
    for gene in genelist:
        writer.writerow([gene, pathways[gene]] + [counts[patient].get(gene, "") for patient in patientlist])
    writer.writerow(["", ""] + patientlist)





#with open("accept.heatmap.tsv", "wt") as f:
    #writer = csv.writer(f, delimiter="\t")
    #writer.writerow(["Patient", "Gene", "Pathway", "Count"])
    #for row in sorted(rows, key=lambda r:(r[2], r[1])):
        #writer.writerow(row)






