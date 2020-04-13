import os, csv, pdb
from collections import defaultdict, Counter


def variant(row):
    return "{}:{} {}".format(row["Chrom"], row["Pos"], row["Change"])



PATH = "/home/ed/Data/accept"
with open(os.path.join(PATH, "accept.artifacts.tsv")) as f_in:
    reader = csv.DictReader(f_in, delimiter="\t")
    artifacts = set(row["Variant"] for row in reader)



subjects = defaultdict(dict)
for fn in os.listdir(os.path.join(PATH, "annotation")):
    subject, _, _, timepoint = fn[:-15].split("-")
    subjects[subject][timepoint] = os.path.join(PATH, "annotation", fn)


mean_vafs = defaultdict(lambda:defaultdict(list))
with open(os.path.join(PATH, "accept_vafs.tsv"), "wt") as f_out:
    writer = csv.writer(f_out, delimiter="\t")
    writer.writerow(["Subject", "Timepoint", "Variant", "VAF"])
    
    for subject, timepoints in subjects.items():
        hethom = Counter()
        vafs = defaultdict(lambda:defaultdict(lambda:"0"))
        for timepoint, fn in timepoints.items():
            with open(fn) as f_in:
                reader = csv.DictReader(f_in, delimiter="\t")
                for row in reader:
                    var = variant(row)
                    if row["Zygosity"] != "unknown" and timepoint != "0":
                        hethom[var] += 1
                    vafs[var][timepoint] = row["VAF"]

        for var, timepoints in vafs.items():
            ref, alt = var.split()[1].split("/")
            indel = "-" in var or len(ref) != len(alt)
            if var not in artifacts and hethom.get(var, 0) != 2:
                for timepoint in ("0", "1", "2"):
                    vaf = timepoints[timepoint]
                    writer.writerow([subject, timepoint, var, vaf])
                    if vaf != "0":
                        mean_vafs[subject][timepoint] += [float(vaf)]
        

with open(os.path.join(PATH, "accept_mean_vafs.tsv"), "wt") as f_out:
    writer = csv.writer(f_out, delimiter="\t")
    writer.writerow(["Subject", "Timepoint", "Mean_VAF"])
    for subject, data in mean_vafs.items():
        for timepoint in ("0", "1", "2"):
            z = data[timepoint]
            writer.writerow([subject, timepoint, sum(z) / len(z) if z else 0])











