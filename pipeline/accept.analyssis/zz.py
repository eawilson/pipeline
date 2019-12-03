import csv
from collections import defaultdict
import pdb



com = {}
with open("accept.annotation.combined.tsv", "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        com[(row["Patient"], row["Variant"])] = row
    
rev = set()
with open("accept.annotation.combined.reviewed.tsv", "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rev.add((row["Patient"], row["Variant"]))

with open("accept.annotation.combined.reviewed.tsv", "at") as f:
    writer = csv.DictWriter(f, reader.fieldnames, delimiter="\t")
    writer.writeheader()
    for key in set(com) - rev:
        writer.writerow(com[key])