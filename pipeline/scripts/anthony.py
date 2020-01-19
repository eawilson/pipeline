import csv
import sys
import pdb
from collections import defaultdict

from pipeline.basespace import BasespaceSession


def main():
    tokens = sys.argv[1:]
    
    with open("/home/ed/Desktop/hncsamples_sequenced.tsv", "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    for index, row in enumerate(rows):
        for token in tokens:
            bs = BasespaceSession(token)
            results = bs.search("samples", experimentname=row["Run"], name=row["Basespace_Sample"])
            if results:
                results.sort(key=lambda r: r["DateCreated"])
                print(results[-1]["Name"], f"{index} of {len(rows)}")
                fastqs = bs.sample_fastqs(results[-1]["Id"])
                for fastq in fastqs:
                    print(fastq["Name"])
                    key = f"projects/head_and_neck/{row["Patient"]}/samples/{fastq["Name"]}"
                    bs.copy_to_s3(fastq["Id"], "omdc-data", key)
                break
                
        else:
            print("MISSING", row["Run"], row["Basespace_Sample"])
            #raise RuntimeError("Missing sample {}.".format(row["Basespace_Sample"]))



if __name__ == "__main__":
    main()



















