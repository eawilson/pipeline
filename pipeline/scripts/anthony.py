import csv
import sys
import pdb

from pipeline.basespace import BasespaceSession
from pipeline.aws_lambda import multipart_upload, invoke
from pipeline import s3_object_exists



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
                sample = results[-1]
                print(sample["Name"], f"{index} of {len(rows)}")
                fastqs = bs.sample_fastqs(results[-1]["Id"])
                for fastq in fastqs:
                    print(fastq["Name"])
                    key = f'fastqs/{row["Run"]}/{fastq["Name"]}'
                    if not s3_object_exists("omdc-data", key):
                        url = bs.file_url(fastq["Id"])
                        multipart_upload({"bucket": "omdc-data", "key": key, "url": url})
                break
        
        else:
            print("MISSING", row["Run"], row["Basespace_Sample"])
            #raise RuntimeError("Missing sample {}.".format(row["Basespace_Sample"]))



if __name__ == "__main__":
    main()



















