
import os
import glob

from pipeline import s3_put, s3_list



for bsrun in ["SCAN002"]:
    stem = f"basespace/Projects/{bsrun}/Samples"
    for sample in os.listdir(stem):
        if sample.startswith(".") or sample.endswith("(2)"):
            continue
        
        if not sample.endswith("Combined"):
            continue
        
        real_sample = sample[:-8].upper()

        print(sample)
        prefix = f"projects/scan/samples/{real_sample}/"
        if not s3_list("omdc-data", prefix):
            for fastq in glob.glob(f"{stem}/{sample}/Files/*fastq.gz"):
                print(fastq)
                s3_put("omdc-data", fastq, prefix=prefix)

