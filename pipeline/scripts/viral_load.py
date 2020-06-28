import os, pdb
from collections import defaultdict

from pipeline import s3_list, s3_get, s3_put, ungzip_and_combine_illumina_fastqs, run, pipe, s3_exists


BUCKET = "omdc-data"
REFERENCE = "reference/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set_plus_hpv_panel.fna"


def main():
    threads = run(["getconf", "_NPROCESSORS_ONLN"]).stdout.strip()
    samples = defaultdict(list)
    for key in s3_list(BUCKET, "projects/head_and_neck/samples", ".fastq.gz"):
        sample = key.split("/")[-1].split("_")[0]
        samples[sample] += [key]
        
    for sample, keys in samples.items():
        prefix = f"projects/head_and_neck/{sample}"
        tsv_file = f"{sample}.reads_by_contig.tsv"
        if s3_exists(BUCKET, f"{prefix}/{tsv_file}"):
            continue
        
        for fastq in os.listdir("."):
            if fastq.endswith(".fastq"):
                raise RuntimeError(f"Unexpected fastq {fastq} in working directory.")
        
        gzipped_fastqs = []
        for n, key in enumerate(sorted(keys)):
            fastq = key.split("/")[-1]
            gzipped_fastqs += [fastq]
            print(f"Downloading {fastq}")
            s3_get(BUCKET, key, fastq)

        fastqs = ungzip_and_combine_illumina_fastqs(gzipped_fastqs, sample)
        for fastq in gzipped_fastqs:
            os.unlink(fastq)

        dude_options = ["-a", "3",
                        "-f", "1",
                        "-t"]
        pipe(["dude"] + fastqs + dude_options)
        deduped_fastqs = []
        for fastq in fastqs:
            deduped_fastqs += ["{}.deduped.fastq".format(fastq[:-6])]
            os.unlink(fastq)
        
        sam_file = f"{sample}.sam"
        with open(sam_file, "wb") as f_out:
            pipe(["bwa", "mem", "-t", threads, REFERENCE] + deduped_fastqs, stdout=f_out)
        for fastq in deduped_fastqs:
            os.unlink(fastq)
            
        pipe(["sam_reads_by_contig", sam_file, tsv_file])
        os.unlink(sam_file)
        s3_put(BUCKET, tsv_file, f"{prefix}/{tsv_file}")
        

if __name__ == "__main__":
    main()

