import os
import sys
import time
import datetime
import argparse
import pdb
from collections import defaultdict
import csv



from pipeline import run, mount_instance_storage, load_panel_from_s3, \
    ungzip_and_combine_illumina_fastqs, s3_put, illumina_readgroup, \
    pipe, create_report, s3_object_exists, s3_open, s3_list_keys, am_i_an_ec2_instance, s3_get
from pipeline.scripts.cfpipeline import cfpipeline


def runner(project, panel, input_csv, genome=None):
    panel_name = panel
    input_csv = os.path.abspath(input_csv)

    if am_i_an_ec2_instance():
        os.chdir(mount_instance_storage())
        panel = load_panel_from_s3(panel)
        genome = panel.properties["reference_fasta"]
        
    if genome is None:
        raise RuntimeError("No reference genome supplied.")

    samples = []
    with open(input_csv, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            samples += [row]
            
    for sample in samples:
        if "-c'" in sample["Sample"]:
            continue
                
        print("Sample {}.".format(sample["Sample"]))
        if s3_object_exists("omdc-data", "projects/{}/{}".format(project, sample["Sample"])):
            print("Exists. skipping.".format(project, sample))
            continue
        
        print("Fetching fastqs...")
        fastqs = []
        s3_fastqs = s3_list_keys("omdc-data", "projects/{}/samples/{}/{}". \
                        format(project, sample["Patient"], sample["Sample"]))
        for s3_fastq in s3_fastqs:
            fastq = s3_fastq.split("/")[-1]
            fastqs += [fastq]
            s3_get("omdc-data", s3_fastq, fastq)

        with open("cfpipeline.txt", "wt") as f:
            stderr = sys.stderr
            sys.stderr = f
            try:
                print("Running pipeline.")
                cfpipeline(fastqs, panel_name, genome, min_family_size=1)
            finally:
                sys.stderr = stderr

        for fastq in fastqs:
            os.unlink(fastq)

        print("Uploading to s3.")
        for filename in os.listdir():
            if os.path.isfile(filename):
                s3_put("omdc-data", filename, prefix="projects/{}/{}".format(project, sample["Sample"]))
                os.unlink(filename)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", help="Sample list.")
    parser.add_argument("--project", help="Project name.")
    parser.add_argument("--panel", help="Panel name.")
    parser.add_argument("--genome", help="Reference genome.", default=None)
    args = parser.parse_args()
    runner(**vars(args))



