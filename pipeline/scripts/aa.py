import os
import sys
import argparse
import pdb
import csv
import subprocess
from pipeline.basespace import Session


from pipeline import run, mount_instance_storage, load_panel_from_s3, \
    ungzip_and_combine_illumina_fastqs, s3_put, illumina_readgroup, \
    pipe, create_report, s3_object_exists, s3_open, s3_list_keys, am_i_an_ec2_instance, s3_get
from pipeline.scripts.cfpipeline import cfpipeline


samples = [
"HNC027-HNC027-c-0",
"HNC035-HNC035-c-0",
"HNC020-HNC020-c-3",
"HNC003-HNC003-c-5",
"HNC024-HNC024-c-2",
"HNC018-HNC018-c-0",
"HNC027-HNC027-c-2",
"HNC021-HNC021-c-1",
"HNC019-HNC019-c-2",
"HNC015-HNC015-c-2",
"HNC021-HNC021-c-2",
"HNC009-HNC009-c-2",
"SCAN403-SCAN403-c-0",
]

token = os.getenv("TOKEN")
if not token:
    sys.exit(1)
bs = Session(token)

#os.chdir("/home/ubuntu/ephemoral")

for name in samples:
    for sample in bs.search("samples", query='(name="{}")'.format(name)):
        print("*****", name)
    
        fastqs = []
        for fastq in bs.sample_fastqs(sample["Id"]):
            print(fastq["Name"])
            with open("temp.fastq", "wb") as f_out:
                bs.download_fileobj(fastq["Id"], f_out)
            os.rename("temp.fastq", fastq["Name"])
            fastqs += [fastq["Name"]]
    
        with open("{}.cfpipeline.txt".format(sample), "wt") as f:
            subprocess.run(["cfpipeline"] + fastqs +
                            ["--panel", "Head_and_Neck",
                            "--genome", "GRCh37",
                            "--min-family-size", "3"],
                            stderr=f, check=True)
        
        break












for sample in samples:
    keys = s3_list_keys("omdc-data", "projects/head_and_neck/samples2")
    print(keys)







sys.exit()


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

    print("Running pipeline.")
    with open("{}.cfpipeline.txt".format(sample["Sample"]), "wt") as f:
        subprocess.run(["cfpipeline"] + fastqs +
                        ["--panel", panel_name,
                        "--genome", genome,
                        "--min-family-size", "1"],
                        stderr=f, check=True)

    for fastq in fastqs:
        os.unlink(fastq)

    print("Uploading to s3.")
    for filename in os.listdir():
        if os.path.isfile(filename):
            s3_put("omdc-data", filename, prefix="projects/{}/{}".format(project, sample["Sample"]))
            os.unlink(filename)













