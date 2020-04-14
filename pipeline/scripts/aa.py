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

os.chdir("/home/ubuntu/ephemoral")

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
    
        with open(f"{name}.cfpipeline.txt", "wt") as f:
            subprocess.run(" ".join(["python3 /home/ubuntu/pipeline/pipeline/scripts/cfpipeline.py"] + fastqs +
                            ["--panel", "Head_and_Neck",
                            "--genome", "GRCh37",
                            "--min-family-size", "3"]),
                            stderr=f, check=True, shell=True)
        
        aws_prefix = f"projects/head_and_neck/{name}"
        for fn in os.listdir("."):
            if not os.path.isdir(fn):
                if not (fn.endswith(".fastq") or fn.endswith(".fastq.gz")):
                    s3_put("omdc-data", fn, aws_prefix)
                os.unlink(fn)





