import os
import sys
from pipeline.aws import s3_put
import pdb

from aireal.basespace.basespace import Session
token = os.getenv("TOKEN")

if not token:
    sys.exit(1)
bs = Session(token)


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

samples = [
"HNC027-HNC027-c-0",
"HNC035-HNC035-c-0",
"HNC024-HNC024-c-2",
"HNC027-HNC027-c-2",
"HNC021-HNC021-c-1",
"HNC021-HNC021-c-2",
"SCAN403-SCAN403-c-0",
]


output_dir = "/home/ubuntu/ephemoral"

for name in samples:
    for sample in bs.search("samples", query='(name="{}")'.format(name)):
        print("*****", name)
        
        for fastq in bs.sample_fastqs(sample["Id"]):
            print(fastq["Name"])
            path = os.path.join(output_dir, "temp.fastq")
            with open(path, "wb") as f_out:
                bs.download_fileobj(fastq["Id"], f_out)
            os.rename(path, os.path.join(output_dir, fastq["Name"]))
    break
            
            
files = os.listdir(output_dir)
for fn in files:
    path = os.path.join(output_dir, fn)
    sample = fn.split("_")[0]
    print(sample)
    s3_put("omdc-data", path, f"projects/head_and_neck/samples2/{sample}")




