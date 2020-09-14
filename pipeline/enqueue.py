from collections import defaultdict
import json
import uuid
import argparse

from pipeline.aws import s3_list
import boto3



BUCKET = "omdc-data"


most_important = [
"HNC006-HNC006-c-3",
"HNC028-HNC028-c-0",
"HNC031-HNC031-c-1",
"HNC024-HNC024-c-2-r",
"HNC038-HNC038-c-0",
"HNC038-HNC038-c-1",
"SCAN024-SCAN024-c-0",
"HNC002-HNC002-c-1",
"HNC033-HNC033-c-0",
"HNC002-HNC002-c-2",
"HNC033-HNC033-c-1",
"HNC024-HNC024-c-0-r",
"HNC022-HNC022-c-0-r",
"HNC025-HNC025-c-1",
"HNC013-HNC013-c-2",
"HNC017-HNC017-c-2",
"HNC019-HNC019-c-3",
"HNC025-HNC025-c-2",
"HNC038-HNC038-c-2",
"HNC018-HNC018-c-1",
"HNC018-HNC018-c-2",
"HNC012-HNC012-c-2",
"HNC014-HNC014-c-2",
"HNC022-HNC022-c-1",
"HNC022-HNC022-c-2",
"HNC027-HNC027-c-1",
"HNC028-HNC028-c-1",
"HNC032-HNC032-c-1"
"HNC034-HNC034-c-0",
]



def enqueue():
    project = "head_and_neck"
    
    sqs = boto3.client("sqs")
    queue_url = sqs.get_queue_url(QueueName="samples")["QueueUrl"]

    complete = set()
    for key in s3_list(BUCKET, f"projects/{project}/analyses", extension=".bam"):
        splitkey = key.split("/")
        if len(splitkey) >= 2:
            complete.add(splitkey[-2])

    fastqs = defaultdict(list)
    for key in s3_list(BUCKET, f"projects/{project}/samples", extension=".fastq.gz"):
        sample = key.split("/")[-1].split("_")[0]
        
        if "-c-" in sample:
            continue
        
        if sample not in complete:
            fastqs[sample] += [f"s3://{BUCKET}/{key}"]
    
    n = 0
    for sample, urls in sorted(fastqs.items()):
        data = {"Script": "cfpipeline",
                "Output": f"s3://{BUCKET}/projects/{project}/analyses/{sample}",
                "Args": urls,
                "Kwargs": {"--sample": sample,
                           "--reference": f"s3://{BUCKET}/reference/37/sequence/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set_plus_hpv_panel.tar.gz",
                           "--panel": f"s3://{BUCKET}/panels/Head_and_Neck.tar.gz",
                           "--vep": f"s3://{BUCKET}/reference/homo_sapiens_refseq_vep_98_GRCh37.tar.gz",
                           "--umi": "thruplex",
                           "--min-family-size": "2"},
                }
        
        message = json.dumps(data)
        print(message, "\n")
        sqs.send_message(QueueUrl=queue_url,
                         MessageBody= message)
        n += 1
        break
    print(f"{n} messages queued.")


if __name__ == "__main__":
    #parser = argparse.ArgumentParser()
    #parser.add_argument('project', help="AWS project name.")
    #parser.add_argument("-p", "--panel", help="Panel name.", required=True)
    #args = parser.parse_args()
    #enqueue(**vars(args))
    enqueue()
            
            
            
