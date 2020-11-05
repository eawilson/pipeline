from collections import defaultdict
import json
import uuid
import argparse

from pipeline.aws import s3_list
import boto3



BUCKET = "omdc-data"


def enqueue(dry_run=False):
    project = "scan"
    panel = "Pan_Cancer"
    analyses = "analyses2"
    
    sqs = boto3.client("sqs")
    queue_url = sqs.get_queue_url(QueueName="samples")["QueueUrl"]

    complete = set()
    for key in s3_list(BUCKET, f"projects/{project}/{analyses}", extension=".bam"):
        splitkey = key.split("/")
        if len(splitkey) >= 2:
            complete.add(splitkey[-2])

    fastqs = defaultdict(list)
    for key in s3_list(BUCKET, f"projects/{project}/samples", extension=".fastq.gz"):
        sample = key.split("/")[-2] or key.split("/")[-3]
        
        if sample not in complete:
            fastqs[sample] += [f"s3://{BUCKET}/{key}"]
    
    n = 0
    for sample, urls in sorted(fastqs.items()):
        data = {"Script": "cfpipeline",
                "Output": f"s3://{BUCKET}/projects/{project}/{analyses}/{sample}",
                "Args": urls,
                "Kwargs": {"--sample": sample,
                           "--reference": f"s3://{BUCKET}/reference/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set.tar.gz",
                           "--panel": f"s3://{BUCKET}/panels/{panel}.tar.gz",
                           "--vep": f"s3://{BUCKET}/reference/vep_101_GRCh37_homo_sapiens_refseq.tar",
                           "--min-family-size": "2",
                           "--umi": "thruplex",
                           }
                }
        
        message = json.dumps(data)
        print(message, "\n")
        if not dry_run:
            sqs.send_message(QueueUrl=queue_url, MessageBody=message)
            n += 1

    print(f"{n} messages queued.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', "--dry-run", action="store_const", const=True, default=argparse.SUPPRESS)
    args = parser.parse_args()
    enqueue(**vars(args))
            
            
            
