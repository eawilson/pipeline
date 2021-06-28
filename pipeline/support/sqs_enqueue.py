from collections import defaultdict
import json
import argparse
import re
import csv
import sys
import pdb

from pipeline.aws import s3_list, s3_exists
import boto3


def pop(arguments, *args, nargs=2, required=False):
    for arg in args:
        if arg in arguments:
            i = arguments.index(arg)
            for n in range(0, nargs):
                try:
                    val = arguments.pop(i)
                except IndexError:
                    sys.exit(f"sqs_enqueue: {arg} is invalid")
            return val
    
    if required:
        sys.exit(f"sqs_enqueue: {arg} is a required argument")



def parse_url(url):
    if url[:5].lower() != ("s3://"):
        sys.exit(f"sqs_enqueue: {url} is not a valid s3 url.")
    parts = url.split("/")
    if len(parts) < 4:
        sys.exit(f"sqs_enqueue: {url} is not a valid s3 url.")
    bucket = parts[2]
    key = "/".join(parts[3:])
    return (bucket, key)



def main():
    if len(sys.argv) < 3:
        sys.exit("sqs_enqueue: Too few arguments")
    args = sys.argv[3:]
    queue_name = sys.argv[1]
    command = sys.argv[2]    
    
    dry_run =  pop(args, "-d", "--dry-run", nargs=1)
    input_url = pop(args, "-i", "--input", required=True)
    output_url = pop(args, "-o", "--output", required=True)
    panel = pop(args, "-p", "--panel", required=True)
    reference = pop(args, "-r", "--reference") or "s3://omdc-data/reference/GRCh37_EBV_HPV_bwa_mem2.tar.gz"
    vep = pop(args, "-v", "--vep") or "s3://omdc-data/reference/vep_104_GRCh37_homo_sapiens_refseq.tar"
    
    if not input_url.endswith("/"):
        input_url = f"{input_url}/"
    if not output_url.endswith("/"):
        output_url = f"{output_url}/"
    if not panel.lower().startswith("s3://"):
        panel = f"s3://omdc-data/panels/{panel}"
    if not panel.endswith(".tar.gz"):
        panel = f"{panel}.tar.gz"
    if not reference.endswith(".tar.gz"):
        reference = f"{reference}.tar.gz"
    
    sqs = boto3.client("sqs")
    try:
        queue_url = sqs.get_queue_url(QueueName=queue_name)["QueueUrl"]
    except sqs.exceptions.QueueDoesNotExist:
        sys.exit(f'sqs_enqueue: "{queue_name}" is not a valid sqs queue')
    
    bams = set()
    bucket, stem = parse_url(output_url)
    for key in s3_list(bucket, stem, extension=".bam"):
        identifier = "/".join(key[len(stem):].split("/")[-3:-1])
        bams.add(identifier)
    
    fastqs = defaultdict(list)
    bucket, stem = parse_url(input_url)
    for key in s3_list(bucket, stem, extension=".fastq.gz"):
        identifier = "/".join(key[len(stem):].split("/")[-3:-1])
        if identifier not in bams:
            fastqs[identifier] += [f"s3://{bucket}/{key}"]
    
    
    n = 0
    for identifier, samples in sorted(fastqs.items()):
        cmd = [command] + sorted(samples) + args + ["--name", identifier.split("/")[-1], "--output", f"{output_url}{identifier}", "--reference", reference, "--vep", vep, "--panel", panel]        
        message = json.dumps(cmd)
        print(message, "\n", file=sys.stderr)
        n += 1
        if not dry_run:
            sqs.send_message(QueueUrl=queue_url, MessageBody=message)

    print("sqs_enqueue: {} messages {}.".format(n, "processed" if dry_run else "queued"), file=sys.stderr)



if __name__ == "__main__":
    main()
