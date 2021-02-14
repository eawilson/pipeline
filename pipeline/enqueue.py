from collections import defaultdict
import json
import argparse
import re
import csv
import sys

from pipeline.aws import s3_list
import boto3



def enqueue(script_args, bucket, project, panel, input="samples", output="analyses", samples=".*", manifest="", ignore_missing=False, dry_run=False):
    sample_regex = re.compile(samples)
    
    sqs = boto3.client("sqs")
    queue_url = sqs.get_queue_url(QueueName="samples")["QueueUrl"]

    complete = set()
    for key in s3_list(bucket, f"projects/{project}/{output}/", extension=".bam"):
        sample = key.split("/")[-1].split("_")[0]
        complete.add(sample)
    
    
    fastqs = defaultdict(list)
    for key in s3_list(bucket, f"projects/{project}/{input}/", extension=".fastq.gz"):
        sample = key.split("/")[-1].split("_")[0]
        if sample not in complete and sample_regex.fullmatch(sample):
            fastqs[sample] += [f"s3://{bucket}/{key}"]
    
    
    sample_args = {}
    if manifest:
        with open(manifest, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sample = row.pop("sample")
                if sample in sample_args:
                    sys.exit(f"Sample {sample} duplicated in manifest")
                args = []
                for key, val in row.items():
                    args.append(f"--{key}")
                    if val:
                        args.append(val)
                sample_args[sample] = args
        
        if not ignore_missing:
            for name in set(fastqs) - set(sample_args):
                print(f"Sample {name} is present in aws but not in the manifest", file=sys.stderr)
            for name in set(sample_args) - set(fastqs):
                print(f"Sample {name} is present in the manifest but not aws", file=sys.stderr)
            if set(sample_args) ^ set(fastqs):
                sys.exit("Missing samples. Aborting")
    
    n = 0
    for sample, urls in sorted(fastqs.items()):
        data = {"Script": "cfpipeline",
                "Output": f"s3://{bucket}/projects/{project}/{output}/{sample}",
                "Input": urls
                "Args": script_args + sample_args.get(sample, []) + [
                        "--sample", sample,
                        "--reference", f"s3://{bucket}/reference/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set_plus_hpv_panel.tar.gz",
                        "--panel", f"s3://{bucket}/panels/{panel}.tar.gz",
                        "--vep", f"s3://{bucket}/reference/vep_101_GRCh37_homo_sapiens_refseq.tar"]
                }
        
        message = json.dumps(data)
        print(message, "\n")
        n += 1
        if not dry_run:
            sqs.send_message(QueueUrl=queue_url, MessageBody=message)

    print("{} messages {}.".format(n, "processed" if dry_run else "queued"))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', "--bucket", help="AWS S3 bucket in which all files are located", required=True)
    parser.add_argument('-j', "--project", help="location of project files. s3://{bucket}/projects/{project}/", required=True)
    parser.add_argument('-p', "--panel", help="location of panel data. s3://{bucket}/panels/{panel}.tar.gz", required=True)
    parser.add_argument('-i', "--input", help="location of input fastqs. s3://{bucket}/projects/{project}/{input}", default=argparse.SUPPRESS)
    parser.add_argument('-o', "--output", help="location to write analysis output files. s3://{bucket}/projects/{project}/{output}", default=argparse.SUPPRESS)
    parser.add_argument('-s', "--samples", help="regular expression to narrow down input fastqs. single quote to prevent shell filename expansion", default=argparse.SUPPRESS)
    parser.add_argument('-m', "--manifest", help="path to tsv file containging additional arguments for individual files", default=argparse.SUPPRESS)
    parser.add_argument('-e', "--ignore-missing", help="proceed even if there are mismatches between the manifest and s3", action="store_const", const=True, default=argparse.SUPPRESS)
    parser.add_argument('-d', "--dry-run", action="store_const", const=True, default=argparse.SUPPRESS)
    args, script_args = parser.parse_known_args()
    enqueue(script_args=script_args, **vars(args))



if __name__ == "__main__":
    main()
