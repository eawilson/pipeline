from collections import defaultdict
import json
import argparse
import re
import csv
import sys
import pdb

from pipeline.aws import s3_list, s3_exists
import boto3



def enqueue(script_args, bucket, project, panel="", input="samples", output="analyses", samples=".*", manifest="", ignore_missing=False, dry_run=False):
    sample_regex = re.compile(samples)

    sqs = boto3.client("sqs")
    queue_url = sqs.get_queue_url(QueueName="samples")["QueueUrl"]


    complete = set()
    for key in s3_list(bucket, f"projects/{project}/{output}/", extension=".bam"):
        run_sample = "/".join(key.split("/")[3:5])
        complete.add(run_sample)
    
    
    fastqs = defaultdict(list)
    for key in s3_list(bucket, f"projects/{project}/{input}", extension=".fastq.gz"):
        run_sample = "/".join(key.split("/")[-3:-1])
        if run_sample not in complete and sample_regex.fullmatch(run_sample):
            fastqs[run_sample] += [f"s3://{bucket}/{key}"]
    
    
    sample_args = {}
    if manifest:
        with open(manifest, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                row = {k.strip(): v.strip() for k, v in row.items()}
                
                run_sample = "/".join([row.pop("run"), row.pop("sample")])
                if run_sample in sample_args:
                    sys.exit(f"Sample {run_sample} duplicated in manifest")
                args = []
                for key, val in row.items():
                    args.append(f"--{key}")
                    if val:
                        args.append(val)
                sample_args[run_sample] = args
        
        if not ignore_missing:
            for name in set(fastqs) - set(sample_args):
                print(f"Sample {name} is present in aws but not in the manifest", file=sys.stderr)
            for name in set(sample_args) - set(fastqs):
                print(f"Sample {name} is present in the manifest but not aws", file=sys.stderr)
            if set(sample_args) ^ set(fastqs):
                sys.exit("Missing samples. Aborting")
            for name in set(fastqs) - set(sample_args):
                fastqs.pop(name)
    
    reference_url = f"s3://{bucket}/reference/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set_plus_hpv_panel.tar.gz"
    vep_url = f"s3://{bucket}/reference/vep_101_GRCh37_homo_sapiens_refseq.tar"
    
    if not s3_exists(bucket, "/".join(reference_url.split("/")[3:])):
        sys.exit(f"{reference_url} does not exist")
    if not s3_exists(bucket, "/".join(vep_url.split("/")[3:])):
        sys.exit(f"{vep_url} does not exist")

    n = 0
    for run_sample, urls in sorted(fastqs.items()):
        try:
            i = sample_args.get(run_sample, []).index("--panel")
            panel = sample_args[run_sample][i+1]
            sample_args[run_sample] = sample_args[run_sample][:i] + sample_args[run_sample][i+2:]
        except ValueError:
            pass
        
        if not panel:
            sys.exit("No panel provided")
        panel_url = f"s3://{bucket}/panels/{panel}.tar.gz"
        if not s3_exists(bucket, "/".join(panel_url.split("/")[3:])):
            sys.exit(f"{panel_url} does not exist")
        
        data = {"Script": "cfpipeline",
                "Output": f"s3://{bucket}/projects/{project}/{output}/{run_sample}",
                "Input": urls,
                "Args": script_args + sample_args.get(run_sample, []) + [
                        "--sample", run_sample.split("/")[1],
                        "--reference", reference_url,
                        "--panel", panel_url,
                        "--vep", vep_url]
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
    parser.add_argument('-p', "--panel", help="location of panel data. s3://{bucket}/panels/{panel}.tar.gz", default=argparse.SUPPRESS)
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
