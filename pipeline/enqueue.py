from collections import defaultdict
import json
import argparse
import re

from pipeline.aws import s3_list
import boto3



def enqueue(bucket, project, panel, input="samples", output="analyses", samples=".*", min_family_size=None, umi=None, cnv=None, sizes=None, run=(), dry_run=False):
    #project = "accept"
    #panel = "Accept"
    #output = "output3"
    
    sample_regex = re.compile(samples)
    
    kwargs = {}
    if min_family_size is not None:
        kwargs["--min-family-size"] = min_family_size
    if umi is not None:
        kwargs["--umi"] = umi
    if cnv is not None:
        kwargs["--cnv"] = cnv
    if sizes is not None:
        kwargs["--sizes"] = sizes
    
    sqs = boto3.client("sqs")
    queue_url = sqs.get_queue_url(QueueName="samples")["QueueUrl"]

    complete = set()
    for key in s3_list(bucket, f"projects/{project}/{output}", extension=".bam"):
        sample = key.split("/")[-1].split("_")[0]
        complete.add(sample)

    fastqs = defaultdict(list)
    for key in s3_list(bucket, f"projects/{project}/{input}", extension=".fastq.gz"):
        sample = key.split("/")[-1].split("_")[0]
        if sample not in complete and sample_regex.fullmatch(sample):
            fastqs[sample] += [f"s3://{bucket}/{key}"]
    
    n = 0
    for sample, urls in sorted(fastqs.items()):
        data = {"Script": "cfpipeline",
                "Output": f"s3://{bucket}/projects/{project}/{output}/{sample}",
                "Args": urls,
                "Kwargs": {"--sample": sample,
                           "--reference": f"s3://{bucket}/reference/GCA_000001405.14_GRCh37.p13_no_alt_analysis_set_plus_hpv_panel.tar.gz",
                           "--panel": f"s3://{bucket}/panels/{panel}.tar.gz",
                           "--vep": f"s3://{bucket}/reference/vep_101_GRCh37_homo_sapiens_refseq.tar",
                           **kwargs,
                           }
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
    parser.add_argument('-m', "--min-family-size", default=argparse.SUPPRESS)
    parser.add_argument('-u', "--umi", default=argparse.SUPPRESS)
    parser.add_argument('-c', "--cnv", default=argparse.SUPPRESS)
    parser.add_argument('-z', "--sizes", default=argparse.SUPPRESS)
    parser.add_argument('-d', "--dry-run", action="store_const", const=True, default=argparse.SUPPRESS)
    args = parser.parse_args()
    enqueue(**vars(args))



if __name__ == "__main__":
    main()
