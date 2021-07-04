from collections import defaultdict
import json
import sys
import pdb
import shlex

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
    
    dry_run =  pop(args, "--dry-run", nargs=1)
    input_extension = pop(args, "-x", "--input-extension") or "fastq.gz"
    output_extension = pop(args, "-X", "--output-extension") or "bam"
    single_sample =  pop(args, "--single-sample", nargs=1)
    input_url = pop(args, "-i", "--input", required=True)
    if not input_url.endswith("/"):
        input_url = f"{input_url}/"
    output_url = pop(args, "-o", "--output", required=True)
    if not output_url.endswith("/"):
        output_url = f"{output_url}/"
    
    if command.startswith("cfpipeline"):
        if "--panel" not in args and "-p" not in args:
            sys.exit(f"sqs_enqueue: --panel is a required argument")
        panel = pop(args, "--panel", "-p")
        if not panel.lower().startswith("s3://"):
            panel = f"s3://omdc-data/panels/{panel}"
        if not panel.endswith(".tar.gz"):
            panel = f"{panel}.tar.gz"
        args.extend(["--panel", panel])
        
        if "--reference" in args or "-r" in args:
            reference = pop(args, "--reference", "-r")
            if not reference.endswith(".tar.gz"):
                reference = f"{reference}.tar.gz"
        else:
            reference = "s3://omdc-data/reference/GRCh37_EBV_HPV_bwa_mem.tar.gz"
        args.extend(["--reference", reference])
        
        if "--vep" not in args and "-v" not in args:
            args.extend(["--vep", "s3://omdc-data/reference/vep_104_GRCh37_homo_sapiens_refseq.tar"])
    
    
    sqs = boto3.client("sqs")
    try:
        queue_url = sqs.get_queue_url(QueueName=queue_name)["QueueUrl"]
    except sqs.exceptions.QueueDoesNotExist:
        sys.exit(f'sqs_enqueue: "{queue_name}" is not a valid sqs queue')
    
    bams = set()
    bucket, stem = parse_url(output_url)
    for key in s3_list(bucket, stem, extension=f".{output_extension}"):
        identifier = "/".join(key[len(stem):].split("/")[-3:-1])
        bams.add(identifier)
    
    fastqs = defaultdict(list)
    bucket, stem = parse_url(input_url)
    for key in s3_list(bucket, stem, extension=f".{input_extension}"):
        identifier = "/".join(key[len(stem):].split("/")[-3:-1])
        if identifier not in bams:
            fastqs[identifier] += [f"s3://{bucket}/{key}"]
    
    
    n = 0
    for identifier, samples in sorted(fastqs.items()):
        cmd = [command] + sorted(samples) + args + ["--name", identifier.split("/")[-1], "--output", f"{output_url}{identifier}"]        
        print(" ".join(shlex.quote(token) for token in cmd), file=sys.stderr)
        message = json.dumps(cmd)
        n += 1
        if not dry_run:
            sqs.send_message(QueueUrl=queue_url, MessageBody=message)
        if single_sample:
            break

    print("sqs_enqueue: {} messages {}.".format(n, "processed" if dry_run else "queued"), file=sys.stderr)



if __name__ == "__main__":
    main()
