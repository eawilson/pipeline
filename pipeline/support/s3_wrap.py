#!/usr/bin/env python3

import os
import subprocess
import sys
import tempfile
import shlex

import boto3



def parse_url(url):
    if url[:5].lower() != ("s3://"):
        sys.exit(f"{url} is not a valid s3 url.")
    parts = url.split("/")
    if len(parts) < 4:
        sys.exit(f"{url} is not a valid s3 url.")
    bucket = parts[2]
    key = "/".join(parts[3:])
    return (bucket, key)



def download_and_untar(client, path, destination=""):
    if path[:5].lower() == ("s3://"):
        bucket, key = parse_url(path)
        downloadname = key.split("/")[-1]
        if destination:
            os.makedirs(destination, exist_ok=True)
            downloadname = os.path.join(destination, downloadname)
        if downloadname.endswith(".tar.gz"):
            path = downloadname[:-7]
        elif downloadname.endswith(".tar"):
            path = downloadname[:-4]
        else:
            path = downloadname

        if not os.path.exists(path):
            print(f"Downloading {key}.", file=sys.stderr)
            client.download_file(bucket, key, downloadname)
            if downloadname.endswith(".tar.gz"):
                tar_args = "-xzf"
            elif downloadname.endswith(".tar"):
                tar_args = "-xf"
            else:
                tar_args = ""
            if tar_args:
                print(f"Unpacking {key}.", file=sys.stderr)
                fn = os.path.basename(downloadname)
                subprocess.run(["tar", tar_args, fn], cwd=destination)
                os.unlink(downloadname)
    return path



def upload(client, fn, url):
    bucket, key = parse_url(url)
    if not key.endswith("/"):
        key = f"{key}/"
    key = "{}{}".format(key, os.path.basename(fn))
    print(f"Uploading {key}.", file=sys.stderr)
    client.upload_file(fn, bucket, key)



def main():
    """ Wrapper around a pipeline (or any other script) to assist
        running on an aws ecs or ec2 instance. All command line args
        are parsed and if any contain a reference to an s3 file then
        that file is downloaded, untarred if needed and the command
        line argument replaced with the path to the the local copy.
        If an --output (-o) argument is present that references an s3
        location then it is replaced with a local temporary directory,
        the contents of which are uploaded to the s3 location after
        script completion. A log of all script output to stderr or
        stdout is also recorded and uploaded with the rest of the
        output files. This log file is named after the --name (-n) 
        argument if present or the first positional argument if not
        preceeded by any keyword arguments. The only cleanup that is
        performed after script completion is the removal of the temp
        output directory and all contained files if they have been
        uploaded to s3.
    """
    original_args = args = sys.argv[1:]
    if len(args) == 0:
        sys.exit("s3_wrap: No command line arguments")
    
    s3_destination = None
    try:
        i = args.index("-o")
    except ValueError:
        try:
            i = args.index("--output")
        except ValueError:
            i = len(args)
    if i + 1 < len(args):
        if args[i + 1].lower().startswith("s3://"):
            s3_destination = args[i + 1]
            args[i + 1] = tempfile.mkdtemp(dir=".")
        output_dir = args[i + 1]
    else:
        output_dir = "."
        
    try:
        i = args.index("-n")
    except ValueError:
        try:
            i = args.index("--name")
        except ValueError:
            i = len(args)
    if i + 1 < len(args):
        name = args[i + 1]
    elif len(args) > 1 and not args[1].startswith("-"):
        name = args[1].split("/")[-1].split(".")[0]
    else:
        name = "script"
    
    s3 = boto3.client("s3")
    for i, arg in enumerate(list(args)):
        if arg.lower().startswith("s3://"):
            args[i] = download_and_untar(s3, arg)
    
    with open(os.path.join(output_dir, "{name}.log.txt"), "wb") as log:
        print(" ".join(shlex.quote(arg) for arg in original_args), file=log, flush=True)
        cp = subprocess.run(["profile"] + args, stderr=subprocess.STDOUT, stdout=log)
        retcode = cp.returncode
        if retcode != 0:
            msg = f"PROCESS EXITED WITH RETURN CODE {retcode}\n"
            log.write(msg.encode())
    
    if s3_destination is not None:
        for fn in os.listdir(output_dir):
            fn = os.path.join(output_dir, fn)
            if os.path.isfile(fn):
                upload(s3, fn, s3_destination)
                os.unlink(fn)
        os.rmdir(output_dir)
    
    sys.exit(retcode)



if __name__ == "__main__":
    main()



