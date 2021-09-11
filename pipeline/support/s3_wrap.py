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



def download_and_unpack(client, path, destination="."):
    if path[:5].lower() == ("s3://"):
        fn = path.split("/")[-1]
        if fn.endswith(".tar.gz"):
            fn = fn[:-7]
            cmd = f"aws s3 cp '{path}' - | tar xzf -"
        elif fn.endswith(".tgz"):
            fn = fn[:-4]
            cmd = f"aws s3 cp '{path}' - | tar xzf -"
        elif fn.endswith(".tar"):
            fn = fn[:-4]
            cmd = f"aws s3 cp '{path}' - | tar xf -"
        #elif fn.endswith(".gz"):
            #fn = fn[:-3]
            #cmd = f"aws s3 cp '{path}' - | gzip -dc >'{fn}'"
        else:
            cmd = f"aws s3 cp '{path}' '{fn}'"
        
        path = fn if destination == "." else os.path.join(destination, fn)
        if not os.path.exists(path):
            print(f"Downloading {fn}", file=sys.stderr)
            subprocess.run(cmd, shell=True, cwd=destination)
        
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
    args = sys.argv[1:]
    if len(args) == 0:
        sys.exit("s3_wrap: No command line arguments")
    command = args[0]
    
    no_log = "--no-log" in args
    if no_log:
        args.remove("--no-log")
    
    profile = ["profile"]
    if "--no-profile" in args:
        args.remove("--no-profile")
        profile = []
    
    
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
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
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
            args[i] = download_and_unpack(s3, arg)
    
    if no_log:
        subprocess.run(profile + args)
    else:
        with open(os.path.join(output_dir, f"{name}.{command}.log.txt"), "wb") as log:
            log.write(" ".join(shlex.quote(arg) for arg in args).encode())
            log.write("\n".encode())
            log.flush()
            cp = subprocess.run(profile + args, stderr=subprocess.STDOUT, stdout=log)
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



