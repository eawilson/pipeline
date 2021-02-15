import os
import subprocess
import json
import pdb
import shlex
import sys

import boto3

from pipeline import mount_instance_storage, am_i_an_ec2_instance



def parse_url(url):
    if url[:5].lower() != ("s3://"):
        raise RuntimeError(f"{url} is not a valid s3 url.")
    parts = url.split("/")
    if len(parts) < 4:
        raise RuntimeError(f"{url} is not a valid s3 url.")
    bucket = parts[2]
    key = "/".join(parts[3:])
    return (bucket, key)



def download(client, path, destination=""):
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
            print(f"Downloading {key}.")
            client.download_file(bucket, key, downloadname)
            if downloadname.endswith(".tar.gz") or downloadname.endswith(".tar"):
                print(f"Unpacking {key}.")
                args = "-xzf" if downloadname.endswith(".gz") else "-xf"
                subprocess.run(["tar", args, os.path.basename(downloadname)], cwd=destination)
                os.unlink(downloadname)
    return path



def upload(client, fn, url):
    bucket, key = parse_url(f"{url}/{fn}")
    client.upload_file(fn, bucket, key)
        
        

def main():
    print("Starting runner...")
    if am_i_an_ec2_instance():
        os.chdir(mount_instance_storage())
    
    s3 = boto3.client("s3")
    sqs = boto3.client("sqs")
    queue_url = sqs.get_queue_url(QueueName="samples")["QueueUrl"]
    while True:
        if any(os.path.isfile(fn) for fn in os.listdir()):
            sys.exit("Working directory is not empty")
        
        response = sqs.receive_message(QueueUrl=queue_url, WaitTimeSeconds=20)
        
        try:
            message = response["Messages"][0]
        except KeyError: # No messages to process
            break
        body = json.loads(message["Body"])
        
        if body["Script"][:5].lower() == ("s3://"):    
            body["Script"] = os.path.join(".", download(s3, body["Script"], destination="downloads"))
        body["Input"] = [download(s3, url) for url in body.get("Input", ())]
        body["Args"] = [download(s3, arg, destination="downloads") for arg in body.get("Args", ())]
        
        try:
            sample = body["Args"][body["Args"].index("--sample")+1]
        except (ValueError, IndexError):
            sample = "pipeline"
        
        command_line = [body["Script"]] + body["Input"] + body["Args"]
        print(" ".join(shlex.quote(token) for token in command_line))
        with open(f"{sample}.log.txt"), "wb") as log:
            completed_process = subprocess.run(command_line, stderr=subprocess.STDOUT, stdout=log)
            if completed_process.returncode != 0:
                log.write(f"PIPELINE EXITED WITH RETURN CODE {completed_process.returncode}\n".encode())
        
        print("Cleaning up.")
        for fn in os.listdir():
            if os.path.isfile(fn):
                if fn not in body["Input"]:
                    upload(s3, fn, body["Output"])
                os.unlink(fn)
        
        try:
            response = sqs.delete_message(QueueUrl=queue_url, ReceiptHandle=message["ReceiptHandle"])
        # I believe this happens if the queue is purged while a job is in progress.
        except Exception: # Actually throws QueueDoesNotExist, but where to import it from?
            pass
    
    
    print("Complete.")



if __name__ == "__main__":
    main()



