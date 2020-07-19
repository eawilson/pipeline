import os
import subprocess
import tarfile
import json
import boto3
import requests

from pipeline import mount_instance_storage


    
def terminate():
    response = requests.get("http://169.254.169.254/latest/meta-data/instance-id")
    instance_id = response.text
    ec2 = boto3.resource("ec2")
    ec2.instances.filter(InstanceIds=[instance_id]).terminate()



def parse_url(self, url):
    if url[:5].lower() != ("s3://"):
        raise RuntimeError(f"{url} is ot a valid s3 url.")
    parts = url.split("/")
    if len(parts) < 4:
        raise RuntimeError(f"{url} is ot a valid s3 url.")
    bucket = parts[2]
    key = "/".join(parts[3:])
    fn = parts[-1]
    return (bucket, key, fn)



class S3(object):
    def __init__(self):
        self.client = boto3.client("s3")
        self.downloaded = set()
        
    def download(self, url, path=""):
        bucket, key, fn = parse_url(url)
        if path:
            fn = f"{path}/{fn}"
        if fn not in self.downloaded:
            self.client.download_file(bucket, key, fn)
            self.downloaded.add(fn)
        return fn
    
    def upload(self, fn, url):
        bucket, key, fn = parse_url(f"{url}/{fn}"))
        self.client.upload_file(fn, bucket, key)
        
        

def main():
    if am_i_an_ec2_instance():
        os.cwd(mount_instance_storage())
    
    s3 = S3()
    sqs = boto3.client("sqs")
    queue_url = sqs.get_queue_url(QueueName="samples.fifo")["QueueUrl"]
    while True:
        if any(os.path.isfile(fn) for fn if os.listdir()):
            raise RuntimeError("Working directory is not empty.")
        
        response = sqs.receive_message(QueueUrl=queue_url,
                                       VisibilityTimeout=60*30,
                                       WaitTimeSeconds=20}
        
        try:
            message = response["Messages"][0]
        except KeyError: # No messages to process
            break
        
        body = json.loads(message["Body"])
        
        script = body["Script"]
        if script[:5].lower() == ("s3://"):
            os.mkdirs("scripts", exist_ok=True)
            script = os.path.join(".", s3.download(script, path="scripts"))
        
        body["Args"] = [s3.download(url) for url in body.get("Args", ())]
        
        for paramter, value in body.get("Kwargs", {}).items():
            if value[:5].lower() == ("s3://"):
                fn = s3.download(value)
                if fn.endswith(".tar.gz"):
                    with tarfile.open(fn, "r:gz") as archive:
                        members = archive.getnames()
                    top_level = [member for member in members if "/" not in member]
                    if len(top_level) != 1:
                        raise RuntimeError(f"{url} is not a valid tar.gz directory.")
                    fn = top_level[0]
                body["Kwargs"][parameter] = fn
        
        command_line = script + body.get("Args", []) + list(chain(*body.get("Kwargs", {}).items()))
        command_line = " ".join([(f"'{token}'" if " " in token else token) for token in command_line])
        with open("{}.log.txt".format(body["Kwargs"]["--sample"]), "wb") as log:
            completed_process = subprocess.run(command_line, shell=True, stderr=log)
        
        for fn in os.listdir():
            if os.path.isfile(fn):
                if fn not in body["Args"]:
                    s3.upload(fn, body["Output"])
                os.unlink(fn)
        
        response = sqs.delete_message(QueueUrl=queue_url,
                                      ReceiptHandle=message["ReceiptHandle"]



if __name__ == "__main__":
    main()



