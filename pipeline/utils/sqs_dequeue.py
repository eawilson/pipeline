#!/usr/bin/env python3

import os
import subprocess
import json
import shlex
import sys

import boto3

from pipeline import mount_instance_storage, am_i_an_ec2_instance



def main():
    """
    """
    print("Starting sqs_dequeue...", file=sys.stderr)
    
    args = sys.argv[1:]
    if len(args) == 0:
        sys.exit("sqs_dequeue: No queue specified")
    sqs = boto3.client("sqs")
    queue_url = sqs.get_queue_url(QueueName=args[0])["QueueUrl"]
    
    if am_i_an_ec2_instance():
        os.chdir(mount_instance_storage())
    
    while True:
        if any(os.path.isfile(fn) for fn in os.listdir()):
            sys.exit("sqs_dequeue: Working directory not empty")
        
        response = sqs.receive_message(QueueUrl=queue_url, WaitTimeSeconds=20)
        if "Messages" not in response:
            break
        
        message = response["Messages"][0]
        command = json.loads(message)
        if not isinstance(command, list):
            print('sqs_dequeue: Malformed message "{command}"', file=sys.stderr)
            continue
        
        print(" ".join(shlex.quote(token) for token in command), file=sys.stderr)
        subprocess.run(["s3_wrap"] + command)
        
        for fn in os.listdir():
            if os.path.isfile(fn):
                 os.unlink(fn)
        
        try:
            response = sqs.delete_message(QueueUrl=queue_url, ReceiptHandle=message["ReceiptHandle"])
        except sqs.exceptions.QueueDoesNotExist:
            # Happens if the queue is purged while a job is in progress.
            pass
    
    print("Complete.", file=sys.stderr)



if __name__ == "__main__":
    main()



