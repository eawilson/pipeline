import os
import sys
import pdb
import re
import csv
import io
import argparse
import pdb
from collections import defaultdict, namedtuple
import requests


from boto3 import client
from botocore.credentials import InstanceMetadataFetcher

from .utils import run



__all__ = ["am_i_an_ec2_instance",
           "s3_get",
           "s3_put",
           "s3_exists",
           "s3_list", 
           "s3_list_samples", 
           "s3_open", 
           "mount_instance_storage"]



def aws_region():
    BASE_URL = "http://169.254.169.254/latest"
    metadata = {}
    session = requests.Session()
    try:
        response = session.put(f"{BASE_URL}/api/token",
                            headers={"X-aws-ec2-metadata-token-ttl-seconds": "60"},
                            timeout=3.0)
    except requests.exceptions.ConnectionError:
        return
    if response.status_code != 200:
        return
    session.headers.update({"X-aws-ec2-metadata-token": response.text})
    
    response = session.get(f"{BASE_URL}/meta-data/placement/region",
                           timeout=2.0)
    if response.status_code == 200:
        return response.text



def spot_interuption():
    BASE_URL = "http://169.254.169.254/latest"
    metadata = {}
    session = requests.Session()
    try:
        response = session.put(f"{BASE_URL}/api/token",
                            headers={"X-aws-ec2-metadata-token-ttl-seconds": "60"},
                            timeout=3.0)
    except requests.exceptions.ConnectionError:
        return {}
    
    if response.status_code == 200:
        session.headers.update({"X-aws-ec2-metadata-token": response.text})        
        response = session.get(f"{BASE_URL}/meta-data/spot/instance-action",
                            timeout=2.0)
        if response.status_code == 200:
            return response.json()
    return {}



def boto3_client(service):
    region_name = aws_region()
    if region_name:
        return client(service, region_name=region_name)
        

def am_i_an_ec2_instance():
    return InstanceMetadataFetcher(timeout=1, num_attempts=1).retrieve_iam_role_credentials()



def s3_get(bucket, key, filename="", prefix=""):
    if not filename:
        filename = key.split("/")[-1]
    if prefix:
        prefix = prefix.rstrip("/")
        filename = f"{prefix}/{filename}"
    s3 = client("s3")
    s3.download_file(bucket, key, filename)



def s3_put(bucket, filename, prefix=""):
    s3 = client("s3")
    basename = os.path.basename(filename)
    print("Uploading {} to S3.".format(basename))
    s3.upload_file(filename, bucket, "{}/{}".format(prefix.rstrip("/"), basename) if prefix else basename)



def s3_exists(bucket, prefix):
    s3 = client("s3")
    response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
    return response["KeyCount"]



def s3_list(bucket, prefix, extension=""):
    """ Returns a dict of all objects in bucket that have the specified prefix and extension.
    """
    s3 = client("s3")
    response = {}
    kwargs = {}
    keys = {}
    while response.get("IsTruncated", True):
        response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix, **kwargs)
        for content in response.get("Contents", ()):
            if content["Key"].endswith(extension):
                keys[content["Key"]] = content
        kwargs = {"ContinuationToken": response.get("NextContinuationToken", None)}
    return keys



def s3_list_samples(bucket, project):
    samples = set()
    for key in s3_list(bucket, "projects/{}".format(project)):
        split_key = key.split("/")
        if len(split_key) > 3:
            samples.add(split_key[2])
    return samples



class s3_open(object):
    def __init__(self, bucket, key, mode="rt"):
        if mode not in ("rt", "rb", "wt", "wb"):
            raise ValueError("Invalid mode {}".format(repr(mode)))
        self.s3 = client("s3")
        self.bucket = bucket
        self.key = key
        self.mode = mode
        self.f_bytes = io.BytesIO()
        if mode.startswith("r"):
            self.s3.download_fileobj(bucket, key, self.f_bytes)
            self.f_bytes.seek(0)
        self.f = io.TextIOWrapper(self.f_bytes) if mode.endswith("t") else self.f_bytes
    
    def close(self):
        if self.mode.startswith("w"):
            self.f_bytes.seek(0)
            self.s3.upload_fileobj(self.f_bytes, self.bucket, self.key)
        self.f.close()
            
    def __enter__(self):
        return self.f
    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.close()



def mount_instance_storage():
    if not am_i_an_ec2_instance():
        sys.exit("Not an EC2 instance")
    
    ephemoral_path = os.path.join(os.path.expanduser("~"), "ephemoral")
    if not os.path.exists(ephemoral_path):
        os.mkdir(ephemoral_path)
    
    # Dont try and mount again if already mounted
    completed = run(["mount"])
    for line in completed.stdout.split("\n"):
        if " on {} type ".format(ephemoral_path) in line:
            print("Instance storage already mounted.")
            return ephemoral_path

    completed = run(["lsblk"])
    devices = completed.stdout.strip("\n").split("\n")[1:]
    unformatted_block_devices = []
    for device_pair in zip(devices, devices[1:] + [""]):
        if not any(device.startswith("└─") or device.startswith("├─") for device in device_pair):
            name, majmin, rm, size, ro, devtype, mountpoint = (device_pair[0].split()  + [""])[:7]
            devname = "/dev/{}".format(name)
            if devtype == "disk" and mountpoint == "":
                completed = run(["sudo", "file", "-s", devname])
                if completed.stdout.strip("\n") == "{}: data".format(devname):
                    unformatted_block_devices += [devname]
                    
    num_devices = len(unformatted_block_devices)
    if num_devices == 0:
        raise RuntimeError("No instance storage devices found.")

    print("Mounting instance storage.")
    if num_devices > 1:
        device = "/dev/md0"
        run(["sudo", "mdadm", "--create", device, "--level=0", f"--raid-devices={num_devices}"] + unformatted_block_devices)
        
        # Alternative approach using lvm to create a jbod volume
        #for device in unformatted_block_devices:
            #run(["sudo", "pvcreate", device])
        #run(["sudo", "vgcreate", "ephemoral_volume_group"] + unformatted_block_devices)
        #device = "ephemoral_logical_volume"
        #run(["sudo", "lvcreate", "-n", device, "-l", "100%FREE", "ephemoral_volume_group"])
    else:
        device = unformatted_block_devices[0]
        
    run(["sudo", "mkfs", "-t", "ext4", device])
    run(["sudo", "mount", device, ephemoral_path])
    run(["sudo", "chmod", "go+rwx", ephemoral_path])
    return ephemoral_path



def s3_resources(*s3_urls):
    for url in s3_urls:
        if url[:5].upper() == "S3://":
            split_url = url.split("/")
            if len(split_url) < 4:
                raise ValueError(f"Invalid S3 url {url}")
            bucket = split_url[2]
            key = "/".join(split_url[3:])
            filename = split_url[-1]
            
        else:
            raise ValueError(f"Invalid S3 url {url}")
        
        s3_get(bucket, key, filename)
        run(["tar", "-xzf", filename])
        run(["rm", filename])
