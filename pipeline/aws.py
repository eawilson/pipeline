import os
import sys
import pdb
import re
import csv
import io
import argparse
import pdb
from collections import defaultdict, namedtuple
import covermi


from boto3 import client
from botocore.credentials import InstanceMetadataFetcher

from .utils import run



__all__ = ["am_i_an_ec2_instance", 
           "s3_put", "s3_object_exists",
           "s3_get_tsv", 
           "s3_list_keys", 
           "s3_list_samples", 
           "s3_open", 
           "mount_instance_storage", 
           "load_panel_from_s3"]



def am_i_an_ec2_instance():
    return InstanceMetadataFetcher(timeout=1, num_attempts=1).retrieve_iam_role_credentials()



def s3_get(bucket, key, filename):
    s3 = client("s3")
    s3.download_file(bucket, key, filename)



def s3_put(bucket, filename, prefix=""):
    s3 = client("s3")
    basename = os.path.basename(filename)
    print("Uploading {} to S3.".format(basename))
    s3.upload_file(filename, bucket, "{}/{}".format(prefix, basename) if prefix else basename)



def s3_object_exists(bucket, prefix):
    s3 = client("s3")
    response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
    return response["KeyCount"]



def s3_list_keys(bucket, prefix, extension=""):
    """ Returns a dict of all objects in bucket that have the specified prefix and extension.
    """
    if extension:
        extension = ".{}".format(extension)
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
    for key in s3_list_keys(bucket, "projects/{}".format(project)):
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
                    
    if len(unformatted_block_devices) == 0:
        raise RuntimeError("No instance storage devices found.")
    elif len(unformatted_block_devices) > 1:
        raise RuntimeError("{} instance storage devices found.".format(len(unformatted_block_devices)))
    else:
        devname = unformatted_block_devices[0]
        run(["sudo", "mkfs", "-t", "ext4", devname])
        print("Mounting instance storage.")
        run(["sudo", "mount", devname, ephemoral_path])
        run(["sudo", "chmod", "go+rwx", ephemoral_path])
        return ephemoral_path



def load_panel_from_s3(panelname):
    s3 = client('s3')

    if not os.path.exists(panelname):
        print("Downloading {} from S3.".format(panelname))
        gziped_panel = "{}.tar.gz".format(panelname)
        s3.download_file("omdc-data", "panels/{}".format(gziped_panel), gziped_panel)
        print("Unpacking {}.".format(panelname))
        run(["tar", "xzf", gziped_panel])
        os.unlink(gziped_panel)
        
    panel = covermi.Panel(panelname)
    assembly = panel.properties.get("assembly", "GRCh37")
    transcript_source = panel.properties.get("transcript_source", "refseq")

    if not os.path.exists(assembly):
        print("Downloading {} from S3.".format(assembly))
        os.mkdir(assembly)
        os.chdir(assembly)
        objects = s3.list_objects(Bucket="omdc-data", Prefix="reference/{}/sequence".format(assembly[-2:])).get("Contents", [])
        if len(objects) != 1:
            raise RuntimeError("Unable to identify reference genome on S3.")
        s3.download_file("omdc-data", objects[0]["Key"], "genome.tar.gz")
        print("Unpacking genome.")
        run(["tar", "xzf", "genome.tar.gz"])
        os.unlink("genome.tar.gz")
        os.chdir("..")

    if not os.path.exists("vep"):
        print("Downloading {} from S3.".format(transcript_source))
        os.mkdir("vep")
        os.chdir("vep")
        objects = s3.list_objects(Bucket="omdc-data", Prefix="reference/{}/{}".format(assembly[-2:], transcript_source)).get("Contents", [])
        if len(objects) != 1:
            raise RuntimeError("Unable to identify vep data on S3.")
        s3.download_file("omdc-data", objects[0]["Key"], "vep.tar.gz")
        print("Unpacking vep data.")
        run(["tar", "xzf", "vep.tar.gz"])
        os.unlink("vep.tar.gz")
        os.chdir("..")
    
    fastas = [fn for fn in os.listdir(assembly) if fn.endswith(".fna")]
    if len(fastas) != 1:
        raise RuntimeError("Must be exactly one fasta in reference genome.")
    panel.properties["reference_fasta"] = os.path.join(assembly, fastas[0])
    return panel

    
