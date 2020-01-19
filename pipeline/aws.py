import subprocess
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



__all__ = ["s3_put", "s3_object_exists", "s3_get_tsv", "s3_list_keys", "s3_list_samples", "s3_open", "run", \
            "mount_basespace", "unmount_basespace", "mount_instance_storage", "list_basespace_fastqs", "ungzip_and_combine_illumina_fastqs", \
            "load_panel_from_s3", "illumina_readgroup", "pipe", "command_line_arguments", "sample_name"]

BUCKET = "omdc-data"
ILLUMINA_FASTQ = re.compile(r"(.+)_S([0-9]{1,2})_L([0-9]{3})_R([12])_001\.fastq(\.gz)?$") # name, s_number, lane, read, gzip


def command_line_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_r1", help="Uncompressed fastq read 1.")
    parser.add_argument("fastq_r2", help="Uncompressed fastq read 2.")
    parser.add_argument("-b", "--bam", help="Matched normal bam to perform tumour/normal subtraction.")
    parser.add_argument("-p", "--panel", help="Directory containing panel data.")
    parser.add_argument("-o", "--output", help="Output directory.", dest="dest_dir", default=argparse.SUPPRESS)
    parser.add_argument("-g", "--genome", help="Directory containing reference genome.", dest="genome_dir", default=argparse.SUPPRESS)
    args = parser.parse_args()
    return vars(args)



def sample_name(*fastqs):
    matches = [ILLUMINA_FASTQ.match(fastq) for fastq in fastqs]
    if not None in matches:
        names = set(match.group(1) for match in matches)
        if len(names) == 1:
            return list(names)[0]
    raise RuntimeError("Inconsistent fastq names.")
    
    

def s3_put(*filenames, prefix=""):
    s3 = client("s3")
    for filename in filenames:
        basename = os.path.basename(filename)
        print("Uploading {} to S3.".format(basename))
        s3.upload_file(filename, BUCKET, "{}/{}".format(prefix, basename) if prefix else basename)



def s3_object_exists(prefix):
    s3 = client("s3")
    response = s3.list_objects_v2(Bucket=BUCKET, Prefix=prefix)
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



#def ungzip_and_combine_illumina_fastqs(*filepaths, destination="", paired_end=True):
    #""" Accepts list of fastqs to be ungzipped and combined. Fastqs will be merged if they only differ by lane number.
        #Output will be written to current working directory.
    #"""
    #fastqs = defaultdict(list)
    #for filepath in filepaths:
        #dirname, basename = os.path.split(filepath)
        #parts = basename.split("_")
        #if parts[-1].endswith(".gz"):
            #parts[-1] = parts[-1][:-3]
        #if not parts[-1] == "001.fastq":
            #raise RuntimeError("Not a fastq {}.".format(filepath))
        #if parts[-2] not in ("R1", "R2") or not parts[-3].startswith("L") or not len(parts[-3]) == 4:
            #raise RuntimeError("Malformed fastq name {}.".format(filepath))
        #parts[-3] = "L000"
        #fastqs[os.path.join(destination, "_".join(parts))] += [filepath]

    #if paired_end and len(fastqs) % 2:
        #raise RuntimeError("Odd number of paired fastqs.")

    #for dest, sources in fastqs.items():
        #with open(dest, "wb") as f:
            #for source in sorted(sources):
                #if source.endswith(".gz"):
                    #pipe(["gzip", "-dc", source], stdout=f)      
                #else:
                    #run(["cat", source], stdout=f)
    #return sorted(fastqs.keys())



def ungzip_and_combine_illumina_fastqs(name=None, source_dir=".", dest_dir="."):
    """ Accepts list of fastqs to be ungzipped and combined. Fastqs will be merged if they only differ by lane number.
        Output will be written to current working directory.
    """
    if name is None:
        name = ".+"
    else:
        name = re.escape(name)
    regex = re.compile(ILLUMINA_FASTQ)
    previous_samples = None
    previous_s_numbers = None
    
    fastqs = defaultdict(list)
    for fn in os.listdir(source_dir):
        match = regex.match(fn)
        if match:
            sample, s_number, lane, read = match.group(1, 2, 3, 4)
            if (previous_samples is not None and sample != previous_samples) or \
                (previous_s_numbers is not None and s_number != previous_s_numbers):
                raise RuntimeError("Mismatched fastqs.")
            previous_samples = sample
            fastqs[read] += [os.path.join(source_dir, fn)]

    if not fastqs:
        raise RuntimeError("No matching fastqs found.")
    if len(fastqs["1"]) != len(fastqs["2"]):
        raise RuntimeError("Unmatched fastqs.")

    ret = []
    for read, fns in sorted(fastqs.items()):
        #if not (len(fns) == 1 and fns[0].endswith(".fastq") and source_dir == dest_dir):
        ret += [os.path.join(dest_dir, "{}_S{}_L000_R{}_001.fastq".format(sample, s_number, read))]
        with open(ret[-1], "wb") as f:
            for fn in sorted(fns):
                if fn.endswith(".gz"):
                    subprocess.run(["gzip", "-dc", fn], stdout=f)      
                else:
                    subprocess.run(["cat", fn], stdout=f)
    return ret



def illumina_readgroup(filepath):
    basename = os.path.basename(filepath)
    sample = "_".join(basename.split("_")[:-4]) # remove the _Sx_Lxxx_Rx_001.fastq from the name
    with open(filepath) as f:
        identifier = f.readline().split(":")
    flowcell = identifier[2]
    return "@RG\\tID:{}\\tSM:{}".format(flowcell, sample)



def run(args):
    """ Run a unix command as a subprocess. Stdout and stderr are captured as a string for review if needed.
        Not to be used for main pipeline comands which should be called with pipe instead.  
    """
    return subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)



def pipe(args, **kwargs):
    """ Runs a main pipeline command. Output is bytes rather than string and is expected to be captured via
        stdout redirection or ignored if not needed. The command is echoed as a bytestring to stderr (if supplied)
        and all stderr diagnostic output from the subprocess is captured via redirection.
    """
    print(" ".join(str(arg) for arg in args), file=sys.stderr)
    return subprocess.run(args, check=True, **kwargs)



def mount_basespace():
    basespace_dir = os.path.join(os.path.expanduser("~"), "basespace")
    if not os.path.exists(basespace_dir):
        os.mkdir(basespace_dir)
        
    completed = run(["mount"])
    for line in completed.stdout.split("\n"):
        if line.startswith("basemount on {}".format(basespace_dir)):
            print("Basespace already mounted.")
            return basespace_dir
        
    print("Mounting basespace.")
    run(["basemount", basespace_dir])
    return basespace_dir



def unmount_basespace():
    basespace_dir = os.path.join(os.path.expanduser("~"), "basespace")
    completed = run(["mount"])
    for line in completed.stdout.split("\n"):
        if line.startswith("basemount on {}".format(basespace_dir)):
            print("Unmounting basespace.")
            run(["basemount", "--unmount", basespace_dir])
            return
        
    print("Basespace not mounted.")



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
    
    
    
    
def list_basespace_fastqs(project="", sample=""):
    basespace_path = os.path.join(os.path.expanduser("~"), "basespace")
    project_regex = re.compile(project)
    sample_regex = re.compile(sample)

    matches = []
    projects_dir = os.path.join(basespace_path, "Projects")
    projects = [project for project in os.listdir(projects_dir) if not project.startswith(".") and project_regex.search(project)]    
    
    for project in projects:
        print("Indexing basespace project {}.".format(project))
        samples_dir = os.path.join(projects_dir, project, "Samples")
        samples = [sample for sample in os.listdir(samples_dir) if not sample.startswith(".") and sample_regex.search(sample)]
        
        for sample in samples:
            files_dir = os.path.join(samples_dir, sample, "Files")
            matches += [os.path.join(files_dir, fn) for fn in os.listdir(files_dir) if not fn.startswith(".") and (fn.endswith(".fastq") or fn.endswith(".fastq.gz"))]
            
    return matches



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

    
