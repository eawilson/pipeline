import subprocess
import os
import sys
import pdb
import re
import gzip
import shutil
from collections import defaultdict, namedtuple

import covermi
import boto3



__all__ = ["run", "mount_basespace", "unmount_basespace", "mount_instance_storage", "list_basespace_fastqs", "ungzip_and_combine_illumina_fastqs", "load_panel_from_s3", "s3_put",
           "illumina_readgroup", "pipe"]



def ungzip_and_combine_illumina_fastqs(*filepaths, destination="", paired_end=True):
    """ Accepts list of fastqs to be ungzipped and combined. Fastqs will be merged if they only differ by lane number.
        Output will be written to current working directory.
    """
    fastqs = defaultdict(list)
    for filepath in filepaths:
        dirname, basename = os.path.split(filepath)
        parts = basename.split("_")
        if parts[-1].endswith(".gz"):
            parts[-1] = parts[-1][:-3]
        if not parts[-1] == "001.fastq":
            raise RuntimeError("Not a fastq {}.".format(filepath))
        if parts[-2] not in ("R1", "R2") or not parts[-3].startswith("L") or not len(parts[-3]) == 4:
            raise RuntimeError("Malformed fastq name {}.".format(filepath))
        parts[-3] = "L000"
        fastqs[os.path.join(destination, "_".join(parts))] += [filepath]

    if paired_end and len(fastqs) % 2:
        raise RuntimeError("Odd number of paired fastqs.")

    for dest, sources in fastqs.items():
        with open(dest, "wb") as f:
            for source in sorted(sources):
                if source.endswith(".gz"):
                    pipe(["gzip", "-dc", source], stdout=f)      
                else:
                    run(["cat", source], stdout=f)
    return sorted(fastqs.keys())



def illumina_readgroup(filepath):
    basename = os.path.basename(filepath)
    sample = "_".join(basename.split("_")[:-4]) # remove the _Sx_Lxxx_Rx_001.fastq from the name
    with open(filepath) as f:
        identifier = f.readline().split(":")
    flowcell = identifier[2]
    return "@RG\\tID:{}\\tSM:{}".format(flowcell, sample)



def run(args):
    return subprocess.run(args, stdout=subprocess.PIPE, stderr=sys.stderr, universal_newlines=True, check=True)



def pipe(args, stdout=None):
    return subprocess.run(args, stdout=stdout, stderr=sys.stderr, check=True)



def mount_basespace():
    basespace_dir = os.path.join(os.path.expanduser("~"), "basespace")
    if not os.path.exists(basespace_dir):
        os.mkdir(basespace_dir)
        
    completed = run(["mount"])
    for line in completed.stdout.split("\n"):
        if line.startswith("basemount on {}".format(basespace_dir)):
            print("Basespace already mounted.")
            return
        
    print("Mounting basespace.")
    run(["basemount", basespace_dir])



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
    


def s3_put(*filenames, prefix=""):
    s3 = boto3.client("s3")
    for filename in filenames:
        basename = os.path.basename(filename)
        print("Uploading {} to S3.".format(basename))
        s3.upload_file(filename, "omdc-data", "{}/{}".format(prefix, basename) if prefix else basename)



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
    s3 = boto3.client('s3')

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
    
    
    
    
    
    
    
    
    
    
    
    
