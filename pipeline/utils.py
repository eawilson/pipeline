import subprocess
import os
import pdb
import re
import gzip
import shutil
from collections import defaultdict, namedtuple

import covermi
import boto3



__all__ = ["run", "mount_basespace", "mount_instance_storage", "unmount_basespace", "list_basespace_fastqs", "copy_fastqs", "reference_genome", "ungzip_and_combine_illumina_fastqs"]



def ungzip_and_combine_illumina_fastqs(*filepaths, destination=""):
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

    for dest, sources in fastqs.items():
        with open(dest, "wb") as f:
            for source in sources:
                if source.endswith(".gz"):
                    run(["gzip", "-dc", source], stdout=f, universal_newlines=False)      
                else:
                    run(["cat", source], stdout=f, universal_newlines=False)
    return fastqs.keys()



def basespace_path():
    return os.path.join(os.path.expanduser("~"), "basespace")



def run(*args, **kwargs):
    return subprocess.run(*args, **{"stdout": subprocess.PIPE, "check": True, "universal_newlines": True, **kwargs})



def mount_basespace():
    basespace_dir = basespace_path()
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
    basespace_dir = basespace_path()
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
            
            name, majmin, rm, size, ro, devtype, mountpoint = (device_pair[0].split()  + [""])[:6]
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
    


def s3_put(name, folder=""):
    s3 = boto3.client("s3")
    s3.upload_file()
    for bucket in s3.buckets.all():
        print(bucket.name)



def list_basespace_fastqs(project="", sample=""):
    project_regex = re.compile(project)
    sample_regex = re.compile(sample)

    matches = defaultdict(list)
    projects_dir = os.path.join(basespace_path(), "Projects")
    projects = [project for project in os.listdir(projects_dir) if not project.startswith(".") and project_regex.search(project)]    
    
    for project in projects:
        print("Indexing basespace project {}.".format(project))
        samples_dir = os.path.join(projects_dir, project, "Samples")
        samples = [sample for sample in os.listdir(samples_dir) if not sample.startswith(".") and sample_regex.search(sample)]
        
        for sample in samples:
            files_dir = os.path.join(samples_dir, sample, "Files")
            matches[sample] += [os.path.join(files_dir, fn) for fn in os.listdir(files_dir) if not fn.startswith(".") and (fn.endswith(".fastq") or fn.endswith(".fastq.gz"))]
            
    return matches
    


def reference_genome(build):
    reference_dir = os.path.join(ngsdata_path(), "bwa_reference", build)
    fastas = [name for name in os.listdir(reference_dir) if name.endswith(".fna")]
    if len(fastas) != 1:
        raise RuntimeError("More than one fasta in {}.".format(reference_dir))
    return os.path.join(reference_dir, fastas[0])



def load_panel_from_s3(panelname):
    if not os.path.exists("reference"):
        os.mkdir("reference")
    
    s3 = boto3.client('s3')

    if not os.path.exists(panelname):
        print("Downloading {}".formar(panel))
        gziped_panel = "{}.tar.gz".format(panelname)
        s3.download_file("omdc-data", "panels/{}".format(gziped_panel), gziped_panel)
        run(["tar", "xzf", gziped_panel])
        os.unlink(gziped_panel)
        
    panel = covermi.Panel(panelname)
    assembly = panel.properties.get("assembly", "GRCh37")
    assembly = assembly[-2:]
    transcript_source = panel.properties.get("transcript_source", "refseq")

    if not os.path.exists("genome"):
        objects = s3.list_objects(Bucket="omdc-data", Prefix="reference/{}/sequence".format(assembly)).get("Contents", [])
        if len(objects) != 1:
            raise RuntimeError("Not able to identify reference genome on S3.")
        s3.download_file("omdc-data", objects[0]["Key"], "genome.tar.gz")
        run(["tar", "xzf", "genome.tar.gz"])
        os.unlink("genome.tar.gz")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    