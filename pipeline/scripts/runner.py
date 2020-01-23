import os
import sys
import time
import datetime
import argparse
import pdb
from collections import defaultdict

from pipeline import run, mount_instance_storage, load_panel_from_s3, \
    ungzip_and_combine_illumina_fastqs, s3_put, illumina_readgroup, \
    pipe, create_report, s3_object_exists, s3_open, s3_list_keys



def runner(project, panel, genome):

    #os.cwd(mount_instance_storage())
    #load_panel_from_s3(panel)
    
    samples = defaultdict(list)
    for key in s3_list_keys("omdc-data", "projects/head_and_neck/samples"):
        sample = key.split("/")[-1].split("_")[0]
        samples[sample] += [key]
    
    pdb.set_trace()
    
    if s3_object_exists("projects/{}/{}".format(s3_project, sample)):
        print("{}/{} already exists, skipping.".format(s3_project, sample))
        return
        
    print("Fetching fastqs...")
    fastqs = ungzip_and_combine_illumina_fastqs(*fastqs)

    panel = load_panel_from_s3(panelname)
    prop = panel.properties

    print("Uploading to s3.")
    for filename in os.listdir():
        if os.path.isfile(filename):
            s3_put(filename, prefix="projects/{}/{}".format(s3_project, sample))

    #except Exception:
        #print("ERROR, SKIPPING")

    #for filename in os.listdir():
        #if os.path.isfile(filename):
            #os.unlink(filename)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("project", help="S3 project name.")
    parser.add_argument("-p", "--panel", help="Panel name.")
    parser.add_argument("-g", "--genome", help="Reference genome.")
    args = parser.parse_args()
    runner(**vars(args))



