import os, pdb, sys
from pipeline import s3_object_exists
import boto3
import csv


projects = ["19-CAPPSeq-5-cf", "19CAPPSEQ007_DLBCL", "19CAPPSEQ008_DLBCL", "19CAPPSEQ009_DLBCL", "ACCEPT-Cf", "ACCEPT1", "Untitled from 190321_NB552023_0004_AH55YVBGXB"]
basespace_path = "/home/ed/basespace/Projects/{}/Samples"


#s3 = boto3.client("s3")
#objects = [content["Key"] for content in s3.list_objects_v2(Bucket="omdc-data", Prefix="projects/accept")["Contents"]]
#objects = [obj for obj in objects if obj.endswith(".annotation.tsv]

#for obj in objects:
    #split_obj = obj.split("/")
    #if len(split_obj) < 4:
        #continue
        
    #folder = split_obj[-2]
    #sample = split_obj[-1]
    #if folder not in sample:
        #split_sample = sample.split(".")
        #split_sample[0] = folder
        #split_obj[-1] = ".".join(split_sample)
        #new_obj = "/".join(split_obj)
        #print(obj, new_obj)
        #s3.copy({'Bucket': 'omdc-data', 'Key': obj}, 'omdc-data', new_obj)
        #s3.delete_object(Bucket='omdc-data', Key=obj)













#s3
s3 = boto3.client("s3")
objects = [content["Key"] for content in s3.list_objects_v2(Bucket="omdc-data", Prefix="projects/accept")["Contents"]]
s3_samples = set([obj.split("/")[2].split(".")[0] for obj in objects if len(obj.split("/")) == 4])


# basespace
bs_samples = {}
for project in projects:
    project_path = basespace_path.format(project)
    for fn in os.listdir(project_path):
        if not fn.startswith(".") and not fn.startswith("B") and not fn.endswith("-c"):
            fn = fn.split(".")[0]
            if fn ==     "10040005-H21161-c-1":
                sample = "10400005-H21161-c-1"
            elif fn ==   "10010028-H32577-c-0":
                sample = "10010027-H32577-c-0"
            elif fn ==   "10010028-H35990-c-1":
                sample = "10010027-H35990-c-1"
            elif fn ==   "10010028-H539-c-2":
                sample = "10010027-H539-c-2"
            else:
                sample = fn
            bs_samples[sample] = [project, fn]
                        

            
            
                
            #files += [[patient_id, timepoint, project, fn]]



with open("accept.samples.tsv", "wt") as f_out:
    writer = csv.writer(f_out, delimiter="\t")
    writer.writerow(["Patient", "Timepoint", "Sample", "Run", "Basespace_Sample"])
    for sample in s3_samples:
        patient_id, sample_id, cee, timepoint = sample.split(".")[0].split("-")
        if cee == "g":
            timepoint = "germline"
        elif patient_id.startswith("S"):
            timepoint = "normal"
        writer.writerow([patient_id, timepoint, sample] + bs_samples.get(sample, []))



            
            
