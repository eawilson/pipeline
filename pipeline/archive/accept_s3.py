import os, sys, pdb
from collections import Counter, defaultdict
import csv
import io
from statistics import median 

from pipeline import s3_get_tsv, s3_list_samples, s3_UploadFileObj
import boto3

BUCKET = "omdc-data"



def clean_annotations(project):
    artifacts = defaultdict(set)
    for row in s3_get_tsv(BUCKET, "projects/{0}/{0}.artifacts.tsv".format(project)):
        #if int(row["Count"]) >= 1:
        artifacts[row["Run"]].add(row["Variant"])

    print("Cleaning...")
    for sample in s3_get_tsv(BUCKET, "projects/{0}/{0}.samples.tsv".format(project)):
        if sample["Timepoint"].isnumeric():
            print(sample["Sample"])
            reader = s3_get_tsv(BUCKET, "projects/{0}/{1}/{1}.annotation.tsv".format(project, sample["Sample"]))
            with s3_UploadFileObj(BUCKET, "projects/{0}/{1}/{1}.annotation.cleaned.tsv".format(project, sample["Sample"])) as f:
                writer = csv.DictWriter(f, reader.fieldnames, delimiter="\t")
                writer.writeheader()
                for row in reader:
                    variant = "{}:{} {}".format(row["Chrom"], row["Pos"], row["Change"])
                    if variant not in artifacts[""] and variant not in artifacts.get(sample["Run"], ()) and row["Impact"] in ("MODERATE", "HIGH"):
                        writer.writerow(row)



def combine_annotations(project, time_points=None):
    count = Counter()
    annotations = []

    print("Combining...")
    for sample in s3_get_tsv(BUCKET, "projects/{0}/{0}.samples.tsv".format(project)):
        if sample["Timepoint"].isnumeric():
            print(sample["Sample"])
            count[sample["Patient"]] += 1
            reader = s3_get_tsv(BUCKET, "projects/{0}/{1}/{1}.annotation.cleaned.tsv".format(project, sample["Sample"]))
            for row in reader:
                variant = "{}:{} {}".format(row["Chrom"], row["Pos"], row["Change"])
                annotations += [{"Variant": variant, "Patient": sample["Patient"], "Timepoint": sample["Timepoint"], **row}]
    
    with s3_UploadFileObj(BUCKET, "projects/{0}/{0}.annotation.combined.tsv".format(project)) as f:
        writer = csv.DictWriter(f, ["Variant", "Patient", "Timepoint"] + reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for row in sorted(annotations, key=lambda a: (a["Variant"], a["Patient"], a["Timepoint"])):
            if time_points is None or count[row["Patient"]] == time_points:
                writer.writerow(row)



def download_annotations(project):
    print("Downloading...")
    for sample in s3_get_tsv(BUCKET, "projects/{0}/{0}.samples.tsv".format(project)):
        print(sample["Sample"])
        reader = s3_get_tsv(BUCKET, "projects/{0}/{1}/{1}.annotation.tsv".format(project, sample["Sample"]))
        with open(sample["Sample"]+".annotation.tsv", "wt") as f:
            writer = csv.DictWriter(f, reader.fieldnames, delimiter="\t")
            writer.writeheader()
            for row in reader:
                writer.writerow(row)



def generate_artifact_list(project):
    samples_key = "projects/{0}/{0}.samples.tsv".format(project)
    samples_tsv = s3_get_tsv(BUCKET, samples_key)
    samples_s3 = s3_list_samples(BUCKET, project)
    
    for sample in samples_s3 - set(row["Sample"] for row in samples_tsv):
        print("WARNING Sample {} is on S3 but is not listed in {}. Skipping.".format(sample, samples_key))
        
    artifacts = defaultdict(lambda:defaultdict(list))
    totals = Counter()

    print("Listing artifacts...")
    for sample in samples_tsv:
        print(sample["Sample"])
        if sample["Sample"] not in samples_s3:
            print("WARNING Sample {} is listed in {} but is not on S3. Skipping".format(sample, samples_key))
            continue
        
        if sample["Timepoint"] in ("germline", "normal"):
            group = ""
        elif sample["Timepoint"].isnumeric():
            group = sample["Run"]
        else:
            continue
        
        totals[group] += 1
        for row in s3_get_tsv(BUCKET, "projects/{0}/{1}/{1}.annotation.tsv".format(project, sample["Sample"])):
            variant = "{}:{} {}".format(row["Chrom"], row["Pos"], row["Change"])
            artifacts[group][variant] += [float("{:.3f}".format(float(row["VAF"])))]
            
    with s3_UploadFileObj(BUCKET, "projects/{0}/{0}.artifacts.tsv".format(project)) as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Variant", "Run", "Count", "Total", "Median_VAF"])
        for group, variants in sorted(artifacts.items()):
            total = totals[group]
            for variant, vafs in sorted(variants.items(), key=lambda i: -len(i[1])):
                # Write all normals/germlines variants but only write sample variants if they occur on every sample in the run and >- 5 samples on run
                if group == "" or (total >= 5 and len(vafs) == total): 
                    writer.writerow([variant, group, len(vafs), total, median(vafs)])



def zz():
    print("????...")
    pathogenic = {}
    comment = {}
    for row in s3_get_tsv(BUCKET, "projects/accept/accept.annotation.combined.reviewed.tsv".format(project)):
        pathogenic[row["Variant"]] = row["RELEVANT"]
        comment[row["Variant"]] = row["COMMENTS"]
    
    annotations = []
    reader = s3_get_tsv(BUCKET, "projects/accept/accept.annotation.combined.tsv".format(project))
    for row in reader:
        if pathogenic.get(row["Variant"], "Y").strip() == "Y":
            annotations += [{"Relevant": pathogenic.get(row["Variant"], "").strip(), "Comment": comment.get(row["Variant"], ""), **row}]
            
    with s3_UploadFileObj(BUCKET, "projects/{0}/{0}.annotation.combined2.tsv".format(project)) as f:
        writer = csv.DictWriter(f, ["Relevant", "Comment"] + reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for row in annotations:
            writer.writerow(row)






if __name__ == "__main__":
    project = "accept"
    zz()
    
    ##download_annotations(project)
    #sys.exit()
    
    #generate_artifact_list(project)
    #clean_annotations(project)
    #combine_annotations(project, time_points=4)
