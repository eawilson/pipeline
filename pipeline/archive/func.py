import pdb
from collections import Counter, defaultdict
import csv
import os
from statistics import median 
from datetime import datetime, timezone

from boto3 import client

from .aws import s3_open, s3_list, s3_list_samples
            
BUCKET = "omdc-data"

__all__ = ["titv", "clean_annotations", "combine_annotations", "download_annotations", "generate_artifact_list"]




def titv(project):
    bases = {"A": "U", "G": "U", "C": "Y", "T": "Y"}
    print("TiTv...")
    with s3_open(BUCKET, "projects/{0}/{0}.titv.tsv".format(project), "wt") as f_out:
        writer = csv.DictWriter(f_out, ["Sample", "Run", "TiTv"], delimiter="\t")
        writer.writeheader()
        
        with s3_open(BUCKET, "projects/{0}/{0}.samples.tsv".format(project)) as f_samples:
            for sample in csv.DictReader(f_samples, delimiter="\t"):
                print(sample["Sample"])
                counts = Counter()
                with s3_open(BUCKET, "projects/{0}/{1}/{1}.annotation.tsv".format(project, sample["Sample"])) as f_annotation:
                    for row in csv.DictReader(f_annotation, delimiter="\t"):
                        ref, alt = row["Change"].split("/")
                        try:
                            counts["ti" if bases[ref] == bases[alt] else "tv"] += 1
                        except KeyError:
                            pass
                    writer.writerow({"Sample": sample["Sample"], "Run": sample["Run"], "TiTv": "{:.2f}".format(float(counts["ti"]) / (counts["tv"] or 1))})



def clean_annotations(project):
    artifacts = defaultdict(set)
    with s3_open(BUCKET, "projects/{0}/{0}.artifacts.tsv".format(project)) as f_artifacts:
        for row in csv.DictReader(f_artifacts, delimiter="\t"):
            #if int(row["Count"]) >= 1:
            artifacts[row["Run"]].add(row["Variant"])

    print("Cleaning...")
    with s3_open(BUCKET, "projects/{0}/{0}.samples.tsv".format(project)) as f_samples:
        for sample in csv.DictReader(f_samples, delimiter="\t"):
            if sample["Timepoint"].isnumeric():
                print(sample["Sample"])
                with s3_open(BUCKET, "projects/{0}/{1}/{1}.annotation.tsv".format(project, sample["Sample"])) as f_annotation:
                    reader = csv.DictReader(f_annotation, delimiter="\t")
                    with s3_open(BUCKET, "projects/{0}/{1}/{1}.annotation.cleaned.tsv".format(project, sample["Sample"]), "wt") as f_out:
                        writer = csv.DictWriter(f_out, reader.fieldnames, delimiter="\t")
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
    
    with s3_open(BUCKET, "projects/{0}/{0}.annotation.combined.tsv".format(project), "wt") as f:
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
            
    with s3_open(BUCKET, "projects/{0}/{0}.artifacts.tsv".format(project), "wt") as f:
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
            
    with s3_open(BUCKET, "projects/{0}/{0}.annotation.combined2.tsv".format(project), "wt") as f:
        writer = csv.DictWriter(f, ["Relevant", "Comment"] + reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for row in annotations:
            writer.writerow(row)



def copy_fastqs_from_basespace_to_s3(s3_project):
    pdb.set_trace()
    s3 = client("s3")
    objects = s3_list(BUCKET, "projects/{0}/samples/".format(s3_project), extension="fastq.gz")
    s3_fastqs = {key.split("/")[-1]: val for key, val in objects.items()}
    basespace_dir = mount_basespace()
    
    with s3_open(BUCKET, "projects/{0}/{0}.samples.tsv".format(s3_project)) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["Run"] and row["Sample"]:
                print(row["Sample"])
                files_dir = os.path.join(basespace_dir, "Projects", row["Run"], "Samples", row.get("Basespace_Sample", row["Sample"]), "Files")
                if os.path.isdir(files_dir):
                    for basespace_fastq in os.listdir(files_dir):
                        basespace_path = os.path.join(files_dir, basespace_fastq)
                        if basespace_fastq.endswith(".fastq.gz"):
                            s3_fastq = basespace_fastq.replace(row["Basespace_Sample"], row["Sample"]) if "Basespace_Sample" in row else basespace_fastq
                            if s3_fastq not in s3_fastqs or os.path.getmtime(basespace_path) > s3_fastqs[s3_fastq]["LastModified"].timestamp():
                                print("Uploading {}".format(basespace_fastq))
                                s3.upload_file(basespace_path, BUCKET,  "projects/{0}/samples/{1}/{2}".format(s3_project, row["Sample"], s3_fastq))



if __name__ == "__main__":
    project = "head_and_neck"
    #generate_artifact_list(project)
    #clean_annotations(project)
    #combine_annotations(project, time_points=4)
    #titv(project)
    copy_fastqs_from_basespace_to_s3(project)
    
    
    
