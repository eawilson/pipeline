import os
import sys
import pdb
from collections import Counter, defaultdict
import csv
import io
from statistics import median 
from boto3 import client    

BUCKET = "omdc-data"

__all__ = ["s3_put", "s3_object_exists", "s3_get_tsv", "s3_list_keys", "s3_list_samples", "s3_UploadFileObj", \
            "clean_annotations", "combine_annotations", "download_annotations", "generate_artifact_list"]





class TsvRows(list):
        pass



def s3_get_tsv(bucket, key):
    s3 = client("s3")

    rows = TsvRows()
    with io.BytesIO() as f:
        try:
            response = s3.download_fileobj(bucket, key, f)
        except s3.exceptions.ClientError: # file does not exist
            raise FileNotFoundError("Unable to retrieve {} from {} bucket".format(key, bucket))

        f.seek(0)
        with io.TextIOWrapper(f) as f2:
            reader = csv.DictReader(f2, delimiter="\t")
            for row in reader:
                rows.append(row)
                
    rows.fieldnames = reader.fieldnames
    return rows



class s3_UploadFileObj(object):
    def __init__(self, bucket, key, text=True):
        self.bucket = bucket
        self.key = key
        self.f_bytes = io.BytesIO()
        self.f_in = io.TextIOWrapper(self.f_bytes, write_through=True) if text == True else self.f_bytes
    
    def upload(self):
        self.f_bytes.seek(0)
        client("s3").upload_fileobj(self.f_bytes, self.bucket, self.key)
        self.close()
        
    def close(self):
        self.f_in.close()
            
    def __enter__(self):
        return self.f_in
    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.upload()



def titv(project):
    base = {"A": "U", "G": "U", "C": "Y", "T": "Y"}
    print("TiTv...")
    with s3_UploadFileObj(BUCKET, "projects/{0}/{0}.titv.tsv".format(project)) as f:
        writer = csv.DictWriter(f, ["Sample", "Run", "TiTv"], delimiter="\t")
        writer.writeheader()
        
        for sample in s3_get_tsv(BUCKET, "projects/{0}/{0}.samples.tsv".format(project)):
            print(sample["Sample"])
            counts = Counter()
            reader = s3_get_tsv(BUCKET, "projects/{0}/{1}/{1}.annotation.tsv".format(project, sample["Sample"]))
            for row in reader:
                ref, alt = row["Change"].split("/")
                print(ref, alt)
                try:
                    counts["ti" if base[ref] == base[alt] else "tv"] += 1
                except KeyError:
                    pass
            print(counts)
            writer.writerow({"Sample": sample["Sample"], "Run": sample["Run"], "TiTv": "{:.2f}".format(float(counts["ti"]) / (counts["tv"] or 1))})



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
    #generate_artifact_list(project)
    #clean_annotations(project)
    #combine_annotations(project, time_points=4)
    titv(project)
    
    
    
    
    
