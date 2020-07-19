import os, sys, pdb
from collections import Counter, defaultdict, OrderedDict
from itertools import chain
import csv
from statistics import median 

import boto3

stem = "/home/ed/Software/data/temp/accept"
HET = "O"
HOM = "X"



class OrderedDefaultdict(OrderedDict):
    def __missing__(self, key):
        self[key] = {}
        return self[key]
    
    

def het_hom_table(outlier_sample_id=None):
    s3 = boto3.client("s3")
    objects = [content["Key"] for content in s3.list_objects_v2(Bucket="omdc-data", Prefix="projects/accept")["Contents"]]

    outlier = OrderedDefaultdict()
    patients = defaultdict(lambda:OrderedDefaultdict())
    samples = [obj for obj in objects if obj.endswith(".annotation.tsv")]
    for obj in samples:
        patient_id, sample_id, cee, timepoint = os.path.basename(obj).split(".")[0].split("-")
        if cee == "g":
            timepoint = "g"
        
        path_in = os.path.join(stem, os.path.basename(obj))
        if not os.path.exists(path_in):
            s3.download_file('omdc-data', obj, path_in)
        
        with open(path_in, "rt") as f_in:
            reader = csv.DictReader(f_in, delimiter="\t")
            fieldnames = reader.fieldnames
            for row in reader:
                variant = "{}:{} {}".format(row["Chrom"], row["Pos"], row["Change"])
                if float(row["VAF"]) >= .35:
                    mark = HET
                    if float(row["VAF"]) > .9:
                        mark = HOM
                    (patients[patient_id] if sample_id != outlier_sample_id else outlier)[(patient_id, sample_id, timepoint)][variant] = mark
    
    if not outlier_sample_id:
        for key in list(patients):
            if len(patients[key]) == 1:
                del patients[key]

    for samples in patients.values():
        for k, v in outlier.items(): # there will only ever be one item
            samples[k] = v
        
        header = sorted(set(v for v in chain(*samples.values())))
        rows = []
        for variants in samples.values():
            rows += [[variants[h] if h in variants else " " for h in header]]
        
        hethom = defaultdict(list)
        for calls in zip(*rows):
            hethom[HOM if HOM in calls else HET] += [calls]
        
        for sample, hom, het in zip(*([samples.keys()] + [["".join(row) for row in zip(*hethom[zyg])] for zyg in (HOM, HET)])):
            print("{:<10s}{:<8s}{}   {}   {}".format(sample[0], sample[1], sample[2], hom, het))
        print(".")
            
            
    
    
    #out_file = os.path.join(stem, "accept.annotation.combined.tsv")
    #with open(out_file, "wt") as f:
        #writer = csv.DictWriter(f, ["Variant", "Patient_Id", "Timepoint"] + fieldnames, delimiter="\t")
        #writer.writeheader()
        #for row in sorted(annotations, key=lambda a: (a["Variant"], a["Patient_Id"], a["Timepoint"])):
            #if patients[row["Patient_Id"]] == 4:
                #row.pop("Quality", None)
                #writer.writerow(row)
            
    #print("Uploading to S3.")
    #s3.upload_file(out_file, "omdc-data", "projects/accept/accept.annotation.combined.tsv")



























if __name__ == "__main__":
    het_hom_table()#outlier_sample_id="H3272")
