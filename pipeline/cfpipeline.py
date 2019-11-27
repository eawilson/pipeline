import os
import sys
import time
import datetime
import contextlib
from pipeline import run, mount_basespace, mount_instance_storage, unmount_basespace, list_basespace_fastqs, ungzip_and_combine_illumina_fastqs, load_panel_from_s3, \
                    s3_put, dedup, illumina_readgroup, pipe, create_report, s3_object_exists
from covermi import Panel, covermimain







def cfpipeline(basespace_project, sample, s3_project, panelname):

    print(basespace_project, sample)
    if s3_object_exists("projects/{}/{}".format(s3_project, sample)):
        print("{}/{} already exists, skipping.".format(s3_project, sample))
        return
    
    threads = "4"
    
    os.chdir(mount_instance_storage())    
    mount_basespace()
    
    fastqs = list_basespace_fastqs(project=basespace_project, sample=sample)
    for fastq in fastqs:
        print(fastq)
        
    print("Fetching fastqs...")
    fastqs = ungzip_and_combine_illumina_fastqs(*fastqs)

    panel = load_panel_from_s3(panelname)
    prop = panel.properties

    for r1_fastq, r2_fastq in zip(fastqs[::2], fastqs[1::2]):

        try:
            if not sample in r1_fastq:
                raise RuntimeError(Wrong fastq.)
            with open("{}.pipeline.txt".format(sample), "wb") as f_report:
                    
                start_time = time.time()
                print("cfPipeline {}, {}.".format(os.path.basename(r1_fastq), os.path.basename(r2_fastq)))
                print("Starting {}.".format(datetime.datetime.now()))
                print("theDuDe allowed = 3, thruplex = {}, min_family_size = {}.".format(prop.get("thruplex", False), prop.get("min_family_size", 1)))
                dedup(r1_fastq, r2_fastq, allowed=3, thruplex=prop.get("thruplex", False), min_family_size=int(prop.get("min_family_size", 1)))
                r1_dedupfastq = "{}.deduped.fastq".format(r1_fastq[:-6])
                r2_dedupfastq = "{}.deduped.fastq".format(r2_fastq[:-6])
                pipe(["wc", r1_fastq, r2_fastq, r1_dedupfastq, r2_dedupfastq], stdout=f_report)
                os.unlink(r1_fastq)
                os.unlink(r2_fastq)
            
                sam_file = "{}.sam".format(sample)
                with open(sam_file, "wb") as f:
                    pipe(["bwa", "mem", "-t", threads, "-R", illumina_readgroup(r1_dedupfastq), prop["reference_fasta"], r1_dedupfastq, r2_dedupfastq], stdout=f, stderr=f_report)
                os.unlink(r1_dedupfastq)
                os.unlink(r2_dedupfastq)

                unsorted_bam_file = "{}.unsorted.bam".format(sample)
                pipe(["samtools", "fixmate", "-O", "bam", sam_file, unsorted_bam_file], stderr=f_report)
                os.unlink(sam_file)
            
                bam_file = "{}.bam".format(sample)
                pipe(["samtools", "sort", "-O", "bam", "-o", bam_file, "-T", "temp", "-@", threads, unsorted_bam_file], stderr=f_report)
                os.unlink(unsorted_bam_file)

                pipe(["samtools", "index", bam_file], stderr=f_report)
                bambai_file = "{}.bai".format(bam_file)

                mpileup_file = "{}.pileup".format(sample)
                pipe(["samtools", "mpileup", "-A", "-d", "10000000", "-o", mpileup_file, "-f", prop["reference_fasta"], bam_file], stderr=f_report)
                
                vcf_file = "{}.vcf".format(sample)
                with open(vcf_file, "wb") as f:
                    pipe(["java", "-jar", "/usr/local/bin/varscan.jar", "mpileup2cns", mpileup_file, "--variants", "--output-vcf", "1", "--min-coverage", "1",
                                                                                                                                    "--min-var-freq", "0", 
                                                                                                                                    "--min-avg-qual", "20",
                                                                                                                                    "--min-reads2", "1"], stdout=f)
                                                                                                                    #"--min-coverage-normal", "1",  
                                                                                                                    #"--min-coverage-tumor", "1",  
                                                                                                                    #"--min-freq-for-hom","0.75", 
                                                                                                                    #"--somatic-p-value", "0.05",
                                                                                                                    #"--strand-filter", "1",
                                                                                                                    #"--validation", "1"
                os.unlink(mpileup_file)

                vepjson_file = "{}.vep".format(sample)
                extra = ["--refseq"] if prop["transcript_source"] == "refseq" else []
                pipe(["perl", "../ensembl-vep/vep", "--verbose", "-i", vcf_file, "-o", vepjson_file, "--no_stats", "--format", "vcf", "--fork", "4", "--json", "--offline", "--everything", 
                        "--assembly", prop["assembly"], "--fasta", prop["reference_fasta"], "--force_overwrite", ] + extra, stderr=f_report)
                annotation_file = create_report(vepjson_file, panel)

                covermi_dir = covermimain(panelname, "", bam_path=bam_file)
                covermi_file = "{}.covermi.tar.gz".format(sample)
                run(["tar", "cvzf", covermi_file, covermi_dir])
                run(["rm", "-r", covermi_dir])
                
                print("Completed {}.".format(datetime.datetime.now()))
                print("Time taken = {} seconds.".format(int(time.time() - start_time)))


            print("Uploading to s3.")
            for filename in os.listdir():
                if os.path.isfile(filename):
                    s3_put(filename, prefix="projects/{}/{}".format(s3_project, sample))

        except Exception:
            print("ERROR, SKIPPING")

        for filename in os.listdir():
            if os.path.isfile(filename):
                os.unlink(filename)

            
            
if __name__ == "__main__":
    project, sample = "", ""
    cfpipeline(basespace_project=project, sample=sample, s3_project="accept", panelname="Accept")



