import os
import sys
import time
import datetime
import contextlib
from pipeline import run, mount_basespace, mount_instance_storage, unmount_basespace, list_basespace_fastqs, ungzip_and_combine_illumina_fastqs, load_panel_from_s3, \
                    s3_put, dedup, illumina_readgroup, pipe, create_report
from covermi import Panel, covermimain







def cfpipeline(basespace_project_regex="", sample_regex="", s3_project=None, panelname=None):
    if not s3_project or not panelname:
        raise RuntimeError("Project and Panel are required arguments.")
    if not basespace_project_regex and not sample_regex:
        raise RuntimeError("Sample or basespace project are required arguments.")
    
    threads = "4"
    
    os.chdir(mount_instance_storage())
    
    mount_basespace()
    fastqs = list_basespace_fastqs(project=basespace_project_regex, sample=sample_regex)
    print("Fetching fastqs.")
    fastqs = ungzip_and_combine_illumina_fastqs(*fastqs)
    #unmount_basespace()

    panel = load_panel_from_s3(panelname)
    prop = panel.properties

    for r1_fastq, r2_fastq in zip(fastqs[::2], fastqs[1::2]):
        
        sample = os.path.splitext(os.path.basename(r1_fastq))[0]
        with open("{}_pipeline.txt".format(sample), "wb") as f_report:
                
            start_time = time.time()
            progress("cfPipeline {}, {}.".format(os.path.basename(r1_fastq), os.path.basename(r2_fastq)), file=f_report)
            progress("Starting {}.".format(datetime.datetime.now()), file=f_report)
            progress("Shaw allowed = 3, thruplex = {}, min_family_size = {}.".format(prop.get("thruplex", False), prop.get("min_family_size", 1)), file=f_report)
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
            run(["perl", "../ensembl-vep/vep", "--verbose", "-i", vcf_file, "-o", vepjson_file, "--no_stats", "--format", "vcf", "--fork", "4", "--json", "--offline", "--everything", 
                    "--assembly", prop["assembly"], "--fasta", prop["reference_fasta"], "--force_overwrite"] + extra, stderr=f_report)
            annotation_file = create_report(vepjson_file, panel)

            covermi_dir = covermimain(panelname, "", bam_path=bam_file)
            covermi_file = "{}.tar.gz".format(covermi_dir)
            run(["tar", "cvzf", covermi_file, covermi_dir])
            run(["rm", "-r", covermi_dir])
            
            progress("Completed {}.".format(datetime.datetime.now()), file=f_report)
            progress("Time taken = {} seconds.".format(int(time.time() - start_time)), file=f_report)


        print("Uploading to s3.")
        for filename in os.listdir():
            if os.path.isfile(filename):
                s3_put(filename, prefix="projects/{}/{}".format(s3_project, sample))
                os.unlink(filename)


if __name__ == "__main__":
    cfpipeline(basespace_project_regex="CAPP", sample_regex="10010014-H3731-c-0", s3_project="accept", panelname="Accept")
    sys.exit()

    for s3_project, panel, sample, basespace_project in [("head_and_neck", "Head_and_Neck", "HNC006-HNC006-c-0", "HNC"),
                                                         ("head_and_neck", "Head_and_Neck", "HNC011-HNC011-c-0", "HNC"),
                                                         ("head_and_neck", "Head_and_Neck", "HNC016-HNC016-c-0", "HNC"),
                                                         ]:
        cfpipeline(basespace_project_regex=basespace_project, sample_regex=sample, s3_project=s3_project, panelname=panel)




