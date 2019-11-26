import os
import sys
import time
import datetime
import contextlib
from pipeline import run, mount_basespace, mount_instance_storage, unmount_basespace, list_basespace_fastqs, ungzip_and_combine_illumina_fastqs, load_panel_from_s3, \
                    s3_put, dedup, illumina_readgroup, pipe, create_report, s3_object_exists
from covermi import Panel, covermimain







def cfpipeline(basespace_project_regex="", sample_regex="", s3_project=None, panelname=None):
    if s3_object_exists("projects/{}/{}".format(s3_project, sample_regex)):
        print "{} already exists, skipping.".format(s3_project, sample_regex)
        return
    
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

        try:
            sample = os.path.splitext(os.path.basename(r1_fastq))[0]
            with open("{}_pipeline.txt".format(sample), "wb") as f_report:
                    
                start_time = time.time()
                print("cfPipeline {}, {}.".format(os.path.basename(r1_fastq), os.path.basename(r2_fastq)))
                print("Starting {}.".format(datetime.datetime.now()))
                print("Shaw allowed = 3, thruplex = {}, min_family_size = {}.".format(prop.get("thruplex", False), prop.get("min_family_size", 1)))
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
                covermi_file = "{}.tar.gz".format(covermi_dir)
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
    cfpipeline(basespace_project_regex="CAPP", sample_regex="10440010-H31479-g-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="CAPP", sample_regex="10460004-H13539-g-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="CAPP", sample_regex="10460004-H16528-g-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10010006-H18838-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10010006-H21802-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10010006-H23564-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10010006-H4615-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10010012-H2236-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10010012-H4128-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10010012-H6615-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10010012-H20580-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10010014-H22658-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10010014-H8033-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10010014-H6036-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10010014-H3731-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10010024-H27176-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10010024-H29436-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10010024-H31762-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10010024-H9061-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ007_DLBCL", sample_regex="10010026-H30277-g-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10010028-H35380-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10010028-H137-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10010028-H2623-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10030011-H18007-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10030011-H2981-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10040020-H21217-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10040020-H24561-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10040020-H26638-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10040020-H2319-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10090001-H13538-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10090001-H16039-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10090001-H17975-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10090001-H25891-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10090001-H13538-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10090001-H16039-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10090001-H25891-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10090001-H17975-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="1009000-H26048-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10090007-H6232-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10090007-H21802-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10090007-H24313-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10090007-H24313-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10090007-H26048-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10090007-H21803-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10090011-H643-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10090011-H32165-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="100900013-H3138-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="100900013-H6158-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="100900013-H8204-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="100900013-H22931-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ007_DLBCL", sample_regex="10090013-H3139-g-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10090023-H25954-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10090023-H28510-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10090023-H30654-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10090023-H8370-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10090025-H27445-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10090025-H30844-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10090025-H33095-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10400005-H18923-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10400005-H23055-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10400005-H4928-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10040005-H21161-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10400005-H23055-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10400005-H18923-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10400005-H21161-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10440010-H33358-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10440010-H31478-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10440010-H1798-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10440010-H15654-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10460004-H21001-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10460004-H32031-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19-CAPPSeq-5-cf", sample_regex="10460004-H16527-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT-Cf", sample_regex="10460004-H19078-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10460004-H21001-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10460004-H19078-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10460004-H32031-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="ACCEPT1", sample_regex="10460004-H16527-c", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10460015-H13440-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10460015-H15640-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10460015-H17828-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10460015-H33712-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10460022-H25388-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10460022-H27835-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10460022-H29964-c-2", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ008_DLBCL", sample_regex="10460022-H34466-c-3", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10010028-H32577-c-0", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10010028-H35990-c-1", s3_project="accept", panelname="Accept")
    cfpipeline(basespace_project_regex="19CAPPSEQ009_DLBCL", sample_regex="10010028-H539-c-2", s3_project="accept", panelname="Accept")    



