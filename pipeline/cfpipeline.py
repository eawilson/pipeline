import os
import sys

from pipeline import run, mount_basespace, mount_instance_storage, unmount_basespace, list_basespace_fastqs, ungzip_and_combine_illumina_fastqs, load_panel_from_s3, \
                    s3_put, dedup, illumina_readgroup, pipe, create_report
from covermi import Panel, covermimain







def cfpipeline(basespace_project="", sample="", s3_project=None, panelname=None):
    if not s3_project or not panelname:
        raise RuntimeError("Project and Panel are required arguments.")
    if not basespace_project and not sample:
        raise RuntimeError("Sample or basespace project are required arguments.")
    
    threads = "4"
    
    os.chdir(mount_instance_storage())
    
    mount_basespace()
    fastqs = list_basespace_fastqs(project=basespace_project, sample=sample)
    print("Fetching fastqs.")
    fastqs = ungzip_and_combine_illumina_fastqs(*fastqs)
    #unmount_basespace()

    panel = load_panel_from_s3(panelname)
    
    for r1_fastq, r2_fastq in zip(fastqs[::2], fastqs[1::2]):
        print("{}, {}".format(os.path.basename(r1_fastq), os.path.basename(r2_fastq)))
        print("Shaw.")
        dedup(r1_fastq, r2_fastq, allowed=3, thruplex=False)
        r1_dedupfastq = "{}.deduped.fastq".format(r1_fastq[:-6])
        r2_dedupfastq = "{}.deduped.fastq".format(r2_fastq[:-6])
    
        print("BWA mem.".format(sample))
        sam_file = "{}.sam".format(sample)
        with open(sam_file, "wb") as f:
            pipe(["bwa", "mem", "-t", threads, "-R", illumina_readgroup(r1_fastq), panel.properties["reference_fasta"], r1_dedupfastq, r2_dedupfastq], stdout=f)
        os.unlink(r1_dedupfastq)
        os.unlink(r2_dedupfastq)

        print("Samtools fixmate.")
        unsorted_bam_file = "{}.unsorted.bam".format(sample)
        run(["samtools", "fixmate", "-O", "bam", sam_file, unsorted_bam_file])
        os.unlink(sam_file)
    
        print("Samtools sort.")
        bam_file = "{}.bam".format(sample)
        run(["samtools", "sort", "-O", "bam", "-o", bam_file, "-T", "temp", "-@", threads, unsorted_bam_file])
        os.unlink(unsorted_bam_file)

        print("Samtools index.")
        run(["samtools", "index", bam_file])
        bambai_file = "{}.bai".format(bam_file)

        print("Mpileup.")
        mpileup_file = "{}.pileup".format(sample)
        run(["samtools", "mpileup", "-A", "-d", "10000000", "-o", mpileup_file, "-f", panel.properties["reference_fasta"], bam_file])
        
        print("Varscan.")
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

        print("VEP.")
        vepjson_file = "{}.vep".format(sample)
        extra = ["--refseq"] if panel.properties["transcript_source"] == "refseq" else []
        run(["perl", "../ensembl-vep/vep", "--verbose", "-i", vcf_file, "-o", vepjson_file, "--no_stats", "--format", "vcf", "--fork", "4", "--json", "--offline", "--everything", 
                  "--assembly", panel.properties["assembly"], "--fasta", panel.properties["reference_fasta"], "--force_overwrite"] + extra)
        annotation_file = create_report(vepjson_file, panel)

        print("Covermi.")
        covermi_dir = covermimain(panelname, "", bam_path=bam_file)
        covermi_file = "{}.tar.gz".format(covermi_dir)
        run(["tar", "cvzf", covermi_file, covermi_dir])
        run(["rm", "-r", covermi_dir])

        filelist = [bam_file, bambai_file, vcf_file, covermi_file, vepjson_file, annotation_file]
        s3_put(*filelist, prefix="{}/{}".format(s3_project, sample))
        for filename in filelist:
            os.unlink(filename)


if __name__ == "__main__":
    for s3_project, panel, sample, basespace_project in [("head_and_neck", "Head_and_Neck", "HNC006", "HNC"),
                                                         ("head_and_neck", "Head_and_Neck", "HNC011", "HNC"),
                                                         ("head_and_neck", "Head_and_Neck", "HNC016", "HNC"),
                                                         ]:
        cfpipeline(basespace_project=basespace_project, sample=sample, s3_project=s3_project, panelname=panel)




