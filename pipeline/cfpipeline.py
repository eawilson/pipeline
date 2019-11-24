import os

from pipeline import run, mount_basespace, mount_instance_storage, unmount_basespace, list_basespace_fastqs, ungzip_and_combine_illumina_fastqs, load_panel_from_s3, \
                    s3_put, dedup, illumina_readgroup
from covermi import Panel







def main():
    s3_project = "EdTest"
    threads = "4"
    project = "CAPP"
    sample = "10010024-H27176-c-0"
    panelname = "Accept"
    
    
    os.chdir(mount_instance_storage())
    
    mount_basespace()
    fastqs = list_basespace_fastqs(project=project, sample=sample)
    print("Fetching fastqs.")
    fastqs = ungzip_and_combine_illumina_fastqs(*fastqs)
    #unmount_basespace()

    panel = load_panel_from_s3(panelname)
    
    for r1_fastq, r2_fastq in zip(fastqs[::2], fastqs[1::2]):
        print("Shaw.")
        dedup(r1_fastq, r2_fastq, allowed=3, thruplex=False)
        #os.unlink(r1_fastq)
        #os.unlink(r2_fastq)
        r1_fastq = "{}.deduped.fastq".format(r1_fastq[:-6])
        r2_fastq = "{}.deduped.fastq".format(r2_fastq[:-6])
    
        print("BWA mem.".format(sample))
        sam_file = "{}.sam".format(sample)
        with open(sam_file, "wb") as f:
            run(["bwa", "mem", "-t", threads, "-R", illumina_readgroup(r1_fastq), panel.properties["reference_fasta"], r1_fastq, r2_fastq], stdout=f, universal_newlines=False)

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
        run(["samtools", "mpileup", "-A", "-d", "10000000", "-o", pileup_file, "-f", panel.properties["reference_fasta"], bam_file])
        
        print("Varscan.")
        vcf_file = "{}.vcf".format(sample)
        with open(vcf_file, "wb") as f:
            run(["java", "-jar", "/usr/local/bin/varscan.jar", "mpileup2cns", mpileup_file, "--variants", "--output-vcf", "1", "--min-coverage", "1",
                                                                                                                            "--min-var-freq", "0", 
                                                                                                                            "--min-avg-qual", "20",
                                                                                                                            "--min-reads2", "1"], stdout=f, universal_newlines=False)
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
                  "--assembly", panel.properties["assembly"], "--fasta", panel.properties["reference_fasta"] + extra]




        print("Covermi.")
        covermi_dir = covermimain(panelname, "", bam_path=bam_file)
        covermi_file - "{}.tar.gz".format(covermi_dir)
        run(["tar", "cvzf", covermi_file, covermi_dir])

        s3_put(bam_file, bambai_file, vcf_file, covermi_file, vep_file, prefix="{}/{}".format(s3_project, sample))
        print("Complete")


if __name__ == "__main__":
    main()




