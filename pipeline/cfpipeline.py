import os

from pipeline import run, mount_basespace, mount_instance_storage, unmount_basespace, list_basespace_fastqs, ungzip_and_combine_illumina_fastqs, load_panel_from_s3, \
                    s3_put, dedup, illumina_readgroup
from covermi import Panel







def main():
    threads = "4"
    project = "CAPP"
    sample = "10010024-H27176-c-0"
    panelname = "Accept"
    
    
    os.chdir(mount_instance_storage())
    
    mount_basespace()
    fastqs = list_basespace_fastqs(project=project, sample=sample)
    print("Fetching fastqs.")
    fastqs = ungzip_and_combine_illumina_fastqs(*fastqs)
    unmount_basespace()

    panel = load_panel_from_s3(panelname)
    
    for r1_fastq, r2_fastq in zip(fastqs[::2], fastqs[1::2]):
        print("Shaw.")
        dedup(r1_fastq, r2_fastq, allowed=3, thruplex=False)
        //os.unlink(r1_fastq)
        //os.unlink(r2_fastq)
        r1_fastq = "{}.deduped.fastq".format(r1_fastq[:-6])
        r2_fastq = "{}.deduped.fastq".format(r2_fastq[:-6])
    
        print("BWA mem.".format(sample))
        sam_file = "{}.sam".format(sample)
        with open(sam_file, "wb") as f:
            completed = run(["bwa", "mem", "-t", threads, "-R", illumina_readgroup(r1_fastq), panel.properties["reference_fasta"], r1_fastq, r2_fastq], stdout=f, universal_newlines=False)

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

        print("Mpileup.")
        pileup_file = "{}.pileup".format(sample)
        run(["bcftools", "mpileup", "-A", "-d", "10000000", "-Ou", "-o", pileup_file, "-f", panel.properties["reference_fasta"], bam_file])
        
        print("Varscan.")java
        vcf_file = "{}.vcf".format(sample)
        with open(vcf_file, "wb") as f:
            run(["java", "-jar", "/usr/local/bin/varscan.jar", "pileup2cns", "--variants", "--output-vcf", "1", "--min-coverage", "1",
                                                                                                                "--min-var-freq", "0", 
                                                                                                                "--min-avg-qual", "20",
                                                                                                                "--min-reads2", "1"], stdout=f, universal_newlines=False)
                                                                                                            #"--min-coverage-normal", "1",  
                                                                                                            #"--min-coverage-tumor", "1",  
                                                                                                            #"--min-freq-for-hom","0.75", 
                                                                                                            #"--somatic-p-value", "0.05",
                                                                                                            #"--strand-filter", "1",
                                                                                                            #"--validation", "1"
        
        print("Complete")






if __name__ == "__main__":
    main()




