import os

from pipeline import run, mount_basespace, mount_instance_storage, unmount_basespace, list_basespace_fastqs, ungzip_and_combine_illumina_fastqs, load_panel_from_s3, s3_put, dedup
from covermi import Panel, 







def main():
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
        dedup(r1_fastq, r2_fastq, allowed=3, thruplex=False)
    
        print("Aligning {} with BWA mem.".format(sample))
        sam_file = "{}.sam".format(sample)
        with open(sam_file, "wb") as f:
            completed = run(["bwa", "mem", "-t", "4", reference_genome("grch38")] + fastqs, stdout=f, universal_newlines=False)

        print("Converting sam to bam.")
        unsorted_bam_file = "{}.unsorted.bam".format(sample)
        with open(unsorted_bam_file, "wb") as f:
            run(["samtools", "view", "-S", "-b", sam_file], stdout=f)
        os.unlink(sam_file)
    
        print("Sorting bam.")
        bam_file = "{}.bam".format(sample)
        run(["samtools", "sort", unsorted_bam_file, "-o" bam_file])
        os.unlink(unsorted_bam_file)



if __name__ == "__main__":
    main()




