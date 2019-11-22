import os

from pipeline import mount_instance_storage, mount_basespace, unmount_basespace, list_basespace_fastqs, copy_fastqs, reference_genome, run, dedup








def main():
    os.chdir(mount_instance_storage())

    mount_basespace()
    matches = list_basespace_fastqs(project="CAPP", sample="10010024-H27176-c-0")
    sample, paths = list(matches.items())[0]
    print("Fetching fastqs.")
    fastqs = copy_fastqs(sample, paths, ".")
    unmount_basespace()

    print("Aligning {} with BWA mem.".format(sample))
    sam_file = "{}.sam".format(sample)
    with open(sam_file, "wb") as f:
        completed = run(["bwa", "mem", "-t", "4", reference_genome("grch38")] + fastqs, stdout=f)

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




