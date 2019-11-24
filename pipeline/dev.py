import pdb, sys, time
from collections import defaultdict, Counter

from pipeline import create_report, dedup
            
            
from covermi import Panel

def main():
    
    create_report("/home/ed/Software/data/temp/z.vep", "/home/ed/Software/data/panels/Accept")
    
    sys.exit()
    
    
    #dedup("/home/ed/Software/data/temp/3reads_R1.fastq", "/home/ed/Software/data/temp/3reads_R2.fastq", 
          #thruplex=False, 
          #allowed=3)
    #sys.exit()
    
    dedup("/home/ed/Software/data/temp/HNC006-HNC006-c-0_S2_L001_R1_001.quarter.fastq", "/home/ed/Software/data/temp/HNC006-HNC006-c-0_S2_L001_R2_001.quarter.fastq", 
          thruplex=True, 
          allowed=3)
    sys.exit()
    
    
    dedup("/home/ed/Software/data/temp/10010014-H3731-c-0_S1_L001_R1_001.fastq", "/home/ed/Software/data/temp/10010014-H3731-c-0_S1_L001_R2_001.fastq",
          thruplex=False,
          allowed=3)
    sys.exit()
    
    
    # fastqs sorted into r1, r2 pairs
    fastqs = sorted(["/home/ed/Software/data/temp/10010014-H3731-c-0_S1_L001_R1_001.fastq.gz", "/home/ed/Software/data/temp/10010014-H3731-c-0_S1_L001_R2_001.fastq.gz"])
    paired_fastqs = ungzip_and_combine_fastqs(*list(zip(fastqs[::2], fastqs[1::2])))
    
    rows = 0
    for read1, read2 in paired_fastqs:
        with open(read1, "rt") as f:
            for row in f:
                rows += 1
    total_reads = rows // 4
    print (total_reads)
    
    
    
    
    




if __name__ == "__main__":
    main()
    #cProfile.run("main()")
