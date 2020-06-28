import csv
import sys
import argparse
from collections import Counter 


def sam_reads_by_contig(f_in, f_out):
    """ Read sam file fn_in and write tsv of total reads by contig to fn_out.
    
    Args:
        fn_in (str): Name or file_obj of input sam file, '-' for stdin.
        fn_out (str): Name or file_obj of output sam file, '-' for stdout.
        
    Raises:
        Usual filesystem exceptions if fn_in or fn_out cannot be read/written.
        RuntimeError if malformed input sam file.
    """
    to_close = []
    try:
        if not hasattr(f_in, "read"):
            if f_in == "-":
                f_in = sys.stdin
            else:
                f_in = open(f_in, "rt", newline="")
                to_close += [f_in]
        reader = csv.reader(f_in, delimiter="\t")

        reads = Counter()
        for row in reader:
            if not row[0].startswith("@"):
                if len(row) < 11:
                    raise RuntimeError("Malformed SAM file, only contains {} columns.".format(len(row)))
                
                reads[row[2]] += 1
                
        if not hasattr(f_out, "write"):
            if f_out == "-":
                f_out = sys.stdout
            else:
                f_out = open(f_out, "wt", newline="")
                to_close += [f_out]
        writer = csv.writer(f_out, delimiter="\t")

        for row in sorted(reads.items()):
            writer.writerow(row)

    finally:
        for f in to_close:
            f.close()
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input SAM file.")
    parser.add_argument("output", help="Output tsv file.")
    args = parser.parse_args()
    cfpipeline(args.input, args.output)



if __name__ == "__main__":
    main()
