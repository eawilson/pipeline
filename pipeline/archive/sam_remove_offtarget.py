from covermi import Entry, Chrom, FileContext
import csv
import sys
import pdb


def sam_remove_offtarget(f_in, f_out, ontarget):
    """ Read sam file fn_in and copy all reads that touch ontarget to fn_out.
    
    Args:
        fn_in (str): Name or file_obj of input sam file, '-' for stdin.
        fn_out (str): Name or file_obj of output sam file, '-' for stdout.
        ontarget (covermi.Gr): Covermi genomic range of ontarget coordinates.
        
    Raises:
        Usual filesystem exceptions if fn_in or fn_out cannot be read/written.
        ValueError, RuntimeError if malformed input sam file.
    """
    to_close = []
    try:
        if not hasattr(f_in, "read"):
            if f_in == "-":
                f_in = sys.stdin
            else:
                f_in = open(f_in, "rt")
                to_close += [f_in]
        if not hasattr(f_out, "write"):
            if f_out == "-":
                f_out = sys.stdout
            else:
                f_out = open(f_out, "wt")
                to_close += [f_out]
        
        reader = csv.reader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")

        total_reads = 0
        offtarget_reads = 0
        for row in reader:
            if not row[0].startswith("@"):
                total_reads += 1
                try:
                    chrom = Chrom(row[2])
                except KeyError:
                    continue
                
                start = stop = int(row[3])
                cigar = row[5]
                if cigar == "*":
                    continue
                
                number = ""
                for char in cigar:
                    if char.isnumeric():
                        number += char
                    else:
                        if char in ("M", "D", "N", "=", "X"):
                            stop += int(number)
                        elif char not in ("I", "S", "H", "P"):
                            raise RuntimeError("Malformed cigar string.")
                        number = ""
                
                if not ontarget.touched_by(Entry(chrom, start, stop)):
                    offtarget_reads += 1
                    continue
                
            writer.writerow(row)
    finally:
        for f in to_close:
            f.close()

    percent_reduction = offtarget_reads * 100 // total_reads
    print("Removed {} offtarget reads ({}% of {})". \
        format(offtarget_reads,percent_reduction,total_reads), file=sys.stderr)
    


def commandline():
    pass


