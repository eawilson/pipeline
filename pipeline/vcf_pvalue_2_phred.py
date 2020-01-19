import csv
import math



def vcf_pvalue_2_phred(fn_in, fn_out):
    """ Varscan annoyingly produces a p-value rather than phred score.
        This is written to the format section of the vcf and not the 
        qual field which is dotted out instead. This function will write
        a new copy of the vcf file with the qual field completed correctly.
    
    Args:
        fn_in (str): Name of input vcf file.
        fn_out (str): Name of output vcf file.
    Returns:
        None
    Raises:
        Will raise the usual filesystem errors if fn_in cannot be read 
        or fn_out cannot be written.
    """

    with open(fn_in, "rt") as f_in:
        reader = csv.reader(f_in, delimiter="\t")
        with open(fn_out, "rt") as f_out:
            writer = csv.writer(f_out, delimiter="\t")

            for row in reader:
                if not row[0].startswith("#"):
                    if row[5] == "." and len(row) > 9:
                        format_keyvals = dict(zip(row[8].split(":"),
                                                  row[9].split(":")))
                        p_value = format_keyvals.get("PVAL", None)
                        if p_value is not None:
                            try:
                                phred = -10 * math.log(float(p_value), 10)
                                row[5] = phred
                            except ValueError:
                                pass
                writer.writerow(row)
