import math
import argparse
import sys



def vcf_pvalue_2_phred(input_vcf, output_vcf):
    """ Varscan annoyingly produces a p-value rather than phred score.
        This is written to the format section of the vcf and not the 
        qual field which is dotted out. This function will write
        a new copy of the vcf file with the qual field completed correctly.
        If the vcf is malformed and the PVAL is not stored in the format
        fields in column 9 then no changes are made any the line is copied
        as is.
    
    Args:
        input_vcf (str): Name of input vcf file.
        output_vcf (str): Name of output vcf file.
    Returns:
        None
    Raises:
        Will raise the usual OSErrors if input_vcf cannot be read
        or output_vcf cannot be written.
    """

    if not input_vcf.endswith(".vcf"):
        sys.exit("Input must be a .vcf file")
    if not output_vcf.endswith(".vcf"):
        sys.exit("Output must be a .vcf file")
    
    with open(input_vcf, "rt") as f_in:
        with open(output_vcf, "wt") as f_out:
            for row in f_in:
                row = row.split("\t")
                if not row[0].startswith("#"):
                    if row[5] == "." and len(row) > 9:
                        format_keyvals = dict(zip(row[8].split(":"),
                                                  row[9].split(":")))
                        try:
                            p_value = float(format_keyvals.get("PVAL", ""))
                        except ValueError:
                            pass
                        else:
                            if p_value == 0:
                                phred = 255
                            else:
                                phred = int(-10 * math.log(float(p_value), 10))
                                if phred > 255:
                                    phred = 255
                            row[5] = str(phred)
                f_out.write("\t".join(row))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_vcf', help="Input vcf file.")
    parser.add_argument("-o", "--output", help="Output vcf file.", dest="output_vcf", default=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        vcf_pvalue_2_phred(**vars(args))
    except OSError as e:
        # File input/output error. This is not an unexpected error so just
        # print and exit rather than displaying a full stack trace.
        sys.exit(str(e))



if __name__ == "__main__":
    main()

