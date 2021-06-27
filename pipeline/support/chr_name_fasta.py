import sys



def main():
    
    if len(sys.argv) == 1:
        fn = "-"
    else:
        fn = sys.argv[1]
    
    if fn == "-":
        f_in = sys.stdin
    else:
        f_in = open(fn, "rt")
    
    try:
        for row in f_in:
            if row.startswith(">"):
                split_row = row.split()
                name = split_row[0][1:]
                row = f">chr{name}\n"
            sys.stdout.write(row)
    finally:
        f_in.close()


if __name__ == "__main__":
    main()
