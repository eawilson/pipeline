import sys



def main():
    for row in sys.stdin:
        if row.startswith(">") and not row.startswith(">chr"):
            row = ">chr" + row[1:]
        sys.stdout.write(row)



if __name__ == "__main__":
    main()
