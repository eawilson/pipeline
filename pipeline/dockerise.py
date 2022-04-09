import argparse
import os



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to input shell script.")
    parser.add_argument("-o", "--output", help="Path to write output Dockerfile to.", default=".")
    args = parser.parse_args()
    
    script = []
    with open(args.input, "rt") as f_in:
        for line in f_in:
            line = line.strip()
            if line and not line.startswith("#"):
                script.append(line)
    
    for i, line in enumerate(script):
        script[i] = ("RUN {}" if i == 0 else "    {}").format(line)
    
    if not os.path.basename(args.output) == "Dockerfile":
        args.output = os.path.join(args.output, "Dockerfile")
    
    with open(args.output, "wt") as f_out:
        f_out.write("FROM public.ecr.aws/lts/ubuntu:18.04_stable\n\n")
        f_out.write(" && \\\n".join(script))
        f_out.write('\n')#\nENTRYPOINT ["cfPipeline"]\n')



if __name__ == "__main__":
    main()
