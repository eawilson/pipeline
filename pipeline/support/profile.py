#!/usr/bin/env python3

import os
import sys
import subprocess
import csv
import time



class Cpu(object):
    def __init__(self):
        table = top_table()
        user_i = table[0].index(" USER")
        self.pre_existing_pids = set(row[:user_i].strip() for row in table[1:])
        
    
    def usage(self):
        names = set()
        cpu_total = 0
        
        table = top_table()
        user_i = table[0].index(" USER")
        cpu_i = table[0].index("%CPU", user_i)
        command_i = table[0].index("COMMAND", cpu_i)
        for row in table[1:]:
            if row[:user_i].strip() not in self.pre_existing_pids:
                name = row[command_i:]
                if name != "top":
                    cpu = float(row[cpu_i+1:row.index(" ", cpu_i+1)])
                    if cpu > 0:
                        names.add(name)
                        cpu_total += cpu
        return [cpu_total, ",".join(sorted(names))]



def run(cmd):
    return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)



def storage():
    table = run("df").stdout.splitlines()
    used_i = table[0].index("1K-blocks") + 10
    avail_i = table[0].index("Available", used_i)
    used = 0
    for row in table[1:]:
        used += int(row[used_i:avail_i-1])
    return used // 1024



def memory():
    return int(run("free").stdout.splitlines()[1].split()[2]) // 1024



def top_table():
    return run(f"top -bn1 -w 512 -u `whoami`").stdout.splitlines()[6:]



def main():
    """ Wrapper around a pipeline (or any other script) to profile
        the running code and record memory, cpu and storage usage
        every 5 seconds and store the output in a tsv file. If an 
        --output (-o) argument is present then the tsv file is stored
        in that directory, otherwise the current working directory is
        used. The tsv file is named after the --name (-n) argument if
        present or the first positional argument if not preceeded by
        any keyword arguments.
    """
    args = sys.argv[1:]
    if len(args) == 0:
        sys.exit("profile: No command line arguments")
    
    try:
        i = args.index("-o")
    except ValueError:
        try:
            i = args.index("--output")
        except ValueError:
            i = len(args)
    if i + 1 < len(args):
        output_dir = args[i + 1]
    else:
        output_dir = "."
        
    try:
        i = args.index("-n")
    except ValueError:
        try:
            i = args.index("--name")
        except ValueError:
            i = len(args)
    if i + 1 < len(args):
        name = args[i + 1]
    elif len(args) > 1 and not args[1].startswith("-"):
        name = args[1].split("/")[-1].split(".")[0]
    else:
        name = "script"
    
    with open(os.path.join(output_dir, f"{name}.profile.tsv"), "wt", newline="") as f:
        out_tsv = csv.writer(f, delimiter="\t")
        out_tsv.writerow(["time(s)", "storage(MB)", "memory(MB)", "cpu(%)", "processes"])

        cpu = Cpu()
        base_storage = storage()
        base_memory = memory()
        base_time = int(time.time())
        process = subprocess.Popen(args)
        retcode = None
        while retcode is None:
            out_tsv.writerow([int(time.time()) - base_time, storage() - base_storage, memory() - base_memory] + cpu.usage())
            time.sleep(5)
            retcode = process.poll()
    
    sys.exit(retcode)
        
        

if __name__ == "__main__":
    main()
