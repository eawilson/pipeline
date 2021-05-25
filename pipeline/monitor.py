import os
import sys
import subprocess


def run(cmd):
    return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)


def used_storage(mount_point)
    return run(f"df | grep '{mount_point}' -").stdout.split()[2]


def used_memory()
    return run("free").stdout.splitlines()[1].split()[2]


def top_table():
    return run(f"top -bn1 -u `whoami`").stdout.splitlines()[7:]


def main():
    pre_existing_pids = set(row.split()[0] for row in top_table)
            
    mount_point = "/home/ubuntu/ephemoral"
    mount_point = "/run/user/1000"
    response = used_memory()
    print(response)

    response = used_storage()
    print(response)





















