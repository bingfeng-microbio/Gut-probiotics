import os
import sys
import argparse
import re
import shutil
from concurrent.futures import ThreadPoolExecutor
import subprocess



def makedir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)

def arg_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq1","-1",help="input reads1")
    parser.add_argument("--fastq2","-2",help="input reads2")
    parser.add_argument("--input","-i",help="input file contains fastq path;reads1\\treads2")
    parser.add_argument("--outdir","-o",help="output dirname")
    parser.add_argument("--threads","-t",help="thread num",default=4)
    parser.add_argument("--fastp_path","-f",help="fastp path")
    args = parser.parse_args()
    return args

def main():
    outdir = os.path.abspath(args.outdir)
    makedir(outdir)
    shell_dir = os.path.join(outdir,"shell")
    result_dir = os.path.join(outdir,"result")
    if os.path.exists(shell_dir):
        shutil.rmtree(shell_dir)
    makedir(shell_dir)
    makedir(result_dir)
    shell_file = os.path.join(shell_dir,"run.sh")
    w=open(shell_file,'w+')
    if args.input:
        f=open(args.input)
        lines = f.readlines()
        for line in lines:
            data = line.strip().split("\t")
            name = os.path.basename(data[0]).split('_')[0]
            if len(data)==2:
                fq1 = data[0]
                fq2 = data[1]
                fq1name = os.path.basename(fq1)
                fq2name = os.path.basename(fq2)
                cmd = f"cd {result_dir} && {fastp} -w {threads} -i {fq1} -I {fq2} -o {result_dir}/{fq1name} -O {result_dir}/{fq2name} -j {result_dir}/{name}.json -h {result_dir}/{name}.html"
            elif len(data)==1:
                fq1 = data[0]
                fq1name = os.path.basename(fq1)
                cmd = f"cd {result_dir} && {fastp} -w {threads} -i {fq1} -o {result_dir}/{fq1name} -j {result_dir}/{name}.json -h {result_dir}/{name}.html"
            final_file = os.path.join(result_dir,name+".html")
            if not os.path.exists(final_file):
                w.write(cmd+'\n')
        f.close()
    else:
        fq1 = args.fastq1
        fq2 = args.fastq2
        if fq1 and fq2:
            name = os.path.basename(fq1).split('_')[0]
            fq1name = os.path.basename(fq1)
            fq2name = os.path.basename(fq2)
            cmd = f"cd {result_dir} && {fastp} -w {threads} -i {fq1} -I {fq2} -o {result_dir}/{fq1name} -O {result_dir}/{fq2name} -j {result_dir}/{name}.json -h {result_dir}/{name}.html"
            final_file = os.path.join(result_dir,name+".html")
            if not os.path.exists(final_file):
                w.write(cmd+'\n')
        elif fq1 and not fq2:
            name = os.path.basename(fq1).split('_')[0]
            fq1name = os.path.basename(fq1)
            cmd = f"cd {result_dir} && {fastp} -w {threads} -i {fq1} -o {result_dir}/{fq1name} -j {result_dir}/{name}.json -h {result_dir}/{name}.html"
            final_file = os.path.join(result_dir,name+".html")
            if not os.path.exists(final_file):
                w.write(cmd+'\n')
    w.close()
    f=open(shell_file)
    lines = f.readlines()
    f.close()
    task_length = len(lines)
    print(f"you have {task_length} tasks not complelted! Please run {shell_file}")


if __name__ == '__main__':
    args = arg_parse()
    threads = args.threads
    if int(threads)>6:
        threads = "4"
    fastp = args.fastp_path
    main()