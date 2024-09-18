import os
import sys
import argparse
import re
import shutil
from concurrent.futures import ThreadPoolExecutor,as_completed
import subprocess
import json
import time
from Bio import SeqIO


def execute_task(file):
    #gene geneLength mapped unmapped
    name = os.path.basename(file).replace(".stat","")
    results = {}
    count_dict = {}
    info_dict = {}
    f=open(file)
    lines = f.readlines()
    for line in lines:
        data = line.split('\t')
        all_data = data[0].split('~')
        acc_name = all_data[0]
        if data[0] != "*":
            gene_length = float(data[1])
            reads_count = float(data[2])
            total_reads = reads_dict[name]
            gene_name = data[0]
            if gene_length>0 and reads_count>0:
                rpkm = reads_count*1000000/float(total_reads)*1000/gene_length
                if gene_name not in info_dict:
                    info_dict[gene_name] = name+'\t'+gene_name
                if gene_name not in count_dict:
                    count_dict[gene_name] = 0
                count_dict[gene_name] += rpkm
    f.close()
    for gene_name,rpkm in count_dict.items():
        if name not in results:
            results[name] = []
        results[name].append(info_dict[gene_name]+'\t'+str(rpkm))
    return results

def parallel_task(task_list,params='',params_list='',func=execute_task,threads_num=20,interval=0):
    executor=ThreadPoolExecutor(max_workers=threads_num)
    result = []
    futures=[]
    for i,task in enumerate(task_list):
        if params_list:
            future=executor.submit(func,task,params_list)
        if params:
            future=executor.submit(func,task,params)
        else:
            future=executor.submit(func,task)
        time.sleep(interval)
        futures.append(future)
    for future in as_completed(futures):
        status = future.done()
        data = future.result()
        result.append(data)
    executor.shutdown(True)
    return result

def makedir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def stat_files(file):
    p = subprocess.run(f"grep -v -P '0\\t0' {file}|wc -l",shell=True,stdout=subprocess.PIPE)
    return p.stdout.decode("utf-8").strip()


#fq1 fq2 reads_count
database = os.path.abspath(sys.argv[1])
fqfile = os.path.abspath(sys.argv[2])
resultdir = os.path.abspath(sys.argv[3])
command = sys.argv[4]
scripts_dir = os.path.dirname(os.path.abspath(__file__))
threads = 5
outdir = os.path.join(resultdir,"map")
shelldir = os.path.join(resultdir,"shell")
makedir(resultdir)
makedir(outdir)
makedir(shelldir)
shell_file = os.path.join(shelldir,"run.sh")

if command == "map":
    w = open(shell_file,'w+')
    f = open(fqfile)
    lines = f.readlines()
    for line in lines:
        data = line.strip().split(" ")
        name = os.path.basename(data[0]).split('_')[0]
        if len(data)==3:
            fq1 = data[0]
            fq2 = data[1]
            cmd = f"bwa mem -t {threads} {database} {fq1} {fq2}|samtools view -t 4 -bS /dev/stdin -o {outdir}/{name}.bam && samtools view -bF 4 {outdir}/{name}.bam >{outdir}/{name}.mapped.bam && samtools sort --threads 4 {outdir}/{name}.mapped.bam -o {outdir}/{name}.sort.mapped.bam && samtools index  {outdir}/{name}.sort.mapped.bam && samtools idxstats {outdir}/{name}.sort.mapped.bam >{outdir}/{name}.stat && rm {outdir}/{name}.bam {outdir}/{name}.mapped.bam"
        elif len(data)==2:
            fq1 = data[0]
            cmd = f"bwa mem -t {threads} {database} {fq1} |samtools view -t 4 -bS /dev/stdin -o {outdir}/{name}.bam && samtools view -bF 4 {outdir}/{name}.bam >{outdir}/{name}.mapped.bam && samtools sort --threads 4 {outdir}/{name}.mapped.bam -o {outdir}/{name}.sort.mapped.bam && samtools index  {outdir}/{name}.sort.mapped.bam && samtools idxstats {outdir}/{name}.sort.mapped.bam >{outdir}/{name}.stat && rm {outdir}/{name}.bam {outdir}/{name}.mapped.bam"
        stat_file = f"{outdir}/{name}.stat"
        if not os.path.exists(stat_file):
            w.write(cmd+'\n')
        else:
            if stat_files(stat_file) =='0':
                w.write(cmd+'\n')
    f.close()
    w.close()
elif command == "rpkm":
    reads_dict = {}
    f = open(fqfile)
    lines = f.readlines()
    for line in lines:
        data = line.strip().split(" ")
        name = os.path.basename(data[0]).split('_')[0]
        if len(data)==3:
            fq1 = data[0]
            fq2 = data[1]
            reads_dict[name] = float(data[2])
        elif len(data)==2:
            fq1 = data[0]
            reads_dict[name] = float(data[1])
    f.close()

    task_list = []
    for file in os.listdir(outdir):
        if re.search("stat$",file):
            task_list.append(os.path.join(outdir,file))
    all_results = parallel_task(task_list,threads_num = 10)
    outfile = os.path.join(resultdir,"map_rpkm.tsv")
    with open(outfile,'w') as w:
        w.write("Sample\tGenename\tRPKM\n")
    for single_result in all_results:
        for name,result_list in single_result.items():
            for results1 in result_list:
                data = results1.split('\t')
                with open(outfile,'a') as w:
                    w.write(results1+'\n')

    