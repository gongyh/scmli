#!/usr/bin/env python3

import os, sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


def pcr_qc(project_name, read1, read2, path_fastqc, path_trim_galore, threads):
    '''
    '''
    #fastqc
    print("fastqc......")
    if path_fastqc:
        os.system(path_fastqc + '/fastqc -o . ' + read1 + ' ' + read2 + "> pcr_pipeline.log 2>&1") 
    else:
        os.system('fastqc -o . ' + read1 + ' ' + read2 + "> pcr_pipeline.log 2>&1")

    #trim
    print("trim_galore......")
    if path_trim_galore:
        os.system(path_trim_galore + 'trim_galore --paired --fastqc --max_n 0 -j ' + threads + ' --gzip ' + read1 + ' ' + read2 + ">> pcr_pipeline.log 2>&1")
    else:
        os.system('trim_galore --paired --fastqc --max_n 0 -j ' + threads + ' --gzip ' + read1 + ' ' + read2 + ">> pcr_pipeline.log 2>&1")

    #jieya hebing
    name1 = read1.split('/')[-1]
    name1 = name1.split(".fq.gz")[0] + "_val_1.fq.gz"
    name2 = read2.split('/')[-1]
    name2 = name2.split(".fq.gz")[0] + "_val_2.fq.gz"
    os.system('gzip -cd ' + name1 + ' ' + name2 + ' > ' + project_name + '.fq')

def pcr_parse_gRNA(lib, fix_seq, number, project_name):
    '''
    XXX
    '''
    print("search......")
    header = True
    gRNA_gene = {}
    with open(lib) as fh:
        for line in fh:
            if header:
                header = False
                continue
            cl = line.split(",")
            gene = cl[0]
            gRNA = cl[1][number[0]:number[1]] 
            gRNA_gene[gRNA] = gene

    fix_seq_len = len(fix_seq)
    gRNAs_dict = {}

    ## extract core function as a function, def core_func(param1,param2): XXX
    ## define thread pool, p = Pool(threads)
    ## map and reduce, p.map(core_func, [inputs])

    with open(project_name + ".fq") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)
            if seq[0:fix_seq_len] == fix_seq: #valid record
                #change left and right equal
                gRNA = str(seq[fix_seq_len+number[0]:fix_seq_len+number[1]])
                if gRNA in gRNAs_dict.keys():
                    gRNAs_dict[gRNA] += 1
                else:
                    gRNAs_dict[gRNA] = 1

    for k, v in gRNAs_dict.items():
        gene = "unknow"
        f = open(project_name + ".counts", "a") #need fix
        if k in gRNA_gene.keys():
            gene = gRNA_gene[k]
        print("%s\t%s\t%d"%(gene,k,v), file = f)


def pcr_count(project_name) :
    '''
    XXX
    '''
    #fixed_seq, e.g. GGTAGAATTGGTCGTTGCCATCGACCAGGC
    print("stats......")
    df1 = pd.read_csv(project_name + ".counts", sep = "\t", header = None, names = ["gene_id", "sequence", "counts", "percent"])
    t = df1["counts"].sum()
    df1["percent"] = df1["counts"]/t
    df1["percent"] = df1["percent"].round(7)    
    df1 = df1.sort_values(by="counts", ascending = False)
    df1.to_csv(project_name + ".percent", sep = "\t", index = False)

    stats = {}
    stats["all_kinds"] = len(df1)
    #stats["lib_kinds"] = 9710
    #stats["un_kinds"] = len(df1[df1["gene_id"]=="unknow"])
    #stats["no_kinds"] = len(df1[df1["gene_id"]!="unknow"])
    #stats["no_sum"] = np.sum(df1[df1["gene_id"]!="unknow"]["counts"])
    #stats["no_average"] = np.average(df1[df1["gene_id"]!="unknow"]["counts"])
    #stats["no_coverage"] = stats["no_kinds"]/9710
    
    return stats
