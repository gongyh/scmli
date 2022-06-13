#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from functools import partial

def pcr_qc(project_name, read1, read2, FASTQC_PATH=None, TRIM_GALORE_PATH=None, threads=8):
    '''
    read1 = args.read1 (file)
    read2 = args.read2 (file)
    '''
    #fastqc
    print("fastqc......")
    if FASTQC_PATH:
        FASTQC_BIN=FASTQC_PATH+"/fastqc"
    else:
        FASTQC_BIN="fastqc"
    os.system(FASTQC_BIN + ' -o . ' + read1 + ' ' + read2 + "> pcr_pipeline.log 2>&1") 

    #trim
    print("trim_galore......")
    if TRIM_GALORE_PATH:
        TRIM_GALORE_BIN=TRIM_GALORE_PATH+"/trim_galore"
    else:
        TRIM_GALORE_BIN="trim_galore"
    os.system(TRIM_GALORE_BIN + ' --paired --fastqc --max_n 0 -j ' + str(threads) + ' --gzip ' + read1 + ' ' + read2 + ">> pcr_pipeline.log 2>&1")
 
   #jieya hebing
    name1 = read1.split('/')[-1]
    name1 = name1.split(".fq.gz")[0] + "_val_1.fq.gz"
    name2 = read2.split('/')[-1]
    name2 = name2.split(".fq.gz")[0] + "_val_2.fq.gz"
    os.system('gzip -cd ' + name1 + ' ' + name2 + ' > ' + project_name + '.fq')

def pcr_parse_gRNA(lib, fix_seq, number=[25,45], project_name="my_project", threads=8):
    
    '''
    lib = args.lib (file.csv)
    fix_seq = args.seq ("GGTAGAATTGGTCGTTGCCATCGACCAGGC")
    
    '''

    print("search......")
    
    pool = Pool(threads)
    gRNA_gene = {}
    
    with open(lib) as fh:
        result = pool.map(partial(search_a,number),[line for line in fh])
        gRNA_gene.update(result)

    fix_seq_len = len(fix_seq)
    gRNAs_dict = {}
    
    with open("../output/my_project.fq") as handle:
        gRNA = Pool(6).map(partial(search_b,fix_seq,fix_seq_len,number),[record for record in SeqIO.parse(handle, "fastq")])
        for key in gRNA:
            gRNAs_dict[key] = gRNAs_dict.get(key, 0) + 1
        del gRNAs_dict[None]

    '''
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
    '''

    with open(project_name+".counts", "w") as f:
        for k, v in gRNAs_dict.items():
            gene = "unknow"
            if k in gRNA_gene.keys():
                gene = gRNA_gene[k]
            f.write("%s\t%s\t%d\n"%(gene,k,v))


def search_a(number,line):
    cl = line.split(",")
    gene = cl[0]
    gRNA = cl[1][number[0]:number[1]]
    result = (gRNA, gene)
    return result

def search_b(fix_seq, fix_seq_len, number, record):
    seq = str(record.seq)
    if seq[0:30] == fix_seq: ##valid record
        gRNA = seq[fix_seq_len+number[0]:fix_seq_len+number[1]]
        return gRNA

def pcr_count(project_name):
    '''
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
