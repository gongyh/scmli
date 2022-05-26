#!/usr/bin/env python3

import os
from parse_gRNA import parse_gRNA
import pandas as pd

def pcr_qc(output_dir, output_name, read1, read2):
    #fastqc
    os.chdir(output_dir)
    os.system('fastqc -o ' + output_dir + ' ' + read1 + ' ' + read2)
    
    #trim
    print("trim_galore")
    os.system('trim_galore --paired --fastqc --max_n 0 -j 4 --gzip ' + read1 + ' ' + read2)

    #jieya hebing
    name1 = read1.split('/')[-1]
    name1 = name1.split(".fq.gz")[0] + "_val_1.fq.gz"
    name2 = read2.split('/')[-1]
    name2 = name2.split(".fq.gz")[0] + "_val_2.fq.gz"
    os.system('gzip -cd ' + name1 + ' ' + name2 + ' > ' + output_name + '.fq')
    

def pcr_count(fqfile, lib, fix_seq) :
    parse_gRNA(fqfile, lib, fix_seq)
    #os.system('python3 parse_gRNA.py my_project.fq NoIMET1_gRNAs.csv GGTAGAATTGGTCGTTGCCATCGACCAGGC > counts.stats')

    df1 = pd.read_csv("counts.stats", sep = "\t", header = None, names = ["gene_id", "sequence", "counts", "percent"])
    t = df1["counts"].sum()
    df1["percent"] = df1["counts"]/t
    df1["percent"] = df1["percent"].round(7)    
    df1 = df1.sort_values(by="counts", ascending = False)
    df1.to_csv("percent.stats", sep = "\t", index = False)


