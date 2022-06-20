#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import matplotlib.pyplot as plt

def pcr_qc(project_name, read1, read2, FASTQC_PATH=None, TRIM_GALORE_PATH=None, threads=8):

    # Fastqc
    print("fastqc......")
    if FASTQC_PATH:
        FASTQC_BIN=FASTQC_PATH+'/fastqc'
    else:
        FASTQC_BIN='fastqc'
    os.system(FASTQC_BIN + ' -o . ' + read1 + ' ' + read2 + '> pcr_pipeline.log 2>&1') 

    # Trim_galore
    print("trim_galore......")
    if TRIM_GALORE_PATH:
        TRIM_GALORE_BIN=TRIM_GALORE_PATH+'/trim_galore'
    else:
        TRIM_GALORE_BIN='trim_galore'
    os.system(TRIM_GALORE_BIN + ' --paired --fastqc --max_n 0 -j ' + str(threads) + ' --gzip ' + read1 + ' ' + read2 + '>> pcr_pipeline.log 2>&1')

    # Unzip and merge files
    name1 = read1.split('/')[-1]
    name1 = name1.split('.fq.gz')[0] + '_val_1.fq.gz'
    name2 = read2.split('/')[-1]
    name2 = name2.split('.fq.gz')[0] + '_val_2.fq.gz'
    os.system('gzip -cd ' + name1 + ' ' + name2 + ' > ' + project_name + '.fq')


def pcr_parse_gRNA(lib, fix_seq, number=[25,45], project_name='my_project', threads=8):
    
    '''
    parse and count target sequence
    lib = args.lib (file.csv)
    fix_seq = args.seq ("GGTAGAATTGGTCGTTGCCATCGACCAGGC")
    '''

    print("search......")

    # Get the gRNA_gene {seq: id, ...} from lib
    # {'GAGTGTGGTGGAATTTGCCG': 'NO01G00240', ...}

    gRNA_gene = {}  
    with open(lib) as fh:
        for line in fh:
            cl = line.split(',')
            gene = cl[0]
            gRNA = cl[1][number[0]:number[1]]
            gRNA_gene[gRNA] = gene

    fix_seq_len = len(fix_seq)
    gRNAs_dict = {}
    
    # Get fixed sequence from file.fastq and count
    # {'GAGTGTGGTGGAATTTGCCG': 3, ...}

    with open(project_name + '.fq') as handle:
        for (title, seq, quality) in FastqGeneralIterator(handle):    
            if seq[0:fix_seq_len] == fix_seq: #valid record
                #change left and right equal
                gRNA = seq[fix_seq_len+number[0]:fix_seq_len+number[1]]
                if gRNA in gRNAs_dict.keys():
                    gRNAs_dict[gRNA] += 1
                else:
                    gRNAs_dict[gRNA] = 1


    # Change count format and add gRNA_gene id
    # NO01G00240  GAGTGTGGTGGAATTTGCCG  3
    # unknow      CCCCCCCCCCGAGTGTGGTG  1

    with open(project_name+'.counts', 'w') as f:
        for k, v in gRNAs_dict.items():
            gene = 'unknow'
            if k in gRNA_gene.keys():
                gene = gRNA_gene[k]
            f.write("%s\t%s\t%d\n"%(gene,k,v))


def pcr_count(project_name):

    # Statistics

    #fixed_seq, e.g. GGTAGAATTGGTCGTTGCCATCGACCAGGC
    print("stats......")
    df1 = pd.read_csv(project_name + '.counts', sep = '\t', header = None, names = ["gene_id", "sequence", "counts", "percent"])
    t = df1['counts'].sum()
    df1['percent'] = df1['counts']/t
    df1['percent'] = df1['percent'].round(7)    
    df1 = df1.sort_values(by='counts', ascending = False)
    df1.to_csv(project_name + '.percent', sep = '\t', index = False)

    stats = {}
    stats['all_kinds'] = len(df1)
    stats['lib_kinds'] = 9709
    stats['un_kinds'] = len(df1[df1['gene_id']=='unknow'])
    stats['no_kinds'] = len(df1[df1['gene_id']!='unknow'])
    stats['all_counts'] = np.sum(df1['counts'])
    stats['un_counts'] = np.sum(df1[df1['gene_id']=='unknow']['counts'])
    no_counts = df1[df1['gene_id']!='unknow']['counts']
    stats['no_counts'] = np.sum(no_counts)
    stats['no_average'] = np.average(df1[df1['gene_id']!='unknow']['counts'])
    stats['no_coverage'] = stats['no_kinds']/9710
    with open(project_name+'.stats','w') as filestats:
        filestats.write(str(stats))

    #plot1 four kinds
    x = ['All_gRNA','Gene_gRNA','Unknow_gRNA','Nano']
    y = [stats['all_kinds'],stats['no_kinds'],stats['un_kinds'],stats['lib_kinds']]
    plt.style.use('seaborn')
    fig, ax = plt.subplots(figsize=(7, 4), facecolor='white', dpi=100)
    plt.bar(x,y,width=0.4)
    plt.ylabel('Number of gRNAs', fontsize=14)
    plt.xticks(fontsize=12)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for a,b in zip(x,y):
        plt.text(a, b, b, ha='center', va='bottom', fontsize=14)
    plt.show()
    fig.savefig('allkinds.png')

    #plot2 four counts
    x = ['All_counts','Gene_counts','Unknow_counts']
    y = [stats['all_counts'],stats['no_counts'],stats['un_counts']]
    plt.style.use('seaborn')
    fig, ax = plt.subplots(figsize=(7.5, 4), facecolor='white', dpi=100)
    plt.bar(x,y,width=0.4)
    plt.ylabel('Number of gRNAs', fontsize=14)
    plt.xticks(fontsize=12)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for a,b in zip(x,y):
        plt.text(a, b, b, ha='center', va='bottom', fontsize=14)
    plt.show()
    fig.savefig('allcounts.png')

    #plot3 frequency
    x = range(len(no_counts))
    y = no_counts
    plt.style.use('seaborn')
    fig, ax = plt.subplots(figsize=(7, 4), facecolor='white', dpi=100)
    plt.scatter(x,y,s=1)
    plt.ylabel('Frequency', fontsize=14)
    plt.xticks(fontsize=12)
    plt.show()
    fig.savefig('Frequency.png')

    return stats
