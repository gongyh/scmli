#!/usr/bin/env python3

import os
import sys
import re
import pandas as pd
import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from lxml import etree

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
    name1 = re.split('/',read1)[-1]
    name1 = re.split('\.fastq|\.fq',name1)[0] 
    name12 = name1 + '_val_1.fq.gz'
    name2 = re.split('/',read2)[-1]
    name2 = re.split('\.fastq|\.fq',name2)[0]
    name22 = name2 + '_val_2.fq.gz'
    os.system('gzip -cd ' + name12 + ' ' + name22 + ' > ' + project_name + '.fq')

    #stats raw_reads
    stats = {}
    html1 = etree.parse(name1+'_fastqc.html',etree.HTMLParser())
    html2 = etree.parse(name2+'_fastqc.html',etree.HTMLParser())
    raw_reads1 = html1.xpath('/html/body/div[3]/div[1]/table/tbody/tr[4]/td[2]')[0].text
    raw_reads2 = html2.xpath('/html/body/div[3]/div[1]/table/tbody/tr[4]/td[2]')[0].text

    stats['raw_reads'] = int((int(raw_reads1)+int(raw_reads2))/2)

    return stats

def pcr_parse_gRNA(stats, lib, fix_seq, number=[25,45], project_name='my_project', threads=8):
    
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
    for i in gRNA_gene.keys():
        gRNAs_dict[i]=0
    gRNAs_dict_original = gRNAs_dict.copy()
    read_counts = 0
    # Get fixed sequence from file.fastq and count
    # {'GAGTGTGGTGGAATTTGCCG': 3, ...}

    with open(project_name + '.fq') as handle:
        file_unknow = open('unknow.seq','w')
        for (title, seq, quality) in FastqGeneralIterator(handle):    
            read_counts += 1
            if seq[0:fix_seq_len] == fix_seq: #valid record
                #change left and right equal
                gRNA = seq[fix_seq_len+number[0]:fix_seq_len+number[1]]
                if gRNA in gRNAs_dict_original.keys(): #gRNAs
                    gRNAs_dict[gRNA] += 1
                elif gRNA in gRNAs_dict.keys():
                    gRNAs_dict[gRNA] += 1
                    file_unknow.write(seq+'\n')
                else:
                    gRNAs_dict[gRNA] = 1
                    file_unknow.write(seq+'\n')
        file_unknow.close()

    #stats all_reads
    stats['all_reads'] = int(read_counts/2) #paired

    '''
    #multi processing
    # from multiprocessing import Pool
    # from functools import partial
    pool = Pool(threads)
    with open(project_name + '.fq') as handle:
        records = list(seq for (title, seq, quality) in FastqGeneralIterator(handle))
    gRNA = pool.map(partial(search,fix_seq,fix_seq_len,number), records)
    for key in gRNA:
        gRNAs_dict[key] = gRNAs_dict.get(key, 0) + 1
    del gRNAs_dict[None]
    '''

    # Change count format and add gRNA_gene id
    # NO01G00240  GAGTGTGGTGGAATTTGCCG  3
    # unknow      CCCCCCCCCCGAGTGTGGTG  1

    with open(project_name+'.counts', 'w') as f:
        for k, v in gRNAs_dict.items():
            gene = 'unknow'
            if k in gRNA_gene.keys():
                gene = gRNA_gene[k]
            f.write("%s\t%s\t%d\n"%(gene,k,v))

    return stats

'''
def search(fix_seq, fix_seq_len, number, record):
    if record[0:30] == fix_seq: ##valid record
        gRNA = record[fix_seq_len+number[0]:fix_seq_len+number[1]]
        return gRNA
'''

def pcr_count(project_name, stats):

    # Statistics

    #fixed_seq, e.g. GGTAGAATTGGTCGTTGCCATCGACCAGGC
    print("stats......")
    df1 = pd.read_csv(project_name + '.counts', sep = '\t', header = None, names = ["gene_id", "sequence", "counts", "percent"])
    t = df1.counts.sum()
    df1.percent = (df1.counts/t)*100
    df1.loc[:,'percent'] = df1.loc[:,'percent'].round(7)    
    df1 = df1.sort_values(by='counts', ascending = False)
    df1.to_csv(project_name + '.percent', sep = '\t', index = False)

    stats['valid_reads'] = np.sum(df1.counts)
    stats['unknow_reads'] = np.sum(df1.loc[df1.gene_id=='unknow','counts'])
    df_gRNAs = df1[(df1.gene_id!='unknow') & (df1.counts!=0)].copy()
    stats['gRNAs_reads'] = np.sum(df_gRNAs.counts)
    stats['all_kinds'] = len(df1)
    stats['lib_kinds'] = 9709
    stats['unknow_kinds'] = len(df1[df1.gene_id=='unknow'])
    stats['gRNAs_kinds'] = len(df_gRNAs)
    stats['all/raw_reads_percent'] = round(stats['all_reads'] / stats['raw_reads'],6)
    stats['valid/all_reads_percent'] = round(stats['valid_reads'] / stats['all_reads'],6)
    stats['gRNAs/valid_reads_percent'] = round(stats['gRNAs_reads'] / stats['valid_reads'],6)
    stats['gRNAs_coverage'] = round(stats['gRNAs_kinds']/9709,6)
    stats['gRNAs_average_all'] = round(np.average(df1.loc[df1.gene_id!='unknow','counts']),6)
    stats['gRNAs_average_present'] = round(np.average(df_gRNAs.counts),6)
    with open(project_name+'.stats','w') as filestats:
        for a,b in stats.items():
            filestats.write('%s\t%s\n'%(a,b))

    return stats


