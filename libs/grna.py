#!/usr/bin/env python3

import os
import sys
import re
import subprocess
import pandas as pd
import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from lxml import etree

# grna pipeline
def grna_pipeline(project_name, read1, read2, FASTQC_PATH, TRIM_GALORE_PATH, threads, lib, fix_seq, number, scmli_dir):

    # quality control, use fastqc
    print("quality control......")
    # check path
    if FASTQC_PATH:
        FASTQC_BIN = os.path.join(FASTQC_PATH, 'fastqc')
    else:
        FASTQC_BIN = 'fastqc'

    # run and save log
    cmd = [FASTQC_BIN, '-o', '.', read1, read2]
    with open('grna_pipeline.log', 'a') as log_file:
        result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        print("FastQC failed with return code:", result.returncode)

    # trim, use trim_galore
    print("trim......")
    # check path
    if TRIM_GALORE_PATH:
        TRIM_GALORE_BIN = os.path.join(TRIM_GALORE_PATH, 'trim_galore')
    else:
        TRIM_GALORE_BIN = 'trim_galore'

    # run and save log
    cmd = [TRIM_GALORE_BIN, '--paired', '--fastqc', '--max_n',
           '0', '-j', str(threads), '--gzip', read1, read2]
    with open('grna_pipeline.log', 'a') as log_file:
        result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        print("Trim galore failed with return code:", result.returncode)

    # Unzip qc/trim output files and merge
    name1_base = re.split(r'/', read1)[-1]
    name1 = re.split(r'\.fastq|\.fq', name1_base)[0]
    name12 = name1 + '_val_1.fq.gz'

    name2_base = re.split(r'/', read2)[-1]
    name2 = re.split(r'\.fastq|\.fq', name2_base)[0]
    name22 = name2 + '_val_2.fq.gz'

    cmd = ['gzip', '-cd', name12, name22]
    with open(f'{project_name}.fq', 'w') as output_file:
        result = subprocess.run(cmd, stdout=output_file)

    # stats, get raw_reads number from fastqc result
    stats = {} 
    html1 = etree.parse(name1+'_fastqc.html', etree.HTMLParser())
    html2 = etree.parse(name2+'_fastqc.html', etree.HTMLParser())
    raw_reads1 = html1.xpath(
        '/html/body/div[3]/div[1]/table/tbody/tr[4]/td[2]')[0].text
    raw_reads2 = html2.xpath(
        '/html/body/div[3]/div[1]/table/tbody/tr[4]/td[2]')[0].text

    stats['raw_reads'] = int((int(raw_reads1)+int(raw_reads2))/2)

    # fix_seq = args.seq("GGTAGAATTGGTCGTTGCCATCGACCAGGC")
    # search gRNAs
    print("search......")

    # Get the gRNA_gene dictionary {seq: id, ...} from lib file
    # {'GAGTGTGGTGGAATTTGCCG': 'NO01G00240', ...}
    gRNA_gene = {}
    lib_kinds = 0
    with open(lib) as fh:
        for line in fh:
            cl = line.split(',')
            gene = cl[0]
            gRNA = cl[1][number[0]:number[1]]
            gRNA_gene[gRNA] = gene
            lib_kinds += 1

    fix_seq_len = len(fix_seq)
    gRNAs_dict = {}
    for i in gRNA_gene.keys():
        gRNAs_dict[i] = 0
    gRNAs_dict_original = gRNAs_dict.copy()
    read_counts = 0
    # Get fixed sequence from file.fastq and count
    # {'GAGTGTGGTGGAATTTGCCG': 3, ...}

    with open(project_name + '.fq') as handle:
        file_unknow = open('unknow.seq', 'w')
        for (title, seq, quality) in FastqGeneralIterator(handle):
            read_counts += 1
            if seq[0:fix_seq_len] == fix_seq:  # valid record
                # change left and right equal
                gRNA = seq[fix_seq_len+number[0]:fix_seq_len+number[1]]
                if gRNA in gRNAs_dict_original.keys():  # gRNAs
                    gRNAs_dict[gRNA] += 1
                elif gRNA in gRNAs_dict.keys():
                    gRNAs_dict[gRNA] += 1
                    file_unknow.write(seq+'\n')
                else:
                    gRNAs_dict[gRNA] = 1
                    file_unknow.write(seq+'\n')
        file_unknow.close()

    # stats all_reads
    stats['all_reads'] = int(read_counts/2)  # paired

    # Change count format and add gRNA_gene id
    # NO01G00240  GAGTGTGGTGGAATTTGCCG  3
    # unknow      CCCCCCCCCCGAGTGTGGTG  1

    with open(project_name+'.counts', 'w') as f:
        for k, v in gRNAs_dict.items():
            gene = 'unknow'
            if k in gRNA_gene.keys():
                gene = gRNA_gene[k]
            f.write("%s\t%s\t%d\n" % (gene, k, v))

    # Statistics
    # fixed_seq, e.g. GGTAGAATTGGTCGTTGCCATCGACCAGGC
    print("stats......")
    df1 = pd.read_csv(project_name + '.counts', sep='\t', header=None,
                      names=["gene_id", "sequence", "counts", "percentage", "percentage_gRNAs", "accumulative_unknow_percentage"])
    # percentage
    t = df1.counts.sum()
    t_gRNAs = df1[df1["gene_id"] != "unknow"].counts.sum()
    df1.percentage = (df1.counts/t)*100
    df1.percentage_gRNAs = (df1[df1["gene_id"] != "unknow"].counts/t_gRNAs)*100
    df1.loc[:, 'percentage'] = df1.loc[:, 'percentage'].round(6)
    df1.loc[:, 'percentage_gRNAs'] = df1.loc[:, 'percentage_gRNAs'].round(6)
    # sort
    df1 = df1.sort_values(by='counts', ascending=False)
    df1 = df1.reset_index(drop=True)
    # accumulative_unknow_percentage
    accumulative_unknow = 0
    accumulative_all = 0
    for i in range(0, len(df1)):
        accumulative_all += df1.counts[i]
        if df1.loc[i, "gene_id"] == "unknow":
            accumulative_unknow += df1.counts[i]
        df1.loc[i, "accumulative_unknow_percentage"] = round(
            (accumulative_unknow/accumulative_all)*100, 3)

    df1.to_csv(project_name + '.percentage', sep='\t', index=False)

    # Mutation target region bed file
    df2 = df1
    df2 = df2[df2['gene_id'] != 'unknow']
    df2 = df2[df2['counts'] != 0]
    df_target = pd.read_csv(scmli_dir+'/test/targets.bed', sep='\t', header=None,
                            names=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'])
    df2 = pd.merge(df2, df_target, left_on='gene_id', right_on='name')
    df2 = df2[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'sequence',
               'counts', 'percentage', 'percentage_gRNAs', 'accumulative_unknow_percentage']]
    df2.to_csv('targets_grna.bed', sep='\t', header=None, index=False)

    print('stats2......')
    stats['valid_reads'] = np.sum(df1.counts)
    stats['unknow_reads'] = np.sum(df1.loc[df1.gene_id == 'unknow', 'counts'])
    df_gRNAs = df1[(df1.gene_id != 'unknow') & (df1.counts != 0)].copy()
    stats['gRNAs_reads'] = np.sum(df_gRNAs.counts)
    stats['all_kinds'] = len(df1)
    stats['lib_kinds'] = lib_kinds
    stats['unknow_kinds'] = len(df1[df1.gene_id == 'unknow'])
    stats['gRNAs_kinds'] = len(df_gRNAs)
    stats['all/raw_reads_percentage'] = round(
        stats['all_reads'] / stats['raw_reads'], 6)
    stats['valid/all_reads_percentage'] = round(
        stats['valid_reads'] / stats['all_reads'], 6)
    stats['gRNAs/valid_reads_percentage'] = round(
        stats['gRNAs_reads'] / stats['valid_reads'], 6)
    stats['gRNAs_coverage'] = round(
        stats['gRNAs_kinds'] / stats['lib_kinds'], 6)
    stats['gRNAs_average_all'] = round(np.average(
        df1.loc[df1.gene_id != 'unknow', 'counts']), 6)
    stats['gRNAs_average_appeared'] = round(np.average(df_gRNAs.counts), 6)
    with open(project_name+'.stats', 'w') as stats_file:
        for a, b in stats.items():
            stats_file.write('%s\t%s\n' % (a, b))

    return stats
