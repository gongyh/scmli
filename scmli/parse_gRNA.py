#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
import sys

def parse_gRNA(fqfile, lib, fix_seq):
    header = True
    gRNA_gene = {}
    with open(lib) as fh:
        for line in fh:
            if header:
                header = False
                continue
            cl = line.split(",")
            gene = cl[0]
            gRNA = cl[1][25:45] #change number
            gRNA_gene[gRNA] = gene

    fix_seq_len = len(fix_seq)
    gRNAs_dict = {}
    
    with open(fqfile) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)
            if seq[0:fix_seq_len] == fix_seq: #valid record
                #change left and right equal
                gRNA = str(seq[fix_seq_len+25:fix_seq_len+45])
                if gRNA in gRNAs_dict.keys():
                    gRNAs_dict[gRNA] += 1
                else:
                    gRNAs_dict[gRNA] = 1
    for k, v in gRNAs_dict.items():
        gene = "unknow"
        f = open("counts.stats", "a")
        if k in gRNA_gene.keys():
            gene = gRNA_gene[k]
        print("%s\t%s\t%d"%(gene,k,v), file = f)
                
















