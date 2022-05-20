#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
import sys

fastq_file = sys.argv[1] # fastq
design_file = sys.argv[2] # NoIMET1_gRNAs.csv
const_seq = sys.argv[3] #GGTAGAATTGGTCGTTGCCATCGACCAGGC or CGCAGGCATAGTTTAGTGGTAGAATTGGTCGTTGCCATCGACCAGGC

header = True
gRNA_gene = {}
with open(design_file) as fh:
    for line in fh:
        if header:
            header = False
            continue
        cl = line.split(",")
        gene = cl[0]
        gRNA = cl[1][25:45]
        gRNA_gene[gRNA] = gene

const_seqlen = len(const_seq)

gRNAs_dict = {}

with open(fastq_file) as handle:
    for record in SeqIO.parse(handle, "fastq"):
        seq = str(record.seq)
        if seq[0:const_seqlen]==const_seq: # valid record
            #confirm left and right seg equal
            gRNA = str(Seq(seq[const_seqlen+25:const_seqlen+45]))
            if gRNA in gRNAs_dict.keys():
                gRNAs_dict[gRNA] += 1
            else:
                gRNAs_dict[gRNA] = 1

for k,v in gRNAs_dict.items():
    gene = "unknown"
    if k in gRNA_gene.keys():
        gene = gRNA_gene[k]
    print("%s\t%s\t%d"%(gene,k,v))
