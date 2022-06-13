from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from functools import partial
import datetime

def search_a(number,line):
    cl = line.split(",")
    gene = cl[0]
    gRNA = cl[1][number[0]:number[1]]
    result = (gRNA, gene)
    return result

pool = Pool(6)
number = [25,45]

gRNA_gene = {}
gRNAs_dict = {}

with open("../test/NoIMET1_gRNAs.csv","r") as fh:
    result = pool.map(partial(search_a,number),[line for line in fh])
    gRNA_gene.update(result)

def search_b(fix_seq, fix_seq_len, number, record):
    seq = str(record.seq)
    if seq[0:30] == fix_seq: ##valid record
        gRNA = seq[fix_seq_len+number[0]:fix_seq_len+number[1]]
        return gRNA

number = [25,45]
fix_seq = "GGTAGAATTGGTCGTTGCCATCGACCAGGC"
fix_seq_len = len(fix_seq)

aa=datetime.datetime.now()

with open("Zhi-li-ku.fq") as handle:
    gRNA = Pool(6).map(partial(search_b,fix_seq,fix_seq_len,number),[record for record in SeqIO.parse(handle, "fastq")])
    for key in gRNA:
        gRNAs_dict[key] = gRNAs_dict.get(key, 0) + 1
    # print(gRNAs_dict)

bb=datetime.datetime.now()
cc=bb-aa
print('multi运行时间是：',cc)


aa=datetime.datetime.now()
with open("Zhi-li-ku.fq") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        seq = str(record.seq)
        if seq[0:fix_seq_len] == fix_seq:  # valid record
            # change left and right equal
            gRNA = str(seq[fix_seq_len+number[0]:fix_seq_len+number[1]])
            if gRNA in gRNAs_dict.keys():
                gRNAs_dict[gRNA] += 1
            else:
                gRNAs_dict[gRNA] = 1

bb=datetime.datetime.now()
cc=bb-aa
print('single运行时间是：',cc)
