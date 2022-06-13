
import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from functools import partial


def search_a(number,line):
    cl = line.split(",")
    gene = cl[0]
    gRNA = cl[1][number[0]:number[1]]
    result = (gRNA, gene)
    return result

def func(x):
    x=int(x)*int(x)


pool = Pool(6)
number = [25,45]


gRNA_gene = {}
gRNAs_dict = {}

# with open("../test/NoIMET1_gRNAs.csv","r") as fh:
#    result = pool.map(partial(search_a,number),[line for line in fh])

aa=datetime.datetime.now()
pool.map(func,[range(10000000)])

bb=datetime.datetime.now()
cc=bb-aa
print('运行时间是：',cc)
