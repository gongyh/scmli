#!/usr/bin/env python3

import os
import sys
import re
import pandas as pd
import numpy as np


def variant_pipeline(args):
    os.system('trim_galore --cores 8 --paired -q 30 --trim-n --max_n 0 --length 70 --gzip %s %s'%(args.read1,args.read2))

    name1 = re.split('/', args.read1)[-1]
    name1 = re.split('\.fastq|\.fq', name1)[0]+'_val_1.fq.gz'
    name2 = re.split('/', args.read2)[-1]
    name2 = re.split('\.fastq|\.fq', name2)[0]+'_val_2.fq.gz'

    os.system('mkdir tmp')
    os.system('snippy --cpus %d --ram 40 --basequal 30 --minqual 0.0 --minfrac 0.0 --report --outdir %s_snippy --tmpdir tmp --ref %s --R1 %s --R2 %s'%(args.threads,args.outname,args.ref,name1,name2))

    os.system('bcftools view -e "(INFO/AO)/(INFO/DP)>0.5" -T %s %s_snippy/snps.vcf.gz > %s_snippy_hq.vcf'%(args.target,args.outname,args.outname))

    os.system("cat %s_snippy_hq.vcf | grep -v '^#' | cut -f8 | awk -F',' '{print $1}' | grep -o 'NO..G.....' | sort | uniq > %s_snippy_hq.gids"%(args.outname,args.outname))
