#!/usr/bin/env python3

import os
import re
import pandas as pd


def variant_pipeline(args):
    print('trim...')
    trim_threads = args.threads if args.threads < 8 else 8
    os.system('trim_galore --cores %d --paired -q 30 --trim-n --max_n 0 --length 70 --gzip %s %s > variant.log 2>&1'%(trim_threads, args.read1, args.read2))

    name1 = re.split('/', args.read1)[-1]
    name1 = re.split('.fastq|.fq', name1)[0]+'_val_1.fq.gz'
    name2 = re.split('/', args.read2)[-1]
    name2 = re.split('.fastq|.fq', name2)[0]+'_val_2.fq.gz'

    print('snippy...')
    os.system('mkdir tmp')
    os.system('snippy --cpus %d --ram 40 --basequal 30 --minqual 0.0 --minfrac 0.0 --report --outdir %s_snippy --tmpdir tmp --ref %s --R1 %s --R2 %s >> variant.log 2>&1'%(args.threads,args.outname,args.ref,name1,name2))
    os.system('cd %s_snippy'%args.outname)
    os.popen('snippy -h').read()
    os.system('cd ..')
    print('search deletion')
    row = re.findall('\d+', os.popen('grep "##" %s_snippy/snps.vcf | wc -l'%args.outname).read())[0]
    df1 = pd.read_csv(args.outname+'_snippy/snps.vcf',sep='\t',skiprows=int(row))
    df2 = df1

    x = len(df1)
    with open('del.bed', 'w') as f:
        x = len(df1)
        for a in range(x):
            if 'del' in df1.iloc[a,7]:
                y = len(df1.iloc[a,3])
                line = df1.iloc[a,0]+'\t'+str(df1.iloc[a,1]-1)+'\t'+str(df1.iloc[a,1]+y-1)+'\n'
                f.write(line)
    os.system('bedtools intersect -a del.bed -b %s -wa -wb > del_target.bed'%args.target)
    os.system('cat %s > snp_del_target.bed'%args.target)
    os.system('cat del_target.bed >> snp_del_target.bed')

    print('filter...')
    os.system('bcftools view -e "(INFO/AO)/(INFO/DP)>0.5" -T %s %s_snippy/snps.vcf.gz > %s_snippy_hq.vcf'%(args.target,args.outname,args.outname))
    os.system('bcftools view -e "(INFO/AO)/(INFO/DP)>0.5" -T snp_del_target.bed %s_snippy/snps.vcf.gz > %s_snippy_target.vcf'%(args.outname,args.outname))
    os.system("cat %s_snippy_hq.vcf | grep -v '^#' | cut -f8 | awk -F',' '{print $1}' | grep -o 'NO..G.....' | sort | uniq > %s_snippy_hq.gids"%(args.outname,args.outname))
    os.system("cat %s_snippy_target.vcf | grep -v '^#' | cut -f8 | awk -F',' '{print $1}' | grep -o 'NO..G.....' | sort | uniq > %s_snippy_target.gids"%(args.outname,args.outname))

    os.system('bedtools intersect -a del.bed -b %s -wa -wb > del_target2.bed'%args.dtarget)
    os.system('cat %s > snp_del_target2.bed'%args.dtarget)
    os.system('cat del_target2.bed >> snp_del_target2.bed')
    os.system('bcftools view -e "(INFO/AO)/(INFO/DP)>0.5" -T snp_del_target2.bed %s_snippy/snps.vcf.gz > %s_snippy_target2.vcf'%(args.outname,args.outname))
    os.system("cat %s_snippy_target2.vcf | grep -v '^#' | cut -f8 | awk -F',' '{print $1}' | grep -o 'NO..G.....' | sort | uniq > %s_snippy_target2.gids"%(args.outname,args.outname))

    df3=pd.read_csv(args.dtarget,sep='\t',header=None,names=['chrom','chromStart','chromEnd','name','score','strand','sequence','counts','percentage','percentage_gRNAs','accumulative_unknow_percentage'])
    df3['variant']=0
    df3['info']=''
    row2 = re.findall('\d+', os.popen('grep "##" %s_snippy_target2.vcf | wc -l'%args.outname).read())[0]
    df4=pd.read_csv(args.outname+'_snippy_target2.vcf',sep='\t',skiprows=int(row2),usecols=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'])
    df4['gene_id']=df4.INFO.str.extract('(NO..G.....)')
    for i in range(len(df4)):
        a=df4.gene_id[i]
        df3.loc[df3['name']==a,'variant'] += 1
        df3.loc[df3['name']==a,'info'] += '%s,%s,%s,%s;'%(df4.iloc[i,0],df4.iloc[i,1],df4.iloc[i,3],df4.iloc[i,4])
    df3.to_csv('target2_variant.txt',sep='\t',index=False)

