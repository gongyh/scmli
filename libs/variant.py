#!/usr/bin/env python3

import os
import re
import subprocess
import pandas as pd


def variant_pipeline(args):
    print('trim......')
    trim_threads = min(args.threads, 8)
    cmd = [
        'trim_galore',
        '--cores', str(trim_threads),
        '--paired',
        '-q', '20',
        '--trim-n',
        '--max_n', '0',
        '--length', '70',
        '--gzip',
        args.read1, args.read2
    ]
    with open('variant.log', 'w') as log_file:
        subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)

    name1 = re.split(r'/', args.read1)[-1]
    name1 = re.split(r'\.fastq|\.fq', name1)[0] + '_val_1.fq.gz'
    name2 = re.split(r'/', args.read2)[-1]
    name2 = re.split(r'\.fastq|\.fq', name2)[0] + '_val_2.fq.gz'

    print('call variants......')
    os.makedirs('tmp', exist_ok=True)
    cmd = [
        'snippy',
        '--cpus', str(args.threads),
        '--ram', str(args.ram),
        '--basequal', '30',
        '--minqual', '0.0',
        '--minfrac', '0.0',
        '--report',
        '--outdir', f'{args.outname}_snippy',
        '--tmpdir', 'tmp',
        '--ref', args.ref,
        '--R1', name1,
        '--R2', name2
    ]

    with open('variant.log', 'a') as log_file:
        subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)

    os.system('snippy --cpus %d --ram %d --basequal 30 --minqual 0.0 --minfrac 0.0 --report --outdir %s_snippy --tmpdir tmp --ref %s --R1 %s --R2 %s >> variant.log 2>&1' %
              (args.threads, args.ram, args.outname, args.ref, name1, name2))
    print('search deletion......')
    row = re.findall(
        r'\d+', os.popen('grep "##" %s_snippy/snps.vcf | wc -l' % args.outname).read())[0]
    df1 = pd.read_csv(args.outname+'_snippy/snps.vcf',
                      sep='\t', skiprows=int(row))
    df2 = df1

    x = len(df1)
    with open('del.bed', 'w') as f:
        f.write('chr1\t1\t2\n')
        x = len(df1)
        for a in range(x):
            if 'del' in df1.iloc[a, 7]:
                y = len(df1.iloc[a, 3])
                line = df1.iloc[a, 0]+'\t' + \
                    str(df1.iloc[a, 1]-1)+'\t'+str(df1.iloc[a, 1]+y-1)+'\n'
                f.write(line)
    os.system(
        'bedtools intersect -a del.bed -b %s -wa -wb > del_target.bed' % args.target)
    with open('del_target.bed', 'a') as f:
        f.write('chr2\t1\t2\n')
    # target1 (9709), concat in target and deletion overlap
    print('filter...')
    os.system('bcftools view -e "QUAL>1" -T %s -O z -o %s_snippy_hq.vcf.gz %s_snippy/snps.vcf.gz' %
              (args.target, args.outname, args.outname))
    os.system('bcftools view -v indels -e "QUAL>1" -T del_target.bed -O z -o %s_snippy_del_target.vcf.gz %s_snippy/snps.vcf.gz' %
              (args.outname, args.outname))
    os.system('bcftools index %s_snippy_hq.vcf.gz' % args.outname)
    os.system('bcftools index %s_snippy_del_target.vcf.gz' % args.outname)
    os.system('bcftools concat -a %s_snippy_hq.vcf.gz %s_snippy_del_target.vcf.gz -o %s_snippy_target.vcf >> variant.log 2>&1' %
              (args.outname, args.outname, args.outname))
    os.system(
        "zcat %s_snippy_hq.vcf.gz | grep -v '^#' | cut -f8 | awk -F'|' '{print $1,$2,$3,$4,$5}' | grep -o 'NO..G.....' | sort | uniq > %s_snippy_hq.gids" % (args.outname, args.outname))
    os.system(
        "cat %s_snippy_target.vcf | grep -v '^#' | cut -f8 | awk -F'|' '{print $1,$2,$3,$4,$5}' | grep -o 'NO..G.....' | sort | uniq > %s_snippy_target.gids" % (args.outname, args.outname))

    # target2 (24), concat
    os.system(
        'bedtools intersect -a del.bed -b %s -wa -wb > del_target2.bed' % args.dtarget)
    with open('del_target2.bed', 'a') as f:
        f.write('chr2\t1\t2\n')
    os.system('bcftools view -e "QUAL>1" -T %s -O z -o %s_snippy_raw_target2.vcf.gz %s_snippy/snps.vcf.gz' %
              (args.dtarget, args.outname, args.outname))
    os.system('bcftools view -v indels -e "QUAL>1" -T del_target2.bed -O z -o %s_snippy_del_target2.vcf.gz %s_snippy/snps.vcf.gz' %
              (args.outname, args.outname))
    os.system('bcftools index %s_snippy_raw_target2.vcf.gz' % args.outname)
    os.system('bcftools index %s_snippy_del_target2.vcf.gz' % args.outname)
    os.system('bcftools concat -a %s_snippy_raw_target2.vcf.gz %s_snippy_del_target2.vcf.gz -o %s_snippy_target2.vcf >> variant.log 2>&1' %
              (args.outname, args.outname, args.outname))
    os.system(
        "cat %s_snippy_target2.vcf | grep -v '^#' | cut -f8 | awk -F'|' '{print $1,$2,$3,$4,$5}' | grep -o 'NO..G.....' | sort | uniq > %s_snippy_target2.gids" % (args.outname, args.outname))

    # stats add info to txt
    df3 = pd.read_csv(args.dtarget, sep='\t', header=None, names=[
                      'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'sequence', 'counts', 'percentage', 'percentage_gRNAs', 'accumulative_unknow_percentage'])
    df3['variant'] = 0
    df3['info'] = ''
    row2 = re.findall(
        r'\d+', os.popen('grep "##" %s_snippy_target2.vcf | wc -l' % args.outname).read())[0]
    df4 = pd.read_csv(args.outname+'_snippy_target2.vcf', sep='\t', skiprows=int(
        row2), usecols=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
    df4['gene_id'] = df4.INFO.str.extract('(NO..G.....)')
    for i in range(len(df4)):
        a = df4.gene_id[i]
        df3.loc[df3['name'] == a, 'variant'] += 1
        df3.loc[df3['name'] == a, 'info'] += '%s,%s,%s,%s;' % (
            df4.iloc[i, 0], df4.iloc[i, 1], df4.iloc[i, 3], df4.iloc[i, 4])
    df3.to_csv('target2_variant.txt', sep='\t', index=False)
