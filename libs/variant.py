#!/usr/bin/env python3

import os
import re
import subprocess
import pandas as pd

# Call variants using the modified Snippy software


def variant_pipeline(args):
    print('trim......')
    # trim_galore recommends using no more than 8 cores
    trim_threads = min(args.threads, 8)
    # Configure parameters
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
    # Run trim_galore and save the log
    with open('variant.log', 'w') as log_file:
        subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)

    # Get the name of the trimmed output file
    name1 = re.split(r'/', args.read1)[-1]
    name1 = re.split(r'\.fastq|\.fq', name1)[0] + '_val_1.fq.gz'
    name2 = re.split(r'/', args.read2)[-1]
    name2 = re.split(r'\.fastq|\.fq', name2)[0] + '_val_2.fq.gz'

    print('call variants......')
    # Snippy requires a pre-created "tmp" folder
    os.makedirs('tmp', exist_ok=True)
    # Configure parameters
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
    # Run snippy and save the log
    with open('variant.log', 'a') as log_file:
        subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)

    # Handle deletion type separately, it may overlap but is not within the target region
    print('search deletion......')
    # Read the snps file, ignoring lines that start with '##'
    cmd = f'grep -c "##" {args.outname}_snippy/snps.vcf'
    header_count = int(subprocess.check_output(cmd, shell=True).strip())
    df_snps = pd.read_csv(args.outname+'_snippy/snps.vcf',
                          sep='\t', skiprows=int(header_count))
    # Determining the position of gene deletion to check if it's in the target
    x = len(df_snps)
    # BED file of deletions
    with open('del.bed', 'w') as f:
        # Write an initial line
        f.write('chr1\t2\t3\n')
        # Iterate through the DataFrame
        for i in range(x):
            # Check if 'del' is in the 8th column 'INFO'
            if 'del' in df_snps.iloc[i, 7]:
                # Get the length of the 4th column 'REF'
                y = len(df_snps.iloc[i, 3])
                # line = 'CHROM' + 'POS'1 + 'POS'2
                line = df_snps.iloc[i, 0]+'\t' + \
                    str(df_snps.iloc[i, 1]-1)+'\t' + \
                    str(df_snps.iloc[i, 1]+y-1)+'\n'
                f.write(line)
    subprocess.run(f'bedtools intersect -a del.bed -b {args.target} -wa -wb > del_target.bed', shell=True)
    with open('del_target.bed', 'a') as f:
        f.write('chr1\t3\t4\n')

    print('filter...')
    # target1: all gRNA target(9709)
    # Extract variants from the original target1 region
    subprocess.run(f'bcftools view -e "QUAL>1" -T {args.target} -O z -o {args.outname}_snippy_hq.vcf.gz {args.outname}_snippy/snps.vcf.gz', shell=True)
    # Extract deletion variants that overlap target1 region
    subprocess.run(f'bcftools view -v indels -e "QUAL>1" -T del_target.bed -O z -o {args.outname}_snippy_del_target.vcf.gz {args.outname}_snippy/snps.vcf.gz', shell=True)
    # Create index files
    subprocess.run(f'bcftools index {args.outname}_snippy_hq.vcf.gz', shell=True)
    subprocess.run(f'bcftools index {args.outname}_snippy_del_target.vcf.gz', shell=True)
    # Merge the two parts of variants
    subprocess.run(f'bcftools concat -a {args.outname}_snippy_hq.vcf.gz {args.outname}_snippy_del_target.vcf.gz -o {args.outname}_snippy_target.vcf >> variant.log 2>&1', shell=True)
    # Extract gene ids of variants from the original target region
    subprocess.run(f"zcat {args.outname}_snippy_hq.vcf.gz | grep -v '^#' | cut -f8 | awk -F'|' '{{print $1,$2,$3,$4,$5}}' | grep -o 'NO..G.....' | sort | uniq > {args.outname}_snippy_hq.gids", shell=True)
    # Extract gene ids of variants from the original and overlap target region
    subprocess.run(f"cat {args.outname}_snippy_target.vcf | grep -v '^#' | cut -f8 | awk -F'|' '{{print $1,$2,$3,$4,$5}}' | grep -o 'NO..G.....' | sort | uniq > {args.outname}_snippy_target.gids", shell=True)

    # target2: screened gRNA target(obtained in gRNA pipeline)
    # Select the deletion and target2 overlap region
    subprocess.run(f'bedtools intersect -a del.bed -b {args.dtarget} -wa -wb > del_target2.bed', shell=True)
    # Write an initial line
    with open('del_target2.bed', 'a') as f:
        f.write('chr2\t1\t2\n')
    # Extract variants from the original target2 region
    subprocess.run(f'bcftools view -e "QUAL>1" -T {args.dtarget} -O z -o {args.outname}_snippy_raw_target2.vcf.gz {args.outname}_snippy/snps.vcf.gz', shell=True)
    # Extract variants from the overlap target2 region
    subprocess.run(f'bcftools view -v indels -e "QUAL>1" -T del_target2.bed -O z -o {args.outname}_snippy_del_target2.vcf.gz {args.outname}_snippy/snps.vcf.gz', shell=True)
    # Create index files
    subprocess.run(f'bcftools index {args.outname}_snippy_raw_target2.vcf.gz', shell=True)
    subprocess.run(f'bcftools index {args.outname}_snippy_del_target2.vcf.gz', shell=True)
    # Merge the two parts of variants
    subprocess.run(f'bcftools concat -a {args.outname}_snippy_raw_target2.vcf.gz {args.outname}_snippy_del_target2.vcf.gz -o {args.outname}_snippy_target2.vcf >> variant.log 2>&1', shell=True)
    # Extract gene ids of variants from the original and overlap target2 region
    subprocess.run(f"cat {args.outname}_snippy_target2.vcf | grep -v '^#' | cut -f8 | awk -F'|' '{{print $1,$2,$3,$4,$5}}' | grep -o 'NO..G.....' | sort | uniq > {args.outname}_snippy_target2.gids", shell=True)

    # Stats and add variant info to target2 dataframe
    df_target = pd.read_csv(args.dtarget, sep='\t', header=None, names=[
                            'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'no'], dtype=str)
    # Read the snps file, ignoring lines that start with '##'
    cmd = f'grep -c "##" {args.outname}_snippy_target2.vcf'
    header_count = int(subprocess.check_output(cmd, shell=True).strip())
    df_vars = pd.read_csv(args.outname+'_snippy_target2.vcf', sep='\t', skiprows=int(header_count),
                          usecols=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
    # Extract variation information
    df_vars['gene_id'] = df_vars.INFO.str.extract('(NO..G.....)')
    df_vars['AODP'] = df_vars.INFO.str.extract('(AO.*DP.*?;)')
    df_vars['type'] = df_vars.INFO.str.extract('(\|.*?\|.*?\|)')
    # Create a dataframe to store the aggregated information
    info_cols = ['#CHROM', 'POS', 'REF', 'ALT', 'AODP', 'type']
    df_vars['info'] = df_vars[info_cols].astype(str).agg(';'.join, axis=1)
    # Merge the variation information into the target dataframe
    variants_grouped = df_vars.groupby('gene_id')['info'].agg(
        [('variant', 'count'), ('info', ' '.join)]).reset_index()
    df_target = df_target.merge(
        variants_grouped, how='left', left_on='name', right_on='gene_id')
    df_target['variant'] = df_target['variant'].fillna(0).astype(int)
    # Save the result
    df_target.to_csv('target2_variant.txt', sep='\t', index=False)
