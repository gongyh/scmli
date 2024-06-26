#!/usr/bin/env python3
import pandas as pd
import numpy as np
import re
import os

'''
Script for identifying gRNA results in experimental protocols.

We used twenty-four 24-well plates for the experiment.
Obtained sequencing data from the same plates and the same wells.
gRNAs results that appear simultaneously in a given plate and well 
will be considered as the result for that specific plate and well.
Then perform whole-genome sequencing for each plate to verify mutations.
'''

# Create a list of 24 numbers since the 24-plate or 24-well
list = ["{:02d}".format(i) for i in range(1, 25)]

def get_score(df1, df2, name1, name2, method):
    # Calculate the logarithm of the counts
    df1['lg_count1'] = np.log2(df1.counts1)
    df1.loc[np.isinf(df1.lg_count1), 'lg_count1'] = 0
    # Convert logarithms to percentages
    df1['lg_count_pct1'] = (df1['lg_count1']/df1['lg_count1'].sum())*100
    df2['lg_count2'] = np.log2(df2.counts2)
    df2.loc[np.isinf(df2.lg_count2), 'lg_count2'] = 0
    df2['lg_count_pct2'] = (df2['lg_count2']/df2['lg_count2'].sum())*100
    # Create a dataFrame to save results
    df3 = pd.DataFrame(
        columns=['seq', 'id', 'plate', 'well', 'a', 'b', 'ab', 'lgpa', 'lgpb', 'lgpab'])
    
    for i in range(len(df1)):
        iseq = df1.sequence1[i]
        # If co-occurrence in the plate and well, calculate score
        if iseq in df2.sequence2.values:
            # gRNA gene id
            iid = df1.gene_id1[i]
            # percentages of gRNA counts
            a = df1.iloc[i, ].percentage1
            b = df2[df2.sequence2 == iseq].percentage2.values[0]
            # Percentages of the logarithms of the counts
            # lgpa lgpb lgpa*lgpb
            lgpa = df1.iloc[i, ].lg_count_pct1
            lgpb = df2[df2.sequence2 == iseq].lg_count_pct2.values[0]
            df3.loc[len(df3)] = [iseq, iid, name1, name2,
                                 a, b, a*b, lgpa, lgpb, lgpa*lgpb]

    # Sort the dataframe by the 'lgpab' score
    df3 = df3.sort_values(by='lgpab', ascending=False)
    df3.a = df3.a.round(3).astype(str)
    df3.b = df3.b.round(3).astype(str)
    df3.ab = df3.ab.round(3).astype(str)
    df3.lgpa = df3.lgpa.round(3).astype(str)
    df3.lgpb = df3.lgpb.round(3).astype(str)
    df3.lgpab = df3.lgpab.round(3).astype(str)
    df3.reset_index(drop=True, inplace=True)

    df3['result'] = df3.apply(lambda row: f"{row['id']},{row['seq']},"
                              f"{row['plate']}_{row['well']}{method}rank,"
                              f"{row['a']},{row['b']},{row['ab']},{row['lgpa']},"
                              f"{row['lgpb']},{row['lgpab']}", axis=1)

    score = ''
    for i in range(3):
        rank = "{:02d}".format(i+1)
        try:
            score += df3.result[i]+'\n'
            score = score.replace('rank', rank)
        except:
            score += 'Na\n'
    with open('location_result.txt', 'a') as f:
        f.write(score)


for n1 in list:
    # Method 01
    # The top three highest values of lgpa*lgpb
    for n2 in list:
        name1 = n1 + 'a'
        name2 = n2 + 'b'
        method = '_01_'
        df1 = pd.read_csv(name1+'.percentage', sep='\t', nrows=100, header=0, names=[
            'gene_id1', 'sequence1', 'counts1', 'percentage1', 'percentage_gRNAs1', 'accumulative_unknow_percentage1'])
        df2 = pd.read_csv(name2+'.percentage', sep='\t', nrows=100, header=0, names=[
            'gene_id2', 'sequence2', 'counts2', 'percentage2', 'percentage_gRNAs2', 'accumulative_unknow_percentage2'])
        get_score(df1, df2, name1, name2, method)

    # Method 02
    # The top three highest values of lgpa*lgpb after removing the maximum values from df1 and df2
    for n2 in list:
        name1 = n1 + 'a'
        name2 = n2 + 'b'
        method = '_02_'
        df1 = pd.read_csv(name1+'.percentage', sep='\t', skiprows=1, nrows=99, header=0, names=[
            'gene_id1', 'sequence1', 'counts1', 'percentage1', 'percentage_gRNAs1', 'accumulative_unknow_percentage1'])
        df2 = pd.read_csv(name2+'.percentage', sep='\t', skiprows=1, nrows=99, header=0, names=[
            'gene_id2', 'sequence2', 'counts2', 'percentage2', 'percentage_gRNAs2', 'accumulative_unknow_percentage2'])
        get_score(df1, df2, name1, name2, method)

    with open('location_result.txt', 'a') as f:
        f.write('\n')

# Select gene targets with gRNA
df_gRNA = pd.read_csv('location_result.txt', header=None,
                  usecols=[0, 2], names=['id', 'no'])
df_target = pd.read_csv('~/app/scmli/test/targets.bed', sep='\t', header=None,
                  names=['chrom', 'chromStart', 'chromEnd', 'id', 'score', 'strand'], dtype=object)
# Retain target information for genes present in gRNA
df_result = pd.merge(df_gRNA, df_target, how='left', on='id')
# Change the format
df_result = df_result[['chrom', 'chromStart', 'chromEnd', 'id', 'score', 'strand', 'no']]
df_result.no = df_result.no.str.replace('a|b', '', regex=True)
df_result = df_result.fillna('Na')

for i in list:
    # each plate
    j = int(i)
    # Select 72 rows from method_01 and 72 rows from method_02, totaling 144 rows 
    df_tmp = df_result[(j-1)*144:j*144]
    # Retain all 72 rows from method_01
    df_c1 = df_tmp[0:72]
    # Remove the rows in method_02 that are present in method_01.
    df_c2 = df_tmp.drop_duplicates('id')
    df_c2 = df_c2[df_c2.no.str.contains('02_..$', regex=True)]
    # Combine the results
    df_tmp = pd.concat([df_c1, df_c2])
    # Save results for each plate
    df_tmp.to_csv('%s.bed' % i, sep='\t', index=False, header=None)
