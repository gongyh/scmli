#!/bin/bash

## 1. read qc, base quality need better than Q30
#fastqc -o . ../clean_data/*.fq.gz

## 2. trim low quality end
trim_galore --paired --fastqc --max_n 0 -j 4 --gzip ../clean_data/PCR-ku_1.clean.fq.gz ../clean_data/PCR-ku_2.clean.fq.gz

## 3. if not strand-specific, merge R1 and R2
gzip -cd PCR-ku_1.clean_val_1.fq.gz PCR-ku_2.clean_val_2.fq.gz > PCR-ku.fq

## 4. parse gRNA
python3 parse_gRNA.py PCR-ku.fq NoIMET1_gRNAs.csv GGTAGAATTGGTCGTTGCCATCGACCAGGC > PCR-ku.stats

## 5. more stats
t=`awk -F'\t' '{total=total+$3};END{print total}' PCR-ku.stats`
awk -F'\t' '{print $0"\t"$3*100.0/'$t'}' PCR-ku.stats | sort -k3,3nr > PCR-ku_pct.stats

