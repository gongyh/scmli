#!/bin/bash

## 1. read qc, base quality need better than Q30
fastqc -o . ../clean_data/*.fq.gz

## 2. trim low quality end
trim_galore --paired --fastqc --max_n 0 -j 4 --gzip ../clean_data/Zhi-li-ku_1.clean.fq.gz ../clean_data/Zhi-li-ku_2.clean.fq.gz

## 3. if not strand-specific, merge R1 and R2
gzip -cd Zhi-li-ku_1.clean_val_1.fq.gz Zhi-li-ku_2.clean_val_2.fq.gz > Zhi-li-ku.fq

## 4. parse gRNA
python3 parse_gRNA.py Zhi-li-ku.fq NoIMET1_gRNAs.csv CGCAGGCATAGTTTAGTGGTAGAATTGGTCGTTGCCATCGACCAGGC > Zhi-li-ku.stats

## 5. more stats
t=`awk -F'\t' '{total=total+$3};END{print total}' Zhi-li-ku.stats`
awk -F'\t' '{print $0"\t"$3*100.0/'$t'}' Zhi-li-ku.stats | sort -k3,3nr > Zhi-li-ku_pct.stats

