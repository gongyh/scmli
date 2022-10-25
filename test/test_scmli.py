import pytest
import sys
import os
import argparse
import pandas as pd
sys.path.append("..")
from scmli import create_arg_parser


def test_grna():
    # parameters
    parser = create_arg_parser()
    
    args = parser.parse_args([
        'gRNA',
        '-l', 'NoIMET1_gRNAs.csv',
        '-s', 'GGTAGAATTGGTCGTTGCCATCGACCAGGC',
        '-r1', 'test_R1.fq.gz',
        '-r2', 'test_R2.fq.gz'])
   
    stats = args.func(args)
    assert os.path.isfile("output/reads.pdf") == True
    assert stats["all_kinds"] == 12649


def test_variant():
    parser = create_arg_parser()
    args = parser.parse_args([
        'variant',
        '-r1', 'test_R1.fq.gz',
        '-r2', 'test_R2.fq.gz',
        '--ref', 'genes.gbk',
        '--target', 'targets.bed',
        '--dtarget', 'filter.bed',
        '-o', 'var',
        '-n', 't1'])
    args.func(args)
    assert os.path.getsize('var/t1_snippy/snps.vcf') > 10
    assert os.path.getsize('var/t1_snippy_hq.vcf') > 10

