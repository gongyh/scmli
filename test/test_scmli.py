import pytest
import sys
import os
import argparse
import pandas as pd
sys.path.append("..")
from scmli import create_arg_parser
from libs.pcr_pipeline import pcr_pipeline
from libs.variant import variant_pipeline
from libs.common import existing_file, existing_dir

def test_pcr():
    '''
    Test PCR pipeline
    '''
    # parameters
    parser = create_arg_parser()
    
    args = parser.parse_args([
        'gRNA',
        '-l','NoIMET1_gRNAs.csv',
        '-s','GGTAGAATTGGTCGTTGCCATCGACCAGGC',
        '-r1','test_R1.fq.gz',
        '-r2','test_R2.fq.gz'])
   
    stats = args.func(args)
    assert os.path.isfile(args.output_name + ".percentage") == True
    assert os.path.isfile("reads.pdf") == True
    assert stats["all_kinds"] == 12649


'''
def test_variant():
    parser = create_arg_parser()
    args = parser.parse_args([
        'variant',
        '-r1','test_R1.fq.gz',
        '-r2','test_R2.fq.gz',
        '--ref','genes.gbk',
        '--target','targets.bed'])
    stats = args.func(args)
    assert os.path.isfile(args.outname + '_snippy_hq.vcf') == True

'''
