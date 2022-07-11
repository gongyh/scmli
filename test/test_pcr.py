import pytest
import sys
import os
import argparse
import pandas as pd
sys.path.append("..")
from libs.pcr_pipeline import pcr_qc, pcr_parse_gRNA, pcr_count
from scmli import create_arg_parser

def test_pcr():
    '''
    Test PCR pipeline
    '''
    # parameters
    parser = create_arg_parser()
    
    args = parser.parse_args([
        '-m','PCR',
        '-l','NoIMET1_gRNAs.csv',
        '-s','GGTAGAATTGGTCGTTGCCATCGACCAGGC',
        '-r1','test_R1.fq.gz',
        '-r2','test_R2.fq.gz'])

    args.read1 = os.path.abspath(args.read1)
    args.read2 = os.path.abspath(args.read2)
    args.lib = os.path.abspath(args.lib)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    args.output_dir = os.path.abspath(args.output_dir)

    # validate each step
    current_dir = os.getcwd()
    scmli_dir = os.path.split(os.path.realpath(__file__))[0]
    os.chdir(args.output_dir)
    pcr_qc(args.output_name, args.read1, args.read2, args.FASTQC_PATH, args.TRIM_GALORE_PATH, args.threads)
    stats = pcr_parse_gRNA(args.lib, args.seq, args.number, args.output_name, args.threads)
    stats = pcr_count(args.output_name, stats)
    os.system('Rscript '+scmli_dir+'/../libs/plot.r '+args.output_name)
    assert os.path.isfile(args.output_name + ".percent") == True
    assert os.path.isfile("frequency.pdf") == True
    assert stats["all_kinds"] == 12649
    os.chdir(current_dir)
