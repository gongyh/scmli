import pytest
import sys
import os
import argparse
import pandas as pd
sys.path.append("..")
from libs.pcr_pipeline import pcr_pipeline
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
   
    stats = args.func(args)
    assert os.path.isfile(args.output_name + ".percentage") == True
    assert os.path.isfile("reads.pdf") == True
    assert stats["all_kinds"] == 12649
