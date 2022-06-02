import pytest
import sys
import os
import argparse
import pandas as pd
sys.path.append("..")
from scmli_lib.pcr_pipeline import pcr_qc, pcr_parse_gRNA, pcr_count

def test_pcr():
    '''
    Test PCR pipeline
    '''
    # parameters
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    args.model='PCR'
    args.lib='NoIMET1_gRNAs.csv'
    args.seq='GGTAGAATTGGTCGTTGCCATCGACCAGGC'
    args.read1='test_R1.fq.gz'
    args.read2='test_R2.fq.gz'
    args.number=[25, 45]
    args.output_name='my_project'
    args.output_dir='output'   

    args.read1 = os.path.abspath(args.read1)
    args.read2 = os.path.abspath(args.read2)
    args.lib = os.path.abspath(args.lib)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    args.output_dir = os.path.abspath(args.output_dir)

    # validate each step

    pcr_qc(args.output_dir, args.output_name, args.read1, args.read2)
    num_gRNAs = pcr_parse_gRNA(args.lib, args.seq, args.number, args.output_name)
    pcr_count(args.output_name)
    assert os.path.isfile(args.output_name + ".percent") == True
    assert os.path.getsize(args.output_name + ".percent") > 100
