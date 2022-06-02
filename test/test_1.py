import pytest
import sys
sys.path.append("..")
from scmli.pcr_pipline import pcr_qc, pcr_parse_gRNA, pcr_count

def test_pcr():
    '''
    Test PCR pipeline
    '''
    # parameters
    output_dir = 'my_project'
    output_name = 'output'
    read1 = 'test/test_R1.fq.gz'
    read2 = 'test/test_R2.fq.gz'
    lib = 'test/NoIMET1_gRNAs.csv'
    seq = 'GGTAGAATTGGTCGTTGCCATCGACCAGGC'
    number = '25 45'
    # validate each step
    pcr_qc(output_dir, output_name, read1, read2)
    num_gRNAs = pcr_parse_gRNA(lib, seq, number, output_name)
    pcr_count(output_name)
    assert num_gRNAs == 10000
