#!/usr/bin/env python3

import argparse
import os
from libs.pcr_pipeline import pcr_qc, pcr_parse_gRNA, pcr_count
import pandas as pd


def create_arg_parser():

    parser = argparse.ArgumentParser(
        description="single cell mutant library inspertion")
    parser.add_argument('-m', '--model', default="TEST",
                        choices=["PCR", "TEST"], help="Choose analysis model")
    parser.add_argument('-l', '--lib', required=True, help="Gene library")
    parser.add_argument('-s', '--seq', required=True,
                        help="The fixed sequence for search")
    parser.add_argument('-r1', '--read1', required=True, help="Read1")
    parser.add_argument('-r2', '--read2', required=True, help="Read2")
    parser.add_argument('-t', '--threads', type=int, default=8, help="Number of threads to use")
    parser.add_argument('-num', '--number', type=int, nargs=2, default=[25, 45], 
                        help="Start and end of the gene position, default='25 45', from the 26-th to the 45-th bases")
    parser.add_argument('-n', '--output_name', default="my_project",
                        help="Prefix of output files, default='my_project'")
    parser.add_argument('-o', '--output_dir', default="output",
                        help="Directory of output files, default='output'")
    parser.add_argument('-p1', '--path_fastqc', help="PATH to fastqc")
    parser.add_argument('-p2', '--path_trim_galore', help="PATH to trim_galore")

    return parser


def check_args(parser):
    #parse args
    args = parser.parse_args()
    #check and modify args
    args.read1 = os.path.abspath(args.read1)
    args.read2 = os.path.abspath(args.read2)
    args.lib = os.path.abspath(args.lib)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    args.output_dir = os.path.abspath(args.output_dir)
    #check software path 
    if args.path_fastqc:
        if os.path.isdir(args.path_fastqc) and (os.path.isfile(args.path_fastqc+"/fastqc") or os.path.isfile(args.path_fastqc+"fastqc")):
            args.path_fastqc = os.path.abspath(args.path_fastqc)
        else:
            print("path to fastqc error")
            os._exit(0)
    else:
        pass

    if args.path_trim_galore:
        if os.path.isdir(args.path_trim_galore) and (os.path.isfile(args.path_trim_galore+"/trim_galore") or os.path.isfile(args.path_trim_galore+"trim_galore")):
            args.path_trim_galore = os.path.abspath(args.path_trim_galore)
        else:
            print("path to trim_galore error")
            os._exit(0)
    else:
        pass

    return args


if __name__ == "__main__":
    try:
        parser = create_arg_parser()
        args = check_args(parser)
    except:
        print("parser error")
        os._exit(0)
    if args.model == "PCR":
        current_dir = os.getcwd()
        os.chdir(args.output_dir)
        pcr_qc(args.output_name, args.read1, args.read2, args.path_fastqc, args.path_trim_galore, args.threads)
        pcr_parse_gRNA(args.lib, args.seq, args.number, args.output_name, args.threads)
        pcr_count(args.output_name)
        os.chdir(current_dir)
        print("Finished!")
    elif args.model == "TEST":
        try:
            print(args)
            print(args.number)
        except:
            print("test pipline error")
    else:
        print("model")
