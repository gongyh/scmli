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
    parser.add_argument('-num', '--number', type=int, nargs=2,
                        default=[25, 45],
                        help="Start and end of the gene position, \
                        '0 10' for the first ten, default='25 45'")
    parser.add_argument('-n', '--output_name', default="my_project",
                        help="Prefix of output files, default='my_project'")
    parser.add_argument('-o', '--output_dir', default="output",
                        help="Directory of output files, default='output'")

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
        pcr_qc(args.output_dir, args.output_name, args.read1, args.read2)
        pcr_parse_gRNA(args.lib, args.seq, args.number, args.output_name)
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
