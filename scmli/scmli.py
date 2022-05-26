#!/usr/bin/env python3

import argparse, os
from pcr_pipline import pcr_qc, pcr_count
from parse_gRNA import parse_gRNA
import pandas as pd

def create_arg_parser():

    parser = argparse.ArgumentParser(description = "single cell mutant library inspertion")
    parser.add_argument('-m', '--model', default = "TEST", choices = ["PCR", "TEST"])
    parser.add_argument('-l', '--lib', required = True)
    parser.add_argument('-s', '--seq', required = True, help = "The fixed sequence for search")
    parser.add_argument('-r1', '--read1', required =True)
    parser.add_argument('-r2', '--read2', required =True)
    parser.add_argument('-n', '--output_name', default = 'my_project')
    parser.add_argument('-o', '--output_dir', default=os.getcwd())

    return parser


def parse_args(parser):
    args = parser.parse_args()
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
        args = parse_args(parser)
    except:
        print("parser error")
        os._exit(0)
    if args.model == "PCR":
        pcr_qc(args.output_dir, args.output_name, args.read1, args.read2)
        pcr_count(args.output_name + ".fq", args.lib, args.seq)
        print(args)
    elif args.model == "TEST":
        try:        
            print(args)
        except:
            print("test pipline error") 
    else:
        print("model")
