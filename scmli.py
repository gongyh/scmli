#!/usr/bin/env python3

import argparse
import os
from libs.grna import grna_pipeline
from libs.variant import variant_pipeline
from libs.common import existing_file, existing_dir


def create_arg_parser():

    parser = argparse.ArgumentParser(
        description="single cell mutant library inspertion")
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_gRNA = subparsers.add_parser('gRNA', help='search gRNAs')
    parser_variant = subparsers.add_parser('variant', help='call variant')

    parser_gRNA.add_argument(
        '-l', dest='lib', required=True, type=existing_file, help="Gene library")
    parser_gRNA.add_argument('-s', dest='seq', required=True,
                             help="The fixed sequence for search")
    parser_gRNA.add_argument(
        '-r1', dest='read1', required=True, type=existing_file, help="Read1")
    parser_gRNA.add_argument(
        '-r2', dest='read2', required=True, type=existing_file, help="Read2")
    parser_gRNA.add_argument('-t', dest='threads',
                             type=int, default=8, help="Number of threads")
    parser_gRNA.add_argument('--number', type=int, nargs=2, default=[25, 45],
                             help="Start and end of the gene position, default='25 45', from the 26-th to the 45-th bases")
    parser_gRNA.add_argument('-n', dest='output_name', default="my_project",
                             help="Prefix of output files, default='my_project'")
    parser_gRNA.add_argument('-o', dest='output_dir', default="output",
                             help="Directory of output files, default='output'")
    parser_gRNA.add_argument('--FASTQC_PATH', help="PATH to fastqc")
    parser_gRNA.add_argument('--TRIM_GALORE_PATH', help="PATH to trim_galore")
    parser_gRNA.set_defaults(func=gRNA)

    parser_variant.add_argument(
        '-r1', dest='read1', required=True, type=existing_file, help='Read1')
    parser_variant.add_argument(
        '-r2', dest='read2', required=True, type=existing_file, help='Read2')
    parser_variant.add_argument(
        '--ref', dest='ref', type=existing_file, help='reference')
    parser_variant.add_argument(
        '--target', dest='target', required=True, type=existing_file, help='target')
    parser_variant.add_argument(
        '--dtarget', dest='dtarget', required=True, type=existing_file, help='dtarget')
    parser_variant.add_argument(
        '-t', dest='threads', type=int, default=8, help='Number of threads')
    parser_variant.add_argument('-n', dest='outname', default='my_project',
                                help="Prefix of output files, default='my_project'")
    parser_variant.add_argument('-o', dest='outdir', default='output',
                                help="Directory of output files, default='output'")
    parser_variant.set_defaults(func=variant)

    return parser


def gRNA(args):
    # check and modify args
    args.seq = args.seq.upper()
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    args.output_dir = os.path.abspath(args.output_dir)
    # check software path
    if args.FASTQC_PATH:
        if os.path.isdir(args.FASTQC_PATH) and (os.path.isfile(args.FASTQC_PATH+"/fastqc") or os.path.isfile(args.FASTQC_PATH+"fastqc")):
            args.FASTQC_PATH = os.path.abspath(args.FASTQC_PATH)
        else:
            print("path to fastqc error")
            os._exit(0)
    else:
        pass

    if args.TRIM_GALORE_PATH:
        if os.path.isdir(args.TRIM_GALORE_PATH) and (os.path.isfile(args.TRIM_GALORE_PATH+"/trim_galore") or os.path.isfile(args.TRIM_GALORE_PATH+"trim_galore")):
            args.TRIM_GALORE_PATH = os.path.abspath(args.TRIM_GALORE_PATH)
        else:
            print("path to trim_galore error")
            os._exit(0)
    else:
        pass

    current_dir = os.getcwd()  # curring working directory
    scmli_dir = os.path.split(os.path.realpath(__file__))[
        0]  # scmli root directory
    os.chdir(args.output_dir)  # change to output directory
    os.system('cp '+scmli_dir+'/doc/result.html ./')
    stats = grna_pipeline(args.output_name, args.read1, args.read2, args.FASTQC_PATH,
                          args.TRIM_GALORE_PATH, args.threads, args.lib, args.seq, args.number, scmli_dir)
    os.system('Rscript '+scmli_dir+'/libs/plot.r '+args.output_name)
    print("grna Finished")
    os.chdir(current_dir)

    return stats


def variant(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    args.outdir = os.path.abspath(args.outdir)
    current_dir = os.getcwd()
    scmli_dir = os.path.split(os.path.realpath(__file__))[0]
    os.chdir(args.outdir)
    variant_pipeline(args)
    print("Variant Finished")
    os.chdir(current_dir)
    print(args)


if __name__ == "__main__":
    try:
        parser = create_arg_parser()
        args = parser.parse_args()
        args.func(args)
    except:
        print("parser error")
        os._exit(0)
    else:
        print("scmli")
