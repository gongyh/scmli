#!/usr/bin/env python3

import argparse
import os
from libs.grna import grna_pipeline
from libs.variant import variant_pipeline
from libs.common import existing_file, existing_dir


def create_arg_parser():

    # Create an argument parser object
    parser = argparse.ArgumentParser(
        description="single cell mutant library inspertion")
    # Create a subparser object to group arguments
    subparsers = parser.add_subparsers(help='sub-command help')

    # Create a parser object for gRNA
    parser_gRNA = subparsers.add_parser('gRNA', help='search gRNAs')
    # Create a parser object for variant
    parser_variant = subparsers.add_parser('variant', help='call variant')

    # Add arguments to the gRNA parser
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

    # Add arguments to the variant parser
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
        '--ram', dest='ram', type=int, default=4, help='ram GB')
    parser_variant.add_argument(
        '-t', dest='threads', type=int, default=8, help='Number of threads')
    parser_variant.add_argument('-n', dest='outname', default='my_project',
                                help="Prefix of output files, default='my_project'")
    parser_variant.add_argument('-o', dest='outdir', default='output',
                                help="Directory of output files, default='output'")
    parser_variant.set_defaults(func=variant)

    # Return the argument parser object
    return parser


def gRNA(args):
    # Convert the sequence to uppercase
    args.seq = args.seq.upper()
    # Create the output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    # Get the absolute path of the output directory
    args.output_dir = os.path.abspath(args.output_dir)
    # Check if the path to the fastqc executable is valid
    if args.FASTQC_PATH:
        fastqc_executable = os.path.join(args.FASTQC_PATH, "fastqc")
        if os.path.isdir(args.FASTQC_PATH) and os.path.isfile(fastqc_executable):
            args.FASTQC_PATH = os.path.abspath(args.FASTQC_PATH)
        else:
            print("path to fastqc error")
            os._exit(0)
    else:
        pass

    # Check if the path to the trim_galore executable is valid
    if args.TRIM_GALORE_PATH:
        trim_galore_executable = os.path.join(args.TRIM_GALORE_PATH, "fastqc")
        if os.path.isdir(args.TRIM_GALORE_PATH) and os.path.isfile(trim_galore_executable):
            args.TRIM_GALORE_PATH = os.path.abspath(args.TRIM_GALORE_PATH)
        else:
            print("path to trim_galore error")
            os._exit(0)
    else:
        pass

    # Get the current working directory
    current_dir = os.getcwd()
    # Get the absolute path of the scmli directory
    scmli_dir = os.path.split(os.path.realpath(__file__))[0]
    # Change the current working directory to the output directory
    os.chdir(args.output_dir)
    # Copy the result.html file from the scmli directory to the output directory
    os.system('cp '+scmli_dir+'/doc/result.html ./')
    # Run the grna pipeline
    stats = grna_pipeline(args, scmli_dir)
    # Run the plot.r to plot
    os.system('Rscript '+scmli_dir+'/libs/plot.r '+args.output_name)
    # Print a message
    print("gRNA Finished")
    # Change the current working directory back to the original
    os.chdir(current_dir)

    # Return the stats
    return stats


def variant(args):
    # Create the output directory if it does not exist
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    # Get the absolute path of the output directory
    args.outdir = os.path.abspath(args.outdir)
    # Get the current working directory
    current_dir = os.getcwd()
    # Get the absolute path of the scmli directory
    scmli_dir = os.path.split(os.path.realpath(__file__))[0]
    # Change the current working directory to the output directory
    os.chdir(args.outdir)
    # Run the variant pipeline
    variant_pipeline(args)
    # Print a message
    print("Variant Finished")
    # Change the current working directory back to the original
    os.chdir(current_dir)


# Define the main function
def main():
    try:
        # Create the argument parser
        parser = create_arg_parser()
        # Parse the arguments
        args = parser.parse_args()
        # Run the function
        args.func(args)
    # Catch any errors
    except:
        print("parser error")
        os._exit(0)
    else:
        print("scmli")


# Run the main function
if __name__ == "__main__":
    main()
