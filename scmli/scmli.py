#!/usr/bin/python3

import argparse, os
from function import pcr_pipline, print_c

def create_arg_parser():
	parser = argparse.ArgumentParser(description = "single cell mutant library inspertion")
	parser.add_argument('-m', '--model', default = 'print_c')
	parser.add_argument('-s', '--seq', default = 'ATGC' )
	parser.add_argument('-r1', '--read1')
	parser.add_argument('-r2', '--read2')
	parser.add_argument('-n', '--project_name', default = 'my_project')
	return parser

def pcr_pipline():

	#fastqc
	print("fastqc")
	os.system('fastqc -o ./ ' + arg.read1 + ' ' + arg.read2)
	
	#trim
	print("trim_galore")
	os.system('trim_galore --paired --fastqc --max_n 0 -j 4 --gzip ' + arg.read1 + ' ' + arg.read2)
	
	#get name
	name1_1 = arg.read1.split("/")[-1]
	name1_2 = name1_1.split(".fq.gz")[0] + "_val_1.fq.gz"
	name2_1 = arg.read2.split("/")[-1]
	name2_2 = name2_1.split(".fq.gz")[0] + "_val_2.fq.gz"
	
	os.system('gzip -cd ' + name1_2 + ' ' + name2_2 + ' > ' + arg.project_name + '.fq')

	os.system('python3 parse_gRNA.py my_project.fq NoIMET1_gRNAs.csv GGTAGAATTGGTCGTTGCCATCGACCAGGC > PCR-ku.stats')


if __name__ == "__main__":
	parser = create_arg_parser()
	arg = parser.parse_args()
	if arg.model == "PCR":
		pcr_pipline()



