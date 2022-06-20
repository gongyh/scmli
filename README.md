# scmli

Single cell mutant library inspection 

## Table of Contents

- [Background](#background)
- [Install](#install)
- [Arguments](#arguments)
- [Usage](#usage)
- [Results](#results)
- [License](#license)

## Background

Identification and description of gRNA mutant library.

## Install

```
git clone https://github.com/gongyh/scmli.git
```

### Software Requirements

Python3 (3.9)<br />
Biopython (1.79) (python package)<br />
pandas (1.4.2) (python package)<br />
argparse (python package, only needed if python<=3.6)<br />
matplotlib (3.5.1)<br />
fastqc (0.11.9)<br />
trim-galore>=0.6.0 (0.6.7)<br />

The tested versions are given in parentheses.


You can install these dependencies using Conda ([Miniconda3](https://docs.conda.io/en/latest/miniconda.html)):
```
conda install -c bioconda pandas biopython matplotlib fastqc trim-galore>=0.6.0
```

## Usage

Sclmi search `reads` which target sequence had been transferred in. It uses `fixed sequence` for filtering successfully transferred plasmid, then aligns 
gene sequence with `gene library file`. Target sequence contains fixed sequence and gene sequence, `number(a b)` is used to locate gene sequence.
```
required: reads(.fq.gz), gene library(.csv), fixed sequence(str), location number(a b)
usage: scmli.py [-h] [-m {PCR,TEST}] -l LIB -s SEQ -r1 READ1 -r2 READ2 [-num NUMBER NUMBER] [-n OUTPUT_NAME] [-o OUTPUT_DIR]
```

## Arguments

```
required arguments:
  -m {PCR,TEST}, --model {PCR,TEST}          Choose analysis model
  -l LIB, --lib LIB                          Gene library file
  -s SEQ, --seq SEQ                          The fixed sequence for search
  -r1 READ1, --read1 READ1                   Read1 file
  -r2 READ2, --read2 READ2                   Read2 file

optional arguments:
  -h, --help                                 Show this help message and exit
  -t NUMBER, --threads NUMBER                Number of threads to use, default = 8
  --number NUMBER NUMBER                     Start and end of the gene position, '0 10' for the first ten, default = "25 45"
  -n OUTPUT_NAME, --output_name OUTPUT_NAME  Prefix of output files, default = "my_project"
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR     Directory of output files, default = "output"
  --FASTQC_PATH                              PATH to fastqc
  --TRIM_GALORE_PATH                         PATH to trim-galore
```

## Test

```
cd scmli
python3 scmli.py -m PCR \
  -l test/NoIMET1_gRNAs.csv \
  -s GGTAGAATTGGTCGTTGCCATCGACCAGGC \
  -r1 test/test_R1.fq.gz \
  -r2 test/test_R2.fq.gz
```

## Results

{OUTPUT_NAME}.fastqc.zip: fastqc results  <br />
{OUTPUT_NAME}.counts:     count file <br />
{OUTPUT_NAME}.percent:    show the counts and percent of target sequences <br />

| gene_id    | sequence             | counts  | percent   |
| ---------- | -------------------- | ------- | --------- |
| NO12G02480 | TCTATCTCAACAGCCACCCG | 17      | 0.0003771 |
| NO03G04750 | ACTTCCTGGTCCTCCCACGA | 17      | 0.0003771 |
| NO08G01490 | TGCCTCAGGAGGGATGATCG | 16      | 0.0003549 |
| NO02G03790 | GAGAACTTTTCATCCTCGCG | 16      | 0.0003549 |
| NO01G05060 | GTTGCCTCTTACCCCACCCA | 15      | 0.0003327 |
| NO14G00420 | TTGATTCGAAGAATGAGTGT | 15      | 0.0003327 |
| .......    | .......              | ....... | ......    |

## License

MIT

