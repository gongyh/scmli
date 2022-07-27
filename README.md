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
lxml (4.9.1) (python package)<br />
argparse (python package, only needed if python<=3.6)<br />
fastqc (0.11.9)<br />
trim-galore>=0.6.0 (0.6.7)<br />
r-base (3.6.1)<br />
r-ggplot2<br />

The tested versions are given in parentheses.


You can install these dependencies using Conda ([Miniconda3](https://docs.conda.io/en/latest/miniconda.html)):
```
conda install pandas biopython lxml r-base r-ggplot2
conda install -c bioconda fastqc trim-galore
```

## Usage

Sclmi search `reads` which target sequence had been transferred in. It uses `fixed sequence` (all sequences before gRNAs) for filtering successfully transferred plasmid, then search
gene sequence with `gene library file`. Target sequence contains fixed sequence and gene sequence, `number(a b)` is used to locate gene sequences in gRNAs.
```
required: reads(.fq.gz), gene library(.csv), fixed sequence(str)
usage: scmli.py [-h] [-m {PCR,TEST}] -l LIB -s SEQ -r1 READ1 -r2 READ2 [-num NUMBER NUMBER] [-n OUTPUT_NAME] [-o OUTPUT_DIR]
```

## Arguments

```
required arguments:
  -m {PCR,TEST}, --model {PCR,TEST}          Choose analysis model
  -l LIB, --lib LIB                          Gene library file
  -s SEQ, --seq SEQ                          Start with the first base, all sequences before gRNAs
  -r1 READ1, --read1 READ1                   Read1 file
  -r2 READ2, --read2 READ2                   Read2 file

optional arguments:
  -h, --help                                 Show this help message and exit
  -t NUMBER, --threads NUMBER                Number of threads to use, default = 8
  --number NUMBER NUMBER                     Start and end of the gene position in gRNAs, default='25 45', from the 26-th to the 45-th bases
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

`file_fastqc.html/zip`: Quality control results(raw data) <br />
`file_val_1/2_fastqc.html/zip`: Quality control results(clean data) <br />
`file_trimming_report.txt`: Trim results
`my_project.counts`:      Raw count result <br />
`my_project.percent`:     Detailed count result <br />

| gene_id    | sequence             | counts  | percent   |
| ---------- | -------------------- | ------- | --------- |
| NO12G02480 | TCTATCTCAACAGCCACCCG | 17      | 0.0003771 |
| NO03G04750 | ACTTCCTGGTCCTCCCACGA | 17      | 0.0003771 |
| NO08G01490 | TGCCTCAGGAGGGATGATCG | 16      | 0.0003549 |
| NO02G03790 | GAGAACTTTTCATCCTCGCG | 16      | 0.0003549 |
| .......    | .......              | ....... | ......    |

`my_project.stats`: Statistical result <br />

| Key                           | Value    |
| -------                       | -------  |  
|raw_reads                      | 50000    |
|all_reads(clean reads)         | 49947    |
|valid_reads                    | 45085    |
|unknow_reads                   | 3393     | 
|gRNAs_reads                    | 41692    |
|all_kinds                      | 12649    |
|lib_kinds                      | 9709     |
|unknow_kinds                   | 2940     |
|gRNAs_kinds                    | 9368     |
|all/raw_reads_percent          | 0.99894  |
|valid/all_reads_percent        | 0.902657 |
|gRNAs/valid_reads_percent      | 0.924742 |
|gRNAs_coverage                 | 0.964878 |
|gRNAs_average_all              | 4.29416  |
|gRNAs_average_present          | 4.45047  |

`unknow.seq`: List of unknow sequences <br />
`my_project.log`: Process log <br />
`reads.plot`: Count of different kinds of reads <br />
`frequency.plot`: Frequency of all gRNAs <br />
`histogram.plot`: Count of different frequency of gRNAs <br />

## License

MIT

