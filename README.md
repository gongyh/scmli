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
r-ggplot2 (3.3.5)<br />
bcftools (1.15)
snippy (4.6.0 modified)

The tested versions are given in parentheses.


You can install these dependencies using Conda ([Miniconda3](https://docs.conda.io/en/latest/miniconda.html)):
```
conda install -c bioconda fastqc trim-galore pandas biopython lxml r-base r-ggplot2 bcftools snippy
```

## Usage
### gRNA model
Sclmi searches `reads` which have target gRNAs sequence. It uses `fixed sequence` (all sequencing bases before gRNAs in forward reads without adapter) for filtering valid reads, then searches
gene-special gRNAs sequence with `gRNAs library file`. The gRNAs library sequence contains universal sequence and gene-special sequence, `number(a b)` is used to locate gene-special sequences in gRNAs. <br />
`gRNAs_library.csv`:  <br />
NO01G00240,ccgggtccgattcccggtgcctgcaGAGTGTGGTGGAATTTGCCGgttttagagctagaaatagcaagttaaaataag <br />
NO01G00250,ccgggtccgattcccggtgcctgcaACACGATAGTCAAGACGCTGgttttagagctagaaatagcaagttaaaataag <br />
  ...... , ...... <br />
```
required: reads(fastq file), fixed sequence(str), gRNAs library(.csv)
```
### variant model <br />

## Arguments
### gRNA model
```
required arguments:
  -l LIB                            gRNAs library file
  -s SEQ                            All sequencing bases before gRNAs in forward reads without adapter
  -r1 READ1                         Read1 fastq file
  -r2 READ2                         Read2 fastq file

optional arguments:
  -h, --help                        Show help message and exit
  -t NUMBER                         Number of threads, default = 8
  --number NUMBER NUMBER            Start and end of the gene-special position in gRNAs,
                                    default='25 45', from the 26-th to the 45-th bases
  -n OUTPUT_NAME                    Prefix of output files, default = "my_project"
  -o OUTPUT_DIR                     Directory of output files, default = "output"
  --FASTQC_PATH                     PATH to fastqc
  --TRIM_GALORE_PATH                PATH to trim-galore
```
### variant model
```
required arguments:
  -r1 READ1                         Read1 fastq file
  -r2 READ2                         Read2 fastq file
  --ref REF                         reference
  --target TARGET                   target

optional arguments:
  -h, --help                        show this help message and exit
  -t THREADS                        Number of threads
  -n OUTNAME                        Prefix of output files, default='my_project'
  -o OUTDIR                         Directory of output files, default='output'
```

## Test

```
cd scmli
python3 scmli.py gRNA \
  -l test/NoIMET1_gRNAs.csv \
  -s GGTAGAATTGGTCGTTGCCATCGACCAGGC \
  -r1 test/test_R1.fq.gz \
  -r2 test/test_R2.fq.gz

python3 scmli.py variant \
  -r1 test/test_R1.fq.gz \
  -r2 test/test_R2.fq.gz \
  --ref test/genes.gbk \
  --target test/targets.bed
```

## Results
### gRNA model
`file_fastqc.html/zip`: Quality control results(raw data) <br />
`file_val_1/2_fastqc.html/zip`: Quality control results(clean data) <br />
`file_trimming_report.txt`: Trim results <br />
`my_project.counts`:      Raw count result <br />
`my_project.percentage`:  Detailed count result <br />

| gene_id    | sequence             | counts  | percentage | percentage_gRNAs | accumulative_unknow_percentage |
| ---------- | -------------------- | ------- | ---------- | ---------------- | ------------------------------ |
| NO12G02480 | TCTATCTCAACAGCCACCCG | 17      | 0.037707   | 0.040775         | 0.0                            |
| NO03G04750 | ACTTCCTGGTCCTCCCACGA | 17      | 0.037707   | 0.040775         | 0.0                            |
| NO08G01490 | TGCCTCAGGAGGGATGATCG | 16      | 0.035489   | 0.040775         | 0.0                            |
| NO02G03790 | GAGAACTTTTCATCCTCGCG | 16      | 0.035489   | 0.040775         | 0.0                            |
| .......    | .......              | ....... | ......     | ......           | ......                         |

`my_project.stats`: Statistical result <br />

| Key                           | Value    |
| -------                       | -------  |  
|raw_reads(paired)              | 50000    |
|all_reads(clean reads,paired)  | 49947    |
|valid_reads                    | 45085    |
|unknow_reads                   | 3393     | 
|gRNAs_reads                    | 41692    |
|all_kinds                      | 12649    |
|lib_kinds                      | 9709     |
|unknow_kinds                   | 2940     |
|gRNAs_kinds                    | 9368     |
|all/raw_reads_percent          | 0.99894  |
|valid/all_reads_percentage     | 0.902657 |
|gRNAs/valid_reads_percentage   | 0.924742 |
|gRNAs_coverage                 | 0.964878 |
|gRNAs_average_all              | 4.29416  |
|gRNAs_average_detected         | 4.45047  |

`unknow.seq`: List of unknow sequences <br />
`my_project.log`: Process log <br />
`reads.plot`: Count of different kinds of reads <br />
`frequency.plot`: Frequency of all gRNAs <br />
`frequency_detected.plot`: Frequency of detected gRNAs <br />
`frequency_distribution.plot`: Count of different frequency of all gRNAs <br />
`frequency_distribution_detected.plot`: Count of different frequency of detected gRNAs <br />
`accumulative_unknow_percentage.plot`: Percentage of accumulative unknow sequences <br />

### variant model

## License

MIT

