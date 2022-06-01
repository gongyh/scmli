# scmli


Single cell mutant library inspection 

## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
- [Arguments](#arguments)
- [Results](#results)
- [License](#license)

## Background
Identification and description of gRNA mutant library 

## Install

```
git clone https://github.com/gongyh/scmli.git
```

### Software Requirements

Miniconda3<br />
Python 3.9 <br />
Biopython 1.78 (python package)<br />
pandas 1.4.2 (python package)<br />
fastqc 0.11.9<br />
trim_galore 0.6.7<br />

Conda install software:
```
conda install -c bioconda pandas biopython fastqc trim-galore
```
You could also install dependencies by other methods.

## Usage
```
required: sequence(.fq.gz), gRNA library(.csv), fixed sequence(str)

usage: scmli.py [-h] [-m {PCR,TEST}] -l LIB -s SEQ -r1 READ1 -r2 READ2 [-num NUMBER NUMBER] [-n OUTPUT_NAME] [-o OUTPUT_DIR]
```

## Test
```
cd scmli #
python3 scmli.py -m PCR \
  -l test/NoIMET1_gRNAs.csv \
  -s GGTAGAATTGGTCGTTGCCATCGACCAGGC \
  -r1 test/test_R1.fq.gz \
  -r2 test/test_R2.fq.gz
```

## Arguments
```
optional arguments:
  -h, --help                                 Show this help message and exit
  -m {PCR,TEST}, --model {PCR,TEST}          Choose analysis model
  -l LIB, --lib LIB                          Gene library file
  -s SEQ, --seq SEQ                          The fixed sequence for search
  -r1 READ1, --read1 READ1                   Read1 file
  -r2 READ2, --read2 READ2                   Read2 file
  -num NUMBER NUMBER, --number NUMBER NUMBER Start and end of the gene position, '0 10' for the first ten
  -n OUTPUT_NAME, --output_name OUTPUT_NAME  Prefix of output files
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR     Directory of output files
```

## Results
percent.stats: show the counts and percent of target sequences

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

