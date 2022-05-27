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

Python 3.x <br />
biopython 1.78 (python package)<br />
pandas 1.4.2 (python package)<br />
argparse 1.1 (python package)<br />
fastqc 0.11.9<br />
trim_galore 0.6.7<br />
Miniconda3<br />

To install python package:
```
pip install --user -r requirements.txt
```
To install software:
```
conda install -c bioconda fastqc
conda install -c bioconda trim-galore
```

## Usage
```
required: sequence(.fq.gz), gRNA library(.csv), fixed sequence(str)

usage: scmli.py [-h] [-m {PCR,TEST}] -l LIB [-s SEQ] [-r1 READ1] [-r2 READ2] [-n PROJECT_NAME] [-o OUTDIR]
                
example: python3 scmli.py -m PCR -l ../test/NoIMET1_gRNAs.csv -s GGTAGAATTGGTCGTTGCCATCGACCAGGC -r1 ../test/test1.fq.gz -r2 ../test/test2.fq.gz 
```

## Test
```
cd {dir of scmli} #
python3 scmli/scmli.py -m PCR -l test/NoIMET1_gRNAs.csv -s GGTAGAATTGGTCGTTGCCATCGACCAGGC -r1 test/test1.fq.gz -r2 test/test2.fq.gz -o {output_dir}
```

## Arguments
```
optional arguments:
  -h, --help            show this help message and exit
  -m {PCR,TEST}, --model {PCR,TEST}
  -l LIB, --lib LIB
  -s SEQ, --seq SEQ     The fixed sequence for search
  -r1 READ1, --read1 READ1
  -r2 READ2, --read2 READ2
  -n OUTPUT_NAME, --output_name OUTPUT_NAME
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
```

## Results
percent.stats: show the counts and percent of target sequences
```
gene_id	 sequence        counts  percent
NO12G02480      TCTATCTCAACAGCCACCCG    17      0.0003771
NO03G04750      ACTTCCTGGTCCTCCCACGA    17      0.0003771
NO08G01490      TGCCTCAGGAGGGATGATCG    16      0.0003549
NO02G03790      GAGAACTTTTCATCCTCGCG    16      0.0003549
NO01G05060      GTTGCCTCTTACCCCACCCA    15      0.0003327
NO14G00420      TTGATTCGAAGAATGAGTGT    15      0.0003327
.......         .......                 ....... ......
```

## License
MIT

