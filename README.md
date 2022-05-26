# scmli


Single cell mutant library inspection 

## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
- [Arguments](#arguments)
- [License](#license)

## Background
Identification and description of gRNA mutant library 
## Install

```
git clone https://github.com/gongyh/scmli.git
```
## Usage
```
required: sequence(.fq.gz), gRNA library(.csv), fixed sequence(str)

useage: scmli.py [-h] [-m {PCR,TEST}] -l LIB [-s SEQ] [-r1 READ1] [-r2 READ2] [-n PROJECT_NAME] [-o OUTDIR]
                
example: scmli.py -m PCR -l NoIMET1_gRNAs.csv -s GGTAGAATTGGTCGTTGCCATCGACCAGGC -r1 test1.fq.gz -r2 test2.fq.gz
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


## License

         
