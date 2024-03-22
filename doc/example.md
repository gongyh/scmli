# Experimental Workflow Demonstration

## Experiment Background
To construct a whole-genome random mutation library in an efficient, reliable and economical manner, we developed a specialized experimental design and corresponding analysis workflow. We utilized 24 plates, each containing 24 wells, as one batch, with each well representing an individual mutation experiment. For the same plate and the same well, we added unique barcodes to the gRNA PCR products and proceeded with next-generation sequencing. In parallel, we performed whole-genome sequencing for samples from each plate. By identifying gRNAs present in both the plate and the well, we could accurately determine the specific gRNA for each well and assess mutations at the targeted locations, validating our results.

## Input Files
1.Sequencing results of gRNA PCR products with barcodes (1 set) <br>
2.Whole-genome sequencing data for each plate (24 sets) <br>
3.Reference genome (.gbk file) <br>
4.gRNA sequences (.csv file) <br>
5.Fixed sequence portion of the gRNA vector tools <br>
6.All gRNAs corresponding to mutation target regions (.bed file) <br>
7.Specific target gRNAs corresponding to mutation target regions (.bed file) <br>

## Analysis Workflow

### 1. Data Demultiplexing
```
#!/usr/bin/bash
seqtk_demultiplex -b seqtk_barcode.tsv -1 gRNA_1.fq.gz -2 gRNA_2.fq.gz -l 10 -d output
```
Demultiplex into 24 plates (01a-24a) and 24 wells (01b-24b) sequences

### 2. Detect gRNA
```
#!/usr/bin/env python3
for i in ["{:02d}".format(i) for i in range(1, 25)]:
    for j in ["a", "b"]:
        r1 = f"{i}{j}_R1.fastq.gz"
        r2 = f"{i}{j}_R2.fastq.gz"
        cmd = f'python scmli.py gRNA -r1 {r1} -r2 {r2} -l gRNAs.csv -sCCCTCCATCCACAGAATCGATATATCTGCCTCGCAGGCATAGTTTAGTGGTAGAATTGGTCGTTGCCATCGACCAGGC -o test -n output'
        os.system(cmd)
```
### 3. Confirm gRNA
```
#!/usr/bin/bash
python3 locate_gRNA.py
```
### 4. Call variants
```
#!/usr/bin/env python3
for i in ["{:02d}".format(i) for i in range(1, 25)]:
    r1 = f"{i}_genome_R1.fastq.gz"
    r2 = f"{i}_genome_R2.fastq.gz"
    cmd = f'python scmli.py variant -r1 {r1} -r2 {r2} --ref test/genes.gbk --target test/targets.bed --dtarget {i}.bed'
	os.system(cmd)
```
