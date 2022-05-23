# Studying complex structural variations in cancer using long reads

Olga Kalinichenko (BI),

Supervisor: Mikhail Kolmogorov (NIH/NCI)

## Introduction

Cancer is caused by genetic changes. There are numerous mutations in tumor cells including 
SNPs, insertions, deletions, inversions, translocations, chromosomal aberrations. 
Small mutations have been studied thoroughly using short read sequencing. However, for big structural variations, it it too hard to determine them using short reads dut to read mapping ambiguities. Recent advances in long read sequencing present an opportunity to solve this problem.
There are many problems to cope with, such as the heterogeneity of tumor cells and loss of information about haplotypes of the cells.

![alt text](https://github.com/madshuttlecock/structural-variations/blob/main/ris.png)


## Aims and objectives

We aim to find complex rearrangements in cancer, including information about haplotypes. The result should be the set of all 'new' chromosomes (which may consist of verious regions pf initial chromosomes) in tumor and a description of structural variations.  

The objectives are as follows:

* Find breakpoints in tumor data
* Develop a script for visualising coverage and breakpoints for haplotypes
* Analyze results and identify possible variations and their origin
* Perform local assembly around breakpoints to sharpen the conclusions


## Methods

We had data of ONT sequencing for tumor and normal cells int he form of an already prepared alignment of tumor and normal reads to GRCh38 reference. Each read was already phased according to its primary alignment (this information was stored in the `HP` tag).

We use Sniffles and Mikhail's own tool (part of HapDup) for determining breakpoints. Sniffles provided a less accurate result, so we focused on the HapDup result.

We found breakpoints with Sniffles using the following commands:

`sniffles --input <tumor alignment in bam format> --vcf tumor_variants_sniffles.vcf --output-rnames --threads 10`  for tumor

`sniffles --input <normal alignment in bam format> --vcf normal_variants_sniffles.vcf --output-rnames --threads 10` for normal

After that, we parsed the results with Python3.

HapDup results were provided by Mikhail.

We found that Sniffles results had a lot of breakpoints present in normal and not in tumor and some breakpoints present in both samples had different orientation. It seemed strange to us, because mutations present in normal shoul very likely be present in tumor also. We believed there were a lot of false positive results, so we decided not to use Snniffles and focused only on HapDup results.

After that, we developed a Python3 script for visualizing read coverage by haplotype and breakpoints. We analysed possible breakage-fusion bridge (a possible way of forming a mutation) events and performed local assembly with Flye assembler (version 2.9-b1768) around the breakpoints to support or contradict our hypothesis. We chose Flye assembler because it was specially designed for ONT reads.

First, we extracted reads covering the region with samtools (version 1.14):

`samtools view -b <tumor alignment in bam format> "chr3:23000000-27500000" > chr3_23-27.5.bam`

Then, we transformed the file to `fastq` format:

`samtools fastq chr3_23-27.5.sam > chr3_23-27.5.new.fastq`

After that, we performed local assembly using Flye:

`flye --nano-raw chr3_23-27.5.new.fastq --out-dir chr3_23-27res`

After that we aligned the resulting contig using minimap2 (version 2.24-r1122):

`minimap2 -ax map-ont  reference.fasta ./chr325-27res/assembly.fasta > chr3_local.sam`

The resulting contig had the same sequence as reference. So, we tracked earlier stages of the alignment (disjointigs and alignment graph) and aligned them to reference.

We examined the alignments in IGV genome browser (version 2.11).


## Results 

We developed a Python3 script for visualizing read coverage by haplotype and breakpoints. It works on a `bam` file and produces pictures in `png` format. It creates visualizations of every chromosome and also a zoom-in of all possible breakage-fusion bridge events.

### Manual


~~~text
usage: script.py [-h] [--chrom_sizes chrom_sizes filename]
                 [--block_size BLOCK_SIZE]
                 [--reference_length REFERENCE_LENGTH] [-o OUT] [--plot]
                 [--breakpoints BREAKPOINTS] [--ylim YLIM] [--zoomin]
                 filename

Create coverage graphics from bam file.

positional arguments:
  filename              Alignment bam filename

options:
  -h, --help            show this help message and exit
  --chrom_sizes chrom_sizes filename
                        Filename of the file with chromosome sizes
  --block_size BLOCK_SIZE
                        Size of the block
  --reference_length REFERENCE_LENGTH
                        Reference length
  -o OUT                Output dir
  --plot                Should we print the graphics? Prints only coverage if
                        breakpoint file is not provided
  --breakpoints BREAKPOINTS
                        Breakpoint csv file
  --ylim YLIM           Y limits on graphs
  --zoomin              Should we print the coverage around interesting
                        breakpoints? Provide breakpoint file
~~~

### Visualization for COLO829

Here, the triangles indicate the breakpoints.

Color indicates the haplotype of the breakpoint.

Triangle direction indicates the following:

* Right: the part of the chromosome on the left of the breakpoint is involved in the rearrangement.
* Left: the part of the chromosome on the right of the breakpoint is involved in the rearrangement.

The transparency indicate read support (the more transparent a triangle is, the less reads support the breakpoint).

Text near breakpoints indicates the chromosome of the second sequence of the breakpoint.

![alt text](https://github.com/madshuttlecock/structural-variations/blob/main/res1.png)

### Visualization for COLO829: zoom-in on chromosome 3

Here is an example of a zoom-in around an interecting position on chromosome 3.


![zoom2](https://user-images.githubusercontent.com/22745262/169700274-430db595-a9f5-4334-80a4-183228331b1e.png)
![zoom1](https://user-images.githubusercontent.com/22745262/169700270-39d4ae0d-d68a-4c26-bb6a-5d47576cf4f1.png)


### Analysis for chromosome 3

There is a potential breakage-fusion-bridge event on chromosome 3 (see above).

We examined it in IGV:
![chr3](https://github.com/madshuttlecock/structural-variations/blob/main/chr3.png)

We performed local assembly in this region and analysed the assembly graph as well as resulting contigs.


<img src="https://github.com/madshuttlecock/structural-variations/blob/main/graph.png" width="100">

We found that indeed in tumor there is a different chromosome structure resembling a breakage-fusion-bridge.

## Future plans

* Polish the script
* Perform local assembly around all possible breakage-fusion-bridge events
* Analyse results and determine the exact compositions of tumor chromosomes around breakpoints
* Find the whole set of rearrangements






