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

We had data of ONT sequencing for tumor and normal cells.

We used an alignment of tumor read to GRCh38 reference, where each read was phased according to its primary alignment.

We used Snniffles and Mikhail's tool based on HapDup for determining breakpoints. Sniffles provided a less accurate result, so we focused on the HapDup result.

After that, we developed a Python3 script for visualizing read coverage by haplotype and breakpoints. We analysed possible breakage-fusion bridge (a possible way of forming a mutation) events and performed local assembly with Flye around the breakpoints to support or contradict our hypothesis.


## Results 

We developed a Python3 script for visualizing read coverage by haplotype and breakpoints. It works on a `bam` file and produces pictures in `png` format. It creates visualizations of every chromosome and also a zoom-in of possible breakage-fusion bridge events.

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








