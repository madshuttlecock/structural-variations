# structural-variations
Studying complex structural variations in cancer using long reads



# Manual


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
