import argparse
from counter import *

parser = argparse.ArgumentParser(description='Create coverage graphics from bam file.')
parser.add_argument('filename', metavar='filename', type=str, 
                    help='Alignment bam filename')

parser.add_argument('--chrom_sizes', metavar='chrom_sizes filename', type=str, action='store', default=None,
                    help='Filename of the file with chromosome sizes')

parser.add_argument('--block_size', dest='block_size', action='store',
                     default=100000,
                    help='Size of the block')


parser.add_argument('--reference_length', dest='reference_length', action='store',
                     default=30000000,
                    help='Reference length')


parser.add_argument('-o', dest='out', action='store',
                     default='.',
                    help='Output dir')


parser.add_argument('--plot', dest='plot', action='store_const',
                     const=True,
                    help='Should we print the graphics? Prints only coverage if breakpoint file is not provided')

parser.add_argument('--breakpoints', dest='breakpoints', action='store',
                     default=None,
                    help='Breakpoint csv file')

parser.add_argument('--ylim', dest='ylim', action='store',
                     default=500,
                    help='Y limits on graphs')

parser.add_argument('--zoomin', dest='zoomin', action='store_const',
                     const=True, help='Should we print the coverage around interesting breakpoints? Provide breakpoint file')

#parser.add_argument('--plot_only', dest='plot_only', action='store_const', TODO
#                     const=True,
#                    help='Only print the graphics for existing coverage')

#parser.add_argument('--plot_only', dest='zoomin_only', action='store_const',
#                     const=True,
#                    help='Only print the graphics for existing coverage')

args = parser.parse_args()

counts = get_counts(args.filename, chr_sizes=args.chrom_sizes, block_size=args.block_size)

for elem in counts:
    try:
        os.mkdir(args.out)
    except:
        pass
    counts[elem].to_csv(str(args.out) + '/' + str(elem)+"_counts.csv")
    
if args.plot == True:
    try:
        os.mkdir(args.out + '/visualization')
    except:
        pass
    print(counts)
    if args.breakpoints is not None:
        somatic = pd.read_csv(args.breakpoints, header=None, sep = '\t')
        somatic.columns = ['chrA', 'posA', 'orA', 'chrB', 'posB', 'orB', 'haplo', 'reads', 'ref_reads']
    else:
        somatic = pd.DataFrame
        somatic.columns = ['chrA', 'posA', 'orA', 'chrB', 'posB', 'orB', 'haplo', 'reads', 'ref_reads']
    visualize(counts, somatic, args.ylim, savefig=True, filename = str(args.out) + '/' + 'visualization/coverage.png')
    
if args.zoomin == True:
    try:
        os.mkdir(args.out + '/visualization')
    except:
        pass
    somatic = pd.read_csv(args.breakpoints, header=None, sep = '\t')
    somatic.columns = ['chrA', 'posA', 'orA', 'chrB', 'posB', 'orB', 'haplo', 'reads', 'ref_reads']
    if args.block_size != 1000:
        counts = get_counts(args.filename, chr_sizes=args.chrom_sizes, block_size=1000)
    plot_nearest(counts, somatic, savefig=False, dirname= args.out + '/' + 'visualization/')
    

#if args.plot_only == True:
#    try:
#        os.mkdir(args.out + '/visualization')
#    except:
#        pass
#    
#    somatic = pd.read_csv(args.breakpoints, header=None, sep = '\t')
#    somatic.columns = ['chrA', 'posA', 'orA', 'chrB', 'posB', 'orB', 'haplo', 'reads', 'ref_reads']
#    visualize(counts, somatic, args.ylim, savefig=True, filename = str(args.out) + '/' + 'visualization/coverage.png')
  