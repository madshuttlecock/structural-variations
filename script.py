import argparse
from counter import *

parser = argparse.ArgumentParser(description='Create coverage graphics from bam file.')

parser.add_argument('--bam', metavar='bam file', type=str, 
                    help='Alignment bam filename. If --count is not provided, this argument is ignored.')

parser.add_argument('--bams', metavar='bam files', type=list, default=None,
                    help='Alignment bam filenames. If --count is not provided, this argument is ignored.')


parser.add_argument('--coverage', metavar='coverage folder', type=str, default='.',
                    help='Counts folder obtained previously.')


parser.add_argument('--chrom_sizes', metavar='chrom_sizes filename', type=str, action='store', default=None,
                    help='Filename of the file with chromosome sizes')

parser.add_argument('--block_size', dest='block_size', action='store',
                     default=100000,
                    help='Size of the block, default = 100000')


parser.add_argument('--reference_length', dest='reference_length', action='store',
                     default=300000000,
                    help='Reference length')


parser.add_argument('-o', dest='out', action='store',
                     default='./script',
                    help='Output dir')


parser.add_argument('--breakpoints', dest='breakpoints', action='store',
                     default=None,
                    help='Breakpoint csv file')

parser.add_argument('--ylim', dest='ylim', action='store',
                     default=500, type=int,
                    help='Y limits on graphs')

parser.add_argument('--count', dest='count', action='store_const',
                     const=True,
                    help='Should we count coverage?')

parser.add_argument('--plot', dest='plot', action='store_const',
                     const=True,
                    help='Plot coverage graphics for chromosomes. Provide --coverage or --count option with bam file. May provide breakpoint ')


parser.add_argument('--potential', dest='potential', action='store_const',
                     const=True, help='Find popential BFB mutations (TODO descr). Should we plot the coverage around interesting breakpoints? Provide breakpoint file.')

parser.add_argument('--variation', dest='variation', action='store_const',
                     const=True, help='Find and plot breakpoints around copy number variations. Provide breakpoint file.')


args = parser.parse_args()


try:
    os.mkdir(args.out)
except:
    pass

print(args.bams)

if args.count:
    if args.chrom_sizes is not None:
        counts = get_counts(args.bam, chr_sizes=args.chrom_sizes, block_size=args.block_size, log_dir = args.out + '/' + 'coverage')   
    else:
        print(args.block_size)
        counts = get_counts(args.bam, reference_length=args.reference_length, block_size=args.block_size, log_dir = args.out + '/' + 'coverage')
else:
    try:
        counts = read_counts(args.coverage)
    except:
        counts = read_counts(args.coverage + '/coverage')
        
if args.plot == True:
    try:
        os.mkdir(args.out + '/visualization')
    except:
        pass
    
    if args.breakpoints is not None:
        somatic = pd.read_csv(args.breakpoints, header=None, sep = '\t')
        somatic.columns = ['chrA', 'posA', 'orA', 'chrB', 'posB', 'orB', 'haplo', 'reads', 'ref_reads', 'a', 'b', 'c']
    else:
        somatic = pd.DataFrame()
        #somatic.columns = ['chrA', 'posA', 'orA', 'chrB', 'posB', 'orB', 'haplo', 'reads', 'ref_reads', 'a', 'b', 'c']
        
    visualize(counts, somatic, int(args.ylim), savefig=True, filename = str(args.out) + '/' + 'visualization/coverage.png')
    
if args.potential == True:
    try:
        os.mkdir(args.out + '/potential')
    except:
        pass
    somatic = pd.read_csv(args.breakpoints, header=None, sep = '\t')
    somatic.columns = ['chrA', 'posA', 'orA', 'chrB', 'posB', 'orB', 'haplo', 'reads', 'ref_reads']
    if args.block_size != 100:
        counts = get_counts(args.filename, chr_sizes=args.chrom_sizes, block_size=100)
    plot_nearest(counts, somatic, savefig=False, dirname= args.out + '/' + 'visualization/')
    
if args.variation == True:
    try:
        os.mkdir(args.out + '/variation')
    except:
        pass
    somatic = pd.read_csv(args.breakpoints, header=None, sep = '\t')
    somatic.columns = ['chrA', 'posA', 'orA', 'chrB', 'posB', 'orB', 'haplo', 'reads', 'ref_reads']
    change, mean = find_changes(counts, somatic, log_dir = args.out + '/variation')
    plot(counts, somatic, change=change, mean=mean, savefig=True, dirname = args.out + '/variation') 
    

    
#if args.plot_only == True:
#    try:
#        os.mkdir(args.out + '/visualization')
#    except:
#        pass
#    
#    somatic = pd.read_csv(args.breakpoints, header=None, sep = '\t')
#    somatic.columns = ['chrA', 'posA', 'orA', 'chrB', 'posB', 'orB', 'haplo', 'reads', 'ref_reads']
#    visualize(counts, somatic, args.ylim, savefig=True, filename = str(args.out) + '/' + 'visualization/coverage.png')
  
