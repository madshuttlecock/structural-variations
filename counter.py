import pandas as pd
import pysam
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import sys
import os



def get_counts(filename, chr_sizes=None, reference_length=None, block_size=100000, unphased=False):
    
    if chr_sizes is not None:
        chr_data = pd.read_csv(chr_sizes, sep='\t', header=None)
        chr_data.columns = ['name', 'siz']
        
        chr_dict = dict(zip(chr_data.name, chr_data.siz))
        
    
    samfile = pysam.AlignmentFile(filename, "rb")
    block_array = dict()
    x = 0
    chr_list = []
    
    for read in tqdm(samfile.fetch()):
        start = read.reference_start
        end = read.reference_end
        
            
        ref = read.reference_name
        
        
        if ref not in block_array:
            if chr_sizes is not None:
                try:
                    reference_length = chr_dict[ref]
                except:
                    print(ref, 'E')
                    continue
            block_array[ref] = np.zeros((2, reference_length//block_size+1))
            chr_list += [ref]
        
        try:
            tag = int(read.get_tag('HP'))-1
        except:
            if unphased:
                tag = 0
            else:
                continue
            
        start_block = start // block_size
        end_block = end // block_size

        for i in range(start_block, end_block+1):
            read_block_begin = max(start, i * block_size)
            read_block_end = min(end, min(reference_length, (i+1)*block_size))
            read_block_len = read_block_end - read_block_begin
            
            block_len = min(reference_length, (i+1)*block_size) - i * block_size
            
            block_array[ref][tag][i] += read_block_len / block_len
            
        
    
    samfile.close()
    res = dict()
    for ref in chr_list:
        if chr_sizes is not None:
            reference_length = chr_dict[ref]
        index = pd.DataFrame(np.arange(0, reference_length + block_size-1, block_size), dtype=int)
        
        
        data = pd.DataFrame(pd.concat([index, pd.DataFrame(block_array[ref], dtype=float).T], axis=1))
        
        
        data.columns = ['beginning', 'haploA', 'haploB']
        res[ref] = data
    
    return res



    
def transform(x):
    return x / 1e6

import matplotlib.pyplot as plt
    
def transform(x):
    return x / 1e6

def visualize(coverage, breakpoints, ylimit=500, savefig=False, filename="coverage.png"):

    c = ['red', 'blue', 'darkgreen']
    m = {"+":'>', '-':'<'}

    fig, axs = plt.subplots(len(coverage), squeeze=False)
    fig.set_figheight(len(coverage) * 6)
    fig.set_figwidth(20)

    
    #if len(coverage) == 1:
    #    axs2 = [axs]
    #    axs = axs2
    
    i=0
    
    #print(pd.DataFrame(coverage['chr15']))

    for elem in coverage:
        #print(elem)
        
        axs[i][0].set(title="Coverage by haplotype " + str(elem))
        axs[i][0].plot(transform(coverage[elem].beginning), coverage[elem].haploA, color='blue', alpha=0.7)
        axs[i][0].plot(transform(coverage[elem].beginning), coverage[elem].haploB, color='darkgreen', alpha=0.7)
        
        axs[i][0].set_ylim((-10, ylimit))
        axs[i][0].set_xlabel("Mbp\n")
        axs[i][0].set_ylabel("Coverage")
        
        upper_general = ylimit - 20 
        j = 0
        for mut in breakpoints.values:
            upper = upper_general - j * 20
            
            
            
            alpha =  (mut[7] / (mut[7]+ mut[8])) * 0.7 + 0.3
            
            if mut[0] == elem:
                axs[i][0].scatter([transform(mut[1])], [upper], marker=m[mut[2]], color= c[mut[6]], s=160, alpha=alpha, label=str(mut[6]))
                axs[i][0].text(transform(mut[1]), upper, s=mut[3])
                j += 1

            if mut[3] == elem:
                axs[i][0].scatter([transform(mut[4])], [upper], marker=m[mut[5]], color=c[mut[6]], s=160, alpha=alpha, label=str(mut[6]))
                axs[i][0].text(transform(mut[4]), upper, s=mut[0])
                j += 1
        
        
        i+=1

    if savefig == False:
        plt.show()
    else:
        fig.savefig(filename)
        plt.show()
    
    


def plot_nearest(coverage, breakpoints, savefig=False, dirname='./coverage_zoom_in'):
    def good(pos, ref_pos):
        return np.abs(transform(pos) - transform(ref_pos)) <= 1
    
    c = ['red', 'blue', 'darkgreen']
    m = {"+":'>', '-':'<'}
    
    for mut in breakpoints.values:
        for shift in [0, 3]:
            ch, pos, ori = mut[0 + shift], mut[1 + shift], mut[2 + shift]
            if ch not in coverage:
                continue
            for mut2 in breakpoints.values:
                for shift2 in [0, 3]:
                    ch2, pos2, ori2 = mut2[0 + shift2], mut2[1 + shift2], mut2[2 + shift2]
                    
                    if (mut == mut2).sum() == len(mut) and shift == shift2: 
                        continue
                    
                    if ch == ch2 and abs(pos - pos2) < 100000: # and ori == ori2:
                        
                        alpha =  (mut[7] / (mut[7]+ mut[8])) * 0.7 + 0.3
                        alpha2 =  (mut2[7] / (mut2[7]+ mut2[8])) * 0.7 + 0.3
                        plt.figure(figsize=(15, 4))
                        plt.ylim(0, 500)
                        upper_general = 480 - 20
                        
                        plt.xlim(transform(pos)-1, transform(pos)+1)
                        
                        
                        plt.plot(transform(coverage[ch].beginning), coverage[ch].haploA, color= 'blue')
                        plt.plot(transform(coverage[ch].beginning), coverage[ch].haploB, color = 'darkgreen')
                        
                        
                        plt.scatter([transform(pos)], [30], marker=m[mut[2]], color= c[mut[6]], s=160, alpha=alpha, label=str(mut[6]))
                        plt.text(transform(pos), 30, s=mut[3-shift] + ' ' + str(mut[4-shift]))
                        plt.text(transform(pos), 0, s=str(mut[1+shift]))
                        
                        plt.scatter([transform(pos2)], [100], marker=m[mut2[2]], color= c[mut2[6]], s=160, alpha=alpha2, label=str(mut2[6]))
                        plt.text(transform(pos2), 100, s=mut2[3-shift2] + ' ' + str(mut2[4-shift2]))
                        plt.text(transform(pos), 70, s=str(mut2[1+shift]))
                        
                        if (mut == mut2).sum() == len(mut):
                            plt.text(transform(pos), 150, s="pair")
                        
                        j = 0
                        for mut3 in breakpoints.values:
                            if mut3[0] != ch and mut3[3] != ch:
                                continue
                            if ((mut == mut3).sum() == len(mut)) or ((mut3 == mut2).sum() == len(mut)):
                                continue
                            
                            upper = upper_general - j * 20
                            alpha =  (mut[7] / (mut[7]+ mut[8])) * 0.7 + 0.3

                            if mut3[0] == ch and good(mut3[1], pos):
                                plt.scatter([transform(mut3[1])], [upper], marker=m[mut3[2]], color= c[mut3[6]], s=160, alpha=alpha, label=str(mut3[6]))
                                plt.text(transform(mut3[1]), upper, s=mut[3])
                                j += 1

                            if mut[3] == ch and good(mut3[4], pos):
                                plt.scatter([transform(mut3[4])], [upper], marker=m[mut3[5]], color=c[mut3[6]], s=160, alpha=alpha, label=str(mut3[6]))
                                plt.text(transform(mut3[4]), upper, s=mut3[0])
                                j += 1
                            
                        plt.xlabel('Mbp')
                        plt.title(ch + ' ' + str(pos) + ' +/- 1 Mbp')
                        if savefig == False:
                            plt.show()
                        else:
                            try:
                                os.mkdir(dirname)
                            except:
                                pass
                            plt.savefig(dirname + '/' + ch + ':' + str(pos) + '+-1Mbp.png')
                            plt.show()