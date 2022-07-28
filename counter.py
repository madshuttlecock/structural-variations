import pandas as pd
import pysam
import numpy as np
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
from sklearn.cluster import \
KMeans, \
AgglomerativeClustering, \
DBSCAN
from scipy.cluster.hierarchy import dendrogram



def get_counts(filename, chr_sizes=None, reference_length=None, block_size=100000, unphased=False, log_dir = '.'):
    try:
        os.mkdir(log_dir)
    except:
        pass
        
    if chr_sizes is not None:
        chr_data = pd.read_csv(chr_sizes, sep='\t', header=None)
        chr_data.columns = ['name', 'siz']
        print(chr_data)
        chr_dict = dict(zip(chr_data.name, chr_data.siz))
        print(chr_dict)
    
    samfile = pysam.AlignmentFile(filename, "rb")
    block_array = dict()# np.zeros((2, reference_length//block_size+1))
    x = 0
    #good = 0
    #good2 = 0
    #summ = 0
    chr_list = []
    
    
    
    for read in tqdm(samfile.fetch()):
        #x += 1
        #if (x == 10):
        #    break
        start = read.reference_start
        end = read.reference_end
        
        #print(j)
        #j += 1
        #if j == 3:
        #    break
            
        ref = read.reference_name
        #print("L", read.reference)
        
        #print(ref)
        
        if ref not in block_array:
            if chr_sizes is not None:
                try:
                    reference_length = chr_dict[ref]
                except:
                    print(ref, 'E')
                    continue
            block_array[ref] = np.zeros((3, reference_length//block_size+1))
            chr_list += [ref]
        
        #print(chr_dict)
        #print(read.query_name, read.reference_start, read.reference_end)
        try:
            tag = int(read.get_tag('HP'))-1
        except:
            tag = 2
            
        #good += 1

        #print("TAG", tag)
        start_block = start // block_size
        end_block = end // block_size

        #s = 0
        for i in range(start_block, end_block+1):
            read_block_begin = max(start, i * block_size)
            read_block_end = min(end, min(reference_length, (i+1)*block_size))
            read_block_len = read_block_end - read_block_begin
            
            block_len = min(reference_length, (i+1)*block_size) - i * block_size
            if block_len == 0:
                print(reference_length, (i+1)*block_size, i * block_size)
            block_array[ref][tag][i] += read_block_len / block_len
            #print(end-start, read_block_len, block_size, block_len)
            #block_array[ref][tag][i][1] += 1
            
        
    
    samfile.close()
    res = dict()
    print(chr_list)
    for ref in chr_list:
        if chr_sizes is not None:
            reference_length = chr_dict[ref]
        index = pd.DataFrame(np.arange(0, reference_length + block_size-1, block_size), dtype=int)
        
        data = pd.DataFrame(pd.concat([index, pd.DataFrame(block_array[ref], dtype=float).T], axis=1))
        
        data.columns = ['beginning', 'haploA', 'haploB', 'unphased']
        res[ref] = data
        
        
        data.to_csv(log_dir + "/" + ref + ".csv")
    
    out = open(log_dir + "/chr_names.txt", "w")
    for elem in chr_list:
        print(elem, file=out)
    return res


def read_counts(log_dir):
    all_cov = dict()
    print(log_dir)
    name_file = open(log_dir + "/chr_names.txt", "r")
    names = name_file.readlines()
    #print(names)
    for name in names:
        name = name.rstrip()
        all_cov[name] = pd.read_csv(log_dir + "/" + name + ".csv")
    return all_cov

def read_counts_raw(log_dir, block_size=0):
    # не учитывает разные возможные названия хромосом  TODO
    all_cov = dict()
    for i in range(1, 23):
        try:
            name = 'chr' + str(i)
            all_cov[name] = pd.read_csv(log_dir + "/" + name + '_' + str(block_size) + "block.csv")
        except:
            try:
                #old format
                name = 'chr' + str(i)
                all_cov[name] = pd.read_csv(log_dir + "/" + name + "_counts.csv")
            except:
                try:
                    name = 'chr' + str(i)
                    all_cov[name] = pd.read_csv(log_dir + "/" + name + ".csv")
                except:
                    pass
    return all_cov

def transform(x):
    return x / 1e6

def visualize(coverage, breakpoints=pd.DataFrame(), ylimit=500, xlimit=None, savefig=False, filename="coverage.png"):

    c = ['red', 'blue', 'darkgreen']
    m = {"+":'>', '-':'<'}

    fig, axs = plt.subplots(len(coverage), squeeze=False)
    fig.set_figheight(len(coverage) * 6)
    fig.set_figwidth(20)

    
    i=0
    
    print(len(coverage))

    for elem in coverage:
        
        axs[i][0].set(title="Coverage by haplotype for chromosome" + str(elem))
        axs[i][0].plot(transform(coverage[elem].beginning), coverage[elem].haploA, color='blue', alpha=0.7, label='haploA')
        axs[i][0].plot(transform(coverage[elem].beginning), coverage[elem].haploB, color='green', alpha=0.7, label= 'haploB')
        try:
            axs[i][0].plot(transform(coverage[elem].beginning), coverage[elem].unphased, color='red', alpha=0.7, label='unphased')
        except:
            pass
        axs[i][0].set_ylim((-10, ylimit))
        axs[i][0].set_xlabel("Mbp\n")
        axs[i][0].set_ylabel("Coverage")
        if xlimit is not None:
            axs[i][0].set_xlim((0, xlimit))
        
        upper_general = ylimit - 20 
        j = 0
        for mut in breakpoints.values:
            upper = upper_general - j * 20
            
            
            
            alpha =  (mut[7] / (mut[7]+ mut[8])) * 0.7 + 0.3
            
            if mut[0] == elem:
                axs[i][0].scatter([transform(mut[1])], [upper], marker=m[mut[2]], color= c[mut[6]], s=160, alpha=alpha)#, label=str(mut[6]))
                axs[i][0].text(transform(mut[1]), upper, s=mut[3])
                j += 1

            if mut[3] == elem:
                axs[i][0].scatter([transform(mut[4])], [upper], marker=m[mut[5]], color=c[mut[6]], s=160, alpha=alpha)#, label=str(mut[6]))
                axs[i][0].text(transform(mut[4]), upper, s=mut[0])
                j += 1
        
        axs[i][0].legend()
        
        i+=1
    
    if savefig == False:
        plt.show()
    else:
        fig.savefig(filename)
        plt.show()

    
import numpy as np
import os

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
                        
                        
                        plt.plot(transform(coverage[ch].beginning), coverage[ch].haploA, color='blue', alpha=0.7, label='haploA')
                        plt.plot(transform(coverage[ch].beginning), coverage[ch].haploB, color='green', alpha=0.7, label= 'haploB')
                        try:
                            plt.plot(transform(coverage[ch].beginning), coverage[ch].unphased, color='red', alpha=0.7, label='unphased')
                        except:
                            pass
                        
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
                            alpha =  (mut3[7] / (mut3[7]+ mut3[8])) * 0.7 + 0.3

                            if mut3[0] == ch and good(mut3[1], pos):
                                plt.scatter([transform(mut3[1])], [upper], marker=m[mut3[2]], color= c[mut3[6]], s=160, alpha=alpha)
                                plt.text(transform(mut3[1]), upper, s=mut3[3])
                                j += 1

                            if mut3[3] == ch and good(mut3[4], pos):
                                plt.scatter([transform(mut3[4])], [upper], marker=m[mut3[5]], color=c[mut3[6]], s=160, alpha=alpha)#, label=str(mut3[6]))
                                plt.text(transform(mut3[4]), upper, s=mut3[0])
                                j += 1
                            
                        plt.xlabel('Mbp')
                        plt.title(ch + ' ' + str(pos) + ' +/- 1 Mbp')
                        plt.legend()
                        if savefig == False:
                            plt.show()
                        else:
                            try:
                                os.mkdir(dirname)
                            except:
                                pass
                            plt.savefig(dirname + '/' + ch + ':' + str(pos) + '+-1Mbp.png')
                            plt.show()
                        
        

        
def extract_chr(data, chr_name):
    first_chr = data[data.chrA==chr_name]
    second_chr = data[data["chrB"]==chr_name]
    dict1 = dict(zip(first_chr.posA, first_chr))
    dict2 = dict(zip(second_chr.posB, second_chr))
    resulting_dict = {**dict1, **dict2}
    
    return pd.concat([first_chr, second_chr], axis=0), resulting_dict

def dump(ans, mean, filename='finder_log.txt'):
    out = open(filename, "w")
    for ch in ans.keys():
        print("CHR", file=out)
        print(ch, file=out)
        for elem in ans[ch]:
            print(elem, file=out)
        print("TRACKS", file=out)
        for track in mean[ch].keys():
            print(track, file=out)
            print(' '.join(map(str, mean[ch][track][0])), file=out)
            print(' '.join(map(str, mean[ch][track][1])), file=out)
    print("END", file=out)
                     
                    
                                
            
def load(filename='finder_log.txt'):
    file = open(filename, 'r')
    s = file.readline()
    #print("A", s)
    ans = dict()
    mean = dict()
    while s.rstrip() != "END":
        s = file.readline()
        #print("B", s)
        ch = s.rstrip()
        ans[ch] = []
        mean[ch] = dict()
        s = file.readline()
        #print("C", s)
        while s.rstrip() != "TRACKS":
            a = int(s.rstrip())
            ans[ch] += [a]
            s = file.readline()
            #print("D", s, s == 'TRACKS')
        s = file.readline()
        while s.rstrip() != "CHR" and s.rstrip() != "END":
            track = s.rstrip()
            #print("T", track)
            mean[ch][track] = dict()
            s = file.readline()
            #print('1', s)
            mean[ch][track][0] = np.array(list(map(float, s.rstrip().split())))
            s = file.readline()
            #print('2', s)
            mean[ch][track][1] = np.array(list(map(float, s.rstrip().split())))
            s = file.readline()
            #print('3', s)
    return ans, mean      
        
        
def find_changes(coverage, breakpoints=None, zoom=100, log_dir = '.'):
    ans = dict()
    mean = dict()
    for chr_name in coverage:
        print(chr_name)
        ans[chr_name] = set()
        mean[chr_name] = dict()
        
        
        for haplo in ["haploA", "haploB", "unphased"]:
            print(haplo)
            result = set()
            try:
                cur_cov = coverage[chr_name][["beginning", haplo]]
            except:
                continue
            cur_cov.columns  = ["beginning", "coverage"]
            block_size = cur_cov.beginning[1] - cur_cov.beginning[0]
            
            if breakpoints is not None:
                pd, mutations = extract_chr(breakpoints, chr_name)
            
            
            buffer = np.zeros(zoom)
            buffer_sum = 0
            buffer_pos = 0
            
            mean[chr_name][haplo] = {1: np.zeros(len(cur_cov))}
            
            for block in range(len(cur_cov)-2, 0, -1):
                buffer_sum -= buffer[buffer_pos]
                buffer[buffer_pos] = cur_cov.coverage[block]
                buffer_sum += buffer[buffer_pos]
                buffer_pos = (buffer_pos + 1) % len(buffer)
                cur_mean = buffer_sum / min(len(cur_cov)-block, len(buffer))
                #print(buffer_sum, cur_mean, min(len(cur_cov)-block, len(buffer)))
                mean[chr_name][haplo][1][block] = cur_mean
            
            
            buffer = np.zeros(zoom)
            buffer_sum = 0
            buffer_pos = 0
            
            mean[chr_name][haplo][0] = np.zeros(len(cur_cov))
            
            for block in range(0, len(cur_cov)-1):
                buffer_sum -= buffer[buffer_pos]
                buffer[buffer_pos] = cur_cov.coverage[block]
                buffer_sum += buffer[buffer_pos]
                buffer_pos = (buffer_pos + 1) % len(buffer)
                
                
                
                cur_mean = buffer_sum / min(block+1, len(buffer))
                mean[chr_name][haplo][0][block] = cur_mean
                
                if abs(mean[chr_name][haplo][1][block+1] - cur_mean) > min(cur_mean, 
                                                                        mean[chr_name][haplo][1][block+1]) * 0.3:
                    i = block
                                 
                    approximate = i * block_size
                    closest = -10000000000
                    if breakpoints is not None:
                        for elem in mutations:
                            if abs(elem - approximate) < abs(closest - approximate):
                                closest = elem
                        if abs(closest - approximate) > block_size*3:
                            continue
                        result.add(closest)
                    else:
                        result.add(approximate)
            ans[chr_name] |= result
    try:
        os.mkdir(log_dir)
    except:
        pass 
    dump(ans, mean, log_dir + '/finder_log.txt') 
    return ans, mean
              
                                
def cluster(chr_positions, threshold=2000000, chr_name=None, dirname='.'):
    agg = AgglomerativeClustering(n_clusters=None, distance_threshold=threshold)
    try:
        clusters = agg.fit_predict(chr_positions)
    except:
        return np.zeros(len(chr_positions))
    #debug
    children = agg.children_
    distance = np.arange(children.shape[0])
    n_dots = np.arange(2, children.shape[0] + 2)
    linkage_matrix = np.column_stack([children, distance, n_dots]).astype(float)
    
    plt.figure(figsize=(40, 10))
    
    names = []
    for i in range(len(chr_positions)):
        names += [str(chr_positions[i]) + ' ' + str(agg.labels_[i])]
    dendrogram(linkage_matrix, labels=names)
    
    plt.title('Дендрограмма для хромосомы' + str(chr_name))
    plt.xlabel('Брейкпоинты');
    try:
        os.mkdir(dirname)
    except:
        pass
    plt.savefig(dirname + '/' + 'dendro_' + chr_name + '.png')
    plt.show()
    return clusters



def plot(coverage, breakpoints=pd.DataFrame(), positions=None, mean=None, load_file=None, savefig=False, dirname='./local'):
    def good(pos, ref_pos_1, ref_pos_2):
        return (ref_pos_1 <= pos) and (pos <= ref_pos_2)
    
    if load is not None:
        positions, mean = load(load_file)
    
    c = ['red', 'blue', 'darkgreen']
    m = {"+":'>', '-':'<'}
    
    for ch in positions:
        
        single_pos = np.array(positions[ch]).reshape(-1, 1)
        #print(single_pos)
        clusters = cluster(single_pos, chr_name=ch, dirname=dirname)
        #return
        pos_by_cluster = dict()

        for (i, el) in enumerate(clusters):
            if el in pos_by_cluster:
                pos_by_cluster[el] += [single_pos[i]]
            else:
                pos_by_cluster[el] = [single_pos[i]]
        
        for cl in pos_by_cluster:
            
            curmin = int(np.min(pos_by_cluster[cl]))
            curmax = int(np.max(pos_by_cluster[cl]))
            delta = (curmax-curmin) // 2
            #print(curmin, delta)
            if delta == 0:
                delta = 1e6
            curmin -= int(delta)
            curmax += int(delta)
            #print(curmin)
            
            plt.figure(figsize=(15, 10))
            
            plt.plot(transform(coverage[ch].beginning), coverage[ch].haploA, color= 'blue')
            plt.plot(transform(coverage[ch].beginning), coverage[ch].haploB, color = 'green')
            plt.plot(transform(coverage[ch].beginning), coverage[ch].unphased, color = 'red')


            if mean is not None:
                plt.plot(transform(coverage[ch].beginning), mean[ch]["haploA"][0], color= 'lightblue')
                plt.plot(transform(coverage[ch].beginning), mean[ch]["haploB"][0], color = 'lightgreen')
                plt.plot(transform(coverage[ch].beginning), mean[ch]["haploA"][1], color= 'darkblue')
                plt.plot(transform(coverage[ch].beginning), mean[ch]["haploB"][1], color = 'darkgreen')
                plt.plot(transform(coverage[ch].beginning), mean[ch]["unphased"][1], color= 'pink')
                plt.plot(transform(coverage[ch].beginning), mean[ch]["unphased"][1], color = 'darkred')

            for pos in pos_by_cluster[cl]:
                #print(cl, pos, curmin, curmax)
                plt.vlines(transform(pos), 0, 20, linewidth=5, color='red')
                
            j = 0
            for mut3 in breakpoints.values:
                if mut3[0] != ch and mut3[3] != ch:
                    continue

                upper_general = 450

                upper = upper_general - j * 20
                alpha =  (mut3[7] / (mut3[7]+ mut3[8])) * 0.7 + 0.3
                f = 0
                
                if mut3[0] == ch and good(mut3[1], curmin, curmax):
                    plt.scatter([transform(mut3[1])], [upper], marker=m[mut3[2]], color= c[mut3[6]], s=160, alpha=alpha, label=str(mut3[6]))
                    plt.text(transform(mut3[1]), upper, s=mut3[3]+' '+str(mut3[4]))
                    f = 1

                if mut3[3] == ch and good(mut3[4], curmin, curmax):
                    plt.scatter([transform(mut3[4])], [upper], marker=m[mut3[5]], color=c[mut3[6]], s=160, alpha=alpha, label=str(mut3[6]))
                    plt.text(transform(mut3[4]), upper, s=mut3[0]+' ' +str(mut3[1]))
                    f = 1
                    
                j += f

            plt.xlabel('Mbp')
            plt.ylim(0, 500)
            plt.xlim(transform(curmin), transform(curmax))
            plt.title(ch + ':' + str(transform(curmin)) + '-' + str(transform(curmax)) + 'Mbp')
            if savefig == False:
                plt.show()
            else:
                try:
                    os.mkdir(dirname)
                except:
                    pass
                plt.savefig(dirname + '/' + ch + ':' + str(transform(curmin)) + '-' + str(transform(curmax)) + 'Mbp' + '.png')
                plt.show()

        
        
