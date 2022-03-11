#!/usr/bin/env python3
#take corr network and genespace as input and create conditional LLS benchmarking data .condLLS
#Author : Junha Cha
#Date : Apr 24th 2021
#usage : ./script.py [edgelist] [genelist.txt] 
#outfile format : [bin_link] [cumLLS] [TP] [FP] [Gene_Coverage] [binLLS]


import sys
import pandas as pd
import numpy as np
from math import * 

#print(sys.version)
#exit()


'''
def calc_pearsonr_fast(expr_data, absort = True):
    genes = expr_data.index
    pcc = expr_data.transpose().corr(method='pearson')
    net = pcc.corr().where(np.triu(np.ones(pcc.shape), k=1).astype(bool)).stack().reset_index()
    net.columns = ['gene1', 'gene2', 'PCC']
    net_cut = net.loc[(net["PCC"] > 0.6) | (net["PCC"] < -0.8)] #cut meaningless links

    #join two columns of gene
    cols = ['gene1', 'gene2']
    net_cut['genepair'] = list(zip(net_cut.gene1, net_cut.gene2))
    
    #rename network rowname
    net_cut.index = net_cut.genepair
    edges_all = net_cut.loc[:,'PCC'] #series

    
    if absort == False:
        edges_all = pd.Series(edges_all).abs().sort_values(ascending=False)
    
    edges_all = pd.Series(edges_all).sort_values(ascending=False)

    print(edges_all)    
    

    return edges_all
'''

def read_edgelist(file_path):
    net = pd.read_csv(file_path, sep='\t')
    net.columns = ['gene1', 'gene2', 'PCC']

    #join two links as tuple
    cols = ['gene1', 'gene2']
    net['genepair'] = list(zip(net.gene1, net.gene2))

    #rename network rowname
    net.index = net.genepair
    edges_all = net.loc[:,'PCC'] #series

    return edges_all

def load_gold_standard_pairs(gsfile):
    network = dict()
    with open(gsfile) as GS:
        for line in GS:
            a, b = line.strip().split("\t")
            #take gs only in the input genelist
            #if a in genelist and b in genelist: 
            if a not in network:
                network[a] = dict()
            if b not in network:
                network[b] = dict()
            network[a][b] = True
            network[b][a] = True
    return network

def run_benchmark_LLS(target_network, gold_standard_pairs, already_sorted = True, ascending = False, max_pair = 1000000,binsize=1000):

    # calculate positive, negative pairs

    pos = 0
    for g in gold_standard_pairs:
        pos+=len(gold_standard_pairs[g])
    
    pos /= 2
    neg = ( len(gold_standard_pairs) * (len(gold_standard_pairs)-1) / 2 ) - pos

    # calculate pr curve

    TP = 0
    FP = 0
    FN = pos
    binTP = 0
    binFP = 0

    count = 0
    print('prior(GS) ratio')
    print("%f"%(float(pos)/float(neg)))

    randomprecision = float(pos)/float(pos+neg)
    priorprob = float(pos)/float(neg)
    #print(priorprob)
    print('Prior PR')
    print(randomprecision)
    already = set()
    benchdata = dict()
    binstats=list()

    if already_sorted == False:
        sorted_keys = sorted(target_network,key=target_network.get,reverse = not ascending)
    else:
        sorted_keys = target_network.keys()


    for pair in sorted_keys:
        if count>=max_pair:
            break

        already.add(pair[0])
        already.add(pair[1])


        binstats.append(target_network[pair])
        count+=1
        if pair[0] in gold_standard_pairs and pair[1] in gold_standard_pairs[pair[0]]:
            TP += 1
            FN -= 1
            binTP+= 1
        else:
            FP += 1
            binFP += 1

        if count%binsize == 0:
            precision = float(TP)/float(TP+FP)
            recall = float(TP)/float(TP+FN)
            oddratio = precision/randomprecision
            if FP!=0:
                try:
                    lls = log((float(TP)/float(FP)) / priorprob) / log(2)
                except ValueError:
                    print('Error: no annotation!')
                    print(TP,FP,priorprob)
                    exit()
            else:
                lls = 20.0
            if binTP!=0:
                if binFP!=0:
                    binlls = log((float(binTP)/float(binFP)) / priorprob) / log(2)
                else:
                    binlls = 20.0
            else:
                binlls = -5.0
            
            # save result
            benchdata[count] = {"cumLLS":lls,"TP":TP,"FP":FP,"MeanBinStatistics":np.mean(binstats),"GeneCoverage":len(already),"BinLLS":binlls}
            #print("%d\t%.3f\t%d\t%d\t%d\t%f"%(count,lls,TP,FP,len(already),binlls))


            # initialize
            binTP=0
            binFP=0
            binstats=list()


    if count%binsize != 0:
            
        precision = float(TP)/float(TP+FP)
        recall = float(TP)/float(TP+FN)
        oddratio = precision/randomprecision
        lls = log((float(TP)/float(FP)) / priorprob) / log(2)

        if binTP!=0:
            if binFP!=0:
                binlls = log((float(binTP)/float(binFP)) / priorprob) / log(2)
            else:
                binlls = 20.0
        else:
            binlls = -5.0
        benchdata[count] = {"cumLLS":lls,"TP":TP,"FP":FP,"MeanBinStatistics":np.mean(binstats),"GeneCoverage":len(already),"BinLLS":binlls}
        #print("%d\t%.3f\t%d\t%d\t%d\t%f"%(count,lls,TP,FP,len(already),binlls))
    return pd.DataFrame().from_dict(benchdata).T


corrnetpath =  sys.argv[1] #"/home3/junhacha/test/network_test/condLLS_test/sorted_network_file.tsv"
#genelistfile = sys.argv[2] #"/home3/junhacha/test/network_test/condLLS_test/genelist.txt" 
gsdatapath = "/home3/junhacha/public_Data/gold_standard_symbol_HNv3" # [name] [des] [gene1] [g2] [g3] ....
#genelist = list() #this is the genespace critical for conLLS

'''
if (method != 'sort' and method != 'absort'):
    print('wrong sorting method, type sort or absort')
    exit()
'''
'''
#genelist into list
with open(genelistfile,'r') as fp:
    for line in fp:
        if line[0] == "#":
            continue
        genelist.append(line.rstrip())
'''

print('loading sorted PCCnet (cut) data')

edges_all = read_edgelist(corrnetpath) # row : gene, column : sample, no NA values, tab seperated

#either sort PPC absolute or regular
'''
print('make CX from single cell data...')

if method == 'absort':
    edges_all = calc_pearsonr_fast(scoretable)
elif method == 'sort':
    edges_all = calc_pearsonr_fast(scoretable, False) 
'''

#load gs data accounting for input genespace
print('loading gold standard list filtered through genelist...')
gspairs = load_gold_standard_pairs(gsdatapath)

                 
#calcuate benchmark and output
print('calculate LLS...')
benchdata = run_benchmark_LLS(edges_all, gspairs)

print('write output...analysis complete')
benchdata.to_csv(f"{corrnetpath}.binlls", sep ='\t')
