#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
parser = ArgumentParser( description = 'Compact Binding Site Analysis')
#
# INPUT FILES
#
parser.add_argument("-d","--dist",dest="dist",action="store",type=str,help="Input Phrases Distance Matrix",required=True,metavar="NPY FILE")
#
# OUTPUT FILES
#
#parser.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=False,help="Output File Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
parser.add_argument("-c","--cutoff",dest="cutoff",action="store",type=float,default=None,help="Cut-off for clustering",metavar="FLOAT")
#parser.add_argument("-t","--threshold",dest="threshold",action="store",type=float,required=False,default=0.25,help="Threshold for Similarity")
#parser.add_argument("-m","--method",dest="method",choices=['jaccard'],default='jaccard',help="Similarity Method")
parser.add_argument("--res0",dest="res0",action="store",type=int,default=1,help="Add this to residue numbering of Protein")
#parser.add_argument("--dt",dest="dt",action="store",type=float,default=1.0,help="Time Step of the saved Phrases")

options = parser.parse_args()

file_in = options.dist
cutoff  = options.cutoff
D = np.load(file_in)['arr_0']
D_0 = (D < cutoff).astype(int)

def similar(D,row):
    nr = len(D)
    val=np.zeros(nr)
    for i in np.arange(nr):
        val[i] = np.linalg.norm(D[row]-D[i])
    return val

def compact(D,idx):
    D_temp = np.copy(D)
    nr = len(D_temp)
    out_idx = np.copy(idx)
    for i in np.arange(nr):
        val = similar(D_temp,row=i)
        ord_i = np.concatenate([np.arange(i+1),np.argsort(val[i+1:])+i+1])
        out_idx = out_idx[ord_i]
        D_new =np.copy(D_temp)
        D_new = D_temp[ord_i][:,ord_i]
        D_temp = np.copy(D_new)
    return D_temp, out_idx

def calc_sim(D):
    S = np.zeros(D.shape)
    nr,nc = D.shape
    for i in np.arange(nr):
        for j in np.arange(i+1,nc):
            S[i,j] = np.linalg.norm(D[i]-D[j])
            S[j,i] = S[i,j]
    return S

def calc_sim_2(D):
    S = np.zeros(D.shape)
    nr,nc = D.shape
    for i in np.arange(nr):
        for j in np.arange(i+1,nc):
            S[i,j] = np.linalg.norm(D[i]-D[j])
            S[j,i] = S[i,j]
    return S

nr,nc = D_0.shape
idx = np.arange(nr)

D_ord = np.copy(D_0)

#nneig = np.sum(D_0,axis=1)
#ord_nneig = np.argsort(nneig)[::-1]
#idx = idx[ord_nneig]
#D_ord = D_ord[nneig][:,nneig]

#
# Calculate Similarities between all residues neigbouring pattern
#
S = calc_sim(D_ord)

#
# Order according to average similarity
#

av_sim = np.mean(S,axis=0)
ord_av_sim = np.argsort(av_sim)
idx = idx[ord_av_sim]

D_ord = D_ord[idx][:,idx]

#
# Compact by re-ordening according to similarity
#

D_final, idx_final = compact(D_ord,idx)

# remove residues without any neigbours

slice = [ i for i in np.arange(len(D_final)) if not np.all(D_final[i] == 0) ]
idx = idx_final[slice]
D_out = D_final[slice][:,slice]

out_plt = plt.matshow(D_out)
out_plt.axes.set_xticks(np.arange(len(slice)))
out_plt.axes.set_yticks(np.arange(len(slice)))
out_plt.axes.set_xticklabels(idx+options.res0)
out_plt.axes.set_yticklabels(idx+options.res0)
plt.grid()
plt.show()






