#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pickle
import phrases
from scipy import stats
from math import *
import MDAnalysis as MD
from argparse import ArgumentParser
parser = ArgumentParser( description = 'Find Binding Sites in Receptor+Ligand MD simulation')
#
# INPUT FILES
#
parser.add_argument("-p","--phrases",dest="phrases",action="store",type=str,default=None,help="Input Phrases list",required=True,metavar="PICKLE FILE")
parser.add_argument("-d","--dist",dest="dist",action="store",type=str,default=None,help="Input Phrases Distance Matrix",required=False,metavar="NPY FILE")
#
# OUTPUT FILES
#
parser.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=False,help="Output File Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
clustering = parser.add_mutually_exclusive_group()
clustering.add_argument("-c","--cutoff",dest="cutoff",action="store",type=float,default=None,help="Cut-off for clustering",metavar="FLOAT")
clustering.add_argument("-l","--clusters",dest="clusters",action="store",type=str,default=None,help="Text file containing the clusters",metavar="CLUST FILE")
parser.add_argument("-t","--threshold",dest="threshold",action="store",type=float,required=False,default=0.25,help="Threshold for Similarity")
#parser.add_argument("-m","--method",dest="method",choices=['jaccard'],default='jaccard',help="Similarity Method")
parser.add_argument("--res0",dest="res0",action="store",type=int,default=1,help="Add this to residue numbering of Protein")
parser.add_argument("--dt",dest="dt",action="store",type=float,default=1.0,help="Time Step of the saved Phrases")


#
options = parser.parse_args()

def read_cluster_file(file_in,r0,nr):
    with open(file_in) as f:
        clusters = []
        labels = np.arange(nr)
        for i,r in enumerate(f.readlines()):
            loc_cluster = []
            for e in r.split():
                resnum = int(e)-int(r0)
                loc_cluster.append(int(e)-r0)
                labels[resnum] = i
            clusters.append(np.array(loc_cluster))
    return labels, clusters

if options.dist != None:
    P = phrases.read_phrases(options.phrases,min_len=3,dist=options.dist)
else:
    P = phrases.read_phrases(options.phrases,min_len=3)
    P.calc_dist()

P.dt = options.dt
p = np.linspace(0,100,1002)
perc = np.percentile(P.D,p)

n_res=len(P.D)

if options.cutoff != None:
    cutoff = options.cutoff
    P.find_cluster(options.cutoff)
else:
    P.labels, P.clusters = read_cluster_file(options.clusters, options.res0,n_res)
P.cluster_phrases(options.threshold)
P.life_time()
P.autocorr_time()

perc_ex=100.*np.sum(P.p_cl>0,axis=0)/float(len(P.p_cl))
life_time=np.zeros(int(max(P.labels)+1))
for i in range(int(max(P.labels)+1)):
    if len(P.LT[i]) > 0:
        life_time[i] = np.average(np.concatenate(P.LT[i]))*P.dt
at = np.sort(P.ac_t[:,:2])*P.dt

if options.out != None:
    with open(options.out+"-binding-site.dat","w") as f:
        if options.cutoff != None:
            f.write("Using Cut-Off                 : {0:10.6f}\n".format(cutoff))
            cutoff_percentile = (p[np.min(np.where(perc>cutoff))]+p[np.max(np.where(perc<cutoff))]) / 2
            f.write("This value corresponds to the   {0:8.4f} percentile\n".format(cutoff_percentile))
        for i in range(1,int(max(P.labels)+1)):
            if len(P.LT[i]) > 0:
                clusters=''
                for e in P.clusters[i].astype(int):
                    clusters=clusters+' {0:3s}'.format(str(e+options.res0))
                f.write("{0:3d}| ({1:80s}) | perc: {2:6.4f} | life-time: {3:8.3f} , {4:8.3f}, {5:8.3f} frames\n".format(i,clusters,perc_ex[i],life_time[i]/1000.0,at[i,0],at[i,1]))
else:
    if options.cutoff != None:
        print("Using Cut-Off                 : {0:10.6f}".format(cutoff))
        cutoff_percentile = (p[np.min(np.where(perc>cutoff))]+p[np.max(np.where(perc<cutoff))]) / 2
        print("This value corresponds to the   {0:8.4f} percentile".format(cutoff_percentile))
    for i in range(1,int(max(P.labels)+1)):
        if len(P.LT[i]) > 0:
            clusters=''
            for e in P.clusters[i].astype(int):
                clusters=clusters+' {0:3s}'.format(str(e+options.res0))
            print("{0:3d}| ({1:80s}) | perc: {2:6.4f} | life-time: {3:8.3f} , {4:8.3f}, {5:8.3f} frames\n".format(i,clusters,perc_ex[i],life_time[i],at[i,0],at[i,1]))


quit()



