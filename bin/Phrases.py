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
subparsers = parser.add_subparsers(help='sub-command help',dest='prog')
#
# TRAJ-READ PROGRAM SUB-PARSER
parser_read = subparsers.add_parser('read-traj',description="read a trajectory, find and cluster binding phrases, basic statistics")
#
# INPUT FILES
#
parser_read.add_argument("-f","--traj",dest="traj",action="store",type=str,default=None,help="Input Trajectory",required=True,metavar="TRAJ FILE")
parser_read.add_argument("-t","--top",dest="top",action="store",type=str,default=None,help="Input Topology",required=True,metavar="TOPOL FILE")
#
# OUTPUT FILES
#
parser_read.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=True,help="Output File Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
parser_read.add_argument("-c","--cutoff",dest="cutoff",action="store",type=float,default=None,help="Cut-off for clustering")
parser_read.add_argument("-b","--begin",dest="begin",action="store",type=int,default=0,help="First frame to read")
parser_read.add_argument("-e","--end",dest="end",action="store",type=int,default=-1,help="Last frame to read")
parser_read.add_argument("-s","--skip",dest="skip",action="store",type=int,default=1,help="number of frame to skip", metavar="INTEGER")
parser_read.add_argument("-j","--jaccard",dest="jac",action="store",default=0.25,help="Jaccard Similarity Threshold")
parser_read.add_argument("-r","--receptor",dest="receptor",action="store",type=str,default="protein",help="Selection string for the Receptor")
parser_read.add_argument("-l","--ligand",dest="ligand",action="store",type=str,default="not protein",help="Selection strin for the Ligand")
parser_read.add_argument("--res0",dest="res0",action="store",type=int,default=1,help="Add this to residue numbering of Protein")
#
# ANALYSIS PROGRAM SUB-PARSER
parser_anal = subparsers.add_parser('anal-phrases',description="read a list of phrases previously generated (possibly a distance matrix between phrases sub-elements), and perform clustering, filtering and basic statistics")
#
# INPUT FILES
#
parser_anal.add_argument("-p","--phrases",dest="phrases",action="store",type=str,default=None,help="Input Phrases list",required=True,metavar="PICKLE FILE")
parser_anal.add_argument("-d","--dist",dest="dist",action="store",type=str,default=None,help="Input Phrases Distance Matrix",required=False,metavar="NPY FILE")
#
# OUTPUT FILES
#
parser_anal.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=False,help="Output File Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
clustering = parser_anal.add_mutually_exclusive_group()
clustering.add_argument("-c","--cutoff",dest="cutoff",action="store",type=float,default=None,help="Cut-off for clustering",metavar="FLOAT")
clustering.add_argument("-l","--clusters",dest="clusters",action="store",type=str,default=None,help="Text file containing the clusters",metavar="CLUST FILE")
parser_anal.add_argument("-t","--threshold",dest="threshold",action="store",type=float,required=False,default=0.25,help="Threshold for Similarity")
parser_anal.add_argument("--dt",dest="dt",action="store",type=float,required=True,default=1.0,help="Time-step of the analyzed trajectory")
#parser.add_argument("-m","--method",dest="method",choices=['jaccard'],default='jaccard',help="Similarity Method")
parser_anal.add_argument("--res0",dest="res0",action="store",type=int,default=1,help="Add this to residue numbering of Protein")
# 
options = parser.parse_args()

# LOCAL FUNCTION DEFINITION
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
            clusters.append(loc_cluster)
    return labels, clusters

if options.prog=='read-traj':
    top             = options.top
    trj             = options.traj
    b               = options.begin
    e               = options.end
    skip            = options.skip
    cutoff          = options.cutoff
    rec_str         = options.receptor
    lig_str         = options.ligand
    res0            = options.res0
    threshold       = options.jac
    lig_dist_cutoff = options.dist_cutoff
    min_phrase_len  = 3
    
    u = MD.Universe(top,trj)
    
    receptor = u.select_atoms(rec_str)
    ligand = u.select_atoms(lig_str)
    
    P = phrases.phrases(u,receptor,ligand,lig_dist_cutoff,min_phrase_len)
    
    P.find_phrases(b,e,skip)
    
    with open(options.out+'-phrases.dat', 'wb') as output:
        pickle.dump(P.phrases, output, pickle.HIGHEST_PROTOCOL)
    
    P.calc_dist()
    np.savez(options.out+"-distance.npz",P.D)
    p = np.linspace(0,100,1002)
    perc = np.percentile(P.D,p)

    if cutoff==None:
        n_val = 551
        c_val = np.linspace(perc[5],perc[60],n_val)
        v_val = np.zeros(n_val)
        
        avg_dist = np.average(P.D)
        
        for i in range(n_val):
            cutoff = c_val[i]
            P.find_cluster(cutoff)
            n_cl = int( np.max(P.labels) + 1 )
            sub_D = P.D[P.centroid[1:]][:,P.centroid[1:]]
            v_val[i] = np.average(sub_D)/avg_dist
        
        buf=20
        
        s1 = np.zeros(n_val)
        s2 = np.zeros(n_val)
        i1 = np.zeros(n_val)
        i2 = np.zeros(n_val)
        for i in range(buf,n_val-buf):
            s1[i], i1[i], r_value1, p_value1, std_err1 = stats.linregress(c_val[buf/2:i],      v_val[buf/2:i]      )
            s2[i], i2[i], r_value2, p_value2, std_err2 = stats.linregress(c_val[i:n_val-buf/2],v_val[i:n_val-buf/2])   
        
        elbow = np.argmax(s1-s2)
        cutoff = c_val[elbow]
        plt.figure()
        plt.plot(c_val,v_val,'.')
        plt.plot(c_val[buf:elbow],i1[elbow]+c_val[buf:elbow]*s1[elbow])
        plt.plot(c_val[elbow:n_val-buf],i2[elbow]+c_val[elbow:n_val-buf]*s2[elbow])
        plt.savefig('{0:s}-elbow-point.png'.format(options.out),fmt="png")

elif options.prog=="anal-phrases":
    threshold = options.threshold
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
    elif options.clusters != None:
        P.labels, P.clusters = read_cluster_file(options.clusters, options.res0,n_res)
    else:
        n_val = 551
        c_val = np.linspace(perc[5],perc[60],n_val)
        v_val = np.zeros(n_val)
        
        avg_dist = np.average(P.D)
        
        for i in range(n_val):
            cutoff = c_val[i]
            P.find_cluster(cutoff)
            n_cl = int( np.max(P.labels) + 1 )
            sub_D = P.D[P.centroid[1:]][:,P.centroid[1:]]
            v_val[i] = np.average(sub_D)/avg_dist
        
        buf=20
        
        s1 = np.zeros(n_val)
        s2 = np.zeros(n_val)
        i1 = np.zeros(n_val)
        i2 = np.zeros(n_val)
        for i in range(buf,n_val-buf):
            s1[i], i1[i], r_value1, p_value1, std_err1 = stats.linregress(c_val[buf/2:i],      v_val[buf/2:i]      )
            s2[i], i2[i], r_value2, p_value2, std_err2 = stats.linregress(c_val[i:n_val-buf/2],v_val[i:n_val-buf/2])   
        
        elbow = np.argmax(s1-s2)
        cutoff = c_val[elbow]
        plt.figure()
        plt.plot(c_val,v_val,'.')
        plt.plot(c_val[buf:elbow],i1[elbow]+c_val[buf:elbow]*s1[elbow])
        plt.plot(c_val[elbow:n_val-buf],i2[elbow]+c_val[elbow:n_val-buf]*s2[elbow])
        plt.savefig('{0:s}-elbow-point.png'.format(options.out),fmt="png")


P.find_cluster(cutoff)
P.cluster_phrases(thresh=threshold)
P.life_time()
P.autocorr_time()

perc_ex=100.*np.sum(P.p_cl>0,axis=0)/float(len(P.p_cl))
life_time=np.zeros(int(max(P.labels)+1))
for i in range(int(max(P.labels)+1)):
    if len(P.LT[i]) > 0:
        life_time[i] = np.average(np.concatenate(P.LT[i]))*P.dt
at = np.sort(P.ac_t[:,:2])*P.dt

#print P.LT
with open(options.out+"-cluster-list.dat","w") as f:
    for i in range(1,int(max(P.labels)+1)):
        clusters=''
        for e in P.clusters[i].astype(int):
            clusters=clusters+' {0:s}'.format(str(e+options.res0))
        f.write("{0:s}\n".format(clusters))

with open(options.out+"-binding-site.dat","w") as f:
    f.write("Using Cut-Off                 : {0:10.6f}\n".format(cutoff))
    cutoff_percentile = (p[np.min(np.where(perc>cutoff))]+p[np.max(np.where(perc<cutoff))]) / 2
    f.write("This value corresponds to the   {0:8.4f} percentile\n".format(cutoff_percentile))
    for i in range(1,int(max(P.labels)+1)):
        if len(P.LT[i]) > 0:
            clusters=''
            for e in P.clusters[i].astype(int):
                clusters=clusters+' {0:3s}'.format(str(e+options.res0))
            f.write("{0:3d}| {1:4d} : ({2:80s}) | perc: {3:6.4f} | life-time: {4:8.3f} , {5:8.3f}, {6:8.3f} ns\n".format(i,P.centroid[i]+1,clusters,perc_ex[i],life_time[i]/1000.,at[i,0]/1000.,at[i,1]/1000.))

quit()
