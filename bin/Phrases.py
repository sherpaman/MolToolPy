#!/usr/bin/env python

from __future__ import print_function
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
parser_read.add_argument("-d","--dist",dest="dist_cutoff",action="store",type=float,default=7.0,help="Distance Cut-off for binding definition")
parser_read.add_argument("-c","--cutoff",dest="cutoff",action="store",type=float,default=None,help="Cut-off for clustering")
parser_read.add_argument("-b","--begin",dest="begin",action="store",type=int,default=0,help="First frame to read")
parser_read.add_argument("-e","--end",dest="end",action="store",type=int,default=-1,help="Last frame to read")
parser_read.add_argument("-s","--skip",dest="skip",action="store",type=int,default=1,help="number of frame to skip", metavar="INTEGER")
parser_read.add_argument("-j","--jaccard",dest="jac",action="store",type=float,default=0.25,help="Jaccard Similarity Threshold")
parser_read.add_argument("-r","--receptor",dest="receptor",action="store",type=str,default="protein",help="Selection string for the Receptor")
parser_read.add_argument("-l","--ligand",dest="ligand",action="store",type=str,default="not protein",help="Selection strin for the Ligand")
parser_read.add_argument("--res0",dest="res0",action="store",type=int,default=1,help="Add this to residue numbering of Protein")
#
# JACCARD ANALYSIS PROGRAM SUB-PARSER
parser_anal = subparsers.add_parser('Jaccard',description="read a list of phrases previously generated (possibly a distance matrix between phrases sub-elements), and perform clustering, filtering and basic statistics")
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
parser_anal.add_argument("-j","--jaccard",dest="jac",action="store",type=float,required=False,default=0.25,help="Threshold for Similarity")
parser_anal.add_argument("--dt",dest="dt",action="store",type=float,required=True,default=1.0,help="Time-step of the analyzed trajectory")
parser_anal.add_argument("--res0",dest="res0",action="store",type=int,default=1,help="Add this to residue numbering of Protein")
# 
#
# LINKAGE ANALYSIS PROGRAM SUB-PARSER
parser_linkage = subparsers.add_parser('Linkage',description="Performe Hierarchical Clustering of Residues based on Phrases Analysis and/or filter the Phrases based on clustering")
#
# INPUT FILES
#
parser_linkage.add_argument("-p","--phrases",dest="phrases",action="store",type=str,default=None,help="Input Phrases list",required=True,metavar="PICKLE FILE")
parser_linkage.add_argument("-d","--dist",dest="dist",action="store",type=str,default=None,help="Input Phrases Distance Matrix",required=False,metavar="NPY FILE")
parser_linkage.add_argument("-s","--topol",dest="top",action="store",type=str,default=None,help="Topology File",required=True,metavar="TOPOL")
parser_linkage.add_argument("-f","--ref",dest="trj",action="store",type=str,default=None,help="Structure File",required=True,metavar="STRUCT")
#
# OUTPUT FILES
#
parser_linkage.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=True,help="Output File Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
parser_linkage.add_argument("-r","--receptor",dest="receptor",action="store",type=str,default=None,help="Selection string for the Receptor",required=True,metavar="STRING")
parser_linkage.add_argument("-t","--threshold",dest="threshold",action="store",type=float,default=None,help="Minimum Probability Threshold for Clustering",required=False,metavar="PROBABILITY")
parser_linkage.add_argument(     "--do_filter",dest="do_filter",action="store_true",default=False,help="Toggle the Cluster Based Filtering of the original Phrases",required=False,metavar="FLOAT")
#
clustering2 = parser_linkage.add_mutually_exclusive_group()
clustering2.add_argument("-e","--do_elbow",dest="do_elbow",action="store_true",default=False,help="Toggle The Elbow Criterion for Cluster Definition")
clustering2.add_argument("-l","--clusters",dest="clusters",action="store",type=str,default=None,help="Text file containing the clusters",metavar="CLUST FILE")
#
#
options = parser.parse_args()

# LOCAL FUNCTION DEFINITION
def read_cluster_file(file_in,r0,nr):
    with open(file_in) as f:
        clusters = [[]]
        labels = np.zeros(nr)
        for i,r in enumerate(f.readlines()):
            loc_cluster = []
            for e in r.split():
                resnum = int(e)-int(r0)
                loc_cluster.append(int(e)-r0)
                labels[resnum] = i+1
            clusters.append(np.array(loc_cluster))
            clusters[0] = [ n for n,i in enumerate(labels) if i == 0 ]
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

    if cutoff != None:
        P.find_cluster(options.cutoff)
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

elif options.prog=="anal-phrases":
    res0            = options.res0
    threshold       = options.jac
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
        P.labels, P.clusters = read_cluster_file(options.clusters,options.res0,n_res)
        P.centroid = [ i[0] for i in P.clusters ]
        print("Will measure Clusters :")
        for i in P.clusters:
            print(i)
        print(" Label:" )
        print(P.labels)
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

elif options.program == "Linkage":
    top      = options.top
    trj      = options.trj
    rec_str  = options.receptor
    base_name = options.out
    
    u = MD.Universe(top,trj)
    receptor = u.select_atoms(rec_str)
    P = phrases.read_phrases(options.phrases,min_len=3)
    if options.dist == None:
        P.calc_dist()
        D = P.D
    else:
        D = np.load(dist)['arr_0']
    
    r0 = receptor.residues[0].resid
    res=np.array([ l.resname+str(l.resid+options.res0) for l in receptor.residues])
    res_idx=np.array([ l.resid for l in receptor.residues])
    d = D[res_idx][:,res_idx]
    idx = np.triu_indices(d.shape[0],1)
    
    #Create The First Figure
    fig = plt.figure(figsize=(10,10))
    
    # Right side Dedrogram
    ax1 = fig.add_axes([0.74,0.1,0.2,0.6])
    Y1 = sch.linkage(d[idx], method='complete')
    Y2 = sch.linkage(d[idx], method='complete')
    Z1 = sch.dendrogram(Y1, orientation='right')
    idx1 = Z1['leaves']
    ax1.set_xticks([])
    ax1.set_yticklabels(res[idx1],size=4)
    
    # Top side Dendrogram
    ax2 = fig.add_axes([0.09,0.75,0.6,0.2])
    Z2 = sch.dendrogram(Y2)
    idx2 = Z2['leaves']
    #ax2.set_xticks([])
    ax2.set_xticklabels(res[idx2],size=4)
    
    # Reordered Large Distance Matrix
    axmatrix = fig.add_axes([0.09,0.1,0.6,0.6])
    im = axmatrix.matshow(d[idx1][:,idx2], aspect='auto', origin='lower', cmap=plt.get_cmap('coolwarm'))
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    
    # Original Distance Matrix (Small)
    axmatrix2 = fig.add_axes([0.74,0.74,0.18,0.18])
    im = axmatrix2.matshow(d, aspect='auto', origin='lower', cmap=plt.get_cmap('coolwarm'))
    #axmatrix.set_xticks([])
    #axmatrix.set_yticks([])
    
    # Colorbar
    cbaxes = fig.add_axes([0.95,0.1, 0.01, 0.6])
    plt.colorbar(im, cax=cbaxes)
    plt.savefig('{0:s}_complete.pdf'.format(base_name),fmt='pdf')
    
    if options.threshold != None:
        if (options.threshold <= 0.0) | (options.threshold >= 1.0):
            print("Threshold ({0:8.3f}) must be a float number between 0.0 and 1.0 (extreemes exluded)".format(options.threshold))
            raise ValueError
        subs=np.where(np.diag(d)< -np.log(options.threshold))[0]
        if len(subs) < 1:
            print("Threshold ({0:8.3f}) Corresponds to 0 residues".format(options.threshold))
            quit()
        d_s=d[subs][:,subs]
        res_s=res[subs]
        idx = np.triu_indices(d_s.shape[0],1)
        #Create The Second Figure
        fig = plt.figure(figsize=(10,10))
        
        # Right side Dendrogram
        ax1 = fig.add_axes([0.74,0.1,0.2,0.6])
        Y1 = sch.linkage(d_s[idx], method='complete')
        Y2 = sch.linkage(d_s[idx], method='complete')
        Z1 = sch.dendrogram(Y1, orientation='right')
        idx1 = Z1['leaves']
        ax1.set_xticks([])
        ax1.set_yticklabels(res_s[idx1],size=7)
        
        # Top side Dendrogram
        ax2 = fig.add_axes([0.09,0.75,0.6,0.2])
        Z2 = sch.dendrogram(Y2)
        idx2 = Z2['leaves']
        #ax2.set_xticks([])
        ax2.set_xticklabels(res_s[idx2],size=7)
        
        # Reordered Large Distance Matrix
        axmatrix = fig.add_axes([0.09,0.1,0.6,0.6])
        im = axmatrix.matshow(d_s[idx1][:,idx2], aspect='auto', origin='lower', cmap=plt.get_cmap('coolwarm'))
        im.set_clim(np.min(d),np.max(d))
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])
        
        # Original Distance Matrix (Small)
        axmatrix2 = fig.add_axes([0.74,0.74,0.18,0.18])
        im = axmatrix2.matshow(d, aspect='auto', origin='lower', cmap=plt.get_cmap('coolwarm'))
        #axmatrix.set_xticks([])
        #axmatrix.set_yticks([])
        
        # Colorbar
        cbaxes = fig.add_axes([0.95,0.1, 0.01, 0.6])
        plt.colorbar(im, cax=cbaxes)
        #plt.show()
        plt.savefig('{0:s}_{1:s}_subset.pdf'.format(base_name,str(options.threshold)),fmt='pdf')
        
        
        
        exit()

# CALCULATE LIFE TIME OF CLUSTERS
P.cluster_phrases(thresh=threshold)
P.life_time()
P.autocorr_time()

perc_ex=100.*np.sum(P.p_cl>0,axis=0)/float(len(P.p_cl))
life_time=np.zeros(int(max(P.labels)+1))
for i in range(int(max(P.labels)+1)):
    if len(P.LT[i]) > 0:
        life_time[i] = np.average(np.concatenate(P.LT[i]))*P.dt
at = np.sort(P.ac_t[:,:2])*P.dt

do_write_clusters=False
if options.prog=="read-traj" :
    do_write_clusters=True
if options.prog == "anal-phrases": 
    if options.clusters != None :
        do_write_clusters=True
        
if do_write_clusters:
    with open(options.out+"-cluster-list.dat","w") as f:
        for i in range(1,int(max(P.labels))+1):
            clusters=''
            for e in P.clusters[i].astype(int):
                clusters=clusters+' {0:s}'.format(str(e+options.res0))
            f.write("{0:s}\n".format(clusters))

        
#print P.LT
with open(options.out+"-binding-site.dat","w") as f:
    if options.prog=="read-traj":
        f.write("Using Cut-Off                 : {0:10.6f}\n".format(cutoff))
        try:
            cutoff_percentile = (p[np.min(np.where(perc>cutoff))]+p[np.max(np.where(perc<cutoff))]) / 2
            f.write("This value corresponds to the   {0:8.4f} percentile\n".format(cutoff_percentile))
        except:
            f.write("This value corresponds to the    0.00 percentile\n")
    else:
        if options.cutoff != None:
            f.write("Using Cut-Off                 : {0:10.6f}\n".format(cutoff))
            try:
                cutoff_percentile = (p[np.min(np.where(perc>cutoff))]+p[np.max(np.where(perc<cutoff))]) / 2
                f.write("This value corresponds to the   {0:8.4f} percentile\n".format(cutoff_percentile))
            except:
                f.write("This value corresponds to the    0.00 percentile\n".format(cutoff_percentile))
        else:
            f.write("Using Manually defined Clusters\n")
    f.write("Using Jaccard Similarity Threshold : {0:10.6f}\n".format(threshold))
    for i in range(1,int(max(P.labels)+1)):
        if len(P.LT[i]) > 0:
            clusters=''
            for e in P.clusters[i].astype(int):
                clusters=clusters+' {0:3s}'.format(str(e+options.res0))
            f.write("{0:3d}| {1:4d} : ({2:80s}) | perc: {3:6.4f} | life-time: {4:8.3f} , {5:8.3f}, {6:8.3f} ns\n".format(i,P.centroid[i]+options.res0,clusters,perc_ex[i],life_time[i]/1000.,at[i,0]/1000.,at[i,1]/1000.))

quit()
