#!/usr/bin/env python

import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.optimize as opt
from scipy.interpolate import interp1d as interp
from sklearn import metrics
import MDAnalysis as MD
import matplotlib.pyplot as plt
from argparse import ArgumentParser
parser = ArgumentParser( description = 'Analyze Binding Sites from Distance Matrix')
#
# INPUT FILES
#
parser.add_argument("-d","--dist",dest="dist",action="store",type=str,default=None,help="Distance Matrix",required=True,metavar="NUMPY")
parser.add_argument("-s","--topol",dest="top",action="store",type=str,default=None,help="Topology File",required=True,metavar="TOPOL")
parser.add_argument("-f","--traj",dest="trj",action="store",type=str,default=None,help="Trajectory",required=True,metavar="TRAJ")
#
# OUTPUT FILES
#
parser.add_argument("-o","--out",dest="out",action="store",type=str,default=None,help="Output base name",required=True,metavar="FILENAME")
#
# OTHER OPTIONS
#
parser.add_argument("-r","--rec",dest="rec",action="store",type=str,default=None,help="Receptor String Selection",required=True,metavar="STRING")
parser.add_argument(     "--list_out",dest="list_out",action="store_true",help="Print List of Residues in Clusters")
parser.add_argument("-t","--threshold",dest="threshold",action="store",type=float,default=None,help="Minimum Probability Threshold",required=False,metavar="PROBABILITY [0.0-1.0]")
parser.add_argument("-c","--clust_cutoff",dest="cluster_cutoff",action="store",type=float,default=None,help="Cutoff for Clustering",required=False,metavar="FLOAT")
parser.add_argument("-m","--clust_min_size",dest="cl_min_sz",action="store",type=int,default=3,help="Minimum size of Cluster to write in the output",required=False,metavar="INTERGER")
parser.add_argument(     "--res0",dest="res0",action="store",type=int,default=1,help="Add this to residue numbering of Protein")
options = parser.parse_args()

def optimal_cutoff(Y,dist_mat,min_size):
    labels = np.array([sch.fcluster(Y,c,criterion='distance') for c in Y[:,2]])
    score = np.array([metrics.silhouette_score(dist_mat,l) for l in labels[:-min_size]])
    c = Y[:-min_size,2]
    f = interp(c,-score,kind='linear')
    opt_c = opt.fmin(f,x0=c[2*min_size])
    return opt_c

dist = options.dist
top  = options.top
trj  = options.trj
base_name = options.out

D = np.load(dist)['arr_0']
u = MD.Universe(top,trj)
receptor = u.select_atoms(options.rec)
r0 = receptor.residues[0].resid
res=np.array([ l.resname+str(l.resid+options.res0) for l in receptor.residues])
res_idx=np.array([ l.resid for l in receptor.residues])
d = D[res_idx][:,res_idx]
idx = np.triu_indices(d.shape[0],1)

fig_size = 20
font_size = 20 * 0.51 * 72 / (2 * len(res))

#Create The First Figure
fig = plt.figure(figsize=(fig_size,fig_size))

# Right side Dedrogram
ax1 = fig.add_axes([0.74,0.1,0.2,0.6])
Y01 = sch.linkage(d[idx], method='complete')
Y02 = sch.linkage(d[idx], method='complete')
Z01 = sch.dendrogram(Y01, orientation='right')
idx1 = Z01['leaves']
ax1.set_xticks([])
ax1.set_yticklabels(res[idx1],size=font_size)

# Top side Dendrogram
ax2 = fig.add_axes([0.09,0.75,0.6,0.2])
Z02 = sch.dendrogram(Y02)
idx2 = Z02['leaves']
#ax2.set_xticks([])
ax2.set_xticklabels(res[idx2],size=font_size)

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
        print("WARNING: Threshold ({0:8.3f}) Corresponds to 0 residues".format(options.threshold))
        quit()
    d_s=d[subs][:,subs]
    res_s=res[subs]
    font_size_s = 20 * 0.51 * 54 / (2 * len(res_s))
    idx = np.triu_indices(d_s.shape[0],1)
    #Create The Second Figure
    fig = plt.figure(figsize=(fig_size,fig_size))
    
    # Right side Dendrogram
    ax1 = fig.add_axes([0.74,0.1,0.2,0.6])
    Y1 = sch.linkage(d_s[idx], method='complete')
    Y2 = sch.linkage(d_s[idx], method='complete')
    Z1 = sch.dendrogram(Y1, orientation='right')
    idx1 = Z1['leaves']
    ax1.set_xticks([])
    ax1.set_yticklabels(res_s[idx1],size=font_size_s)
    
    # Top side Dendrogram
    ax2 = fig.add_axes([0.09,0.75,0.6,0.2])
    Z2 = sch.dendrogram(Y2)
    idx2 = Z2['leaves']
    #ax2.set_xticks([])
    ax2.set_xticklabels(res_s[idx2],size=font_size_s)
    
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

if options.list_out:
    if options.cluster_cutoff == None:
        if options.threshold != None:
            clust_c = optimal_cutoff(Y1,d_s,options.cl_min_sz)
        else:
            clust_c = optimal_cutoff(Y01,d,options.cl_min_sz)
    else:
        clust_c = options.cluster_cutoff
    if options.threshold != None:
        lab_opt = sch.fcluster(Y1,t=clust_c,criterion='distance') # CLUSTERING USING OPTIMIZED OR CHOOSEN CUTOFF
        C = [ subs[np.where(lab_opt==i)[0]] for i in range(1,max(lab_opt)+1) if len(np.where(lab_opt==i)[0])>=options.cl_min_sz ]
    else:
        lab_opt = sch.fcluster(Y01,t=clust_c,criterion='distance') # CLUSTERING USING OPTIMIZED OR CHOOSEN CUTOFF
        C = [ np.where(lab_opt==i)[0] for i in range(1,max(lab_opt)+1) if len(np.where(lab_opt==i)[0])>=options.cl_min_sz ]
    print("Writing  Residue List")
    fo=open(options.out+".resname.dat",'w+')
    fo.write("# CLUSTERING WITH CUTOFF : {0:.3f}\n".format(clust_c))
    for n,i in enumerate(C):
        fo.write("{:30s} {:2d}".format(options.dist,n+1))
        for j in i:
            fo.write(" {:6s}".format(res[j]))
        fo.write("\n")
    fo.close()
    fo=open(options.out+".resnum.dat",'w+')
    fo.write("# CLUSTERING WITH CUTOFF : {0:.3f}\n".format(clust_c))
    for n,i in enumerate(C):
        fo.write("{:30s} {:2d}".format(options.dist,n+1))
        for j in i:
            fo.write(" {:3d}".format(j))
        fo.write("\n")
    fo.close()

quit()
