#!/usr/bin/env python

import numpy as np
import scipy.cluster.hierarchy as sch
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
parser.add_argument("-r","--rec",dest="rec",action="store",type=str,default=None,help="Receptor String Selection",required=True,metavar="FILENAME")
parser.add_argument("-t","--threshold",dest="threshold",action="store",type=float,default=None,help="Distance Threshold",required=False,metavar="FILENAME")
parser.add_argument("--res0",dest="res0",action="store",type=int,default=1,help="Add this to residue numbering of Protein")
options = parser.parse_args()

dist = options.dist
top  = options.top
trj  = options.trj
base_name = options.out


D = np.load(dist)['arr_0']
u = MD.Universe(top,trj)
receptor = u.select_atoms(options.rec)
r0 = receptor.residues[0].resid
res=np.array([ l.resname+str(l.resid+1) for l in receptor.residues])
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

quit()
