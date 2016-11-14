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
parser.add_argument("-f","--traj",dest="traj",action="store",type=str,default=None,help="Input Trajectory",required=True,metavar="TRAJ FILE")
parser.add_argument("-t","--top",dest="top",action="store",type=str,default=None,help="Input Topology",required=True,metavar="TOPOL FILE")
#
# OUTPUT FILES
#
parser.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=True,help="Output File Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
parser.add_argument("-c","--cutoff",dest="cutoff",action="store",type=float,default=0.0,help="Cut-off for clustering")
parser.add_argument("-b","--begin",dest="begin",action="store",type=int,default=0,help="First frame to read")
parser.add_argument("-e","--end",dest="end",action="store",type=int,default=-1,help="Last frame to read")
parser.add_argument("-s","--skip",dest="skip",action="store",type=int,default=1,help="number of frame to skip", metavar="INTEGER")
#parser.add_argument("-p","--plot",dest="plot",action="store_true",default=False,help="toggle auto-saving matrix plot")
parser.add_argument("-r","--receptor",dest="receptor",action="store",type=str,default="protein",help="Selection string for the Receptor")
parser.add_argument("-l","--ligand",dest="ligand",action="store",type=str,default="not protein",help="Selection strin for the Ligand")
#
options = parser.parse_args()


top = options.top
trj = options.traj
b    = options.begin
e    = options.end
skip = options.skip
cutoff = options.cutoff
rec_str = options.receptor
lig_str = options.ligand

min_phrase_len = 3
threshold = 0.25
lig_dist_cutoff = 6.0

u = MD.Universe(top,trj)

receptor = u.select_atoms(rec_str)
ligand = u.select_atoms(lig_str)

P = phrases.phrases(u,receptor,ligand,lig_dist_cutoff,min_phrase_len)

P.find_phrases(b,e,skip)


with open(options.out+'-phrases.dat', 'wb') as output:
    pickle.dump(P, output, pickle.HIGHEST_PROTOCOL)

P.calc_dist()
np.savez(options.out+"-distance.npz",P.D)

perc = np.percentile(P.D,np.linspace(0,100,1002))

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
    
    
    for i in range(buf,n_val-buf):
        slope1[i]  intercept1, r_value1, p_value1, std_err1 = stats.linregress(c_val[buf/2:i],      v_val[buf/2:i]      )
        slope2[i], intercept2, r_value2, p_value2, std_err2 = stats.linregress(c_val[i:n_val-buf/2],v_val[i:n_val-buf/2])   
    
    elbow = np.max(slope1-slope2)
    cutoff = c_val[elbow]

cutoff_percentile = (perc[np.min(np.where(perc>cutoff))]+perc[np.max(np.where(perc<cutoff))]) * 50

print("Using Cut-Off                 : {0:10.6f}".format(cutoff))
print("This value corresponds to the   {0:8.4f} percentile".format(cutoff_percentile))

P.find_cluster(cutoff)
P.cluster_phrases()
P.life_time()

perc=np.sum(P.p_cl>0,axis=0)/float(len(P.p_cl))
life_time=np.zeros(int(max(P.labels)+1))

#print P.LT
with open(options.out+"-binding-site.dat","w") as f:
    for i in range(1,int(max(P.labels)+1)):
        life_time[i] = np.average(P.LT[i])*P.dt
        clusters=''
        for e in P.clusters[i]:
            clusters=clusters+' {0:3s}'.format(e)
        f.write("{0:3d}|{1:4d} : ({2:s}) | perc: {3:6.4f} | life-time: {4:10.6f} ns\n".format(i,P.centroid[i],clusters,perc[i],life_time[i]/1000.0))


quit()
