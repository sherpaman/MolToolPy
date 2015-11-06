#!/usr/bin/env python

import ts
import matplotlib.pyplot as plt
import numpy as np
import mdtraj as md
from argparse import ArgumentParser
parser = ArgumentParser( description = 'Calculate Mutual Information')
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
parser.add_argument("-s","--stride",dest="stride",action="store",type=int,default=1,help="time", metavar="INTEGER")
#parser.add_argument("-d","--ndim",dest="ndim",action="store",type=int,default=1,help="nuber of dimensions", metavar="INTEGER")
parser.add_argument("-n","--nbins",dest="nbins",action="store",type=int ,default=10,help="number of bins", metavar="INTEGER")
parser.add_argument("-b","--opt",dest="opt",action="store_true",default=False,help="toggle bins optimization")
parser.add_argument("-p","--plot",dest="plot",action="store_true",default=False,help="toggle auto-saving matrix plot")
#
options = parser.parse_args()

f_traj = options.traj
f_top = options.top
f_out = options.out
stride = options.stride

t = md.load(f_traj,top=f_top,stride=stride)
n_fr = len(t)
Ca = t.top.select('name CA')

aver_str = np.average(t.xyz[:,Ca], axis=0)

dat  = np.swapaxes(t.xyz[:,Ca] - aver_str,1,2)

RMSF = np.sqrt(np.average(np.sum(dat,axis=1)**2,axis=0))
del(t)
DATA= ts.TimeSer(dat,n_fr,dim=3,nbins=options.nbins,reshape=False)
DATA.calc_bins(opt=options.opt)

M, E = DATA.mutual_info_omp()

fig0=plt.figure()
plt.plot(range(1,DATA.rep+1),RMSF)

fig = plt.figure()
ax  = fig.add_subplot(111)
mat = ax.matshow(M)
fig.colorbar(mat)
plt.show()
#tx = ax.get_xticks().astype(int)
#ty = ax.get_yticks().astype(int)
#ax.set_xticklabels(ticks[tx])
#ax.set_yticklabels(ticks[ty])
if options.plot:
        fig.savefig(f_out.split('.')[0]+".svg",format='svg')

np.savetxt(f_out.split('.')[0]+".dat",M)

quit()
