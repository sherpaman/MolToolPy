#!/usr/bin/env python

import ts
import matplotlib.pyplot as plt
import numpy as np
from argparse import ArgumentParser
parser = ArgumentParser( description = 'Calculate Mutual Information')
#
# INPUT FILES
#
parser.add_argument("-i","--inp",dest="inp",action="store",type=str,default=None,help="input Data",required=True,metavar="DAT FILE")
#
# OUTPUT FILES
#
parser.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=True,help="Output File Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
parser.add_argument("-s","--stride",dest="stride",action="store",type=int,default=1,help="time", metavar="INTEGER")
parser.add_argument("-d","--ndim",dest="ndim",action="store",type=int,default=1,help="nuber of dimensions", metavar="INTEGER")
parser.add_argument("-n","--nbins",dest="nbins",action="store",type=int ,default=10,help="number of bins", metavar="INTEGER")
parser.add_argument("-b","--opt",dest="opt",action="store_true",default=False,help="toggle bins optimization")
parser.add_argument("-p","--plot",dest="plot",action="store_true",default=False,help="toggle auto-saving matrix plot")
parser.add_argument("-x","--interleave",dest="interleave",action="store_true",default=False,help="toggle interleaving of data")
#
#
options = parser.parse_args()

def interleave(data,ndim):
	nfr, nrep = data.shape	
	out = np.zeros(data.shape)
	for i in range(nrep/ndim):
                for j in range(ndim):
		        out[:,ndim*i+j]   = data[:,j*(nrep/ndim)+i]
	return out

f_dat = options.inp
f_out = options.out
stride = options.stride

dat   = np.loadtxt(f_dat)
dat   = dat[::stride]

if (options.interleave) & (options.ndim != 1):
        dat = interleave(dat,options.ndim)

DATA= ts.TimeSer(dat,len(dat),dim=options.ndim,nbins=options.nbins)
DATA.calc_bins(opt=options.opt)

M, E = DATA.mutual_info_omp()

fig = plt.figure()
ax  = fig.add_subplot(111)
mat = ax.matshow(M,cmap=plt.get_cmap('RdBu'))
fig.colorbar(mat)
plt.show()
#tx = ax.get_xticks().astype(int)
#ty = ax.get_yticks().astype(int)
#ax.set_xticklabels(ticks[tx])
#ax.set_yticklabels(ticks[ty])
if options.plot:
        fig.savefig(f_out.split('.')[0]+".svg",format='svg')

np.savetxt(f_out,M)

quit()
