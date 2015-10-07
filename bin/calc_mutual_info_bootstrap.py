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
parser.add_argument("-r","--resample",dest="resample",action="store",type=int,default=None,help="Bootstrap Resampling", metavar="INTEGER")
parser.add_argument("-s","--skip",dest="skip",action="store",type=int,default=1,help="skipping rows", metavar="INTEGER")
parser.add_argument("-n","--nbins",dest="nbins",action="store",type=int ,default=10,help="number of bins", metavar="INTEGER")
parser.add_argument("-b","--opt",dest="opt",action="store_true",default=False,help="toggle bins optimization")
parser.add_argument("-p","--plot",dest="plot",action="store_true",default=False,help="toggle mutual info matrix plot")
#
#
options = parser.parse_args()

f_in  = options.inp 
f_out = options.out
skip  = options.skip
f_out_mi = f_out.split('.')[0]+'_MI.svg'
f_out_dev = f_out.split('.')[0]+'_DEV.svg'

data   = np.loadtxt(f_in)
p=np.where(data>120)
data[p] = data[p] - 360

d_temp = data[::skip]

d = ts.TimeSer(d_temp,n_data=len(d_temp),dim=1,nbins=options.nbins)
d.calc_bins(opt=options.opt)

M , E = d.mutual_info_for()
M_B = d.mutual_info_bootstrap(resample=options.resample)

MI_aver=np.average(M_B,axis=2)
MI_var=np.var(M_B,axis=2)

DEV=np.abs(MI_aver - M)/np.sqrt(MI_var)

ticks = np.array([ 4,5,6,10,11,12,16,17,18,22,23,24,0 ],dtype=int)

fig1 = plt.figure()
ax1  = fig1.add_subplot(111)
mat1 = ax1.matshow(M)
fig1.colorbar(mat1)
tx1 = ax1.get_xticks().astype(int)
ty1 = ax1.get_yticks().astype(int)
ax1.set_xticklabels(ticks[tx1])
ax1.set_yticklabels(ticks[ty1])

fig2 = plt.figure()
ax2  = fig2.add_subplot(111)
mat2 = ax2.matshow(DEV)
fig2.colorbar(mat2)
tx2 = ax2.get_xticks().astype(int)
ty2 = ax2.get_yticks().astype(int)
ax2.set_xticklabels(ticks[tx2])
ax2.set_yticklabels(ticks[ty2])

if options.plot:
	plt.show()

fig1.savefig(f_out_mi,format='svg')
fig2.savefig(f_out_dev,format='svg')

