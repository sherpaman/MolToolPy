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
parser.add_argument("-o","--out",dest="out",action="store",type=str,default=True,required=False,help="Output File Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
parser.add_argument("-s","--skip",dest="skip",action="store",type=int,default=1,help="skipping rows", metavar="INTEGER")
parser.add_argument("-n","--nbins",dest="nbins",action="store",type=int ,default=10,help="number of bins", metavar="INTEGER")
parser.add_argument("-b","--opt",dest="opt",action="store_true",default=False,help="toggle bins optimization")
parser.add_argument("-p","--plot",dest="plot",action="store_true",default=False,help="toggle mutual info matrix plot")
#
#
options = parser.parse_args()
skip  = options.skip

data   = np.loadtxt(options.inp)
#p=np.where(data>120)
#data[p] = data[p] - 360

d_temp = data[::skip]



d = ts.TimeSer(d_temp,n_data=len(d_temp),dim=1,nbins=options.nbins)
d.calc_bins(opt=options.opt)

M , E = d.mutual_info_for()

#ticks = np.array([ 4,5,6,10,11,12,16,17,18,22,23,24,0 ],dtype=int)

fig = plt.figure()
ax  = fig.add_subplot(111)
mat = ax.matshow(M)
fig.colorbar(mat)
#tx = ax.get_xticks().astype(int)
#ty = ax.get_yticks().astype(int)
#ax.set_xticklabels(ticks[tx])
#ax.set_yticklabels(ticks[ty])
if options.plot:
	plt.show()
	if options.out != None:
		fig.savefig(options.out,format='svg')

quit()

