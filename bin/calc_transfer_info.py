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
parser.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=False,help="Output File Name",metavar="DAT FILE")
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

f_in  = options.inp 
f_out = options.out
skip  = options.skip

data   = np.loadtxt(f_in)
d_temp = data[::skip]

d = ts.TimeSer(d_temp,n_data=len(d_temp),dim=1,nbins=options.nbins)
d.calc_bins(opt=options.opt)

T,D = d.transfer_entropy_for(time=3)

#ticks = np.array([4,5,6,10,11,12,16,17,18,22,23,24,0],dtype=int)

if options.plot:
	fig = plt.figure()
	ax  = fig.add_subplot(111)
	mat = ax.matshow(D)
	fig.colorbar(mat)
#	tx = ax.get_xticks().astype(int)
#	ty = ax.get_yticks().astype(int)
#	ax.set_xticklabels(ticks[tx])
#	ax.set_yticklabels(ticks[ty])
	plt.show()

quit()

