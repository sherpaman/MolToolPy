#!/usr/bin/env python

import numpy as np
import ts

#
# SET THE INTERFACE MODULES
#
from argparse import ArgumentParser
parser = ArgumentParser( description = 'Calculate Mutual Information')
#
# INPUT FILES
#
parser.add_argument("-f","--ref",dest="ref",action="store",type=str,default=None,help="Reference Data",required=True,metavar="DAT FILE")
parser.add_argument("-d","--dat",dest="dat",action="store",type=str,default=None,help="Data to be analysed vs. reference",required=True,metavar="DAT FILE")
#
# OUTPUT FILES
#
parser.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=False,help="Output File Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
parser.add_argument("-s","--skip",dest="skip",action="store",type=int,default=1,help="skipping time", metavar="INTEGER")
parser.add_argument("-b","--begin", dest="begin", action="store",type=int,default=1,   help="first column to use",              metavar="INTEGER")
parser.add_argument("-e","--end", dest="end", action="store",type=int,default=0,help="last columb to use",               metavar="INTEGER")
#
#
options = parser.parse_args()

f_ref = options.ref
f_ss  = options.dat
f_out = options.out
skip  = options.skip

ref   = np.loadtxt(f_ref)
ss    = np.loadtxt(f_ss,dtype=int)
nframes, nres = ss.shape
res_a=options.begin
if options.end == 0:
        res_b=nres
else:
        res_b=options.end

time = ref[:,0]
order_base=np.argsort(time)
CHI_skip=[]
ss_skip=[]	
for i in range(0,skip):
	order=order_base[i::skip]
	chi=ref[order,1]
	ss_skip.append(ss[order])
	CHI=ts.TimeSer(chi,n_data=len(chi),dim=1,nbins=12)
	CHI.calc_bins(opt=True)
	CHI_skip.append(CHI)

M = np.zeros((res_b-res_a+1,skip))
f=open(f_out,'w')
for j in np.arange(res_b-res_a+1):
		f.write("{0:4d} ".format(res_a+j))
		for i in range(0,skip):
			SS = ts.TimeSer(ss_skip[i][:,res_a+j-1],n_data=len(ss_skip[i]),dim=1,nbins=8,bins=[np.linspace(0,7,9)],dtype=int)
			M[j,i] , E = CHI_skip[i].mutual_info_other_for(SS)
			f.write("{0:8.3f} ".format(M[j,i]))
		f.write("{0:8.3f} \n".format(np.mean(M[j,:])))

f.close()
quit()

