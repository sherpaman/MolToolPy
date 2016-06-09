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
parser.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=True,help="Output File Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
parser.add_argument("-s","--stride",dest="stride",action="store",type=int,default=1,help="time", metavar="INTEGER")
parser.add_argument("-b","--opt",dest="opt",action="store_true",default=False,help="toggle bins optimization")
parser.add_argument("-p","--plot",dest="plot",action="store_true",default=False,help="toggle auto-saving matrix plot")
#
parser.add_argument("--x_ref",dest="x_ref",action="store_true",default=False,help="toggle interleaving of reference data")
parser.add_argument("--x_tar",dest="x_tar",action="store_true",default=False,help="toggle interleaving of target data")
parser.add_argument("--ndim1",dest="ndim1",action="store",type=int,default=1,help="nuber of dimensions of ref", metavar="INTEGER")
parser.add_argument("--ndim2",dest="ndim2",action="store",type=int ,default=1,help="number of dimensions of target", metavar="INTEGER")
parser.add_argument("--nbins1",dest="nbins1",action="store",type=int,default=10,help="nuber of bins of ref", metavar="INTEGER")
parser.add_argument("--nbins2",dest="nbins2",action="store",type=int,default=10,help="nuber of bins of target", metavar="INTEGER")
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

f_ref = options.ref
f_tar = options.dat
f_out = options.out

stride = options.stride
ndim1  = options.ndim1
ndim2  = options.ndim2
nbins1 = options.nbins1
nbins2 = options.nbins2

ref   = np.loadtxt(f_ref)
tar   = np.loadtxt(f_tar)

if (options.x_ref) & (options.ndim1 != 1):
    ref = interleave(ref,options.ndim1)

if (options.x_tar) & (options.ndim2 != 1):
    tar = interleave(tar,options.ndim2)

REF=[]
TAR=[]

try:
    nfr, ref_rep = ref.shape
except:
    nfr=len(ref)
    ref_rep=1
nfr , tar_rep = tar.shape

M = np.zeros( ( ref_rep/ndim1, tar_rep/ndim2 , stride ) )

for i in range(stride):
    t_ref = ref[i::stride]
    t_tar = tar[i::stride]
    REF=ts.TimeSer(t_ref,len(t_ref),dim=ndim1,nbins=nbins1)
    TAR=ts.TimeSer(t_tar,len(t_tar),dim=ndim2,nbins=nbins2)
    M[:,:,i], e = REF.mutual_info_other_omp(TAR)

if (ref_rep/ndim1 == 1)&(tar_rep/ndim2 == 1):
    out = open(f_out,'w')
    out.write("{0:4d} {1:10.6f} {1:10.6f}\n".format(i,np.average(M[0,0]), np.sqrt(np.var(M[0,0]))))
    out.close()
elif (ref_rep/ndim1 == 1)|(tar_rep/ndim2 == 1):
    out = open(f_out,'w')
    if (ref_rep/ndim1 == 1):
        for i in range(tar_rep/ndim2):
            out.write("{0:4d} {1:10.6f} {1:10.6f}\n".format(i,np.average(M[0,i]), np.sqrt(np.var(M[0,i]))))
    elif (tar_rep/ndim2 == 1):
        for i in range(ref_rep/ndim1):
            out.write("{0:4d} {1:10.6f} {1:10.6f}\n".format(i,np.average(M[i,0]), np.sqrt(np.var(M[i,0]))))
    out.close()
else:
    np.savetxt(f_out.split('.')[0]+"MI_AVERAGE.dat",np.average(M,axis=2))
    np.savetxt(f_out.split('.')[0]+"MI_STDEV.dat",np.sqrt(np.var(M,axis=2)))
quit()

