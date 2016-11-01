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
parser.add_argument("-n","--nbins",dest="nbins",action="store",type=int ,default=10,help="number of bins", metavar="INTEGER")
parser.add_argument("-b","--opt",dest="opt",action="store_true",default=False,help="toggle bins optimization")
#
options = parser.parse_args()

f_traj = options.traj
f_top = options.top
f_out = options.out
stride = options.stride

t = md.load(f_traj,top=f_top,stride=stride)

Ca = t.top.select('name CA')
aver_str = np.average(t.xyz[:,Ca], axis=0)
dat1  = np.swapaxes(t.xyz[:,Ca] - aver_str,1,2)

phi_ndx , phi = md.compute_phi(t)
psi_ndx , psi = md.compute_psi(t)
nres,n_fr = phi.transpose()[:-1].shape
dat2 = np.zeros([nres,2,n_fr]) 

dat2[:,0,:] = phi.transpose()[:-1]
dat2[:,1,:] = psi.transpose()[1:]

del(t)
DATA1= ts.TimeSer(dat1,n_fr,dim=3,nbins=options.nbins,reshape=False)
DATA2= ts.TimeSer(dat2,n_fr,dim=2,nbins=options.nbins,frame_row=False,reshape=False)

DATA1.calc_bins(opt=options.opt)
DATA2.calc_bins(opt=options.opt)

M1,E1 = DATA1.mutual_info_omp()
M2,E2 = DATA2.mutual_info_omp()
T1,D1 = DATA1.transfer_entropy_omp(time=3)
T2,D2 = DATA2.transfer_entropy_omp(time=3)

np.savetxt(f_out+".MI.RMSF.dat",M1)
np.savetxt(f_out+".TI.RMSF.dat",T1)
np.savetxt(f_out+".DI.RMSF.dat",D1)

np.savetxt(f_out+".MI.DIH.dat",M2)
np.savetxt(f_out+".TI.DIH.dat",T2)
np.savetxt(f_out+".DI.DIH.dat",D2)

quit()
