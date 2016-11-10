#!/usr/bin/env python

import ts
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as md
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
parser.add_argument("-b","--begin",dest="begin",action="store",type=int,default=0,help="Start reading the trajectory at this frame", metavar="INTEGER")
parser.add_argument("-e","--end",dest="end",action="store",type=int,default=-1,help="Stop reading the trajectory at this frame", metavar="INTEGER")
parser.add_argument("-s","--stride",dest="stride",action="store",type=int,default=1,help="time", metavar="INTEGER")
parser.add_argument("-n","--nbins",dest="nbins",action="store",type=int ,default=10,help="number of bins", metavar="INTEGER")
parser.add_argument("--opt",dest="opt",action="store_true",default=False,help="toggle bins optimization")
#
options = parser.parse_args()

f_traj = options.traj
f_top = options.top
f_out = options.out
begin = options.begin
end = options.end
stride = options.stride

t = md.Universe(f_top,f_traj)

sel = 'name CA'

CA = t.select_atoms(sel)
aver_str = np.zeros((len(CA),3))
n_fr = (end - begin)/int(stride)+1

print ("Calculating average structure over {0:d} frames".format(n_fr))
for i in t.trajectory[begin:end:stride]:
    CA = t.select_atoms(sel)
    aver_str = aver_str + CA.atoms.coordinates()/n_fr
print ("Done!")

dat1 = np.zeros((len(CA),3,n_fr))

print("Start Calculating Fluctuations")
for i in t.trajectory[begin:end:stride]:
     CA = t.select_atoms(sel)
     fr=(i.frame - begin) / stride
     dat1[:,:,fr] = CA.atoms.coordinates() - aver_str     
print("Done!")

del(t)
DATA1= ts.TimeSer(dat1,n_fr,dim=3,nbins=options.nbins,reshape=False)

DATA1.calc_bins(opt=options.opt)

M1,E1 = DATA1.mutual_info_omp()
T1,D1 = DATA1.transfer_entropy_omp(time=2)

np.savetxt(f_out+".MI.RMSF.dat",M1)
np.savetxt(f_out+".TI.RMSF.dat",T1)
np.savetxt(f_out+".DI.RMSF.dat",D1)

quit()
