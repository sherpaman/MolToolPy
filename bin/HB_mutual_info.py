#!/usr/bin/env python

import numpy as np
from ts import TimeSer
import hb
import xpm
import gro
import sys 
import os
import math
import gzip
import json
import matplotlib.pyplot as plt

#
# SET THE INTERFACE MODULES
#
from argparse import ArgumentParser
parser = ArgumentParser( description = 'Calculate Mutual Information')
#
# INPUT FILES
#
parser.add_argument("-f","--gro",dest="file_gro",action="store",type=str,default=None,help="Input GRO File",required=True,metavar="GRO FILE")
parser.add_argument("-l","--log",dest="file_log",action="store",type=str,default=None,help="Input LOG File",required=True,metavar="LOG FILE")
parser.add_argument("-x","--xpm",dest="file_xpm",action="store",type=str,default=None,help="Input XPM File",required=True,metavar="XPM FILE")
parser.add_argument("-m","--mat",dest="file_mut",action="store",type=str,default=None,help="Mutual Info File",required=False,metavar="DAT FILE")
#
# OUTPUT FILES
#
parser.add_argument("-o","--out",dest="file_out",action="store",type=str,default=None,required=False,help="Output Base Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
#parser.add_argument("--nx",  dest="nx",  action="store",type=int,default=12,  help="nbin 1st VAR",                    metavar="INTEGER")
#parser.add_argument("--ny",  dest="ny",  action="store",type=int,default=12,  help="nbin other VAR",                  metavar="INTEGER")
parser.add_argument("--col1",dest="col1",action="store",type=int,default=0,   help="first column to read",            metavar="INTEGER")
parser.add_argument("--col2",dest="col2",action="store",type=int,default=-1,help="last column to read",             metavar="INTEGER")
parser.add_argument("--fr1", dest="fr1", action="store",type=int,default=0,   help="first line to read",              metavar="INTEGER")
parser.add_argument("--fr2", dest="fr2", action="store",type=int,default=-1,help="last line to read",               metavar="INTEGER")
parser.add_argument("--skip",dest="skip",action="store",type=int,default=1,   help="skip every this number of frames",metavar="INTEGER")
parser.add_argument("-t","--transfer",dest="transfer",action="store_true",default=False,help="Toggle calculation of Transfer Entropy")
options = parser.parse_args()

save_objects = []
out_array    = {}

#nbin_x   = options.nx
#nbin_y   = options.ny

print "Reading File..."

mol = gro.read_gro(options.file_gro)
hbonds  = hb.HBonds(mol=mol,xpm=options.file_xpm,log=options.file_log)
DATA = np.transpose(hbonds.nr_xpm.array)

if options.col1 < 0 :
        print " You must select a positive starting column (--col1 : {0:d})".format(options.col1)
        exit(2)
if options.col2 != -1:
        if options.col2 < options.col1:
                print " Last colum (--col2 : {1:4d}) must be higher than first (--col1 {0:4d})".format(options.col1,options.col2)
                exit(2)
        if options.col2 > ncol :
                print " Last colum (--col2 : {0:4d})  will be lowered to the actual number of column in the file ({1:4d}".format(options.col2,ncol-1)
                options.col2 = ncol
if options.fr1 < 0:
        options.fr1 = 0
if options.fr2 != -1:
        if options.fr2 <= options.fr1:
                print " Last colum (--fr2 : {1:6d}) must be higher than first (--fr1 {0:6d})".format(options.fr1,options.fr2)
                exit(2)
if options.skip < 1:
        print " You must select a positive skip tome (--skip : {0:4d})".format(options.skip)
        exit(2)


DATA = DATA[options.fr1:options.fr2:options.skip,options.col1:options.col2]

nframes, ncol = DATA.shape

hb_ts = TimeSer(data=DATA,n_data=nframes,dim=1,nbins=2,dtype=int)

if options.file_mut != None:
    M = np.loadtxt(options.file_mut)
else:
    M,E,P = hb_ts.mutual_info_for()

if options.transfer :
    Transfer = hb_ts.transfer_entropy_for(time=4)
    D = hb_ts.calc_direction(Transfer)
    T, [ [M_t,  E_t,  P_t ], [M_t1, E_t1, P_t1 ] ] = Transfer
    plt.matshow(D,origin='lower')
    plt.colorbar()
    plt.figure()
    plt.plot(np.sum(D,axis=0))
    plt.show()

n_hb = M.shape[0]

percentile = 1. - ( 10. / ( n_hb ** 2 ) ) # Will show at least the top-5
quantile_m = np.sort(np.ravel(M))[int(percentile*n_hb**2)-1] 

pm = np.where(M >= quantile_m)

print " Mutual Info {0:8.6f} percentile (val={1:10.6f};n={2:4d})".format(percentile,quantile_m,len(pm))

#print "{0:>10s}{1:>30s}{2:>30s}".format("Mut.Info","HB1","HB2")
#for i in range(len(pm[0])):
#    if pm[0][i] != pm[1][i]:
#        print "{0:>10.6f}{1:>30s}{2:>30s}".format(M_o[pm[0][i],pm[1][i]], str(hbonds.hblist[pm[0][i]])[:30], str(hbonds.hblist[pm[1][i]])[:30])

if options.file_out != None: 
    out_array["Mutual Info"]         = M.tolist()
    out_array["Joint Entropies"]     = E.tolist()
    out_array["Joint Probabilities"] = P.tolist()
    if options.transfer:
        out_array["Transfer Entropy"] = T.tolist() 


for i in out_array.keys():
    print " Out : {0:s}".format(i)

save_objects.append(out_array)

json_out=options.file_out.split('.')[0]+'.json.gz'

if options.transfer :
    tran_out=options.file_out.split('.')[0]+'_transfer-entropy.dat'
    np.savetxt(tran_out,T)

if save_objects != []:
    with gzip.GzipFile(json_out, 'w') as outfile:
        for obj in save_objects:
            outfile.write(json.dumps(obj) + '\n')



# reading
#with gzip.GzipFile(options.file_in, 'r') as isfile:
#    for line in infile:
#        obj = json.loads(line)
#


