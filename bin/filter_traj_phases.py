#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import phrases
#import MDAnalysis as MD
from argparse import ArgumentParser
parser = ArgumentParser( description = 'Filter Tarjectories Based on Phrases')
#
#
# INPUT FILES
#
#parser.add_argument("-p","--phrases",dest="phrases",action="store",type=str,default=None,help="Phrases File",required=True,metavar="PICKLE FILE")
#parser.add_argument("-f","--traj"   ,dest="traj"   ,action="store",type=str,default=None,help="Input Trajectory",required=True,metavar="TRAJ FILE")
#parser.add_argument("-t","--top"    ,dest="top"    ,action="store",type=str,default=None,help="Input Topology",required=True,metavar="TOPOL FILE")
parser.add_argument("-l","--list"   ,dest="bs_list",action="store",type=str,default=None,help="File containing phrases to filter",required=True,metavar="TXT FILE")
#
# OUTPUT FILE
#
parser.add_argument("-o","--out"   ,dest="out",action="store",type=str,default=None,help="Frames ndx file",required=True,metavar="NDX FILE")
#
options = parser.parse_args()

bs_list = open(options.bs_list,'r')

trj_list     = []
name_list    = []
phrases_list = []

for r in bs_list.readlines():
    row = r.split()
    trj_list.append(row[0])
    name_list.append(row[1])
    phrases_list.append([ int(i) for i in [ row[n] for n in range(2,len(row)) ] ])

open_trj_list = list(set(trj_list))
open_trj_list.sort()

for t in trj_list:
    print("{0:s} is element {1:d}".format(t,open_trj_list.index(t)))
bs_list.close()

P=[]

print("Starting Reading Phrases:")
for i in open_trj_list:
    print("{0:s}".format(i))
    P.append(phrases.read_phrases('{0:s}-phrases.dat'.format(i)))
    print("DONE!")

ndx=open(options.out,'w')

print("Starting Writing NDX")
for n,p in enumerate(trj_list):
    name_long="{0:s}".format(name_list[n])
    for i in phrases_list[n]:
        name_long=name_long+"_{0:d}".format(i)
    print("[ {0:s}.{1:s} ]\n".format(trj_list[n],name_long))
    idx_p = open_trj_list.index(p)
    c = phrases_list[n]
    frames = np.where([ np.any(np.array([ set(c) == set(c) & set(i) for i in j])) for j in P[idx_p].phrases ])[0]
    ndx.write("[ {0:s}.{1:s} ]\n".format(trj_list[n],name_long))
    for f,i in enumerate(frames):
        ndx.write(" {0:6d}".format(i))
        if (f%15 == 0) & (f != 0):
            ndx.write("\n")
    ndx.write("\n")

num_trj=len(trj_list)
for i in range(num_trj):
    for j in range(i+1,num_trj):
        a = trj_list[i].split('.')[0]
        b = trj_list[j].split('.')[0]
        if a == b:
            idx_a = open_trj_list.index(trj_list[i])
            idx_b = open_trj_list.index(trj_list[j])
            print( "{0:s} vs. {1:s} [{2:d}-{3:d}]".format( trj_list[i], trj_list[j], idx_a, idx_b ) )
            name_long_a="{0:s}".format( name_list[i] )
            for ia in phrases_list[i]:
                name_long_a=name_long_a+"_{0:d}".format(ia)
            name_long_b="{0:s}".format( name_list[j] )
            for jb in phrases_list[j]:
                name_long_b=name_long_b+"_{0:d}".format(jb)
            print( "[ {0:s}.{1:s}_&_{2:s}.{3:s} ]\n".format( trj_list[i], name_long_a, trj_list[j], name_long_b ) )
            c_a = phrases_list[i]
            c_b = phrases_list[j] 
            f_a = np.array([ np.any(np.array([ set(c_a) == set(c_a) & set(p) for p in p_f])) for p_f in P[idx_a].phrases ])
            f_b = np.array([ np.any(np.array([ set(c_b) == set(c_b) & set(p) for p in p_f])) for p_f in P[idx_b].phrases ])
            frames = np.where(f_a & f_b)[0]
            ndx.write( "[ {0:s}.{1:s}_&_{2:s}.{3:s} ]\n".format( trj_list[i], name_long_a, trj_list[j], name_long_b ) )
            for f,o in enumerate(frames):
                ndx.write(" {0:6d}".format(o))
                if (f%15 == 0) & (f != 0):
                    ndx.write("\n")
            ndx.write("\n")

ndx.close()
