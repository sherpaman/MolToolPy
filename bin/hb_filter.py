#!/usr/bin/env python

import hb 

from argparse import ArgumentParser
parser = ArgumentParser( description = 'Filter HB')
#
# INPUT FILES
#
parser.add_argument("-d","--dat",dest="dat",action="store",type=str,default=None,help="HB Perc. Data",required=True,metavar="DAT FILE")
parser.add_argument("-f","--gro",dest="gro",action="store",type=str,default=None,help=".gro file",required=True,metavar="GRO FILE")
#
# OUTPUT FILES
#
parser.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=False,help="Output File Name",metavar="DAT FILE")
#
# VAR ARGUMENTS
#
parser.add_argument("--a0",dest="a_0",action="store",type=int,required=True,help="Interval 1 start", metavar="INTEGER")
parser.add_argument("--a1",dest="a_1",action="store",type=int,required=True,help="Interval 1 end", metavar="INTEGER")
parser.add_argument("--b0",dest="b_0",action="store",type=int,required=True,help="Interval 2 start", metavar="INTEGER")
parser.add_argument("--b1",dest="b_1",action="store",type=int,required=True,help="Interval 2 end", metavar="INTEGER")



parser.add_argument("-t","--threshold",dest="p_thresh",action="store",type=float,default=45,help="Percentage Threshold", metavar="FLOAT")
#
#
options = parser.parse_args()

def betw(i,a,b):
	return ( a <= i ) & ( i <= b )

def outs(i,a,b):
	return not betw(i,a,b)

dat = options.dat
gro = options.gro

a_0 = options.a_0
a_1 = options.a_1

b_0 = options.b_0
b_1 = options.b_1

p_t = options.p_thresh

d = hb.HBonds()
d.read_file_perc(dat,gro)
hb_sub=[i for i in d.hblist if (( betw(i.don[1],a_0,a_1) & betw(i.acc[1],b_0,b_1) ) | ( betw(i.acc[1],a_0,a_1) & betw(i.don[1],b_0,b_1) )) & (i.perc > p_t) ]

for i in hb_sub:
    print i

quit()

