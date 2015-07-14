#!/usr/bin/env python

import math, sys, time
import numpy as np
import scipy.stats
import hb, xpm

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-l","--hb" , dest="log" , action="store", type="string", default=None, help="hbond.log",metavar="LOGFILE")
#parser.add_option("-x","--xpm", dest="xpm", action="store", type="string", default=None, help="hbmax.xpm",metavar="XPMFILE")

(options, args) = parser.parse_args()

hbonds=hb.hbonds(name='Hbonds')

l_hb, red_list_hb, ref_hb=hbonds.read_log(options.log)

mediated=[]
med_to_list=[]
list_to_med=list(np.zeros(len(l_hb)))

for n in np.arange(len(l_hb)):
	# IF DONOR IS A WATER
	if hb.is_protein(l_hb[n][0]):
		for m in np.arange(n+1,len(l_hb)):
			# FIND THE OTHER HBONDS HAVING AS ACCEPTOR THE SAME WATER DONOR "i"
			if (l_hb[m][1] == l_hb[n][0]):
				temp=[ l_hb[n][0], l_hb[m][1]]
				if temp in mediated:
					med_to_list[mediated.index(temp)].append(n)
					med_to_list[mediated.index(temp)].append(m)
					list_to_med[m].append(temp)
					list_to_med[n].append(temp)
				else:
					mediated.append(temp)
					med_to_list.append([n,m])
					list_to_med[n]=[temp]
					list_to_med[m]=[temp]
	elif  hb.is_protein(l_hb[n][1]): #IF 
		for m in np.arange(n+1,len(l_hb)):
			# FIND THE OTHER HBONDS HAVING AS DONOR THE SAME WATER DONOR "i"
			if (l_hb[m][0] == l_hb[n][1]):
				temp=[ l_hb[m][0], l_hb[n][1]]
				if temp in mediated:
					med_to_list[mediated.index(temp)].append(n)
					med_to_list[mediated.index(temp)].append(m)
					list_to_med[m].append(temp)
					list_to_med[n].append(temp)
				else:
					mediated.append(temp)
					med_to_list.append([n,m])
					list_to_med[n]=[temp]
					list_to_med[m]=[temp]

for i,m in enumerate(mediated):
	print m
	print med_to_list[i]

 
#data_raw=xpm.read_xpm(options.xpm)
