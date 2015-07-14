#!/usr/bin/env python

import re
import numpy
from optparse import OptionParser
import hb
import xpm

from optparse import OptionParser

parser = OptionParser()

parser.add_option("-l", "--hb", dest="hb",  action="store", type="string", default=None, help="DNA hbond.log",metavar="LOGFILE")
parser.add_option("-g","--gro", dest="gro", action="store", type="string", default=None, help="gro file",metavar="GROFILE")

(options, args) = parser.parse_args()

hbonds=hb.hbonds(name='DNA')

hbonds.read_file_perc_red(options.hb,options.gro)
wc_list=hb.hbonds(name='Watson-Crick pairs')

for h in hbonds.hblist:
	if (h.don[0]=='DA') | (h.don[0]=='DA5') | (h.don[0]=='DA3') \
	 & (h.acc[0]=='DT') | (h.acc[0]=='DT5') | (h.acc[0]=='DT3'):
		if (h.don_atom=='N6') & (h.acc_atom=='O4'):
			wc_list.hblist.append(h)
	elif (h.don[0]=='DT') | (h.don[0]=='DT5') | (h.don[0]=='DT3') \
	   & (h.acc[0]=='DA') | (h.acc[0]=='DA5') | (h.acc[0]=='DA3'):
		if (h.don_atom=='N3') & (h.acc_atom=='N1'):
			wc_list.hblist.append(h)
	elif (h.don[0]=='DG') | (h.don[0]=='DG5') | (h.don[0]=='DG3') \
	   & (h.acc[0]=='DC') | (h.acc[0]=='DC5') | (h.acc[0]=='DC3'):
		if (h.don_atom=='N1') & (h.acc_atom=='N3'):
			wc_list.hblist.append(h)
		elif (h.don_atom=='N2') & (h.acc_atom=='O2'):
			wc_list.hblist.append(h)
	elif (h.don[0]=='DC') | (h.don[0]=='DC5') | (h.don[0]=='DC3') \
	   & (h.acc[0]=='DG') | (h.acc[0]=='DG5') | (h.acc[0]=='DG3'):
		if (h.don_atom=='N4') & (h.acc_atom=='O6'):
			wc_list.hblist.append(h)

wc_list.refresh()
#for h in wc_list.hblist:
#	print " (%d) %s --%s--H--%s-- %s (%d)"  %(h.don[1],h.don[0],h.don_atom,h.acc_atom,h.acc[0],h.acc[1])

for i in range(1,hbonds.mol.nres/2+1):
	for h in wc_list.hblist:
		#print h.don, h.don_atom, h.acc, h.acc_atom
		if (h.don[1]==i) & (h.perc > 50.0):
			print " (%4d) %4s --%s-H---%s-- %-4s (%4d) : %6.2f"  %(h.don[1],h.don[0],h.don_atom,h.acc_atom,h.acc[0],h.acc[1],h.perc)
		elif (h.acc[1]==i) & (h.perc > 50.0):
			print " (%4d) %4s --%s---H-%s-- %-4s (%4d) : %6.2f"  %(h.acc[1],h.acc[0],h.acc_atom,h.don_atom,h.don[0],h.don[1],h.perc)

