#!/usr/bin/env python

import re
import os
import numpy
from optparse import OptionParser
import hb
import xpm

parser = OptionParser()

parser.add_option("-l","--log", dest="log", action="store", type="string", default=None, help="hbond.log",metavar="LOGFILE")
parser.add_option("-x","--xpm", dest="xpm", action="store", type="string", default=None, help="xpm file",metavar="XPMFILE")
parser.add_option("-g","--gro", dest="gro", action="store", type="string", default=None, help="gro file",metavar="GROFILE")
parser.add_option("-o","--out", dest="out", action="store", type="string", default=None, help="output file",metavar="DATFILE")
#parser.add_option("-r"   , dest="redundant", action="store_true", help="Write out also redundat h-hbonds")
parser.add_option("-v"   , dest="verbose", action="store_true", help="Be verbose")

(options, args) = parser.parse_args()

out_base=os.path.basename(options.out)
out_dir=os.path.dirname(options.out)
if out_dir=="":
	out_dir="."


united_hbonds=hb.hbonds(log=options.log,mol=options.gro,file_xpm=options.xpm)

united_hbonds.write_file(out_dir+'/'+out_base+'.dat')

exit()
