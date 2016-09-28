#!/usr/bin/env python

import re
import os
import numpy
from optparse import OptionParser
from hb import HBonds
import xpm

parser = OptionParser()

parser.add_option("-l","--log", dest="log", action="store", type="string", default=None, help="hbond.log",metavar="LOGFILE")
parser.add_option("-x","--xpm", dest="xpm", action="store", type="string", default=None, help="xpm file",metavar="XPMFILE")
parser.add_option("-g","--gro", dest="gro", action="store", type="string", default=None, help="gro file",metavar="GROFILE")
parser.add_option("-o","--out", dest="out", action="store", type="string", default=None, help="output file",metavar="DATFILE")
parser.add_option("-r","--red", dest="red", action="store_true",                         help="Write out redundat h-hbonds")

(options, args) = parser.parse_args()

out_base=os.path.basename(options.out)
out_dir=os.path.dirname(options.out)
if out_dir=="":
    out_dir="."

united_hbonds = HBonds(log=options.log,mol=options.gro,xpm=options.xpm,red=options.red)

united_hbonds.write_file(out_dir+'/'+out_base+'.dat')

exit()
