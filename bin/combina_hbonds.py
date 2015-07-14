#!/usr/bin/env python

import math, sys, time
import numpy as np
import scipy.stats
import hb

from optparse import OptionParser

parser = OptionParser()

parser.add_option("--hb1", dest="hb1", action="store", type="string", default=None, help="First percent.dat",metavar="LOGFILE")
parser.add_option("--hb2", dest="hb2", action="store", type="string", default=None, help="Second percent.dat",metavar="LOGFILE")
parser.add_option("--gro", dest="gro", action="store", type="string", default=None, help=".gro file",metavar="GROFILE")

(options, args) = parser.parse_args()

hb1=hb.hbonds(name='First')
hb2=hb.hbonds(name='Second')

hb1.read_file_perc(options.hb1,options.gro)
hb2.read_file_perc(options.hb2,otpions.gro)

hb3=hb.hbonds(name='Combination')
hb3=hb.merge_two(hb1,hb2)

print hb3






