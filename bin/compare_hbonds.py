#!/usr/bin/env python

from hb import Hbonds, HBondsCompare
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-a","--hb1" , dest="hb1",       action="store", type="string", default=None, help="First hbond.log",     metavar="LOGFILE")
parser.add_option("-b","--hb2" , dest="hb2",       action="store", type="string", default=None, help="Second hbond.log",    metavar="LOGFILE")
parser.add_option("-g","--gro1", dest="gro1",      action="store", type="string", default=None, help="gro file",            metavar="GROFILE")
parser.add_option("-f","--gro2", dest="gro2",      action="store", type="string", default=None, help="gro file",            metavar="GROFILE")
parser.add_option("-p","--perc", dest="perc",      action="store", type="float",  default=0.0,  help="Percentage Threshold",metavar="REAL")
parser.add_option("-r","--red" , dest="redundant", action="store_true",                         help="Write out also redundat h-hbonds")

(options, args) = parser.parse_args()

hb1 = HBonds(name='First')
hb2 = HBonds(name='Second')

hb1.read_file_perc(options.hb1,options.gro1)
hb2.read_file_perc(options.hb2,options.gro2)

compare = HBondsCompare(first=hb1,second=hb2)

compare.do_comparison(options.perc)

print compare






