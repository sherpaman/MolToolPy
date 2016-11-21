#!/usr/bin/env python

#import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
parser = ArgumentParser( description = 'Calculate X-CORRELATION Matrix from g_covar ascii outpoot')
#
# INPUT FILES
#
parser.add_argument("-i","--inp",dest="inp",action="store",type=str,default=None,help="g_covar ASCII output",required=True,metavar="ASCII DAT FILE")
#
# OUTPUT FILES
#
parser.add_argument("-o","--out",dest="out",action="store",type=str,default=None,required=True,help="Output File Name",metavar="DAT FILE")


options = parser.parse_args()

raw_dat = np.loadtxt(options.inp)

N,dim = raw_dat.shape


D = raw_dat.reshape(np.sqrt([N*dim,N*dim]).astype(int))

n_res = int(np.sqrt(N*dim)/3)

v = np.zeros(n_res)
c = np.ones([n_res,n_res])

for i in np.arange(n_res):
    v[i] = np.sqrt(np.trace(D[3*i:3*(i+1),3*i:3*(i+1)]))
    for j in np.arange(i):
        c[i,j] = np.trace(D[3*i:3*(i+1),3*j:3*(j+1)]) / (v[i]*v[j])
        c[j,i] = c[i,j]

np.savez(options.out,c)
#np.savetxt(options.out+'.txt',c)

plt.matshow(c,cmap=plt.get_cmap('coolwarm'))
plt.colorbar()
plt.savefig(options.out+".svg",fmt='svg')
#plt.show()

quit()
