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
ix = np.triu_indices(len(c),1)
p = np.linspace(0,100,1002)
perc = np.percentile(c[ix],p)
SIGN = (c > p[981]).astype(int) - (c < p[21]).astype(int)
plt.matshow(c,cmap=plt.get_cmap('coolwarm'),vmin=-1.,vmax=+1)
plt.colorbar()
plt.savefig(options.out+".svg")
plt.matshow(SIGN,cmap=plt.get_cmap('coolwarm'))
plt.title("Significance regions")
plt.savefig(options.out+".significance.svg")
#plt.show()

quit()
