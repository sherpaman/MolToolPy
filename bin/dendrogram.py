#!/usr/bin/env python

import os, sys
import numpy as np
import scipy.cluster.hierarchy as hc
import matplotlib.pyplot as plt

file_in = sys.argv[1]
D = np.load(file_in)['arr_0']
idx=np.triu_indices(D.shape[0],1)
Z = hc.linkage(D[idx],method='complete')
plt.subplot(2,2,1)
hc.dendrogram(Z,distance_sort='ascending',truncate_mode='level',p=10)
plt.subplot(2,2,2)
hc.dendrogram(Z,distance_sort='ascending',truncate_mode='level',p=2)
plt.subplot(2,1,2)
hc.dendrogram(Z,distance_sort='ascending')
plt.show()
