#!/usr/bin/env python

###########################################################################################################
#Clustering script no. 3                                                                                  #
#Modifications according to "Sammalkorpi, Karttunen and Haataja, J. Phys. Chem. B 2007, 111, 11722-11733" #
# Suggested cutoff values: 0.45 0.50 0.60 all in nm, the distances are converted into angstrom            #
#Cut off values can be increased with 25 %                                                                #
###########################################################################################################

###########################################################################################################
# Since the script do not take periodic boundaries into considerations, the input trajectory files        #
# needs to be treated in gmx trjconv first with -pbc cluster                                              #
###########################################################################################################


import argparse
import numpy as np
from numpy import array
import MDAnalysis
import MDAnalysis.analysis
import MDAnalysis.analysis.distances
import MDAnalysis.tests.datafiles
from MDAnalysis.tests.datafiles import PSF,DCD, PDB, XTC, GRO
import copy

parser = argparse.ArgumentParser(description='SDS MICELLE SCRIPT USING THREE CUTOFF VALUES. [ You can use the suggested cutoff values and/or perhaps increase those with 25% ]')
parser.add_argument('-f', dest='xtc_in' , action='store', type=str,  help='-f for input trajectory file (xtc)')
parser.add_argument('-s', dest='gro_in' , action='store', type=str,  help='-s for input topology file (gro, pdb)')
parser.add_argument('-o', dest='file_out', action='store', type=str,  help='-o for output file name')
parser.add_argument('-b', dest='begin' , action='store', type=int, help='-b for frame to start from in trajectory')
parser.add_argument('--dt', dest='stride' , action='store', type=int, help='--dt stride in ps thorough the trajectory')
parser.add_argument('--c1', dest='r1' , action='store', type=float, help='--c1 for the first cut_off (Angstrom) (5.6 suggested)')
parser.add_argument('--c2', dest='r2' , action='store', type=float, help='--c2 for the second cut_off (6.25 suggested)')
parser.add_argument('--c3', dest='r3' , action='store', type=float, help='--c3 for the third cut_off (7.5 suggested)')

args = parser.parse_args()
#indx1 = 195 #The first SDS molecule in the gro-file is numbered 195, in the tpr file it is 193
begin = args.begin
stride = args.stride
cut_off = (args.r1, args.r2, args.r3)

#for calculating the distance matrix between the center of mass of all the SDS molecules
def cal_mat_dist_com(u,res1,nres):
    dist = np.zeros((nres,nres))
    com=np.zeros([nres,3])
    for i in range(0,nres):
            i_res = u.select_atoms('resid {0:d}'.format(i+res1))
            com[i]=i_res.center_of_mass()
    for i in range(0,nres):
        for j in range(i+1,nres):
            dist[i,j] = np.linalg.norm(com[i]-com[j]) #Take the lenght of the vector
            dist[j,i] = dist[i,j]
    return dist

#A different function to calculate the distance matrix of the C5 atoms.
#Takes around 3 sec
def cal_mat_dist_C5_norm(u, res1, nres):
    dist = np.zeros((nres,nres))
    res_coord = np.zeros((nres,3))
    for i in range(0,nres):
        res_coord[i]= u.select_atoms('resid {0:d} and name C5'.format(i+res1)).positions
    for i in range(0,nres):
        for j in range(i+1,nres):
            dist[i,j] = np.linalg.norm(res_coord[i] - res_coord[j])
            dist[j,i] = dist[i,j]
    return dist

#Function for calculating the distance matrix between the C12 atoms in the SDS molecules
def cal_mat_dist_C12_norm(u, res1, nres):
    dist = np.zeros((nres,nres))
    res_coord = np.zeros((nres,3))
    for i in range(0,nres):
        res_coord[i]= u.select_atoms('resid {0:d} and name C12'.format(i+res1)).positions
    for i in range(0,nres):
        for j in range(i+1,nres):
            dist[i,j] = np.linalg.norm(res_coord[i] - res_coord[j])
            dist[j,i] = dist[i,j]
    return dist

#This function returns the starting stems
#c is the contact matrix of the type boolian
def stems(c,minsz=3):
    d = copy.deepcopy(c) #Make a dopy of the contact matrix
    st = []
    while (np.amax(np.sum(d, axis=1)) >= minsz ): #Keep finding stems, until the min size is reached
        nneighb=np.sum(d, axis=1) #The number of neighbors is the sum of the rows
        most_nneighb = np.argmax(nneighb) #The element (SDS) with the most neighbors
        stbuild = np.where(d[most_nneighb,:]==True)[0]  #find the neighbours in the stem
        st.append(stbuild) #Append to a list
        d[:,stbuild] = False #Set the coulums and rows just used to false
        d[stbuild,:] = False
    return st


def merge_stem (st):
    ns = len(st)
    if ns == 0: #If no stems are found, e.g. if no clusters is above minsz
        return st
    check_list = np.zeros(ns,dtype=int) #This list keep track of the stems, which is merged.
    merge_st=[]
    for i in range(ns-1):
        if check_list[i] == 1: #If Check_list = 1, stems have allready been merged. Therefore continue with the next i.
            continue
        l_merge=copy.deepcopy(st[i]) #Make a copy of the stem
        for j in range(i+1,ns):
            for e in np.nditer(l_merge):
                if np.any(e==st[j]): #If any element in the stem, is the same as any element in any of the other stems, merge. 
                    check_list[i] = 1
                    check_list[j] = 1
                    l_merge = np.unique(np.concatenate([l_merge,st[j]]))
                    break
        merge_st.append(l_merge)
    if check_list[ns-1]==0: 
        merge_st.append(st[-1])
    # return the list of micelles ordered by size
    for idx in np.argsort([ len(i) for i in merge_st ]):
        out_ms.append(merge_st[idx+1])
    return out_ms

#Function for checking three different cutoff values for clustering. 
#The vector cut_off contain these three values.
# The distances matrices are given, and checked if any element is below or equal to each cut off value
# The all in all 9 matrices are then sum, and at least one of the cut off value should be valid, for the SDS molecule
# to continue being clustered. 
# The function will then return a boolian matrix (contact matrix) from which stems can be found, merged and grown. 
def check_contact (cut_off, dist_com, dist_C5, dist_C12):
    test = [ 1 ,2 ,3 ]
    nres = dist_com.shape[0]
    D = np.zeros([nres,nres],dtype=int)
    for i,c in enumerate(cut_off):
        Dcom = (dist_com <= c).astype(int)
        DC5 = (dist_C5 <= c).astype(int)
        DC12 = (dist_C12 <= c).astype(int)
        D = D + (Dcom + DC5 + DC12 >= test[i] ).astype(int)
    return D >= 1

#For growing of the stems
def grow_stem(ms, check_list,contact):
    if len(ms) == 0:
        return ms
    c_ms=copy.deepcopy(ms)
    for i,s in enumerate(ms):
        if check_list[i] == 0:
            continue
        list_new=[]
        for j in c_ms[i]: #enter the i.th list containing the values j
            neighbor = np.where(contact[:,j])[0] # go into the contact matrix to look for neighbor for each j
            for k in neighbor:
                if not np.any(k == c_ms[i]):      #if the neighbor is not already on the list, append it to the list
                    #now append k to the cluster(i) = ms[i]
                    list_new.append(k) 
        if (list_new != []):
            check_list[i]=1
            c_ms[i]=np.insert(c_ms[i],len(c_ms[i]),np.array(list_new))
        else:
            check_list[i]=0
    return c_ms, check_list

#Loading of the universe.
u=MDAnalysis.Universe(args.gro_in, args.xtc_in)
print "Universe loaded"
sds = u.select_atoms('resname SDS').residues
nres = len(sds)
print 'Number of SDS molecules found in system {0:d}'.format(nres)
indx1 = sds[0].resid

#The loop over the trajectory. Calculating the three distance matrices for each timestep, merge and grow them.  
out=np.zeros([(len(u.trajectory)-begin)/stride+1,nres],dtype=int)
for ts in u.trajectory[begin::stride]:
    f = (u.trajectory.frame-begin)/stride
    dist_com = cal_mat_dist_com(u,indx1,nres)
    dist_C5 = cal_mat_dist_C5_norm(u,indx1,nres)
    dist_C12 = cal_mat_dist_C12_norm(u,indx1,nres)
    contact = check_contact(cut_off, dist_com, dist_C5, dist_C12)
    s = stems(contact)
    ms = merge_stem(s)
    check_list = np.ones(len(ms))
    while (np.any(check_list==1)):        
        ms, check_list = grow_stem(ms, check_list, contact)
    ms = merge_stem(ms)
    for r,m in enumerate(ms):
        out[f,m] = r+1
    print "Frame {0:d} done".format(f)
np.savez(args.file_out,out)
#np.savez(args.file_out+'.micelle_cluster.npz',out)


quit ()
