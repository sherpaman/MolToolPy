#!/usr/bin/env python
import os
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import MDAnalysis as MD
import MDAnalysis.analysis.distances as dist

from argparse import ArgumentParser
parser = ArgumentParser( description = 'Perform Binding Sites Analysis')
#
# INPUT FILES
#
parser.add_argument("-l","--list_sites",dest="list_sites",action="store",type=str,default=None,help="List of binding sites to analyze Matrix",required=False,metavar="TXT")
parser.add_argument("-t","--topol",dest="top",action="store",type=str,default=None,help="Topology Path",required=False,metavar="TOPOL")
parser.add_argument("-f","--traj",dest="trj",action="store",type=str,default=None,help="Trajectory Path",required=False,metavar="TRAJ")
parser.add_argument("--npz",dest="npz",action="store_true",default=False,help="Toggle use of Numpy files instead of trajectory")
#
# OUTPUT FILES
#
parser.add_argument("-o","--out",dest="out",action="store",type=str,default=None,help="Output base name",required=True,metavar="FILENAME")
#
# OTHER OPTIONS
#
parser.add_argument("-b","--begin",dest="begin",action="store",type=int,default=0,help="First frame to read")
parser.add_argument("-e","--end",dest="end",action="store",type=int,default=-1,help="Last frame to read")
parser.add_argument("-s","--skip",dest="skip",action="store",type=int,default=1,help="number of frame to skip", metavar="INTEGER")
parser.add_argument("-c","--dist_cutoff",dest="cutoff",action="store",type=float,default=7.0,help="Distance Cutoff")

options = parser.parse_args()

top       = options.top
trj       = options.trj
base_out  = options.out
b         = options.begin
e         = options.end
skip      = options.skip
l_sites   = options.list_sites
npz       = options.npz
dist_cutoff = options.cutoff # Distance Cutoff in Aangstrom

smooth = 2
conv = np.ones(smooth)

def read_list_sites(file_in):
    with open(file_in) as f:
        raw = []
        for rl in f.readlines():
            if rl[0] != "#":
                raw.append(rl.split())
        list_a = [ i for i in raw if i[1]=="A" ]
        list_b = [ i for i in raw if i[1]=="B" ]
    return list_a, list_b

def calc_acfs(data):
    acfs = np.zeros(len(data))
    conv = np.ones(2, dtype = np.int32)
    # compute the ACFs for the chunk of binding sites passed in argument
    
    # extract the data for this site and this lipid to avoid dimensionality nightmares
    red = np.copy(data)
    # first value is the sum of all single-frame occupancy
    acfs[0] = float(np.sum(red))
    # if the ACF reached 0 already, no need to compute the rest
    for k in range(1,len(data)):
        if not acfs[k-1]:
            acfs[k:] = 0
            break
        # ACF at this point
        red = np.convolve(red, conv, "valid") == 2
        acfs[k] = float(np.sum(red))
        # normalization
    acfs = acfs/acfs[0]
    return acfs

def _exp2(x,l1,l2,a):
    return a * np.exp(-l1*x) + (1-a) * np.exp(-l2*x)

def scan_tau_2(a,n=50):
    nz = np.where(a==0)[0][0]
    z_range = np.linspace(10,nz,n).astype(int)
    l = np.zeros([len(z_range),3])
    l0 = a[0] - a[1]
    for i,z in enumerate(z_range):
        l[i], o = opt.curve_fit(_exp2,np.arange(z),a[:z],p0=[l0,l0/100.,0.5],bounds=( [0.,0.,0.],[np.inf, np.inf,1.0] ))
    max_l=np.argmax(1./l[:,0])
    #return l[max_l], [z_range ,l]
    return l[max_l]

def opt_tau_2(a,n=50):
    nz = np.where(a==0)[0][0]
    z_range = np.linspace(10,nz,n).astype(int)
    l0 = a[0] - a[1]
    l, o = opt.curve_fit(_exp2,np.arange(nz),a[:nz],p0=[l0,l0/100.,0.5],bounds=( [0.,0.,0.],[np.inf, np.inf,1.0] ))
    return l

if not npz:
    l_a , l_b = read_list_sites(l_sites)
    
    l_traj_a = np.unique([i[0] for i in l_a ])
    l_traj_b = np.unique([i[0] for i in l_b ])
    l_traj = list(set(l_traj_a) | set(l_traj_b))
    l_traj.sort()
    
    l_u = [ MD.Universe("{p_top:s}/{n_top:s}.tpr".format(p_top=top,n_top=i),"{p_trj:s}/{n_trj:s}.xtc".format(p_trj=trj,n_trj=i)) for i in l_traj ] 
    
    text_out = open('{0:s}_ANAL.txt'.format(base_out),'w')
    
    # The analysis is performed one trajectory at a time
    for ntrj,trj in enumerate(l_traj):
        u = l_u[ntrj]
        if e == -1:
            e = len(u.trajectory)
        # the list of bindibg site to be analyzed for each trajectory is created
        sub_l_a = [ i for i in l_a if i[0] == trj]
        sub_l_b = [ i for i in l_b if i[0] == trj]
        num_a = len(sub_l_a)
        num_b = len(sub_l_b)
        s1 = [ [] for i in range(num_a)]
        s2 = [ [] for j in range(num_b)]
        name_a = [ '' for i in range(num_a)]
        name_b = [ '' for j in range(num_b)]
        # per each couple of binding site "A" - "B" (on the receptor and ligand respectively
        # the atom selection is stored in a list 
        for n_i,i in enumerate(sub_l_a):
            n_a = i[2]
            bs_a = [ int(i[l]) for l in range(3,len(i)) ]
            str1 = "resnum"
            for r in bs_a:
                str1 = "{0:s} {1:d}".format(str1,r)
            s1[n_i] = u.select_atoms(str1)
            name_a[n_i]=""
            for ra in s1[n_i].residues:
                    name_a[n_i] = "{0:s} {1:s}_{2:s}".format(name_a[n_i],ra.resname,str(ra.resid))
            for n_j,j in enumerate(sub_l_b):
                n_b = j[2]
                bs_b = [ int(j[k]) for k in range(3,len(j)) ]
                str2 = "resnum"
                for p in bs_b:
                    str2 = "{0:s} {1:d}".format(str2,p+25)
                s2[n_j] = u.select_atoms(str2)
                name_b[n_j]=""
                for rb in s2[n_j].residues:
                    name_b[n_j] = "{0:s} {1:s}_{2:s}".format(name_b[n_j],rb.resname,str(rb.resid))
                print name_a[n_i], name_b[n_j]
        # A 2-fold nested list is used to store an array of minimal residue-residue distances per each couple "A"-"B" 
        M = [ [ [] for j in np.arange(num_b) ] for i in np.arange(num_a) ]
        d = np.zeros(2) 
        w = np.zeros(2) 
        o = np.zeros([num_a,num_b,(e-b)/skip])
        time = np.zeros([num_a,num_b,(e-b)/skip])
        # Per each frame ...
        for ts in u.trajectory[b:e:skip]: 
            # ... each site in the Receptor (Site A) ...
            for i in np.arange(num_a):
                # ... and each site in the Ligand (Site B) ...
                for j in np.arange(num_b):
                    # ... the minimum distance between each residue in Site A and each residue in Site B is calculated ...
                    M[i][j] = [ [ np.min(dist.distance_array(s_1.positions,s_2.positions,box=u.dimensions)) for s_2 in s2[j].residues] for s_1 in s1[i].residues ]
                    d[0] = np.mean(np.min(M[i][j],axis=0)) # Average minumum distance per each residue in Site A
                    d[1] = np.mean(np.min(M[i][j],axis=1)) # Average minumum distance per each residue in Site B
                    w[0] = np.var(np.min(M[i][j],axis=0))  
                    w[1] = np.var(np.min(M[i][j],axis=1))
                    time[i,j,(ts.frame-b)/skip-1] = ts.time
                    # ... STD.DEV Weighted Average ...
                    o[i,j,(ts.frame-b)/skip-1] = np.average(d,weights=w)
        text_out.write("TRAJ {0:s} ANAL\n".format(trj))
        np.savez('{0:s}_{1:s}.npz'.format(base_out,trj),occupancy=o,time=time,name_s1=name_a,name_s2=name_b)
        for i in range(num_a):
            for j in range(num_b):
                if smooth > 1:
                    running = np.convolve(o[i,j],conv,"valid")/smooth
                else:
                    running = o[i,j]
                occ=100.0 * float(np.sum(running<dist_cutoff))/len(running)
                if occ > 5.0:
                    acsf = calc_acfs(running<dist_cutoff)
                    try:
                        l = opt_tau_2(acsf)
                    except:
                        print ("Error analysis Traj : {0:s} , [{1:s} - {2:s}]".format(trj,name_a[i],name_b[j]))
                        text_out.close()
                        raise
                    dt = u.trajectory[skip].time - u.trajectory[0].time 
                    t1 = (1./l[0]) * dt / 1000.0
                    t2 = (1./l[1]) * dt / 1000.0
                    p1 = 100*l[2]
                    p2 = 100-p1
                    text_out.write( "BINDING : [{0:s}]-[{1:s}] {2:6.2f} , t1={3:.2f} ({4:6.2f}), t2={5:.2f} ({6:6.2f})\n".format(name_a[i],name_b[j],occ,t1,p1,t2,p2))
                else:
                    text_out.write( "BINDING : [{0:s}]-[{1:s}] {2:6.2f}\n".format(name_a[i],name_b[j],occ))
else:
    list_npz = [ i for i in os.listdir('./') if ( i[:len(base_out)] == base_out ) & ( i[-3:]=='npz' ) ]
    list_npz.sort()
    text_out = open('{0:s}_ANAL.txt'.format(base_out),'w')
    for trj in list_npz:
        print("opening : {0:s}\n".format(trj))
        text_out.write("{0:s} ANALYSIS\n".format(trj))
        data = np.load(trj)
        o = data['occupancy']
        time = data['time']
        name_a = data['name_s1']
        name_b = data['name_s2']
        num_a = len(name_a)
        num_b = len(name_b)
        for i in range(num_a):
            for j in range(num_b):
                if smooth > 1:
                    running = np.convolve(o[i,j],conv,"valid")/smooth
                else:
                    running = o[i,j]
                occ=100.0 * float(np.sum(running<dist_cutoff))/len(running)
                if occ > 5.0:
                    b = (running<dist_cutoff).astype(int)
                    bind = 1 - np.trim_zeros(1-b) # Remove Bound satates frome the start and end
                    acsf_off = calc_acfs(bind)
                    acsf_on  = calc_acfs(1-b)
                    try:
                        l_off = opt_tau_2(acsf_off)
                        l_on  = opt_tau_2(acsf_on)
                    except:
                        print ("Error analysis Traj : {0:s} , [{1:s} - {2:s}]".format(trj,name_a[i],name_b[j]))
                        text_out.close()
                        raise
                    dt = time[i,j,1] - time[i,j,0] 
                    #
                    t1_on = (1./l_on[0]) * dt / 1000.0
                    t2_on = (1./l_on[1]) * dt / 1000.0
                    p1_on = 100*l_on[2]
                    p2_on = 100-p1_on
                    #
                    t1_off = (1./l_off[0]) * dt / 1000.0
                    t2_off = (1./l_off[1]) * dt / 1000.0
                    p1_off = 100*l_off[2]
                    p2_off = 100-p1_off
                    text_out.write( " BINDING :[{0:s}]-[{1:s}]:{2:6.2f}\n".format(name_a[i],name_b[j],occ))
                    text_out.write( " K_on    :t1={0:7.2f}({1:6.2f}):t2={2:7.2f}({3:6.2f})\n".format(t1_on ,p1_on ,t2_on ,p2_on ))
                    text_out.write( " K_off   :t1={0:7.2f}({1:6.2f}):t2={2:7.2f}({3:6.2f})\n".format(t1_off,p1_off,t2_off,p2_off))
                else:
                    text_out.write( " BINDING :[{0:s}]-[{1:s}]:{2:6.2f}\n".format(name_a[i],name_b[j],occ))
text_out.close()
quit()
            
