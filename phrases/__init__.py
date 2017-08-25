import numpy as np
import scipy.optimize as opt
import scipy.cluster.hierarchy as sch
import pickle
import itertools

def _exp2(x,l1,l2,a):
 return a * np.exp(-l1*x) + (1-a) * np.exp(-l2*x)


def read_phrases(f_dat,*args,**kwargs):
    with open(f_dat,'rb') as f:
        L=pickle.load(f)
        out = phrases(min_len=kwargs.get('min_len',2))
        if kwargs.has_key("dist"):
            out.D = np.load(kwargs["dist"])['arr_0']
        out.phrases = L
    return out

def _jaccard(a,b):
    return float(len(set(a) & set(b))) /  float(len(set(a) | set(b)))

def _mj(a,c,t):
    if len(a) == 0:
        return 0
    list_j = [ _jaccard(a,i) for i in c ] # Jaccad Similarity
    if max(list_j) < t :
        return 0
    else:
        return np.argmax(list_j)
        
def _similarity(a,c,t):
    if len(a) == 0:
        return 0
    list_s = [ float(len(set(a) & set(i))) / len(i) for i in c ] 
    if max(list_s) < t :
        return 0
    else:
        return np.argmax(list_s)

def _mj_subset(a,c,t):
    if len(a) == 0:
        return 0
    list_j = [ _j_subset(a,i) for i in c ] # Jaccad Similarity
    if max(list_j) < t :
        return 0
    else:
        return np.argmax(list_j)

def _j_subset(a,c):
    if len(a) == 0:
        return 0
    if len(a) <= len(c):
        min_list = a
        min_len = len(a)
        max_list = c
        max_len = len(c)
    elif len(a) > len(c):
        min_list = c
        min_len = len(c)
        max_list = a
        max_len = len(a)
    list_j = [ _jaccard(i,min_list) for i in itertools.combinations(max_list,min_len) ]
    return np.max(list_j)
        
def _percent_of_cluster(a,c,t): 
    list_j = [ float(len(set(a) & set(i))) /  float(len(set(i))) for i in c ] 
    if max(list_j) < t :
        return 0
    else:
        return np.argmax(list_j)

def _autocorr(data, trim_val=1):
    n_fr = len(data)
    data_t = _trim(data, trim_val)
    mu = np.average(data_t)
    sigma = np.var(data_t)
    result = np.correlate(data_t-mu, data_t-mu, mode='full')/(sigma**2)
    norm = np.arange(result.size/2+1,0,-1)
    result = result[result.size/2:]/norm
    result = result / result[0]
    n_r = len(result)
    if n_r == n_fr:
        return result
    else:
        out = np.zeros(n_fr)
        out[:n_r] = result
        return out

def _autocorr2(data, trim_val=1):
    n_fr = len(data)
    data_t = _trim(data, trim_val)
    result = np.correlate(data_t, data_t, mode='full')
    norm = np.arange(result.size/2+1,0,-1)
    result = result[result.size/2:]/norm
    result = result / result[0]
    n_r = len(result)
    if n_r == n_fr:
        return result
    else:
        out = np.zeros(n_fr)
        out[:n_r] = result
        return out

def _acfs(occupancies):
    # Code from Arnarez.
    acfs = np.zeros(len(occupancies))
    convolve_array = np.ones(2, dtype = np.int32)
    # compute the ACFs for the chunk of binding sites passed in argument
    
    # extract the data for this site and this lipid to avoid dimensionality nightmares
    reduced_occupancies = np.copy(occupancies)
    # first value is the sum of all single-frame occupancy
    acfs[0] = float(np.sum(reduced_occupancies))
    # if the ACF reached 0 already, no need to compute the rest
    for k in range(1,len(occupancies)):
        if not acfs[k]:
            acfs[k] = 0
            break
        # ACF at this point
        reduced_occupancies = np.convolve(reduced_occupancies, convolve_array, "valid") == 2
        acfs[k] = float(np.sum(reduced_occupancies))
        # normalization
    acfs = acfs/acfs[0]
    return acfs

def _autocorr_time(data):
    if len(data.shape) == 2:
        n_fr, n_lig = data.shape
        ac = np.array([ _autocorr( data[:,i], val=1 )  for i in range(n_lig) ] )
        ac_t = np.zeros((n_lig,3))
        for i in range(n_lig):
            ac_t[i], o = opt.curve_fit(_exp2,np.arange(n_fr),ac[i],p0=(0.1,0.9,0.9),bounds=( [0.,0.,0.],[np.inf, np.inf,1.0] ) )
            ac_t[i,0] = 1. / ac_t[i,0]
            ac_t[i,1] = 1. / ac_t[i,1]
    else:
        n_fr = len(data)
        ac = _autocorr( data )
        ac_t, o = opt.curve_fit(_exp2,np.arange(n_fr),ac,p0=(0.1,0.9,0.9),bounds=( [0.,0.,0.],[np.inf, np.inf,1.0] ) )
    return ac, ac_t

def _trim(data,val):
    data_out = np.copy(data)
    return np.trim_zeros(data_out - val) + val

def _gromos(D,cutoff,min_sz):
    M=np.copy(D)
    cent=[len(M)]
    elem = np.zeros(len(M))
    idx  = np.arange(len(M))
    cl_s = min_sz 
    n_cl = 1
    while (cl_s >= min_sz):
        loc_elem = np.zeros(len(M))
        nngb  = np.sum(M<cutoff, axis=1)
        if len(nngb) < min_sz :
            break
        loc_cent = np.argmax(nngb)
        clust = np.where(M[loc_cent]<cutoff)[0]
        cent.append(idx[loc_cent])
        cl_s = len(clust)
        elem[idx[clust]] = n_cl
        loc_elem[clust] = n_cl
        M=M[loc_elem==0][:,loc_elem==0]
        idx = idx[loc_elem==0]  
        n_cl += 1
    return elem,cent

    
def _calc_lifetime_old(data):
    n_fr , n_cl = data.shape
    LT = []
    for n in np.arange(n_cl):
        d0 = np.copy(data[:,n])
        m = np.max(d0)
        loc_t=[]
        for i in range(int(m)-1):
            d = (np.concatenate([[0],d0,[0]])>0).astype(int)
            up   = np.array([ i   for i in np.arange(1,n_fr+1) if (d[i] - d[i-1] ==  1) ])
            down = np.array([ i   for i in np.arange(2,n_fr+2) if (d[i] - d[i-1] == -1) ])
            loc_t= np.concatenate([loc_t,down-up])
            d0[d0!=0] -= 1
        LT.append(loc_t)
    return LT

def _calc_lifetime(data,n_cl):
    # 
    # This version assumes that you have one entry per ligand
    #
    if len(data.shape) > 1:
        n_fr , n_lig = data.shape
    else:
        n_fr = len(data)
        n_lig = 1
    LT = [ [ ] for j in range(int(n_cl)+1) ]
    for n in range(n_lig):
        for c in range(int(n_cl)+1):
            d0 = (data[:,n] == c).astype(int)
            d = (np.concatenate([[0],d0,[0]])).astype(int)
            up   = np.array([ i   for i in np.arange(1,n_fr+1) if (d[i] - d[i-1] ==  1) ])
            down = np.array([ i   for i in np.arange(2,n_fr+2) if (d[i] - d[i-1] == -1) ])
            LT[c].append(down-up)
    return LT

def _calc_lifetime_per_ligand_fast(data,n_cl):
    LT = [ [ ] for j in range(int(n_cl)+1) ]
    if len(data.shape) > 1:
        n_fr , n_lig = data.shape
        for n in range(n_lig):
            for c in range(int(n_cl)+1):
                d0 = (data[:,n] == c).astype(int)
                d = (np.concatenate([[0],d0,[0]])).astype(int)
                jump = d[1:]-d[:-1]
                up   = np.where(jump== 1)[0]
                down = np.where(jump==-1)[0]
                LT[c].append(down-up)
    else:
        n_fr = len(data)
        for c in range(int(n_cl)+1):
            d0 = (data == c).astype(int)
            d = (np.concatenate([[0],d0,[0]])).astype(int)
            jump = d[1:]-d[:-1]
            up   = np.where(jump== 1)[0]
            down = np.where(jump==-1)[0]
            LT[c].append(down-up)
    return LT

class phrases:
    def __init__(self,universe=None,receptor=None,ligand=None,cutoff=None,min_len=2):
        self.universe = universe
        self.receptor = receptor
        self.ligand = ligand
        self.cutoff = cutoff
        self.min_len = min_len
        self.phrases = []
        self.D = None

    def find_phrases_old(self,b,e,skip):
        self.phrases = []
        print ("Start reading trajectory")
        old_t=0.0
        for fr in self.universe.trajectory[b:e:skip]:
            self.dt = fr.time - old_t
            old_t = fr.time
            phrase=[]
            for lr in self.ligand.residues:
                p = list(np.unique(self.receptor.select_atoms("around %f global group grp" %(self.cutoff), grp=lr).resids).astype(int))
                if len(p) >= self.min_len:
                    phrase.append(p)
            self.phrases.append(phrase)
        print("Done!")

    def find_phrases(self,b=0,e=-1,skip=1):
        self.phrases = []
        print ("Start reading trajectory")
        old_t=0.0
        for fr in self.universe.trajectory[b:e:skip]:
            self.dt = fr.time - old_t
            old_t = fr.time
            phrase=[]
            # Could be used List comprehension:
            # p = [ list(np.unique(receptor.select_atoms("around %f global group grp" %(self.cutoff), grp=lr)) if len(list(np.unique(receptor.select_atoms("around %f global group grp" %(self.cutoff), grp=lr))) > self.min_len else [] for lr in self.ligand.residues ]
            #
            l_by_res = [ self.ligand.atoms[np.where(self.ligand.resids == r)[0]] for r in np.unique(self.ligand.resids) ] # create a list of AtomGroups divided by Ligand residues 
            for lr in l_by_res:
                p = list(np.unique(self.receptor.select_atoms("around %f global group grp" %(self.cutoff), grp=lr).resids).astype(int))
                if len(p) >= self.min_len:
                    phrase.append(p)
                else:
                    phrase.append([])
            self.phrases.append(phrase)
        print("Done!")

    def calc_dist(self):
        if self.phrases == [] :
            self.find_phrases()
        p_a = []
        self.nobs = 0
        for t in self.phrases:
            for f in t:
                self.nobs = self.nobs + 1 
                for a in f:
                    p_a.append(a)
        C = np.ones((max(p_a)+1,max(p_a)+1))
        for t in self.phrases:
            for f in t:
                for i in range(len(f)):
                    C[f[i],f[i]] += 1
                    for j in range(i+1,len(f)):
                        C[f[i],f[j]] += 1
                        C[f[j],f[i]] += 1
        self.D = -np.log(C/self.nobs)

    def find_cluster(self,cutoff):
        if type(self.D).__module__ != np.__name__:
            self.calc_dist()
        e, c = _gromos(self.D,cutoff,self.min_len)
        self.centroid = c
        self.labels   = e
        self.clusters = [ np.where(e==i)[0] for i in range(int(max(e)+1)) ]
        return
    
    def complete_linkage(self,cutoff=0.0,threshold=np.inf):
        if threshold < np.inf:
            t_idx = np.where(np.diagonal(self.D)<threshold)[0]
            D_t = self.D[t_idx][:,t_idx]
        else:
            D_t = self.D
        idx = np.triu_indices(D_t.shape[0],1)
        Z = sch.linkage(D_t[idx],method='complete')
        L = sch.fcluster(Z,t=cutoff,criterion='distance')
        if threshold < np.inf:
            self.clusters = [ t_idx[np.where(L == i )[0]] for i in np.arange(max(L))+1 if len(np.where(L==i)[0]) >= self.min_len ]
            self.labels = np.zeros(self.D.shape[0],dtype=int)
            self.linkage = Z
            for n,i in enumerate(t_idx):
                self.labels[i] = L[n]
        else:
            self.clusters = [ np.where(L == i )[0] for i in np.arange(max(L))+1 if len(np.where(L==i)[0]) >= self.min_len ]
            self.labels   = L
            self.linkage  = Z
        return

    def cluster_phrases(self,thresh=0.2):
        self.phrases_cl = [ [ _mj(p,self.clusters,thresh) for p in t ] for t in self.phrases]
        n_fr = len(self.phrases_cl)
        n_cl = int(max(self.labels))
        self.p_cl = np.zeros((n_fr,n_cl+1))
        for i in range(n_fr):
            for j in self.phrases_cl[i]:
                self.p_cl[i,j] += 1
        return
    
    def filter_phrases(self):
        frames_list = [ [] ]
        for n,c in enumerate(self.clusters):
            if n > 0:
                frames_list.append( np.where([ np.any(np.array([ set(c) == set(c) & set(i) for i in j])) for j in P.phrases ])[0] )
        return frames_list
        
    def life_time_old(self):
        self.LT=_calc_lifetime_old(self.p_cl,)
        return
    
    def life_time(self):
        self.phrases_cl = np.array(self.phrases_cl)
        self.LT=_calc_lifetime(self.phrases_cl,max(self.labels))
        return
    
    def autocorr_time(self):
        data = np.array(self.phrases_cl)
        clusters=self.clusters
        n_fr, n_lig = data.shape
        n_cl = len(clusters)
        self.ac = np.array([ [ _autocorr( (data[:,i]==j).astype(int) )  for i in range(n_lig) ] for j in range(n_cl) ])
        self.ac_t = np.zeros((n_cl,n_lig,3))
        for i in range(n_cl):
            for j in range(n_lig):
                self.ac_t[i,j], o = opt.curve_fit(_exp2,np.arange(n_fr),self.ac[i,j],p0=(0.1,0.9,0.9),bounds=( [0.,0.,0.],[np.inf, np.inf,1.0] ) )
                self.ac_t[i,j,0] = 1. / self.ac_t[i,j,0]
                self.ac_t[i,j,1] = 1. / self.ac_t[i,j,1]
        return
    
