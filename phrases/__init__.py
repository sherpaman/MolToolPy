
class phrases(universe,receptor,ligand,cutoff,min_len):
    self.universe = universe
    self.receptor = receptor
    self.ligand = ligand
    self.cutoff = cutoff
    self.min_len = min_len
    self.phrases = []
    self.D = None

    def find_phrases(self,b,e,skip):
        self.phrases = []
        print ("Start reading trajectory")
        for fr in u.trajectory[b:e:skip]:
            dt = fr.time - old_t
            old_t = fr.time
            phrase=[]
            for lr in self.ligand.residues:
                p = list(np.unique(self.receptor.select_atoms("around %f global group grp" %(cutoff), grp=lr).resids).astype(int))
                if len(p) > min_phrase_len:
                    phrase.append(p)
            self.phrases.append(phrase)
        print("Done!")

    def calc_dist(self)
        p_a = []
        self.nobs = 0
        for t in self.phrases:
            for f in t:
                self.nobs = nobs + 1 
                for a in f:
                    p_a.append(a)
        
        C = np.zeros((max(p_a),max(p_a)))
        for t in self.phrases:
            for f in t:
                for i in range(len(f)):
                    C[f[i]-1,f[i]-1] += 1
                    for j in range(i+1,len(f)):
                        C[f[i]-1,f[j]-1] += 1
                        C[f[j]-1,f[i]-1] += 1
        
        self.D = - np.log2(C/nobs)

    def find_cluster(self):
        e, c = self._gromos(self.D,cutoff,min_sz=self.min_len)
        self.labels   = e
        self.clusters = [ np.where(e==i)[0] for i in range(int(max(e)+1)) ]

    def cluster_phrases(self,threshold=0.5):
        self.phrases_cl = [ [ self._mj(p,clusters,threshold) for p in t ] for t in self.phrases]
        n_fr = len(self.phrases_cl)
        n_cl = int(max(self.labels))
        self.p_cl = np.zeros((n_fr,n_cl+1))
        for i in range(n_fr):
            for j in self.phrases_cl[i]:
                self.p_cl[i,j] += 1

    def life_time(self):
        self.life_time=self._calc_lifetime(self.p_cl)

    def _mj(a,c,t): 
        list_j = [ float(len(set(a) & set(i))) /  float(len(set(a) | set(i))) for i in c ] # Jaccad Similarity
        if max(list_j) < t :
            return 0
        else:
            return np.argmax(list_j)

    def _gromos(D,cutoff,min_sz=self.min_len):
        M=np.copy(D)
        cent=[len(M)]
        elem = np.zeros(len(M))
        idx  = np.arange(len(M))
        cl_s = min_sz 
        n_cl = 1
        while (cl_s >= min_sz):
            loc_elem = np.zeros(len(M))
            nngb  = np.sum(M<c, axis=1)
            if len(nngb) < min_sz :
                break
            loc_cent = np.argmax(nngb)
            clust = np.where(M[loc_cent]<c)[0]
            cent.append(idx[loc_cent])
            cl_s = len(clust)
            elem[idx[clust]] = n_cl
            loc_elem[clust] = n_cl
            M=M[loc_elem==0][:,loc_elem==0]
            idx = idx[loc_elem==0]  
            n_cl += 1
        return elem,cent

        
    def _calc_lifetime(data):
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
    
