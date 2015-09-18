import numpy as np
import math as m
import scipy.special as sp
import copy
from mi import mi


class TimeSer:
        def __init__(self,data,n_data,dim,nbins=12,bins=None,prob=None,reshape=True,frame_row=True,dtype=float):
                self.n_data = int(n_data)
                self.dtype  = dtype
                self        = np.array(data,dtype=dtype)
                self.dim    = int(dim)
                self.rep    = self._calc_rep()
                self.nbins  = nbins
                self.min    = np.min(data)
                self.max    = np.max(data)
                #self.bins   = [ np.linspace(np.min(self.data),np.max(self.data),nbins+1) ] # IT MUST BE A *LIST* OF *1-DIM ARRAYS*
                self.bins   = bins
                self.reshape= reshape
                self.prod=np.ones(dim+1,dtype=int)
                self.frame_row=frame_row
                self._shape(force=self.reshape)
                self.prob=prob
                if self.prob == None:
                        self.prob_av=False
                else:
                        if self._check_prob_shape():
                                self.prob_av=True
                self.entropy=None
                self.entropy_av=False
        
        def _calc_rep(self):
                p=1
                for i in range(len(self.data.shape)):
                        p*=self.data.shape[i]
                #print " Found {0:4d} replica(s)".format(p / (self.n_data * self.dim))
                return int(p / (self.n_data * self.dim))

        def _shape(self,force=False):
                #if self.reshape == False:
                if  ( self._check_shape() & (not force) ) :
                        #print "Data already reshaped"
                        return
                self.reshape = False
                proper_shape_data = np.zeros((self.rep,self.dim,self.n_data))
                #row, lines = self.data.shape
                if ( (self.dim == 1)  & (self.rep == 1) ):
                        for k in np.arange(self.n_data):
                                # IN_DATA [ N_sample, N_rep ]
                                # DATA [ N_rep, DIM, N_sample ]
                                proper_shape_data [0,0,k] = self.data[ k ]
                else:
                        if self.frame_row:
                                # IN_DATA [ N_sample, N_rep ]
                                for i in np.arange(self.rep):
                                        for j in np.arange(self.dim):
                                                for k in np.arange(self.n_data):
                                                        # OUT_DATA [ N_rep, DIM, N_sample ]
                                                        proper_shape_data [i,j,k] = self.data[ k , i * self.dim + j ]
                        else:
                                # IN_DATA [ N_rep, N_sample ]
                                for i in np.arange(self.rep):
                                        for j in np.arange(self.dim):
                                                for k in np.arange(self.n_data):
                                                        # OUT_DATA [ N_rep, DIM, N_sample ]
                                                        proper_shape_data [i,j,k] = self.data[ i * self.dim + j, k ]
                self.data = proper_shape_data
                return

        def _check_shape(self):
                return ( self.data.shape == ( self.rep, self.dim, self.n_data ) )
                
        def _check_prob_shape(self):
                return ( self.prob.shape == tuple([self.rep]) + tuple( self.nbins ) )
        
        def calc_bins(self):
                self.bins = []
                for d in np.arange(self.dim):
                        self.bins.append(np.linspace(np.min(self.data[:,d,:]),np.max(self.data[:,d,:]),self.nbins+1))
                return
                                
        def calc_prob(self):
                if self.reshape == True:
                        self._proper_shape()
                # PROB ( N_rep, DIM_1, DIM_2 ... )
                self.prob = np.zeros( np.hstack((self.rep, self.nbins)) )
                for i in np.arange(self.rep):
                        DATA = self.data[i]
                        if np.any(self.bins == None):
                                self.calc_bins()
                        if self.dim == 1:
                                histo , null = np.histogram(DATA,bins=self.bins[0])
                                self.prob[i] = 1.0 * histo / self.n_data
                        else:
                                histo , null = np.histogramdd(np.transpose(DATA),bins=self.bins)
                                self.prob[i] = 1.0 * histo / self.n_data
                self.prob_av=True
                return
        
        def calc_entropy(self,force=False):
                if ( self.entropy_av & (not force) ):
                        print "Entropy already available"
                        return
                if not (self.prob_av & (not force)) :
                        self.calc_prob()
                self.entropy = np.zeros(self.rep)
                for i in np.arange(self.rep):
                        self.entropy[i] = -np.sum(sp.xlogy(self.prob[i],self.prob[i]))/m.log(2)
                self.entropy_av=True
                return

        def mutual_info_bootstrap(self,resample=100,joint=False):
                # Can be used to calculate a confidence interval for Mutual Information
                mi = np.zeros((self.rep,self.rep,resample))
                if joint:
                        for n in np.arange(resample):
                                c = np.random.choice(np.arange(self.n_data),self.n_data,replace=True)
                                sample = TimeSer(self.data[...,c],self.n_data,self.dim,self.nbins,reshape=False,dtype=self.dtype)
                                mi[:,:,n] = sample.mutual_info_for()[0]
                else:
                        for n in np.arange(resample):
                                d = []
                                for r in np.arange(self.rep):
                                        c = np.random.choice(np.arange(self.n_data),self.n_data,replace=True)
                                        d.append(self.data[r,:,c].reshape((self.dim,self.n_data)))
                                data = np.dstack(d).reshape((self.rep,self.dim,self.n_data))
                                sample = TimeSer(data,self.n_data,self.dim,self.nbins,reshape=False,dtype=self.dtype)
                                mi[:,:,n] = sample.mutual_info_for()[0]
                return mi
                                        
        
        def mutual_info(self):
                if not self.entropy_av:
                        self.calc_entropy()
                P_joint = np.zeros( np.hstack( (self.rep, self.rep, self.nbins, self.nbins) ) )
                E_joint = np.zeros((self.rep,self.rep))
                M       = np.zeros((self.rep,self.rep))
                total_step = ( self.rep * ( self.rep - 1 ) ) / 2
                #n=0
                for s in np.arange(self.rep - 1):
                        E_joint[s,s] = self.entropy[s] + self.entropy[s]
                        M[s,s]       = 0
                        for i in np.arange(self.nbins):
                                P_joint[s,s,i,i] = self.prob[s,i]
                        for o in np.arange(s+1,self.rep):
                                #n = n + 1
                                #print " Completed {0:6.3%}\r".format(float(n)/total_step)
                                DATA = np.transpose(np.vstack((self.data[s],self.data[o])))
                                histo , null = np.histogramdd(DATA,bins=(self.nbins,self.nbins)) # concatenate bins list with "+"
                                P_joint[s,o] = 1.*histo / self.n_data
                                E_joint[s,o] = -np.sum(sp.xlogy(P_joint[s,o],P_joint[s,o]))
                                E_joint[o,s] = E_joint[s,o]
                                M[s,o] = self.entropy[s] + self.entropy[o] - E_joint[s,o]
                                M[o,s] = M[s,o]
                E_joint[self.rep-1,self.rep-1] = self.entropy[self.rep-1] + self.entropy[self.rep-1]
                M[self.rep-1,self.rep-1]       = 0
                return M, E_joint, P_joint
        
        def mutual_info_for(self):
                # CALCULATE MUTUAL INFO OF EACH PAIR OF REPLICAS
                if not self.entropy_av:
                        self.calc_entropy()
                P_joint = np.zeros( np.hstack( (self.nbins,self.nbins,self.rep,self.rep) ) )
                E_joint = np.zeros((self.rep,self.rep))
                M       = np.zeros((self.rep,self.rep))
                if self.dtype == float :
                        M, E_joint, P_joint = mi.r_mutualinfo(np.transpose(self.data),self.entropy,self.bins,self.n_data,self.rep,self.nbins)
                elif self.dtype == int :
                        M, E_joint, P_joint = mi.i_mutualinfo(np.transpose(self.data),self.entropy,self.bins,self.n_data,self.rep,self.nbins)
                return np.transpose(M), np.transpose(E_joint), np.transpose(P_joint)
        
        def mutual_info_simple(self,nrep1=0):
                # CALCULATE MUTUAL INFO OF A SINGLE REPLICA AGAINST ALL THE OTHERS
                if nrep1 not in np.arange(self.rep):
                        print "you should indicate a nrep between 0 and %d" %(self.rep)
                        return None
                if not self.entropy_av:
                        self.calc_entropy()
                P_joint = np.zeros( np.hstack( (self.rep,self.nbins,self.nbins) ) )
                E_joint = np.zeros(self.rep)
                M       = np.zeros(self.rep)
                E_joint[nrep1] = 2 * self.entropy[nrep1]
                M[nrep1] = 0
                for s in np.arange(self.rep):
                        if s != nrep1:
                                DATA = np.transpose(np.vstack((self.data[nrep1],self.data[s])))
                                histo , null = np.histogramdd(DATA,bins=(self.nbins,self.nbins)) # concatenate bins list with "+"
                                P_joint[s] = 1.*histo / self.n_data
                                E_joint[s] = -np.sum(sp.xlogy(P_joint[s],P_joint[s]))
                                M[s] = self.entropy[nrep1] + self.entropy[s] - E_joint[s]
                return M, E_joint, P_joint

        def mutual_info_simple_for(self,nrep1=0):
                # CALCULATE MUTUAL INFO OF A SINGLE REPLICA AGAINST ALL THE OTHERS
                if nrep1 not in np.arange(self.rep):
                        print "you should indicate a nrep between 0 and %d" %(self.rep)
                        return None
                if not self.entropy_av:
                        self.calc_entropy()
                P_joint = np.zeros( np.hstack( (self.nbins,self.nbins,self.rep) ) )
                E_joint = np.zeros(self.rep)
                M       = np.zeros(self.rep)
                if self.dtype == float :
                        M, E_joint, P_joint = mi.r_mutualinfo_simp\
                                (np.transpose(self.data),self.entropy,self.bins,nrep1,self.n_data,self.rep,self.nbins) 
                elif self.dtype == int :
                        M, E_joint, P_joint = mi.i_mutualinfo_simp\
                                (np.transpose(self.data),self.entropy,self.bins,nrep1,self.n_data,self.rep,self.nbins)
                return np.transpose(M), np.transpose(E_joint), np.transpose(P_joint)
                
        def mutual_info_other(self,other):
                if self.n_data != other.n_data:
                        print "The Number of Observation in the two time series are different ({0:d} != {0:d})".format(self.n_data,other.n_data)
                        return None, None
                if not self.entropy_av:
                        self.calc_entropy()
                if not other.entropy_av:
                        other.calc_entropy()
                P_joint = np.zeros( np.hstack( (self.rep, other.rep, self.nbins, other.nbins) ) )
                E_joint = np.zeros((self.rep,other.rep))
                M       = np.zeros((self.rep,other.rep))
                total_step=self.rep*other.rep
                for s in np.arange(self.rep):
                        for o in np.arange(other.rep):
                                # DATA* -> Transposed[ Self.DIM+Other.DIM, N_sample ]
                                #print " Completed {0:6.3%}\r".format(float(s*o)/total_step)
                                DATA = np.transpose(np.vstack((self.data[s],other.data[o])))
                                histo , null = np.histogramdd(DATA,bins=self.bins+other.bins) # concatenate bins list with "+"
                                P_joint[s,o] = 1.*histo / self.n_data
                                E_joint[s,o] = -np.sum(sp.xlogy(P_joint[s,o],P_joint[s,o]))
                                M[s,o] = self.entropy[s] + other.entropy[o] - E_joint[s,o]
                return M, E_joint, P_joint
        
        def mutual_info_other_for(self,other):
                if self.n_data != other.n_data:
                        print "The Number of Observation in the two time series are different ({0:d} != {0:d})".format(self.n_data,other.n_data)
                        return None, None
                if not self.entropy_av:
                        self.calc_entropy()
                if not other.entropy_av:
                        other.calc_entropy()
                P_joint = np.zeros(np.hstack( (self.rep, other.rep, self.nbins, other.nbins) ) )
                E_joint = np.zeros((self.rep,other.rep))
                M       = np.zeros((self.rep,other.rep))
                d1 = np.transpose(self.data)
                d2 = np.transpose(other.data)

                if (self.dtype == int):
                        if (other.dtype == int):
                                M, E_joint, P_joint =  mi.i_mutualinfo_other(d1,d2,    \
                                                        self.entropy, other.entropy,   \
                                                        self.bins, other.bins,   \
                                                        self.n_data,                   \
                                                        self.rep,     other.rep,       \
                                                        self.nbins,   other.nbins)
                        elif (other.dtype == float):
                                M, E_joint, P_joint = mi.ir_mutualinfo_other(d1,d2,     \
                                                        self.entropy, other.entropy,   \
                                                        self.bins[0], other.bins[0],   \
                                                        self.n_data,                   \
                                                        self.rep,     other.rep,       \
                                                        self.nbins,   other.nbins)
                elif (self.dtype == float):
                        if (other.dtype == int):
                                M, E_joint, P_joint = mi.ri_mutualinfo_other(d1,d2,     \
                                                        self.entropy, other.entropy,   \
                                                        self.bins[0], other.bins[0],   \
                                                        self.n_data,                   \
                                                        self.rep,     other.rep,       \
                                                        self.nbins,   other.nbins)
                        elif (other.dtype == float):
                                M, E_joint, P_joint =  mi.r_mutualinfo_otherd(1,d2,     \
                                                        self.entropy, other.entropy,   \
                                                        self.bins[0], other.bins[0],   \
                                                        self.n_data,                   \
                                                        self.rep,     other.rep,       \
                                                        self.nbins,   other.nbins)
                return np.transpose(M), np.transpose(E_joint), np.transpose(P_joint)
        
        def digitized(self):
                other = copy.deepcopy(self)
                other.data = np.zeros([self.rep,self.dim,self.n_data],dtype=int)
                #values = np.zeros((self.dim,self.nbins[0]))
                #if self.dim == 1 :
                for r in np.arange(self.rep):
                        for d in np.arange(self.dim):
                                for k in np.arange(self.n_data):
                                        #print r,d,k, self[r,d,k]
                                        other.data[r,d,k] = int(np.max(np.where( self.bins[d][:-1] <= self[r,d,k] )))
                #else:
                        #for r in np.arange(self.rep):
                                #for d in np.arange(self.dim):
                                        #for k in np.arange(self.n_data):
                                                #other.data[r,d,k] = np.max(np.where( self.bins[d][:-1] <= self[r,d,k] ))
                return other
        
        def traj(self,time=2,nbins=None):
                if nbins == None:
                         #
                         # THIS IS THE OPTIMAL NUMBER OF BINS THAT SHOULD BE USED 
                         # TO RESOLVE ALL THE POSSIBLE TIME SERIES (TO BE CHECKED)
                         #
                         # SINCE THE PROCESS INVOLVES THE CALCULATION OF A JOINT
                         # PROBABILITY, A MATRIX OF DIM nbins*(nbins*(nbin^dim)) IS CREATED:
                         # BE CAREFUL SINCE THIS WILL QUICKLY SATURATE THE MEMORY!!
                         #
                         # EXAMPLE:
                         #  nbin =10; dim =1 ; time =2 ->  nbins = 10^2 
                         #                              -> nbins*nbins*nbin = 10^5 
                         #
                         #  nbin =10, dim =1 ; time =3 -> nbins = 10^3
                         #                             -> nbins*nbins*nbin = 10^7 !!!!
                         # 
                         #
                        nbins= int(self.nbins) ** int( self.dim * time )
                
                prod_t = np.ones(time+1,dtype=int)
                for i in np.arange(1,self.dim+1):
                        self.prod[i] = self.prod[i-1] * self.nbins
                for i in np.arange(1,time+1):
                                prod_t[i] =  prod_t[i-1] * self.prod[self.dim]
                otherk_nbins   = nbins
                otherk1_nbins  = nbins * ( self.nbins ** self.dim )
                otherk_bins    = [ np.linspace(0, otherk_nbins-1, otherk_nbins+1) ]
                otherk1_bins   = [ np.linspace(0,otherk1_nbins-1,otherk1_nbins+1) ]
                other_nd       = self.n_data-time
                other_k  = TimeSer(np.zeros((self.rep,1,other_nd)),other_nd,1,nbins=otherk_nbins ,bins=otherk_bins ,prob=None,reshape=False,dtype=int)
                other_k1 = TimeSer(np.zeros((self.rep,1,other_nd)),other_nd,1,nbins=otherk1_nbins,bins=otherk1_bins,prob=None,reshape=False,dtype=int)
                digit = self.digitized()
                for r in np.arange(self.rep):
                        for k in np.arange(self.n_data-time):
                                hash_num = 0
                                for t in np.arange(time):
                                        #
                                        # IN THIS IMPLEMENTATION THE MOST SIGNIFICATIVE "DIGIT" OF THE HASH_NUMBER
                                        # IS DETERMINED BY THE LAST (IN TIME) VALUE OF TIME SERIE. THIS REQUIRE
                                        # THAT THE NUMBER OF BINS IS DIFFERENT FOR THE O(K) AND O(K+1) PROBABILITY
                                        #
                                        # N.B. : THE DATA HAVE TO BE DIGITIZED
                                        #
                                        hash_num = hash_num + int(np.vdot(digit.data[r,:,k+t],self.prod[:-1])) * prod_t[t]
                                hash_num1 = hash_num + int(np.vdot(digit.data[r,:,k+time],self.prod[:-1])) * prod_t[time]
                                other_k.data[r,0,k]  = hash_num
                                other_k1.data[r,0,k] = hash_num1
                return other_k, other_k1
                
        
        def transfer_entropy_for(self,time=2,nbins=None):
                # CALCULATE THE TRANSFER ENTROPIES BETWEEN THE REPLICAS
                ok, ok1 = self.traj(time,nbins)
                #print "Calculating Mutual Info of traj(K)s"
                M,  E,  P  = ok.mutual_info_for()
                #print "Calculating Mutual Info of traj(K+1)s with respect of traj(K)s"
                # WHICH ONE IS CORRECT ?!
                M1, E1, P1 = ok.mutual_info_other_for(ok1)
                #M1, E1, P1 = ok1.mutual_info_other(ok)
                T = np.zeros((self.rep,self.rep))
                # T(I,J) = H(I) - H(J) - H(I_k+1, J_k) + H(I_k,J_k)
                #print "Calculating Transfer Entropy"
                for s in np.arange(self.rep):
                        for o in np.arange(self.rep):
                                T[s,o] = ok1.entropy[s] + ok.entropy[o] - M1[s,o] + M[s,o]
                #print "Done"
                return T, [ [M,  E,  P ], [M1, E1, P1 ] ] 
        
        
        
        def transfer_entropy(self,time=2,nbins=None):
                # CALCULATE THE TRANSFER ENTROPIES BETWEEN THE REPLICAS
                if nbins == None:
                         #
                         # THIS IS THE OPTIMAL NUMBER OF BINS THAT SHOULD BE USED 
                         # TO RESOLVE ALL THE POSSIBLE TIME SERIES (TO BE CHECKED)
                         #
                         # SINCE THE PROCESS INVOLVES THE CALCULATION OF A JOINT
                         # PROBABILITY, A MATRIX OF DIM nbins*(nbins*(nbin^dim)) IS CREATED:
                         # BE CAREFUL SINCE THIS WILL QUICKLY SATURATE THE MEMORY!!
                         #
                         # EXAMPLE:
                         #  nbin =10; dim =1 ; time =2 ->  nbins = 10^2 
                         #                              -> nbins*nbins*nbin = 10^5 
                         #
                         #  nbin =10, dim =1 ; time =3 -> nbins = 10^3
                         #                             -> nbins*nbins*nbin = 10^7 !!!!
                         # 
                         #
                        nbins= int(self.nbins) ** int( self.dim * time )

                ok, ok1 = self.traj(time,nbins)
                #print "Calculating Mutual Info of traj(K)s"
                M,  E,  P  = ok.mutual_info()
                #print "Calculating Mutual Info of traj(K+1)s with respect of traj(K)s"
                # WHICH ONE IS CORRECT ?!
                M1, E1, P1 = ok.mutual_info_other(ok1)
                #M1, E1, P1 = ok1.mutual_info_other(ok)
                T = np.zeros((self.rep,self.rep))
                # T(I,J) = H(I) + H(J) - H(I_k+1, J_k) + H(I_k,J_k)
                #print "Calculating Transfer Entropy"
                for s in np.arange(self.rep):
                        for o in np.arange(self.rep):
                                T[s,o] = ok1.entropy[s] + ok.entropy[o] - M1[s,o] + M[s,o]
                #print "Done"
                return T, [ [M,  E,  P ], [M1, E1, P1 ] ]
        
        def calc_direction(self,T):
                D = np.zeros((self.rep,self.rep))
                Tran = T[0]
                Rate = np.array([(T[1][0][1]-T[1][1][1])[i,i] for i in np.arange(self.rep) ])
                for i in np.arange(self.rep):
                        for j in np.arange(self.rep):
                                 D[i,j] = Tran[i,j]/Rate[i] - Tran[j,i]/Rate[j]
                return D
                                        
                                        
                
                
                
                                        
                                        
                        
                
