from scipy.interpolate import interp2d, interp1d, RectBivariateSpline, LinearNDInterpolator
from scipy.ndimage.filters import gaussian_filter as gauss
from skimage.feature import peak_local_max
import matplotlib.pyplot as plt
import numpy as np
import math as m
import sys

def _read_dim(f,rd_names=True):
    """
        Reads a PLUMED GRID or SUMHILLS file, counting the number of 
        dimensions and nthe numebr of bins on each dimension.
        
        Arguments
        ---------
        f : str
            A Plumed GRID or SUMHILLS file name. 
        
        Returns
        -------
        out : a list of integers, the number of bins on each dimension of the GRID
        
    """
    out=[]
    if rd_names:
        names=[]
    with open(f,'r') as fi:
        for line in fi:
            if line[:2] != "#!":
                break
            if line[7:12] == "nbins":
                out.append(int(line.split()[-1]))
                if rd_names:
                    names.append(line.split()[-2][6:])
    names.reverse()
    out.reverse()
    return out,names

def read_energy_field(file_in,g_grad=True,**kargs):
    shape,names=_read_dim(file_in)
    data =np.loadtxt(file_in)
    ndim = len(shape)
    FIELD = [ data[:,ndim-i-1].reshape(shape) for i in range(ndim) ]
    E = data[:,ndim].reshape(shape)
    return energy_field(energy=E,field=FIELD,notes=str(file_in),names=names,**kargs)

class hills:
    def __init__(self, time, coord, sigmas, heights, cv_names, biasf=None, ranges=None):
        self.time = time
        self.coord = coord
        self.sigmas = sigmas
        self.heights = heights
        self.biasf = biasf
        self.ranges = ranges
        self.cv_names=cv_names
    
    
class energy_field:
    def __init__ (self, energy=None, gradient=None, field=None, g_grad=True, notes=None, names=None, interpolate=False):
        self.energy   = energy
        self.energy_g = gauss(self.energy,1.0,mode='constant')
        self.shape    = self.energy.shape
        self.ndim     = len(self.shape)
        self.idx  = np.meshgrid(*[range(self.shape[d]) for d in range(self.ndim)],indexing='ij' )
        if field == None:
            self.field = self.idx
        else:
            self.field = field
        if np.all(energy != None):
            if np.all(gradient == None):
                if g_grad:
                    self.gradient= np.array(np.gradient(self.energy_g))
                else:
                    self.gradient= np.array(np.gradient(self.energy))
            else:
                self.gradient = gradient
            self.idxs = np.vstack([i.ravel() for i in self.idx]).transpose()
            if interpolate:
                self._do_interpolate(kind="Energy")
                self._do_interpolate(kind="Gradient")
            else:
                self.interp_E = None
                self.interp_grad = None
            self.slc =[ [slice(0,1,1)] * len(self.shape) for i in range(self.ndim) ]
            self.field_interp = [] # FUNC[axis](index) -> value
            self.field_invert = [] # FUNC[axis](value) -> index
            self.axis=[]
            for i in range(self.ndim):
                self.slc[i][i] = slice(None)
                self.axis.append(self.field[i][self.slc[i]].ravel())
                self.field_interp.append(interp1d(range(self.shape[i]),self.axis[i],bounds_error=False,fill_value='extrapolate'))
                self.field_invert.append(interp1d(self.axis[i],range(self.shape[i]),bounds_error=False,fill_value='extrapolate'))
            self.dR = np.array([ np.average(np.diff(self.axis[i])) for i in range(self.ndim) ])
        self.notes=notes
        self.names=names

    def __setattr__(self,name,value):
        if name=='energy':
            self.__dict__[name]=np.array(value,dtype=float)
        elif name=='gradient':
            self.__dict__[name]=np.array(value,dtype=float)
        elif name=='field':
            if len(value)!=len(self.energy.shape):
                print "field has a number of elements [%d] not consistent with energy shape [%d]" %(len(value),len(self.energy.shape))
                raise AttributeError
            else:
                temp_list=[]
                for  n,f in enumerate(value):
                    if f.shape!=self.energy.shape:
                        print "One field [%d] has a dimension [%s] not consistent with energy [%s]" %(n,f.shape,self.energy.shape)
                        raise AttributeError
                    else:
                        temp_list.append(f)
                self.__dict__[name]=temp_list
        else:
            self.__dict__[name]=value

    def __str__(self):
        if self.isdefined:
            out = "Energy Field [%s] \n" %(self.notes)
        return out
    
    def _do_interpolate(self,kind=None,force=False):
        if force == True:
            if kind == "Energy":
                self.interp_E = LinearNDInterpolator(self.idxs,self.energy.ravel())
            elif kind == "Gradient":
                self.interp_grad = LinearNDInterpolator(self.idxs,self.gradient.reshape([self.ndim,len(self.idxs)]).swapaxes(0,1))
        else:
            if kind == "Energy":
                if self.interp_E != None:
                    print(" Energy interpolator already present.")
                    print(" If you want to force the calculation use the force = True option")
                else:
                    self.interp_E = LinearNDInterpolator(self.idxs,self.energy.ravel())
            if kind == "Gradient":
                if self.interp_grad != None:
                    print(" Gradient interpolator already present.")
                    print(" If you want to force the calculation use the force = True option")
                else:
                    self.interp_grad = LinearNDInterpolator(self.idxs,self.gradient.reshape([self.ndim,len(self.idxs)]).swapaxes(0,1))
        return
        
    def idx2point(self,idx):
        """
            Takes as a input an array of integer indexes and 
            returns the corresponding point in the CV space. 
        """
        if len(idx) != len(self.shape):
            print ("Number of dimension inconsistent")
            raise ValueError
            return
        p = tuple( np.int(np.rint(i)) for i in idx )
        return np.array([self.field[d][p] for d in range(self.ndim) ])
    
    def interp_idx2point(self,idx):
        """
            Convert coordinates from a continous coordinate in the indexes space
            to coordinates in the CV.
            It is used to convert interpolated index values to CV values.
             e.g. : 
                For semplicity immagine a 1D system.
                if the argument is idx = 5.8 the result is a corresponding 
                linear inerpolation of CV vales self.field[5] and self.field[6]
                
                  FUNC(5.8) = self.field[5] + 0.8 * (self.field[6] - self.field[5])
                  
        """
        if len(idx) != len(self.shape):
            print ("Number of dimension inconsistent")
            raise ValueError
            return
        p = tuple( np.int(np.rint(i)) for i in idx )
        return np.array([ self.field_interp[i](p) for i in range(self.ndim) ])
    
    def point_to_idx(self,p):
        """
            Takes as a input a point in the CV space and returns an array 
            containing the indexes correspoinding to the neares point on the GRID. 
        """
        if len(p) != len(self.shape):
            print ("Number of dimension inconsistent")
            raise ValueError
            return
        return np.array([ min(max(np.rint(self.field_invert[i](p[i])),0),self.shape[i]-1) for i in range(self.ndim) ],dtype=int)
    
    def search_min(self,E_tr=0.0,min_dist=5):
        coordinates = peak_local_max(-self.energy, min_distance=10,threshold_abs=abs(E_tr),exclude_border=False)
        P_min=np.array([ self.idx2point(i) for i in coordinates ])
        return coordinates, P_min
            
class neb:
    def __init__( self, npoint=3, start=np.zeros(2), end=np.zeros(2), energy=energy_field()):
        if np.array(start).shape[0]!=np.array(end).shape[0]:
            print "start and end point have inconsistent dimensions [%d,%d]" %(start.shape[0],end.shape[0])
            raise AttributeError
        else:
            self.start = start
            self.end   = end
        self.ndim   = len(self.start)
        self.energy = energy
        self.npoint = npoint
        self.tau_p  = np.zeros([self.npoint,self.ndim])
        self.tng    = np.zeros([self.npoint,self.ndim])
        self.step   = (self.end-self.start)/(self.npoint-1)
        self.neb    = np.array([ self.start + i * self.step for i in range(self.npoint) ])
        self.neb_E  = np.zeros([self.npoint])
        self.neb_G  = np.zeros([self.npoint,self.ndim])
        self.neb_F  = np.zeros([self.npoint,self.ndim])
        self.calc_grad()
        self.do_elastic()
        

    def __str__(self):
        return " Nudged-elastic band\n Emax   : %20.10f\n Eref   : %20.10f\n Points : %20d\n Kmax   : %20.10f\n dK     : %20.10f" %(self.Emax, self.E_ref, self.npoint, self.kmax, self.dk)

    def __iter__(self):
        return self.neb

    def __len__(self):
        return self.npoint

    def __setattr__(self,name,value):
        if name=='energy':
            if self!=None:
                if not isinstance(value,energy_field):
                    print "Energy should be passed as an energy_field instance!"
                    raise AttributeError
                elif len(value.energy.shape)!=self.ndim:
                    print "The energy field has a dimensionality [%d] different from the nab [%d]" %(value.energy.shape[0,self.dim])
                    raise AttributeError
                else:
                    self.__dict__[name]=value
            else:
                self.__dict__[name]=value
        elif ((name=='start')|(name=='end')):
            self.__dict__[name]=np.array(value,dtype=float)
        else:
            self.__dict__[name]=value

    def do_neb(self):
        self.neb = np.array([ self.start + i * self.step for i in range(self.npoint) ])

    def do_elastic(self):
        self.E_ref  = max(self.neb_E[0],self.neb_E[-1])
        self.Emin   = min(self.neb_E[0],self.neb_E[-1])
        self.spring = np.zeros(self.npoint)
        self.ts     = np.argmax(self.neb_E)
        self.Emax   = self.neb_E[self.ts]
        self.kmax   = (self.Emax - self.E_ref)
        self.dk     = self.kmax / 2.0
        return

    def check_energy(self):
        if self.energy.isdefined :
            print "Energy Defined"
        else:
            print "Energy not Defined"
        return self.energy.isdefined
        
    def F_E(self,*args,**kargs):
        return self.energy.interp_E(*args,**kargs)

    def F_G(self,*args,**kargs):
        return self.energy.interp_grad(*args,**kargs)[0]

    def calc_grad(self):
        for n in range( self.npoint ) :
            try:
                self.neb_E[n] = self.F_E(self.neb[n])
                self.neb_G[n] = self.F_G(self.neb[n])
            except:
                print "\n Error at point %d " %(n)
                print self.neb[n], self.F_E(self.neb[n]), self.F_G(self.neb[n])[0]
        return

    def max_curv(self,eps=0.5):
        c=np.zeros(self.npoint)
        for n in range(1,self.npoint-1):
            fE = self.F_E
            fG = self.F_G
            d  = self.neb_G[n]/np.linalg.norm(self.neb_G[n])
            ep = self.neb[n]+eps*d
            em = self.neb[n]-eps*d
            e1  = self.neb_G[n]
            f2  = fG(ep)-fG(em) / (2*eps)
            try:
                e2  = f2 - np.dot(f2, e1) * e1
            except:
                print f2
                print e1
                raise
            c[n] = np.linalg.norm(self.neb_G[n]) / np.dot(self.neb_G[n],f2)
        return np.amax(c)

    def calc_force(self):
        k=1.0
        n_max=np.argmax(self.neb_E)
        if (len(self.neb_E)==len(self.neb)):
            self.neb_F  = np.zeros([self.npoint,self.ndim])
            l=+1
            self.tau_p[0]  = self.neb[1] - self.neb[0]
            self.tau_p[-1] = self.neb[-1] - self.neb[-2]
            for i in range(1,self.npoint-1):
                self.tau_p[i] = ( self.neb[i+1] - self.neb[i] ) / 2.
            self.tau = np.linalg.norm(self.tau_p, axis=1)
            self.tng = np.array([ l/d for l,d in zip(self.tau_p,self.tau)])
            # F_p = -self.neb_G + np.transpose(np.multiply(np.transpose(self.tng),np.einsum('...i,...i',self.neb_G,self.tng)))
            F_p           = np.array( [ -self.neb_G[i] + np.dot(self.neb_G[i],self.tng[i]) * self.tng[i] for i in range(1,self.npoint-1) ])
            # F_t = k * np.transpose(np.multiply(np.transpose(self.tng),self.tau))
            F_t           = np.array( [ k * self.tau[i] * self.tng[i] for i in range(1,self.npoint-1) ])
            self.neb_F[1:-1] = F_t + F_p
            self.neb_F[n_max] = self.neb_F[n_max] - 2 * np.dot(self.neb_F[n_max],self.tng[n_max]) * self.tng[n_max]
        return
        
    def calc_force_2(self):
        k=1.0
        n_max=np.argmax(self.neb_E)
        if (len(self.neb_E)==len(self.neb)):
            self.neb_F  = np.zeros([self.npoint,self.ndim])
            l=+1
            self.tau_p[0]  = self.neb[1] - self.neb[0]
            self.tau_p[-1] = self.neb[-1] - self.neb[-2]
            for i in range(1,self.npoint-1):
                self.tau_p[i] = ( self.neb[i+1] - self.neb[i] ) / 2.
            self.tau = np.linalg.norm(self.tau_p, axis=1)
            self.tng = np.array([ l/d for l,d in zip(self.tau_p,self.tau)])
            F_p = -self.neb_G + np.transpose(np.multiply(np.transpose(self.tng),np.einsum('...i,...i',self.neb_G,self.tng)))
            #F_p           = np.array( [ -self.neb_G[i] + np.dot(self.neb_G[i],self.tng[i]) * self.tng[i] for i in range(1,self.npoint-1) ])
            F_t = k * np.transpose(np.multiply(np.transpose(self.tng),self.tau))
            #F_t           = np.array( [ k * self.tau[i] * self.tng[i] for i in range(1,self.npoint-1) ])
            self.neb_F[1:-1] = F_t[1:-1] + F_p[1:-1]
            self.neb_F[n_max] = self.neb_F[n_max] - 2 * np.dot(self.neb_F[n_max],self.tng[n_max]) * self.tng[n_max]
        return
    
