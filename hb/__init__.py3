import re
import sys
import gro
import xpm
import copy
import numpy
import scipy.stats
import scipy.optimize

#REGULAR EXPRESSION DEFINITIONS
re_nucleotide=re.compile('D[ATCG][0-9]+')
re_protein=re.compile('[A-Z]{3}[0-9]+')

re_N=re.compile("(D[ACTG])([1-9]{1}[0-9]*)([HNOP][A-Z,0-9]*)") # MATCHES NUCLEOTIDES
re_P=re.compile("([ACGHILMNPSTV][A-Z]{2})([1-9]{1}[0-9]*)([HNOP][A-Z,0-9]*)") # MATCHES PROTEIN RESIDUES
re_R=re.compile("(D[ACTG]|[ACGHILMNPSTV][A-Z]{2})([1-9]{1}[0-9]*)([HNOP][A-Z,0-9]*)") # MATCHES BOTH
re_R_sim=re.compile("(D[ACTG]|[ACGHILMNPSTV][A-Z]{2})([1-9]{1}[0-9]*)")

def _autocorr(data):
    mu     = numpy.average(data)
    sigma  = numpy.var(data)
    result = numpy.correlate(data-mu, data-mu, mode='full')/(sigma**2)
    norm   = numpy.arange(result.size/2+1,0,-1)
    return result[result.size/2:]/norm

def _exp2(x,l1,l2,a):
    return a * numpy.exp(-l1*x) + (1-a) * numpy.exp(-l2*x)

def uniq_ref(l):
    """
    Auxiliary Function that accept as a signle argument a list and returns: 
        
        - A list containing the unique elements of the input list,
        - A list of index pointing to the unique elements in the input list,
        - A list of references from the unique elements to all their occurences in the input list.
        
    """
    out=[]; out_id = []; ref=[];
    for n,i in enumerate(l):
        if i not in out:
            out.append(i)
            out_id.append(n)
            ref.append([n])
        else:
            ref[out.index(i)].append(n)
    return out, out_id, ref

def is_water(name):
    return (name[:3]=="WAT")|(name[:3]=="SOL")|(name[:3]=="HOH")

def is_nucleotide(name):
    return (name[0]=="D")

def residue_split(name,m):
    return m.find_atom_res(name)

def uniq(list):
    checked=[]
    for i in list:
        if i not in checked:
            checked.append(i)
    return checked

def fscore(p1,p2,t1,t2):
    n1=int(p1*t1/100.0)+1
    n2=int(p2*t2/100.0)+1
    odds, pvalue = scipy.stats.fisher_exact([ [ n1 , (t1 - n1)+1 ], [ n2, (t2 - n2)+1 ] ])
    #print "F-SCORE [ [ %d , %d ], [ %d, %d] ] = %f" %(n1 , (t1 - n1)+1 ,  n2 , (t2 - n2)+1, pvalue )
    return pvalue

def bootstrap(data, num_samples, statistic, alpha):
    n       = len(data)
    idx     = numpy.random.randint(0, n, (num_samples, n))
    samples = data[idx]
    stat    = numpy.sort(statistic(samples, 1))
    return numpy.array([stat[int((alpha/2.0)*num_samples)], stat[int((1-alpha/2.0)*num_samples)]])

class HBond:
    def __init__(self,don='res1',don_atom='',acc='res2',acc_atom='',h='',perc=0.0,perc_ci=numpy.zeros(2),perc_c=0.0,mediated=0.0,both=0.0,nfr=0):
        self.acc      = acc
        self.acc_atom = acc_atom
        self.don      = don
        self.don_atom = don_atom
        self.perc     = perc     # PERCENTAGE
        self.perc_ci  = perc_ci  # PERCENTAGE CONFIDENCE INTERVAL
        self.perc_c   = perc_c   # CONFIDENCE
        self.mediated = mediated
        self.both     = both
        self.nfr      = nfr
    
    def __str__(self):
        return "%11s %11s %6.2f %6.2f %6.2f %d" %(self.don+self.don_atom,self.acc+self.acc_atom,self.perc,0.0,0.0,self.nfr)
    
    def __repr__(self):
        return str(self)+"\n"
    
    def renum(self,func):
        self.don[1] = func(self.don[1])
        self.acc[1] = func(self.acc[1])
    
    def is_same_couple(self,other):
        return (self.acc==other.acc)&(self.don==other.don)
    
    def fscore(self,other):
        n1=int(self.perc*self.nfr/100.0)
        n2=int(other.perc*other.nfr/100.0)
        odds, pvalue = scipy.stats.fisher_exact([ [ n1 , perc*self.nfr - n1 ],[ n2 , perc*other.nfr - n2 ] ])
        return pvalue
    

class HBonds: 
    def __init__(self,name='New List',mol=None,red=False,log=None,xpm=None,perc=None,conf=None):
        self.name   = name
        self.log    = log
        self.xpm    = xpm
        #self.mol    = mol
        self.hblist = []
        self.red_hb = []
        self.ref_hb = []
        self.nbonds = 0 # number of HBonds (non-redundant hb if merge is performed)
        self.nrhb   = 0 # number of redundand HBonds
        self.nfr    = 0
        self.red    = red
        self.conf   = conf # confidence interval for bootstrap analysis
        if self.log != None:
            self.read_log(self.log,self.red)
            if self.xpm != None:
                check_nhb, self.nfr = numpy.shape(self.xpm.array)
                if ( check_nhb != self.nrhb ):
                    print("XPM and Log file do not correspond (%d /= %d)" %(check_nhb,self.nrhb))
                    raise ValueError
                print("FOUND %4d     REDUNDANT HB" % self.nrhb)
                if self.red == False :
                    print("      %4d NON-REDUNDAT HB" % self.nbonds)
                    self.nr_xpm = self.merge_xpm()
                if self.conf != None:
                    self.perc_bootstrap(conf=self.conf)
                else:
                    self.calc_perc()
        else:
            if perc!=None:
                self.read_perc(perc,self.red)

    def __setattr__(self,name,value):
        if (name=='mol') :
            if isinstance(value,gro.Molecule):
                print("mol Object found!")
                self.__dict__[name] = value
            elif type(value) == str :
                print("Start reading mol file")
                self.__dict__[name] = gro.read_gro(value)
                print("Done!")
            else:
                if value != None :
                    print("Wrong Type mol: {0:s}".format(type(value)))
                self.__dict__[name]=None
        elif (name=='xpm') :
            if isinstance(value,xpm.Xpm):
                print("xpm Object found!")
                self.__dict__[name] = value
            elif type(value) == str :
                print("Start reading xpm file")
                self.__dict__[name] = xpm.read_xpm(value)
                print("Done!")
            else:
                if value != None :
                    print("Wrong Type xpm: {0:s}".format(type(value)))
                self.__dict__[name]=None
        else:
            self.__dict__[name]=value
          
    def __iter__(self):
        return iter(self.hblist)
    
    def __getitem__(self,key):
        return self.hblist[key]
    
    def __len__(self):
        return len(self.hblist)
        
    def __str__(self):
        out=self.name+'\n'
        for i in self.hblist:
            out = out + str(i)+'\n'
        return out
    
    def __repr__(self):
        return str(self)
    
    def merge_xpm(self):
        nr_xpm       = copy.deepcopy(self.xpm)
        nr_xpm.rows  = self.nbonds
        nr_xpm.cols  = self.nfr
        nr_xpm.yaxis = numpy.arange(nr_xpm.rows)
        nr_xpm.ylabel= 'Non-Redundant Hydrogen Bond Index'
        nr_xpm.title = 'Non-Redundant Hydrogen Bond Existence Map'
        nr_xpm.array = numpy.zeros((nr_xpm.rows,nr_xpm.cols)) # COLS AND ROWS IN XPM HAVE DIFFERENT MEANING THAN IN THE ACTUAL ARRAY!!!
        for i in numpy.arange(self.nbonds):
            nr_xpm.array[i,:] = numpy.max(self.xpm.array[self.ref_hb[i],:],axis=0)
        return nr_xpm
    
    def renum(self,func):
        for hb in self.hblist:
            hb.renum(func)
    
    def refresh(self):
        self.nbonds=len(self.hblist)
        self.create_names()
        
    def write_file(self,filename):
        out=open(filename,'w')
        out.write(self.name+'\n')
        for i in self.hblist:
            out.write(str(i)+'\n')
        out.close()
    
    def read_log(self,file_in,red=False):
        """
        Populate the HBonds object with the HBond found in a
        file .log generated by g_HBond [-g] merging the bonds
        for the same residue couple 
            
                        self.red_hb        = list of redundant h-bond   
                        self.hblist        = list of non redundant h-bonds
                        self.ref_hb        = list of redundand h-bonds, stored as elements of self.red_hb, 
                                             that has been merged as element of self.ref_hb
        
        """
        fi=open(file_in,'r')
        raw=fi.readlines()
        fi.close()
        self.hblist=[]
        self.red_hb=[]
        self.ref_hb=[]
        self.nrhb = len(raw)-1
        self.red = red
        if self.red==False:
            l_sim=[]
            hblist = []
            for l in raw[1:]:
                l_sim.append([ list(i)[0]+list(i)[1] for n,i in enumerate(re_R.findall(l)) if (n != 1) ])
                hblist.append([ list(i) for i in re_R.findall(l)])
                d_r , a_r = l_sim[-1][0] , l_sim[-1][1]
                d_a , a_a = hblist[-1][0][2], hblist[-1][2][2]
                self.red_hb.append(HBond(don=d_r, don_atom=d_a, acc=a_r, acc_atom=a_a, perc=0.0, mediated=0.0, both=0.0, nfr=self.nfr))
            l_hb,l_hb_id,l_ref = uniq_ref(l_sim)
            for i in l_hb_id:
                don_res , acc_res = l_sim[i][0],l_sim[i][1]
                self.hblist.append(HBond(don=don_res, acc=acc_res, perc=0.0, mediated=0.0, both=0.0, nfr=self.nfr))
            self.ref_hb = l_ref
            self.nbonds = len(l_hb)
        else:
            for l in raw[1:]:
                hb = re_R.findall(l)
                d_r , a_r = hb[0][0]+hb[0][1] , hb[2][0]+hb[2][1]
                d_a , a_a = hb[0][2], hb[2][2]
                self.hblist.append(HBond(don=d_r, don_atom=d_a, acc=a_r, acc_atom=a_a, perc=0.0, mediated=0.0, both=0.0, nfr=self.nfr))
            self.nbonds = self.nrhb
        self.red = red
                
        
    def read_log_old(self,file_in,file_gro,red=False):
        """
        Populate the HBonds object with the HBond found in a
        file .log generated by g_HBond [-g] merging the bonds
        for the same residue couple 
            
                        self.red_hb        = list of redundant h-bond   
                        self.hblist        = list of non redundant h-bonds
                        self.ref_hb        = list of redundand h-bonds, stored as elements of self.red_hb, 
                                             that has been merged as element of self.ref_hb
        
        """
        fi=open(file_in,'r')
        self.mol = file_gro
        print(self.mol)
        self.name   = file_in
        list_hb     = []
        self.nbonds = 0
        self.nrhb   = 0
        raw_fi      = fi.readline()
        while raw_fi:
            #
            # Exclude the comment lines and those with less than 3 field
            #
            if (raw_fi.split()[0] != '#')&(len(raw_fi.split())>2):
                donor   = raw_fi.split()[0]
                hydrogen= raw_fi.split()[1]
                acceptor= raw_fi.split()[2]
                #
                # Extract the residue name, including the number
                #
                don_res=self.mol.find_atom_res(donor)
                acc_res=self.mol.find_atom_res(acceptor)
                                #
                #if is_nucleotide(donor):
                #   don_res=re_nucleotide.findall(donor)[0]
                #else:
                #   don_res=re_protein.findall(donor)[0]
                #if is_nucleotide(acceptor):
                #   acc_res=re_nucleotide.findall(acceptor)[0]
                #else:
                #   acc_res=re_protein.findall(acceptor)[0]
                                #
                don_atom = self.mol.find_atom_name(donor)
                acc_atom = self.mol.find_atom_name(acceptor)
                
                if not red:
                    self.red_hb.append(HBond(don=don_res, acc=acc_res, perc=0.0, mediated=0.0, both=0.0,nfr=self.nfr))
                    new_hb   = [don_res[0]+"_"+str(don_res[1]),acc_res[0]+"_"+str(acc_res[1])]
                    if new_hb not in list_hb:
                        list_hb.append(new_hb)
                        self.ref_hb.append([self.nrhb])
                        self.hblist.append(HBond(don=don_res, acc=acc_res, perc=0.0, mediated=0.0, both=0.0,nfr=self.nfr))
                        self.nbonds += 1
                    else:
                        self.ref_hb[list_hb.index(new_hb)].append(self.nrhb)
                    self.nrhb += 1
                else:
                    self.hblist.append(HBond(don=don_res,don_atom=don_atom, acc=acc_res,acc_atom=acc_atom, perc=0.0, mediated=0.0, both=0.0,nfr=self.nfr))   
            raw_fi=fi.readline()
        self.red = red
        fi.close()

    def calc_perc(self):
        self.nrhb, self.nfr = self.xpm.array.shape
        if self.red == False:
            for n,i in enumerate(self.hblist):
                ref0 = self.xpm.array[self.ref_hb[n][0],:]
                self.red_hb[self.ref_hb[n][0]].perc = 100.0 * ref0.sum() / len(ref0)
                for r in range(1,len(self.ref_hb[n])):
                    ref1 = self.xpm.array[self.ref_hb[n][r],:]
                    self.red_hb[self.ref_hb[n][r]].perc = 100.0 * ref1.sum() / len(ref1)
                    ref0 = numpy.logical_or(ref0,ref1)
                self.hblist[n].perc  = 100.0 * ref0.sum() / len(ref0)
                self.hblist[n].nfr = self.nfr
        else:
            for n,i in enumerate(self.hblist):
                self.hblist[n].perc  = 100.0 * self.xpm.array[n,:].sum() / self.nfr
                self.hblist[n].nfr = self.nfr
    
    def autocorr_time(self,b=0,e=-1):
        self.ac_t = numpy.zeros([self.nbonds,3])
        if b != 0:
            if b >= len(self.hblist):
                print("b (%d) should be lower than %d" %(b, len(self.hblist)))
                raise ValueError
        if e > 0:
            if e < b:
                print("e (%d) should be higher than b (%d)" %(e, b))
                raise ValueError
            if e >= len(self.hblist):
                 print("e (%d) should be lower than %d" %(e, len(self.hblist)))
                 raise ValueError
        if e < 0:
            e = len(self.hblist) + e
        for n in numpy.arange(b,e,dtype=int):
            d0   = self.xpm.array[self.ref_hb[n][0],:]
            for r in range(1,len(self.ref_hb[n])):
                d1 = self.xpm.array[self.ref_hb[n][r],:]
                d0 = numpy.logical_or(d0,d1)
            self.ac_t[n], o = scipy.optimize.curve_fit(_exp2,numpy.arange(self.nfr),_autocorr(d0),p0=(0.01,1.0,0.9),bounds=( [0.,0.,0.],[numpy.inf, numpy.inf,1.0] ) )
    
    def calc_lifetime(self,b=0,e=-1):
        LT = []
        if b != 0:
            if b >= len(self.hblist):
                print("b (%d) should be lower than %d" %(b, len(self.hblist)))
                raise ValueError
        if e > 0:
            if e < b:
                print("e (%d) should be higher than b (%d)" %(e, b))
                raise ValueError
            if e >= len(self.hblist):
                 print("e (%d) should be lower than %d" %(e, len(self.hblist)))
                 raise ValueError
        if e < 0:
            e = len(self.hblist) + e
        for n in numpy.arange(b,e,dtype=int):
            d0   = self.xpm.array[self.ref_hb[n][0],:]
            for r in range(1,len(self.ref_hb[n])):
                d1 = self.xpm.array[self.ref_hb[n][r],:]
                d0 = numpy.logical_or(d0,d1)
            d    = numpy.concatenate([[0],d0,[0]])
            up   = numpy.array([ i+1 for i in numpy.arange(len(d)-1) if (d[i+1] - d[i] == 1) ])
            down = numpy.array([ i for i in numpy.arange(1,len(d)) if (d[i] - d[i-1] == -1) ])
            try: 
                LT.append(down-up)
            except:
                print("Problem calculating lifetime of HB %d [%s - %s ] " %(n, self.hblist[n].acc, self.hblist[n].don))
                raise
        return LT
    
    def calc_perc_red(self):
        for n,i in enumerate(self.hblist):
            ref  = self.xpm.array[n,:]
            self.hblist[n].perc = 100.0 * ref.sum() /len(ref0)
            self.hblist[n].nfr = self.nfr

    def perc_bootstrap(self,conf=0.95,nsample=1000):
        import time
        import sys
        print("Starting Bootstrap")
        t0 = time.time()
        if not self.red:
            for n,i in enumerate(self.hblist):
                ref0 = self.xpm.array[self.ref_hb[n][0],:]
                self.red_hb[self.ref_hb[n][0]].perc_ci = bootstrap(ref0,nsample,numpy.sum,1.0-conf) * 100.0 / len(ref0)
                self.red_hb[self.ref_hb[n][0]].perc_c  = conf
                self.red_hb[self.ref_hb[n][0]].perc    = numpy.mean(self.red_hb[self.ref_hb[n][0]].perc_ci)
                for r in range(1,len(self.ref_hb[n])):
                    ref1 = self.xpm.array[self.ref_hb[n][r],:]
                    self.red_hb[self.ref_hb[n][r]].perc_ci = bootstrap(ref1,nsample,numpy.sum,1.0-conf) * 100.0 / len(ref1)
                    self.red_hb[self.ref_hb[n][r]].perc_c  = conf
                    self.red_hb[self.ref_hb[n][r]].perc    = numpy.mean(self.red_hb[self.ref_hb[n][r]].perc_ci)
                    ref0 = numpy.logical_or(ref0,ref1)
                self.hblist[n].perc_ci = bootstrap(ref0,nsample,numpy.sum,1.0-conf) * 100.0 / len(ref0)
                self.hblist[n].perc_c  = conf
                self.hblist[n].perc = numpy.mean(self.hblist[n].perc_ci)
                self.hblist[n].nfr = self.nfr
                t = time.time()
                ETA =  ( t - t0 ) * ( float(self.nbonds) / ( n +1 ) - 1. )
                sys.stdout.write(" %4d/%4d calculations ETA: %6.1f sec.\r" %(n,self.nbonds,ETA)) 
        else:
            for n,i in enumerate(self.hblist):
                ref0 = self.xpm.array[n,:]
                self.hblist[n].perc_ci = bootstrap(ref0,nsample,numpy.sum,1.0-conf) * 100.0 / len(ref0)
                self.hblist[n].perc_c  = conf
                self.hblist[n].perc = numpy.mean(self.hblist[n].perc_ci)
                self.hblist[n].nfr = self.nfr
                t = time.time()
                ETA =  ( t - t0 ) * ( float(self.nbonds) / ( n +1 ) - 1. )
                sys.stdout.write(" %4d/%4d calculations ETA: %6.1f sec.\r" %(n,self.nbonds,ETA))  
        sys.stdout.write (" Completed in %6.1f sec.                \n" %(t-t0))

    def read_perc(self,filein,red=False):
        fi          = open(filein,'r')
        raw         = fi.readlines()
        self.name   = raw[0]
        self.hblist = []
        if not red:
            l_sim=[]
            hblist = []
            for l in raw[1:]:
                line = l.split()
                l_sim.append([ list(i)[0]+list(i)[1] for i in re_R_sim.findall(l) ])
                d_r , a_r = l_sim[-1][0] , l_sim[-1][1]
                p = float(line[2])
                nfr = int(line[5]) 
                self.hblist.append(HBond(don=d_r, acc=a_r, perc=p, mediated=0.0, both=0.0, nfr=nfr))
        else:
            for l in raw[1:]:
                line = l.split()
                hb = re_R.findall(l)
                d_r , a_r = hb[0][0]+hb[0][1] , hb[1][0]+hb[1][1]
                d_a , a_a = hb[0][2], hb[1][2]
                p = float(line[2])
                nfr = int(line[5])
                self.hblist.append(HBond(don=d_r, don_atom=d_a, acc=a_r, acc_atom=a_a, perc=p, mediated=0.0, both=0.0, nfr=nfr))
        self.nfr = nfr
        self.red = red
    
    def read_file_perc(self,filein,filegro):
        """
            Populates the HBonds object using a percentage file already
            generated.
            
            Only the uncommented lines including 6 fields will be read
            
            The format of the file should be the following:
            
            Donor     Acceptor      Direct   Mediated       Both      Total Frames
            RES999    RES000         99.99      33.33      66.66       9999
            
        """
        fi        = open(filein,'r')
        raw_fi    = fi.readline()
        self.name = raw_fi[:-1]
        self.mol  = filegro
        
        n = 0
        while raw_fi:
            line = raw_fi.split()
            if ( line[0] == '#' ):
                raw_fi=fi.readline()
                continue
            if ( len(line) == 6 ):
                n = n + 1
                t_don  = self.mol.split_res(line[0])
                t_acc  = self.mol.split_res(line[1])
                t_perc = float(line[2])
                t_nfr  = int(line[5])
                self.hblist.append(HBond(acc=t_acc,don=t_don,perc=t_perc,nfr=t_nfr))
                raw_fi=fi.readline()
                continue
            raw_fi=fi.readline()
        print("read %d lines from file %s %d frames" %(n,filein,t_nfr))
        self.nbonds = len(self.hblist)
        fi.close()
        self.refresh()
    
    def read_file_perc_red(self,filein,filegro):
        """
            Populates the HBonds object using a percentage file already
            generated.
            
            Only the uncommented lines including 6 fields will be read
            
            The format of the file should be the following:
            
            Donor       Acceptor          Direct   Mediated       Both      Total Frames
            RES999NAME  RES000NAME         99.99      33.33      66.66       9999
            
        """
        fi        = open(filein,'r')
        raw_fi    = fi.readline()
        self.name = raw_fi[:-1]
        self.mol  = filegro
        
        n = 0
        while raw_fi:
            line=raw_fi.split()
            if ( len(line) == 6):
                if ( line[0] != '#' ):
                    n = n + 1
                    t_don    = self.mol.find_atom_res(line[0])
                    t_acc    = self.mol.find_atom_res(line[1])
                    don_atom = self.mol.find_atom_name(line[0])
                    acc_atom = self.mol.find_atom_name(line[1])
                    t_perc   = float(line[2])
                    t_nfr    = int(line[5])
                    self.hblist.append(HBond(acc=t_acc,acc_atom=acc_atom,don=t_don,don_atom=don_atom,perc=t_perc,nfr=t_nfr))
            raw_fi=fi.readline()
        print("read %d lines from file %s %d frames" %(n,filein,t_nfr))
        self.nbonds = len(self.hblist)
        fi.close()
        self.refresh()
    
    def list_acc(self):
        return uniq( [ hb.acc for hb in self ] )
    
    def list_acc_n(self):
        return uniq( [ int(re_R_sim.findall(hb.acc)[0][1]) for hb in self ] )
    
    def list_don(self):
        return uniq( [ hb.don for hb in self ] )
    
    def list_don_n(self):
        return uniq( [ int(re_R_sim.findall(hb.don)[0][1]) for hb in self ] )
    
    def list_res(self):
        return uniq(self.list_acc()+self.list_don())
    
    def list_res_n(self):
        return uniq(self.list_acc_n()+self.list_don_n())
    
    def merge_acc(self, other):
        return uniq(self.list_acc()+other.list_acc())
    
    def merge_acc_n(self, other):
        return uniq(self.list_acc_n()+other.list_acc_n())
    
    def merge_don(self, other):
        return uniq(self.list_don()+other.list_don())
    
    def merge_don_n(self, other):
        return uniq(self.list_don_n()+other.list_don_n())
    
    def merge_res(self, other):
        return uniq(self.list_res()+other.list_res())
    
    def merge_res_n(self, other):
        return uniq(self.list_res_n()+other.list_res_n())

def merge_two(first, second):
    """
        Merge two HBonds object by performing the weightd averages
        based on the number of number of frames (self.nfr)
    """
    merged=HBonds()
    check2=numpy.zeros(len(second.hblist),dtype=int)
    for hb1 in first.hblist:
        check_is_same_couple=0
        for n2,hb2 in enumerate(second.hblist):
            if hb1.is_same_couple(hb2):
                check_is_same_couple=1
                check2[n2]=1
                merged.hblist.append(HBond(acc=hb1.acc,don=hb1.don,perc=(hb1.perc*hb1.nfr+hb2.perc*hb2.nfr)/float(hb1.nfr+hb2.nfr),nfr=hb1.nfr+hb2.nfr))
                break
        if check_is_same_couple==0:
            merged.hblist.append(HBond(acc=hb1.acc,don=hb1.don,perc=float(hb1.perc*hb1.nfr)/float(hb1.nfr+hb2.nfr),nfr=hb1.nfr+hb2.nfr))
    for n2, hb2 in enumerate(second.hblist):
        if check2[n2]==0:
            check2[n2]=1
            merged.hblist.append(HBond(acc=hb2.acc,don=hb2.don,perc=float(hb2.perc*hb2.nfr)/float(hb1.nfr+hb2.nfr),nfr=hb1.nfr+hb2.nfr))
    merged.refresh()
    return merged

class HBComp:
    def __init__(self,don=['res1',0],acc=['res2',0],perc1=0.0,perc1_2=0.0,perc2=0.0,perc2_2=0.0,nfr1=0,nfr2=0,fisher=None,conf=None):
        self.acc=acc
        self.don=don
        self.perc1=perc1
        self.perc1_2=perc1_2
        self.perc2=perc2
        self.perc2_2=perc2_2
        self.nfr1=nfr1
        self.nfr2=nfr2
        self.fisher=fisher
        self.conf=conf
    
    def __str__(self):
        if self.fisher != None:
            return "%3s | %4d | %3s | %4d |    %6.2f    |    %6.2f    | %12.10f\n" %(self.don[0],self.don[1],self.acc[0],self.acc[1],self.perc1,self.perc2,self.fisher)
        if self.conf != None:
            return "%3s | %4d | %3s | %4d | %5.1f-%5.1f | %5.1f-%5.1f | %6.4f\n" %(self.don[0],self.don[1],self.acc[0],self.acc[1],self.perc1,self.perc1_2,self.perc2,self.perc2_2,1.0-self.conf)

    def __repr__(self):
        print(str(self))
    
    def renum(self,func):
        self.don[1]=func(self.don[1])
        self.acc[1]=func(self.acc[1])

class HBondsCompare:
    """
        This Class can be used to do comparison between two HBonds object
        
    """
    def __init__(self, first=HBonds(name='First' ,mol=None,red=False,log=None,xpm=None,conf=None),\
                          second=HBonds(name='Second',mol=None,red=False,log=None,xpm=None,conf=None)):
        self.first=first
        self.second=second
        self.don=self.first.merge_don(self.second)
        self.acc=self.first.merge_acc(self.second)
        #self.res=self.first.merge_res(self.second)
        self.acc_n=self.first.merge_acc_n(self.second)
        self.don_n=self.first.merge_don_n(self.second)
        #self.res_n=self.first.merge_res_n(self.second)
        #
        # Matrix M1 and M2 will contain the data of percentage of existance
        # and the nfr of the thwo HBonds objects
        #
        self.M1=numpy.zeros([len(self.don_n),len(self.acc_n),3])
        self.M2=numpy.zeros([len(self.don_n),len(self.acc_n),3])
        #
        # Matrix E is a binary matrix for internal use:
        # E[i,j] == 1  if a HBond between i an j exists.
        #
        self.E =numpy.zeros([len(self.don_n),len(self.acc_n)],dtype=int)
        #
        # Matrix F contains the Fisher scores of the difference of percentage
        # F[i,j] = Fisher_Score( M1[i,j], M2[i,j] )
        #
        self.F =numpy.zeros([len(self.don_n),len(self.acc_n)])
        self.compare=[]
    
    def __str__(self):
        out="%3s | %4s | %3s | %4s | %10s  | %10s  | %10s\n" %('Don','#','Acc','#','Perc.1','Perc.2','P-Val')
        for hb in self.compare:
            out=out+str(hb)
        return out
    
    def __repr__(self):
        print("%3s | %4s | %3s | %4s | %10s  | %10s  | %10s" %('Don','#','Acc','#','Perc.1','Perc.2','P-val'))
        for hb in self.compare:
            print(hb)
    
    def renum(self,func):
        """
            Permits to renumber the residues using a user-defined function
        """
        for hb in self.compare:
            hb.renum(func)
    
    def do_comparison(self,perc=0.0,verbose=False):
        """
            After creating the object, this function permits to actually
            perform a comparison.
            
            In the current version only the Exact-Fisher score for the 
            evaluation of significance of the differences is implemented. 
        """
        for i in range(self.first.nbonds):
            self.M1[self.don_n.index(self.first.hblist[i].don[1]),self.acc_n.index(self.first.hblist[i].acc[1]),0]  = self.first.hblist[i].perc
            self.M1[self.don_n.index(self.first.hblist[i].don[1]),self.acc_n.index(self.first.hblist[i].acc[1]),1]  = self.first.hblist[i].nfr
            self.E [self.don_n.index(self.first.hblist[i].don[1]),self.acc_n.index(self.first.hblist[i].acc[1])] += 1
        for i in range(self.second.nbonds):
            self.M2[self.don_n.index(self.second.hblist[i].don[1]),self.acc_n.index(self.second.hblist[i].acc[1]),0]  = self.second.hblist[i].perc
            self.M2[self.don_n.index(self.second.hblist[i].don[1]),self.acc_n.index(self.second.hblist[i].acc[1]),1]  = self.second.hblist[i].nfr
            self.E [self.don_n.index(self.second.hblist[i].don[1]),self.acc_n.index(self.second.hblist[i].acc[1])] += 1
        for d in numpy.arange(len(self.don_n)):
            for a in numpy.arange(len(self.acc_n)):
                if (self.E[d,a]!=0)&((self.M1[d,a,0]>perc)|(self.M2[d,a,0]>perc)):
                    if verbose:
                        print(" [%d,%d] F-score: fscore(%5.2f,%5.2f,%d,%d)" %(d,a,self.M1[d,a,0],self.M2[d,a,0],self.first.hblist[0].nfr,self.second.hblist[0].nfr))
                    self.F[d,a] = fscore(self.M1[d,a,0],self.M2[d,a,0],self.first.hblist[0].nfr,self.second.hblist[0].nfr)
                    #print a,d,self.F[a,d]
                    self.compare.append(hbcomp(acc=self.acc[a],don=self.don[d],perc1=self.M1[d,a,0],nfr1=self.M1[d,a,1],perc2=self.M2[d,a,0],nfr2=self.M2[d,a,1],fisher=self.F[d,a]))


    def compare_bootstrap(self,confidence=0.95,perc=0.0,verbose=False):
        """
            After creating the object, this function permits to actually
            perform a comparison.
            
            In the current version only the Exact-Fisher score for the 
            evaluation of significance of the differences is implemented. 
        """
        test_conf_first=0
        test_conf_second=0
        for i in self.first.hblist:
            if i.perc_c != confidence :
                test_conf_first += 1
        for i in self.second.hblist:
            if i.perc_c != confidence :
                test_conf_second += 1
        
        if (test_conf_first+test_conf_second) != 0 :
            self.first.perc_bootstrap(confidence)
            self.second.perc_bootstrap(confidence)
                    
        for i in range(self.first.nbonds):
            self.M1[self.don_n.index(self.first.hblist[i].don[1]),self.acc_n.index(self.first.hblist[i].acc[1]),0]  = self.first.hblist[i].perc_ci[0]
            self.M1[self.don_n.index(self.first.hblist[i].don[1]),self.acc_n.index(self.first.hblist[i].acc[1]),1]  = self.first.hblist[i].perc_ci[1]
            self.M1[self.don_n.index(self.first.hblist[i].don[1]),self.acc_n.index(self.first.hblist[i].acc[1]),2]  = self.first.hblist[i].perc
            self.E [self.don_n.index(self.first.hblist[i].don[1]),self.acc_n.index(self.first.hblist[i].acc[1])]   += 1
        for i in range(self.second.nbonds):
            self.M2[self.don_n.index(self.second.hblist[i].don[1]),self.acc_n.index(self.second.hblist[i].acc[1]),0]  = self.second.hblist[i].perc_ci[0]
            self.M2[self.don_n.index(self.second.hblist[i].don[1]),self.acc_n.index(self.second.hblist[i].acc[1]),1]  = self.second.hblist[i].perc_ci[1]
            self.M2[self.don_n.index(self.second.hblist[i].don[1]),self.acc_n.index(self.second.hblist[i].acc[1]),2]  = self.second.hblist[i].perc
            self.E [self.don_n.index(self.second.hblist[i].don[1]),self.acc_n.index(self.second.hblist[i].acc[1])]   += 1
        for d in numpy.arange(len(self.don_n)):
            for a in numpy.arange(len(self.acc_n)):
                if (self.E[d,a]!=0) & ( (self.M1[d,a,2] > perc          ) | (self.M2[d,a,2] > perc          ) ) & ( (self.M1[d,a,0] > self.M2[d,a,1]) | (self.M1[d,a,1] < self.M2[d,a,0]) ):
                    self.compare.append(hbcomp(acc=self.acc[a],don=self.don[d],perc1=self.M1[d,a,0],perc1_2=self.M1[d,a,1],nfr1=self.first.nfr,perc2=self.M2[d,a,0],perc2_2=self.M2[d,a,1],nfr2=self.second.nfr,conf=confidence))
