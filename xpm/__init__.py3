import sys
import numpy
import re

class CScale:
    def __init__(self, c=0, list_ch=[], list_v=[], list_c=[]):
        self.list_ch = list_ch
        self.list_v  = list_v
        self.list_c  = list_c
        self.c       = c

class Xpm:
  
    def __init__(self,title='',legend='',xlabel='',ylabel='',typ='Continuous',name='', cols=0, rows=0, colors=0, cop=0, array=numpy.array(0),scal=CScale(),xaxis=numpy.array(0), yaxis=numpy.array(0)):
        self.title  = title
        self.name   = name
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.type   = typ
        self.cols   = cols
        self.rows   = rows
        self.colors = colors
        self.cop    = cop
        self.array  = array
        self.scal   = scal
        self.xaxis  = xaxis
        self.yaxis  = yaxis

    def __setattr__(self,name,value):
        if (name=='cols')|(name=='rows')|(name=='colors')|(name=='cop'):
            self.__dict__[name]=int(value)
        elif (name=='array'):
            self.__dict__[name]=numpy.array(value,dtype=int)
        elif (name=='xaxis')|(name=='yaxis'):
            self.__dict__[name]=numpy.array(value)
        else:
            self.__dict__[name]=value

    #def write_xpm(self,file_out):
        #fo=open(file_out,'w')
        #return  

def read_xpm(file_in):
    import _pyxpm
    d = _pyxpm.read(file_in)
    cTab = CScale(c=len(d[4]),list_ch=[ i[0] for i in d[4] ], list_c = [ i[1] for i in d[4] ],list_v = numpy.arange(len(d[4])))
    xpm = Xpm(title=file_in,rows=d[0],cols=d[1],cop=d[2],colors=d[3],scal=cTab,array=d[5][:, ::-1],xaxis=numpy.arange(d[1]),yaxis=numpy.arange(d[0]))
    return xpm

def main(argv=None):
    import matplotlib.pyplot
    if argv==None:
        argv=sys.argv
    file_in=argv[1]
    data_xpm=read_xpm(file_in);
    matplotlib.pyplot.matshow(data_xpm.array)
    matplotlib.pyplot.show()
    return 0

if __name__ == "__main__":
    sys.exit(main())
    
