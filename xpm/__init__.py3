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

    def write_xpm(self,file_out):
        #fo=open(file_out,'w')
                return  

def read_xpm(file_in,verbose=False):
    fi=open(file_in,'r')
    if fi.readline()!='/* XPM */\n':
        raise IOError
        fi.close()
        exit(2)
    xpm=Xpm()
    header=5
    while (header!=0):
        newline=fi.readline()
        if newline.split()[1]=='title:':
            if verbose:
                print("Found a Title")
            xpm.title=newline.split('"')[1]
            header -= 1 
        elif newline.split()[1]=='legend:':
            if verbose:
                print("Found a Legend")
            xpm.legend=newline.split('"')[1]
            header -= 1
        elif newline.split()[1]=='x-label:':
            if verbose:
                print("Found a x-label")
            xpm.xlabel=newline.split('"')[1]
            header -= 1
        elif newline.split()[1]=='y-label:':
            if verbose:
                print("Found a y-label")
            xpm.ylabel=newline.split('"')[1]
            header -= 1
        elif newline.split()[1]=='type:':
            if verbose:
                print("Found a type")
            xpm.typ=newline.split('"')[1]
            header -= 1
        else:
            if verbose:
                print("Trash")
    #XPM reading starts
    start_xpm=1
    while(start_xpm!=0):
        newline=fi.readline()
        if (newline.split()[0]=='static')&(newline.split()[1]=='char'):
            xpm.name=newline.split()[2][1:-2]
            start_xpm -= 1
        else:
            raise IOError
            fi.close()
            exit(2)
    newline=fi.readline()
    xpm.cols, xpm.rows, xpm.colors, xpm.cop = newline[1:-3].split()
    left=xpm.colors
    xpm.scal.c=xpm.colors
    while(left!=0):
        newline=fi.readline()
        xpm.scal.list_ch.append(newline[1:xpm.cop+1])
        xpm.scal.list_c.append(newline[5+xpm.cop:12+xpm.cop])
        xpm.scal.list_v.append(newline.split('"')[3])
        left -= 1
    newline=fi.readline()
    xaxis=[]
    yaxis=[]
    while(newline.split()[0]=='/*'):
        line_split=newline.split()
        if (line_split[1] == 'x-axis:'):
            for i in range(2,len(line_split)-2):
                xaxis.append(float(line_split[i]))
        elif (line_split[1] == 'y-axis:'):
            for i in range(2,len(line_split)-2):
                yaxis.append(float(line_split[i]))
        else:
            if verbose:
                print("Very Strange")
        newline=fi.readline()
    xpm.xaxis=numpy.array(xaxis)
    xpm.yaxis=numpy.array(yaxis)
    xpm.array=numpy.zeros([xpm.rows,xpm.cols])
    row_left=xpm.rows
    while(row_left!=1):
        row=newline
        r=xpm.rows-row_left
        for c in range(xpm.cols):
            xpm.array[r,c]=xpm.scal.list_ch.index(row[1+c*xpm.cop:1+(c+1)*xpm.cop])
        row_left -= 1
        newline=fi.readline()
    row=newline
    r=xpm.rows-1
    for c in range(xpm.cols):
        xpm.array[r,c]=xpm.scal.list_ch.index(row[1+c*xpm.cop:1+(c+1)*xpm.cop])
    xpm.array=numpy.flipud(xpm.array)
    return xpm

#
# 
# Some functions are COPIED FROM 
#
# GromacsWrapper: xpm.py
# Copyright (c) 2012 Oliver Beckstein <orbeckst@gmail.com>
# Copyright (c) 2010 Tsjerk Wassenaar <tsjerkw@gmail.com>
#: compiled regular expression to parse the colors in the xpm file::
#:
#:   static char *gromacs_xpm[] = {
#:   "14327 9   2 1",
#:   "   c #FFFFFF " /* "None" */,
#:   "o  c #FF0000 " /* "Present" */,
#:
#: Matches are named "symbol", "color" (hex string), and "value". "value"
#: is typically autoconverted to appropriate values with
#: :class:`gromacs.fileformats.convert.Autoconverter`.
#: The symbol is matched as a `printable ASCII character`_ in the range
#: 0x20 (space) to 0x7E (~).
#:
#: .. _`printable ASCII character`: http://www.danshort.com/ASCIImap/indexhex.htm
#
COLOUR = re.compile("""\
        ^.*"                   # start with quotation mark
        (?P<symbol>[\x20-\x7E])# printable ASCII symbol used in the actual pixmap: 'space' to '~'
        \s+                    # white-space separated
        c\s+                   # 'c' to prefix colour??
        (?P<color>\#[0-9A-F]+) # colour as hex string (always??)
        \s*"                   # close with quotes
        \s*/\*\s*              # white space then opening C-comment /*
        "                      # start new string
        (?P<value>.*)          # description/value as free form string
        "                      # ... terminated by quotes
        """, re.VERBOSE)
META = re.compile("""^.*/\*\s+(?P<name>.*):\s*"(?P<value>.*)" """,re.VERBOSE)

def col(c):
    """Parse colour specification"""
    m = COLOUR.search(c)
    if m == None:
        return False
    else:
        value = m.group('value')
        color = m.group('symbol')
        return color, value

def meta(c):
    """Parse meta specification"""
    m = META.search(c)
    if m == None:
        return False
    else:
        value = m.group('value')
        name = m.group('name')
        return name, value

def unquote(s):
    """Return string *s* with quotes ``"`` removed."""
    return s[1+s.find('"'):s.rfind('"')]

def uncomment(s):
    """Return string *s* with C-style comments ``/*`` ... ``*/`` removed."""
    return s[2+s.find('/*'):s.rfind('*/')]

def parse(file_in):
    """Parse the xpm file and populate :attr:`XPM.array`."""   
    with open(file_in) as fi:
        m_dat = [fi.readline()]
        while not m_dat[-1].startswith("static char *gromacs_xpm[]"):
            m_dat.append(fi.readline())
        m_dat = dict([meta(i) for i in m_dat if meta(i) != False ])
        dim = fi.readline()
        nx, ny, nc, nb = [int(i) for i in unquote(dim).split()]
        xpm=Xpm(title=m_dat["title"],name=file_in,xlabel=m_dat["x-label"],ylabel=m_dat["y-label"],typ=m_dat["type"],cols=nx,rows=ny,colors=nc,cop=nb)
        # The next dim[2] lines contain the color definitions
        # Each pixel is encoded by dim[3] bytes, and a comment
        # at the end of the line contains the corresponding value
        colors = dict([col(fi.readline()) for i in range(nc)])
        convert={}
        v=0
        for i in colors:
            convert[i] = v
            v +=1
        # make an array containing all possible values and let numpy figure out the dtype
        dtype = numpy.array(convert.values()).dtype
        # pre-allocate array
        data = numpy.zeros((nx/nb, ny), dtype=dtype)
        iy = 0
        xval = []
        yval = []
        for line in fi:
            if line.startswith("/*"):
                # lines '/* x-axis:' ... and '/* y-axis:' contain the
                # values of x and y coordinates
                s = uncomment(line).strip()
                if s.startswith('x-axis:'):
                    xval.extend([x for x in s[7:].split()])
                elif s.startswith('y-axis:'):
                    yval.extend([y for y in s[7:].split()])
                continue
            s = unquote(line)
            # Joao M. Damas <jmdamas@itqb.unl.pt> suggests on gmx-users (24 Oct 2014)
            # that the next line should read:
            #
            #  data[:, iy]  =  [colors[j[k:k+nb]] for k in range(0,nx*nb,nb)]
            #
            # "if one is using higher -nlevels for the .xpm construction (in g_rms, for example)"
            # However, without a test case I am not eager to change it right away so in
            # case some faulty behavior is discovered with the XPM reader then this comment
            # might be helpful. --- Oliver 2014-10-25
            data[:, iy] = [convert[s[k:k+nb]] for k in range(0,nx,nb)]
            iy += 1  # for next row
    xpm.xaxis = numpy.array(xval)
    xpm.array = data[:, ::-1]
    xpm.yaxis = numpy.array(yval)
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
    
