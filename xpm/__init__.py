import sys
import numpy

class c_scale:
	def __init__(self, c=0, list_ch=[], list_v=[], list_c=[]):
		self.list_ch = list_ch
		self.list_v  = list_v
		self.list_c  = list_c
		self.c       = c
	
class xpm_data:
	def __init__(self,title='',legend='',xlabel='',ylabel='',typ='Continuous',name='',cols=0,rows=0,colors=0,cop=0,array=numpy.array(0),scal=c_scale(),xaxis=numpy.array(0),yaxis=numpy.array(0)):
		self.title=title
		self.name=name
		self.xlabel=xlabel
		self.ylabel=ylabel
		self.type=typ
		self.cols=cols
		self.rows=rows
		self.colors=colors
		self.cop=cop
		self.array=array
		self.scal=scal
		self.xaxis=xaxis
		self.yaxis=yaxis

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
		fo=open(file_out,'w')

def read_xpm(file_in,verbose=False):
	fi=open(file_in,'r')
	if fi.readline()!='/* XPM */\n':
		raise IOError
		fi.close()
		exit(2)
	xpm=xpm_data()
	header=5
	while (header!=0):
		newline=fi.readline()
		if newline.split()[1]=='title:':
			if verbose:
				print "Found a Title"
			xpm.title=newline.split('"')[1]
			header -= 1	
		elif newline.split()[1]=='legend:':
			if verbose:
				print "Found a Legend"
			xpm.legend=newline.split('"')[1]
			header -= 1
		elif newline.split()[1]=='x-label:':
			if verbose:
				print "Found a x-label"
			xpm.xlabel=newline.split('"')[1]
			header -= 1
		elif newline.split()[1]=='y-label:':
			if verbose:
				print "Found a y-label"
			xpm.ylabel=newline.split('"')[1]
			header -= 1
		elif newline.split()[1]=='type:':
			if verbose:
				print "Found a type"
			xpm.typ=newline.split('"')[1]
			header -= 1
		else:
			if verbose:
				print "Trash"
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
				print "Very Strange"
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
	
