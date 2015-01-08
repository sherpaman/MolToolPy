#!/usr/bin/env python

import numpy as np
from scipy.special import erfc

class ResContrib():
	def __init__(self,files=None,n_rep=0,n_res=0,protein="Protein Name",ligand="Ligand Name",site="Site",resid_file=None):
		self.files=files
		self.n_rep=n_rep
		self.protein=protein
		self.ligand=ligand
		self.site=site
		self.initialize()
		if resid_file:
			self.name_residues(resid_file)


	def __setattr__(self,name,value):
		if (name=='n_replicas')|(name=='n_res'):
			self.__dict__[name]=int(value)
		elif (name=="protein")|(name=="ligand")|(name=="site"):
			self.__dict__[name]=str(value)
		else:
			self.__dict__[name]=value
			
	def initialize(self):
		ene=[]
		std=[]
		self.n_rep=len(self.files)
		for n,f_in in enumerate(self.files):
			data = self.read_file(f_in)
			if n==0:
				self.n_res=len(data)
			ene.append(data[0])
			std.append(data[1])
		self.energy=np.array(ene)
		self.stdev=np.array(std)
		self.aver_energy()

	def read_file(self,filename):
		read_data=np.loadtxt(filename)
		return read_data[:,0], read_data[:,1]
		
	def aver_energy(self):
		self.average = np.average(self.energy,axis=0,weights=(1.0/self.stdev)**2)
		self.error   = np.sqrt(1.0/np.sum(self.stdev**-2,axis=0))
		self.signif=erfc(abs(self.average)/self.error)
	
	def name_residues(self,filename):
		f=open(filename,"r")
		self.residue=[]
		for r in f.readlines():
			self.residue.append(r[:-1])
		f.close()

	def print_all(self):
		self.n_res=len(self.average)
		string_out=""
		for i in np.arange(self.n_res):
			if self.residue!=None:
				string_out = string_out + "%8s %6s %4s %8.4f +/- %8.4f\n" %(self.residue[i],self.protein,self.site,self.average[i], self.error [i])
			else:
				string_out = string_out + "%8s %6s %4s %8.4f +/- %8.4f\n" %(i+1,self.protein,self.site,self.average[i], self.error [i])
		self.print_out=string_out
		print string_out

	def print_selection(self,selection):
		string_out   = ""
		print_average= self.average[selection]
		print_error  = self.error  [selection]
		if self.residue!=None:
			print_res=[]
			for i,r in enumerate(self.residue):
				if selection[i]:
					print_res.append(r)
		else:
			print_res=(np.arange(n_res)+1)[selection]
		for i in np.arange(len(print_average)):
			string_out = string_out + "%8s %6s %4s %8.4f +/- %8.4f\n" %(print_res[i],self.protein,self.site,print_average[i], print_error [i])
		return string_out

	def print_comparison(self,other,selection):
		string_out   = ""
		print_average_1= self.average[selection]
		print_error_1  = self.error  [selection]
		print_average_2= other.average[selection]
		print_error_2  = other.error  [selection]

		#IT ASSUME THE RESIDUE LIST IS THE SAME IN THE TWO SYSTEMS
		if self.residue!=None:
			print_res=[]
			for i,r in enumerate(self.residue):
				if selection[i]:
					print_res.append(r)
		else:
			print_res=(np.arange(n_res)+1)[selection]
		for i in np.arange(len(print_average_1)):
			diff_signif=erfc(np.abs(print_average_1[i]-print_average_2[i])/(print_error_1[i]+print_error_2[i]))
			string_out = string_out + "%8s %6s %4s %8.4f +/- %8.4f | %6s %4s %8.4f +/- %8.4f : %8.4f \n" %(print_res[i], self.protein, self.site, print_average_1[i], print_error_1[i], other.protein, other.site, print_average_2[i], print_error_2[i],diff_signif)
		return string_out

	def print_comparison_three(self,other,another,selection):
		string_out   = ""
		print_average_1= self.average[selection]
		print_error_1  = self.error  [selection]
		print_average_2= other.average[selection]
		print_error_2  = other.error  [selection]
		print_average_3= another.average[selection]
		print_error_3  = another.error  [selection]

		#IT ASSUME THE RESIDUE LIST IS THE SAME IN THE THREE SYSTEMS
		if self.residue!=None:
			print_res=[]
			for i,r in enumerate(self.residue):
				if selection[i]:
					print_res.append(r)
		else:
			print_res=(np.arange(n_res)+1)[selection]
		for i in np.arange(len(print_average_1)):
			diff_signif_2=erfc(np.abs(print_average_1[i]-print_average_2[i])/(print_error_1[i]+print_error_2[i]))
			diff_signif_3=erfc(np.abs(print_average_1[i]-print_average_3[i])/(print_error_1[i]+print_error_3[i]))
			string_out = string_out + "%8s %6s %4s %8.4f +/- %8.4f | %6s %4s %8.4f +/- %8.4f : %8.4f | %6s %4s %8.4f +/- %8.4f : %8.4f\n" %(print_res[i], self.protein, self.site, print_average_1[i], print_error_1[i], other.protein, other.site, print_average_2[i], print_error_2[i],diff_signif_2, another.protein, another.site, print_average_3[i], print_error_3[i],diff_signif_3)
		return string_out

	
	def print_significant(self,pval,e_thresh):
		significant_sel=np.logical_and(self.signif<pval,np.abs(self.average)>e_thresh)
		self.print_out=self.print_selection(significant_sel)
		print self.print_out
		
	def compare(self,other,pval,e_thresh):
		interesting_1 = np.logical_and( self.signif<pval,np.abs( self.average)>e_thresh)
		interesting_2 = np.logical_and(other.signif<pval,np.abs(other.average)>e_thresh)
		list_compare=np.logical_or(interesting_1,interesting_2)
		print self.print_comparison(other,list_compare)
	
	def compare_three(self,other,another,pval,e_thresh):
		interesting_1 = np.logical_and( self.signif<pval,np.abs( self.average)>e_thresh)
		interesting_2 = np.logical_and(other.signif<pval,np.abs(other.average)>e_thresh)
		interesting_3 = np.logical_and(another.signif<pval,np.abs(another.average)>e_thresh)
		list_compare=np.logical_or(np.logical_or(interesting_1,interesting_2),interesting_3)
		print self.print_comparison_three(other,another,list_compare)

if __name__ == "__main__":
	
	proteins=["WT","R218H","R218P"]
	sites=["TR1","TR2","TR4"]
	
	num_replicas=3
	
	DATA=[]
	
	for p in proteins:
		for s in sites:
			files=[]
			for n in range(1,num_replicas+1):
				f="%s_%s_%d_contrib.dat" %(p,s,n)
				files.append(f)
			data=ResContrib(files=files,protein="HSA_%s" % p,ligand="T44",site=s,resid_file="resid_list.dat")
			DATA.append(data)
			#data.print_significant(0.001,1.0)
	DATA[0].compare_three(DATA[3],DATA[6],0.001,1.0)
	DATA[1].compare_three(DATA[4],DATA[7],0.001,1.0)
	DATA[2].compare_three(DATA[5],DATA[8],0.001,1.0)

