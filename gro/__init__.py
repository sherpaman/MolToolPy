import re

def read_gro(gro):
	f=open(gro)
	data=f.readlines()
	f.close()
	natoms=int(data[1])
	res_list=[]
	for i in range(2,natoms+2):
		line=data[i]
		num_res   = int(line[0:5])
		name_res  = line[5:9].rstrip().lstrip()
		name_atom = line[10:15].rstrip().lstrip()
		num_atom  = int(line[15:20])
		A=atom(name=name_atom,num=num_atom,resname=name_res,resnum=num_res)
		if i==2:
			active_res=res(name=name_res,num=num_res,atoms=[])
			active_res.atoms.append(A)
		else:
			if (active_res.name != name_res) | (active_res.num != num_res):
				#print "NEW RESIDUE"
				res_list.append(active_res)
				active_res=res(name=name_res,num=num_res,atoms=[])
			#print "ADDING ATOM"
			active_res.atoms.append(A)
		#print "RES %5d : %4s | ATOM %5d : %5s" %(num_res,name_res,num_atom,name_atom)
	res_list.append(active_res)
	out_mol=molecule(name=data[0],residues=res_list,nres=len(res_list))
	for r in out_mol.residues:
		r.mol=out_mol
		r.natoms=len(r.atoms)
		for a in r.atoms:
			a.res=r
	return out_mol

class molecule:
	def __init__(self,name='molecule',residues=[],nres=0):
		self.name=name
		self.residues=residues
		self.nres=nres
	
	def __str__(self):
		return "Molecule: %s (%d residues)" %(self.name,self.nres)
	
	def atom(self,num):
		for r in self.residues:
			for a in r.atoms:
				if a.num==num:
					return a
	
	def find_atom_in_res(self,res,atom):
		for n,r in enumerate(self.residues):
			if (r == str(res) ) | ( n == int(res) ):
				for a in r.atoms:
					if atom==str(a):
						return a.num
	
	def find_atom(self,atom):
		for r in self.residues:
			for a in r.atoms:
				if atom==str(a):
					return a.num
	
	def find_atom_name(self,atom):
		for r in self.residues:
			for a in r.atoms:
				if atom==str(a):
					return a.name
	
	def find_atom_res(self,atom):
		for r in self.residues:
			for a in r.atoms:
				if atom==str(a):
					return [r.name, r.num]
	
	def split_name(self,atom):
		for r in self.residues:
			for a in r.atoms:
				if atom==str(a):
					return [r.name, r.num, a.name]
	
	def split_res(self,resid):
		for r in self.residues:
			if (r.name+str(r.num))==resid:
				return [r.name, r.num]
	
class res:
	def __init__(self,name='residue',num=0,atoms=[],natoms=0,mol=molecule()):
		self.name=name
		self.num=num
		self.atoms=atoms
		self.natoms=natoms
		self.mol=mol
	
	def __str__(self):
		return "Residue: %s (%d atoms)" %(self.name,self.natoms)
	
	def __repr__(self):
		return self.__str__()
	
	def atom_num_name(self,num):
		for a in self.atoms:
			if a.num==num:
				return a
				
	def atom_name_num(self,name):
		for a in self.atoms:
			if a.name==name:
				return a.num

class atom:
	def __init__(self,name='',num=0,resname=0,resnum=0,resid=res()):
		self.name=name
		self.num=num
		self.resname=resname
		self.resnum=resnum
		self.resid=res
	
	def __str__(self):
		return "%s%d%s" %(self.resname,self.resnum,self.name)
	
	def __repr__(self):
		self.__str__()
