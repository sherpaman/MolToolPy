import re

def read_gro(gro):
	f        = open(gro)
	data     = f.readlines()
        natoms   = int(data[1])
	res_list = []
        f.close()
	for i in range(2,natoms+2):
		line      = data[i]
		num_res   = int(line[0:5])
		name_res  = line[5:9].rstrip().lstrip()
		name_atom = line[10:15].rstrip().lstrip()
		num_atom  = int(line[15:20])
		if i==2:
			active_res =  Res(name=name_res,  num=num_res, atoms=[])
                        A          = Atom(name=name_atom, num=num_atom, resname=name_res, resnum=num_res,resid=active_res)
			active_res.atoms.append(A)
		else:
			if (active_res.name != name_res) | (active_res.num != num_res):
				res_list.append(active_res)
				active_res = Res(name=name_res,  num=num_res, atoms=[])
                        A = Atom(name=name_atom, num=num_atom, resname=name_res, resnum=num_res,resid=active_res)
			active_res.atoms.append(A)
	res_list.append(active_res)
	out_mol = Molecule(name=data[0][:-1], residues=res_list, nres=len(res_list))
	for r in out_mol.residues:
		r.mol    = out_mol
		r.natoms = len(r.atoms)
		for a in r.atoms:
			a.res = r
	return out_mol

class Molecule:
	def __init__(self, name='molecule', residues=[], nres=0):
		self.name     = name
		self.residues = residues
		self.nres     = nres
	
        def __iter__(self):
                return iter(self.residues)
        
        def __getitem__(self,key):
                return self.residues[key]
        
        def __len__(self):
                return len(self.residues)
        
	def __str__(self):
		return "Molecule: %s (%d residues)" %(self.name, self.nres)
        
        def __repr__(self):
                return "Molecule: %s (%d residues)" %(self.name, self.nres)
	
        def atom_iterator(self):
                return iter([ atom for res in self for atom in res ])
        
	def atom(self,num):
		return [a for a in self.atom_iterator() if a.num==num ][0]
	
	def find_atom_in_res(self,res,atom):
                return [ a.num for a in self.atom_iterator() if ( (a.name == atom) & ( (r.num == int(res)) | (r.name == str(res) ) ) ) ][0]
	
	def find_atom(self,atom):
                return [ a.num for a in self.atom_iterator() if ( str(a) == atom ) ][0]
	
	def find_atom_name(self,atom):
                return [ a.name for a in self.atom_iterator() if ( str(a) == atom ) ][0]
        
	def find_atom_res(self,atom):
                return [ [r.name, r.num] for r in self for a in r if ( str(a) == atom ) ][0]
	
	def split_name(self,atom):
		return [ [r.name, r.num, a.name] for r in self for a in r if ( str(a) == atom ) ][0]
	
	def split_res(self, res):
                return [ [r.name, r.num] for r in self if ( str(r) == res ) ][0]
		
	
class Res:
	def __init__(self, name='residue', num=0, atoms=[], natoms=0, mol=Molecule()):
		self.name   = name
		self.num    = num
		self.atoms  = atoms
		self.natoms = natoms
		self.mol    = mol
        
        def __iter__(self):
                return iter(self.atoms)
        
        def __getitem__(self,key):
                return self.atoms[key]
        
        def __len__(self):
                return len(self.atoms)
	
	def __str__(self):
		return "%s%d" %(self.name, self.num)
        
        def __repr__(self):
                return "%s%d" %(self.name, self.num)
	
	def atom_num_name(self, num):
                return [ a for a in self if a.num == num ][0]
				
	def atom_name_num(self,name):
                return [ a.num for a in self if a.name == name ][0]

class Atom: 
	def __init__(self,name='', num=0, resname=0, resnum=0, resid=Res()):
		self.name    = name
		self.num     = num
		self.resname = resname
		self.resnum  = resnum
		self.resid   = resid
	
	def __str__(self):
		return "%s%d%s" %(self.resname, self.resnum, self.name)
        
        def __repr__(self):
                return "%s%d%s" %(self.resname, self.resnum, self.name)
        
