import chimera
import numpy as np
from math import *

def geocentre(mol):
	atoms = mol.atoms
	c=[]
	for a in atoms:
		c.append(list(a.coord().data()))
	coord=np.array(c)
	return coord.mean(axis=0)
	
def max_radius(mol,gc):
	atoms = mol.atoms
	c=[]
	rmax=0.0
	for a in atoms:
		r=0.
		v=np.array(a.coord().data())-gc
		for i in v:
			r=r+i**2.
		r=sqrt(r)
		if r>rmax:
			rmax=r
	return rmax

def mol_translate(mol, p):
	for a in mol.atoms:
		x = a.coord().x + p[0]
		y = a.coord().y + p[1]
		z = a.coord().z + p[2]
		new_coord=chimera.Point(x,y,z)
		a.setCoord(new_coord)
	
		
