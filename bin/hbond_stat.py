#!/usr/bin/env python
from sys import argv,stderr
#Prende in input il nome di un file che contiene, i dati di coppie di residui per ogni frame.
#Ogni riga ha il seguente formato:
#frame atom1_id res1_name res1_id atom1_name  atom2_id res2_name res2_id atom2_name   ...........
#0       8661    T          273       N3        8577      T        271       O2P     0.287049 4.688220
#L'output è un dizionario
#diz[(res1,res2)=frequenza

def group_values(filename):
  hbond={}
  local={}
  resname={}
  prev_frame=-1
  tot_frame=0
  for line in  f:  
    flags=line.split()
    
    frame=int(flags[0])
    res1 =int(flags[3])
    res2 =int(flags[7])
    
    resname[res1]=flags[2]
    resname[res2]=flags[6]
    
    if frame<>prev_frame:
      prev_frame=frame
      tot_frame+=1
      for k in local.keys():
	try:
	  hbond[k]+=1
	except KeyError:
	  hbond[k]=1	
      local={}
      stderr.write("\rframe    %d  " %(frame))
      
      
    if res1<=res2:
      local[res1,res2]=1
    else:
      local[res1,res2]=1
      
  stderr.write("\n")   
  return hbond
