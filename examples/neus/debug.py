import lammpsWalker
import random 
import numpy as np
import ctypes


wlkr = lammpsWalker.lammpsWalker("input.diala", "log", index=0)
        
wlkr.setTimestep(1.0)
        
# set the dynamics to sample by langevin
wlkr.command("fix 1 all nve")
# the langevin fix here sets the temperature to 310K and the friction
# coefficient to 30 ps-1
wlkr.command("fix 2 all langevin 310.0 310.0 30.0 20874")    
        
# set colvars for this walker (currently, alanine dipeptide dihedrals)
wlkr.colvars.append(['dihedral', 5, 7, 9, 15])
wlkr.colvars.append(['dihedral', 7, 9, 15, 17]) 

# now specify colvars to dynamics routines
wlkr.setColvars()

# now we initialize the starting coordinates from the entry points library
data = np.load("entryPoints/in_1_w0.entryPoints.npy")

temp_indx = random.randint(0, data.shape[0]-1)
wlkr.setConfig(data[temp_indx].ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        
cvs = wlkr.getColvars()

print cvs

wlkr.command("run 1000 post no")

cvs = wlkr.getColvars()

print cvs