# #!/software/python-2.7-el6-x86_64/bin/python
"""
Created on Mon Mar 24 12:42:46 2014

The main script file for a neus implementation. 

author: Jeremy Tempkin
"""

import sys
from time import time
import fileIO
import copy
import numpy as np
import errors
import basisFunctions
import lammpsWalker


def main():
    """
    definition of main scripting file for debugging purposes.
    """
    # time the execution
    starttime = time()
    
    # in the serial version, lets just set rank=0 
    rank = 0
    
    #---------------SET PARAMETERS -----------
    if rank == 0: print "Setting up umbrella sampling parameters."
    
    params = {}
    
    # set the scratch directory for the calculation files
    params['scratchdir'] = "/Users/jtempkin/enhanced_sampling_toolkit/umbrella_sampling/debug_US"
    params['inputFilename'] = "/Users/jtempkin/enhanced_sampling_toolkit/umbrella_sampling/input.diala"
    
    # here we will set the umbrella sampling parameters in the params dict
    params['ncells'] = 144
    params['cellWidth'] = 60.0
    params['nwalkers'] = 1
    params['walkerSteps'] = 10000
    params['stepLength'] = 10
    params['Ftype'] = 'transition'
    
    # lets set the dynamics parameters that are needed to specify the walker 
    params['temperature'] = 310.0

    #--------------- INITIALIZATION--------
    # only allow root rank to build files
    if rank == 0: 
        print "Building scratch directory."
        # construct the wkdir. This is where the intermediate dynamcs files will 
        # be written. 
        fileIO.makeWkdir(params['scratchdir'])

    # construct the umbrella data structure
    if rank == 0: print "Initializing the umbrella structure."
    
    # create the partition object 
    system = basisFunctions.partition(params['ncells'])


    # now we construct the umbrella windows
    if rank == 0: print "Building the umbrellas."
    
    # specify the list of boundaries, nboxes, 1/2 width of the windows
    umbParams = {}
    
    # right now, we will hardcore a 12x12 array using Erik's routine for gridding a space
    umbParams['cvrange'] = np.array([map(float,entry.split(",")) for entry in "-180,180,12,30;-180,180,12,30".split(";")])
    umbParams['wrapping'] = np.array(map(float, "1,1".split(',')))
    
    system.umbrellas = fileIO.createUmbrellas(umbParams)
    
    #------------ MAIN LOOP --------------
    # sample umbrellas and construct F
    if rank == 0: print "Entering main loop."
    
    for i in range(len(system.umbrellas)):
        print "Rank", rank, "sampling umbrella", i, "."
        
        # lets instantiate a walker object to sample this window. 
        wlkr = lammpsWalker.lammpsWalker(params['inputFilename'])
        
        # minimize structure prior to dynamics
        wlkr.minimize()
        
        # set the dynamics to sample by langevin
        wlkr.command("fix 1 all nve")
        wlkr.command("fix 2 all langevin 310.0 310.0 30.0 20874")
    
        
        # set colvars for this walker (currently, alanine dipeptide dihedrals)
        wlkr.colvars.append(['dihedral', 5, 7, 9, 15])
        wlkr.colvars.append(['dihedral', 7, 9, 15, 17]) 
        
        # set an array of starting/stoping restraints for equilibration
        restraint = [[0.0, 100.0], [0.0, 100.0]]
        
        # now specify colvars to dynamics routines
        wlkr.setColvars()
        
        # equilibrate the walker to the target point in CV space
        wlkr.equilibrate(system.umbrellas[i].center, restraint, 100000)
        

        # enter the sampling routine 
        try: 
            system.sample(wlkr, params['walkerSteps'], i, 0, params, rank)
        except errors.DynamicsError:
            print "Rank", rank, "sampling error occured in umbrella", i, "."
            continue

        
        # now we are done populating the samples array, close the walker
        wlkr.close()
    

    #----------------WRITE OUT DATA-------------
    # now allow rank 0 to process data. 
    if rank == 0:
        print system.F
        fileIO.writeMat(system.F, params['wkdir'] + "/F.out")
    
        print "Entering eigenvalue routine."
        # solve eigenvalue problem for F
        system.getZ()
        fileIO.writeMat(system.z, params['wkdir'] + "/z.out")
    
        print "Computing the sensitivities."
        bounds = system.getlogBound(system.F)
        fileIO.writeMat(bounds, params['wkdir'] + "/bounds.out")


    # now we will perform an analysis of the data and increase sampling of 
    # windows with high variance

    if rank == 0:
        print "Done!"
        print "Total wallclock time was: " + str(time() - starttime) + " seconds."

    
    return 0

# run as a main file if called. 
if __name__ == '__main__':
    main()
