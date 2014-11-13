#!/software/python-2.7-el6-x86_64/bin/python
"""
Created on Mon Mar 24 12:42:46 2014

The main script file for the umbrella sampling code. 

author: Jeremy Tempkin
"""

import sys
from time import time
import fileIO
import copy
from mpi4py import MPI
import errors
import basisFunctions
import lammpsWalker


def main():
    """
    definition of main scripting file for debugging purposes.
    """
    # time the execution
    starttime = time()
    
    #--------------MPI INIT---------------------
    # init communicator
    comm = MPI.COMM_WORLD
    # get proc rank 
    rank = comm.Get_rank()
    nprocs = comm.Get_size()
    
    #---------------SET PARAMETERS -----------
    if rank == 0: print "Setting up umbrella sampling parameters."
    
    params = {}
    
    # set the scratch directory for the calculation files
    params['scratchdir'] = "/Users/jeremytempkin/Documents/enhanced_sampling_toolkit/umbrella_sampling/debug_US"
    
    # here we will set the umbrella sampling parameters in the params dict
    params['ncells'] = 3
    params['cellWidth'] = 60.0
    params['nwalkers'] = 1
    params['walkerSteps'] = 1000
    params['stepLength'] = 10
    params['Ftype'] = 'transition'
    
    # lets set the dynamics parameters that are needed to specify the walker 
    
    params['inputFilename'] = 'syn.in'
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
    umbs = []
    # these are for the two-dimensional landscape of pairwise distances in the 
    # synuclein protein based off PRE data 
    umbs.append([0, 350, 70, 5])
    umbs.append([0, 70, 14, 5])
    
    system.umbrellas = system.createUmbrellas(umbs)
    
    #------------ MAIN LOOP --------------
    # sample umbrellas and construct F
    if rank == 0: print "Entering main loop."
    
    for i in range(rank, len(system.umbrellas), nprocs):
        print "Rank", rank, "sampling umbrella", i, "."
        
        # lets instantiate a walker object to sample this window. 
        wlkr = lammpsWalker.lammpsWalker(params['inputFilename'])
        
        # equilibrate the walker to the target point in CV space
        
        # set colvars for recording samples
        
        try: 
            system.sample(wlkr, params['walkerSteps'], i, 0, params, rank)
        except errors.DynamicsError:
            print "Rank", rank, "sampling error occured in umbrella", i, "."
            continue
        
        # now we are done populating the samples array, close the walker
        wlkr.close()
    
    #-----------------MPI COMMUNICATION------------
    """
    Here we communicate the rows of F to each processor to set up and solve the 
    eigenvalue problem. 
    """
    
    # now reduce the F matrix at root, first making a send buffer, 
    system.Fbuff = copy.deepcopy(system.F)
    comm.Reduce([system.Fbuff, MPI.DOUBLE], [system.F, MPI.DOUBLE], op=MPI.SUM, root=0)
        
    # also share the F_error matrix
    system.F_error_buff = copy.deepcopy(system.F_error)
    comm.Reduce([system.F_error_buff, MPI.DOUBLE], [system.F_error, MPI.DOUBLE], op=MPI.SUM, root=0)
    
    # issue a global MPI barrier to ensure data has been recieved. 
    comm.Barrier()
    
    
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
