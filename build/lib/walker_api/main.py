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
import basisFunctions
from mpi4py import MPI
import errors

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
    
    #---------------READ PARAMETERS -----------
    if rank == 0: print "Reading input."
    sysParams = {}

    fileIO.readSysInput(sys.argv[1], sysParams, trustMe = True)

    if rank == 0: print "Read from ", sys.argv[1]

    #--------------- INITIALIZATION--------
    # only allow root rank to build files
    if rank == 0: 
        print "Building working directory."
        # construct the wkdir. This is where the intermediate dynamcs files will 
        # be written. 
        fileIO.makeWkdir(sysParams)

    # construct the umbrella data structure
    if rank == 0: print "Initializing the umbrella structure."
    system = basisFunctions.partition(sysParams['ncells'])

    if rank == 0: print "Building the umbrellas."
    system.umbrellas = fileIO.createUmbrellas(sysParams)


    #------------ MAIN LOOP --------------
    # sample umbrellas and construct F
    if rank == 0: print "Entering main loop."
    
    for i in range(rank, len(system.umbrellas), nprocs):
        print "Rank", rank, "sampling umbrella", i, "."
        try: 
            system.sample(sysParams['walkerSteps'], i, 0, sysParams, rank)
        except errors.DynamicsError:
            print "Rank", rank, "sampling error occured in umbrella", i, "."
            continue

        # now analyze the time series generated using the acor function
        """
        for j in range(system.F.shape[1]):
            tau, mean, sigma = system.computeAcor(system.umbrellas[i].basisFnxTimeSeries[j])
            system.F_error[i][j] = sigma
        """
    
    #-----------------MPI COMMUNICATION------------
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
        fileIO.writeMat(system.F, sysParams['wkdir'] + "/F.out")
    
        print "Entering eigenvalue routine."
        # solve eigenvalue problem for F
        system.getZ()
        fileIO.writeMat(system.z, sysParams['wkdir'] + "/z.out")
    
        print "Computing the sensitivities."
        bounds = system.getlogBound(system.F)
        fileIO.writeMat(bounds, sysParams['wkdir'] + "/bounds.out")

    
    # now we will perform an analysis of the data and increase sampling of 
    # windows with high variance

    if rank == 0:
        print "Done!"
        print "Total wallclock time was: " + str(time() - starttime) + " seconds."

    
    return 0

# run as a main file if called. 
if __name__ == '__main__':
    main()
