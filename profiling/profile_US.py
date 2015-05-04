# -*- coding: utf-8 -*-
"""
This is the profiling application to test the usage of the partition module for an overlap US calculation. 
"""

import sys
import pickle
import random
import os
from time import time
import fileIO
import copy
import numpy as np
import lammpsWalker
import partition_mod
import h5py
import observables
import dynamics

def setParameters(params):
    """
    For organization of main(), I've set up a  separate routine to initialize the dictionary of parameters.
    """
    # set the scratch directory for the calculation files
    params['scratchdir'] = os.path.abspath("data/profile_US.scratch")
    params['inputFilename'] = os.path.abspath("data/input.enk")
    params['systemFile'] = os.path.abspath("data/system.24x24")
    params['startingPoints']  = os.path.abspath("data/24x24.startingPoints")
    params['logFilename'] = params['scratchdir'] + "/log"
    
    # here we will set the umbrella sampling parameters in the params dict which passes these variables
    # NOTE THAT SOME OTHER LAMMPS VARIABLES ARE SPECIFIED IN THE INPUT FILE
    params['ncells'] = 24**2
    params['walkerSteps'] = 1000
    params['iterations'] = 1
    params['stepLength'] = 5
    params['basisType'] = "Pyramid"
    params['transitionMatrixType'] = 'overlap'
    params['cvrange'] = np.array([map(float,entry.split(",")) for entry in "-180.0,180.0,24,15;-180.0,180.0,24,15".split(";")])
    params['wrapping'] = np.array(map(float, "1,1".split(',')))

    params['temperature'] = 310.0
    params['damping_coefficient'] = 30.0 # in LAMMPS UNITS HERE
    
    return params
    
def dragWalker(wlkr, umbrella, nSteps=25):
    """
    This routine drags a walker toward a point in collective variable space.
    """
    # now let's prepare the walker in the target window by dragging the 
    restraint = [[1000.0, 1000.0],[1000.0, 1000.0]]
    config = wlkr.getColvars()    
    interp = []
    for i in range(len(config)):
        interp.append(np.linspace(config[i], umbrella.center[i], nSteps))
        
    v = np.asarray(interp)
    
    for i in range(nSteps):
        wlkr.equilibrate(v[:,i], restraint, "10000")
        
    # once we have the system here let's equilibrate it at the target, tightening the restraints here. 
    restraint = [[1000.0, 10000.0],[1000.0, 10000.0]]
    wlkr.equilibrate(umbrella.center, restraint, "100000")
    
    

"""
---------------------------------------
                MAIN
---------------------------------------
"""
def main(debug=False):
    """
    definition of main scripting file for debugging purposes.
    """
    # time the execution
    startTime = time()
    
    if debug: print "WARNING: Debug passed as True. Will dump maximum set of internal values."

    #---------------SET SIMULATION PARAMETERS -----------

    params = {}
    
    params = setParameters(params)
    
    print "parameters initialized as: ", params
        
    print "Building scratch directory."
    # construct the wkdir. This is where the intermediate dynamics files will be written.
    fileIO.makeWkdir(params['scratchdir'])

    # make sure we've built the scratchdir since we'll assume it exists for the rest of the code
    assert os.path.exists(params['scratchdir']), "The scratch directory was not created."
    
    # we're going to move the wokring directory into the data directory
    os.chdir(os.path.dirname(params['inputFilename']))
 

    #--------------- INITIALIZATION OF PARTITION INSTANCE--------
    # construct the umbrella data structure
    # commented out below this block is the origional construction lines for the partition module. 
    # here we're simply loading this from a binary pickle written to file. This obviates the need to 
    # reconstruct the neighborlist every call
    if os.path.exists(params['systemFile']):
        with open(params['systemFile'], "r") as f_handle:
            system = pickle.load(f_handle)
    else:

        system = partition_mod.partition(params['ncells'])

        if debug: print system

        # now we construct the umbrella windows
        system.umbrellas = system.createUmbrellas(params['cvrange'], params['wrapping'], basisType=params['basisType'])
        
        with open(params['systemFile'], "w") as f_handle:
            pickle.dump(system, f_handle, protocol=2)

    # double check we've created as many windows as we said in ncells
    assert len(system.umbrellas) == params['ncells'], "The number of cells specified does not equal the number of umbrellas created."

    # debug prints
    if debug:
        for i in range(len(system.umbrellas)):
            print "umbrella", i, system.umbrellas[i].center, system.umbrellas[i].width


    # ------------------ SETTING THE DYNAMICS PARAMETERS FOR SAMPLING THE WINDOWS ----------------
    
    # now set up partition with a list of associated ranks
    system.rank_window_index = range(len(system.umbrellas))
    
    # -------------- INITIALIZATION OF SYSTEM OBSERVABLES -------
    # here we add a list of observables to compute
    pmf_prototype = observables.pmf("pmf", np.zeros((180,180)))
    system.addObservable(pmf_prototype)
    """
    #  now we check to make sure we've correctly set up the local observables
    for window in system.umbrellas:
        assert hasattr(window, "local_observables"), str(rank) + ": The local observables in window " + str(system.umbrellas.index(window)) + " was not created."

    for window in system.umbrellas:
        assert len(window.local_observables) == len(system.observables), str(rank) + ": The number of window observables does not match the number of partition observables."
    """

    #----------------MAIN LOOP----------------
    print "Sampling via equilibrium US."
    for k in range(params['iterations']):

        # ---------- RESET DATA USED FOR CURRENT ITERATION -------------
        # reset the M matrix to zeros for this iteration of the algorithm
        # this matrix will get populated from the inner loop
        system.M = np.zeros((len(system.umbrellas), len(system.umbrellas)))
        system.nsamples_M = np.zeros((len(system.umbrellas), len(system.umbrellas)))
        
        # start a timer on root rank
        st = time()

        # this loop does one iteration of the update
        for i in system.rank_window_index:
            print "sampling window", i
            # lets instantiate a walker object to sample this window.

            wlkr = lammpsWalker.lammpsWalker(params['inputFilename'], params['logFilename'], index=i, debug=debug)
            
            # set langevin dynamics object and pass it to the walker.
            dyn = dynamics.langevin(params['temperature'], params['damping_coefficient'], shake = True)            
            wlkr.setDynamics(dyn)
            
            wlkr.setTimestep(1.0)
            
            wlkr.drawVel()

            # set colvars for this walker
            wlkr.addColvars('c1', "dihedral", [9, 30, 37, 44])
            wlkr.addColvars('c2','dihedral', [30, 37, 44, 64])
            
            #wlkr.setOutput("traj.out", "dcd", params['scratchdir'] + "/" + str(i) + ".dcd", nSteps=1000)
            #let's add an output for writing out the configurations in hdf5 format so we have an estimate for the local distribution
            if os.path.exists(params['scratchdir'] + "/ep." + str(i) + ".h5py"):
                f_handle = h5py.File(params['scratchdir'] + "/ep." + str(i) + ".h5py", "r")
                config = random.choice(f_handle.keys())
                wlkr.setConfig(f_handle[config][:])
                wlkr.drawVel()
                f_handle.close()
            
            else:
                try:
                    if os.path.exists(params['startingPoints'] + "/start." + str(i) + ".npy"):
                        
                        config = np.load(params['startingPoints'] + "/start." + str(i) + ".npy")
                        
                        wlkr.setConfig(config)
                        
                        wlkr.drawVel()
                    else:
                        
                        dragWalker(wlkr, system.umbrellas[i], nSteps=5)
                
                        assert system.umbrellas[i](wlkr.getColvars(), system.umbrellas) > 0.0, "The walker in " + str(i) + " did not initialize in the support of the window."  
                    
                        np.save(params['startingPoints'] + "/start." + str(i), wlkr.getConfig())
                        
                except AssertionError:  
                    pass
            
            # enter the sampling routine. 
            system.sample_US(wlkr, params['walkerSteps'], i, k, params, debug=True)
            
            # now we are done populating the samples array, close the walker
            wlkr.close()
            del wlkr

            if debug: print system.umbrellas[i].local_observables[0].data

        # report time it too for root to move out of the barrier (ie all ranks have sampled)
        print "Elapsed time in sampling routine on root rank after barrier: ", time() - st

        # now let's time the communication section
        st = time()

        # ------------ EIGENVALUE PROBLEM --------------        

        """
        Step 2) solution of the eigenvalue problem at each rank
        """
        # at rank, compute G and z
        system.updateF('overlap')

        system.k += 1
        
        system.computeZ()
        
        if debug: print system.z
        
        """
        Step 3) estimation of observables on each rank
        """
        system.computeObservables()
        

        """
        Step 4) all reduction of averaged observables to each processor
        """
        if debug: print system.observables[0].data

        """
        Step 5) writing local and global estimates of the observables on root rank only
        """

        f_handle = h5py.File(params['scratchdir'] + "/F.out", "a")            
        f_handle.create_dataset("F." + str(k), data=system.F)
        f_handle.close()
        
        f_handle = h5py.File(params['scratchdir'] + "/z.out", "a")            
        f_handle.create_dataset("z." + str(k), data=system.z)
        f_handle.close()
        
        for obs in system.observables:
            f_handle = h5py.File(params['scratchdir'] + "/" + obs.name, "a")
            f_handle.create_dataset(obs.name + "." + str(k), data=obs.data)
            f_handle.close()
        

        print "Elapsed time in communication routine on root rank after barrier: ", time() - st

    print "Done!"
    print "Total wall clock time was: " + str(time() - startTime) + " seconds."


    return 0

# run as a main file if called.
if __name__ == '__main__':
    if len(sys.argv) == 2:
        if sys.argv[1] in ["debug","Debug"]:
            main(debug=True)
        else:
            main()
    else: main()

