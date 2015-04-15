# -*- coding: utf-8 -*-
"""
Head unit testing module for the Enhanced Sampling Toolkit Package.

This file executes a series of Unit Tests written for the lammpsWalker package. 
"""

# import header
import pytest
import time
import subprocess
import numpy as np
# now we will proceed through our various module tests. 
"""
Let's outline the tests for the lammpsWalker module:

We'll test the initialization of a lammpsWalker object, verify that we got an object that is a lammpsWalker that inherits from the walker base. We'll move through the API for the lammpsWalker and test each routine using a particle on a flat surface system and also for an alanine dipeptide system. This should give a good unit test on the functionality in the walker API level of the code.
"""

# write out the current Python environment to a file for a record. 
subprocess.call("pip list > pyenv.test", shell=True)

@pytest.fixture
def dialaSetup():
    """
    This routine sets up the data structures needed for the following unit tests. This includes a copy of a LAMMPS walker with the alanine dipeptide and a 2D particle walker. 
    """
    import lammpsWalker
    wlkr = lammpsWalker.lammpsWalker("input.diala","log.diala.test")
    ref_config = np.loadtxt("diala.structure").flatten()
    ref_vel = np.loadtxt("vel.ref").flatten()
    
    return wlkr, ref_config, ref_vel    

def test_initialization(dialaSetup):
    """
    This test initiallizes a LAMMPS walker and tests if the object returned by the initialization routine is of the correct type and inhereted from the lammpsWalker class. 
    """
    import lammpsWalker
    
    wlkr, ref_config, ref_vel = dialaSetup
    
    assert type(wlkr) == lammpsWalker.lammpsWalker, "walker object type not a lammps walker."
    
    assert isinstance(wlkr, lammpsWalker.lammpsWalker), "Lammps walker did not inheret from lammpsWalker class."

def test_getConfig(dialaSetup):
    """
    1) check that the getConfig returns a numpy array that is equivalent to the reference array from the data file. 
    """
    wlkr, ref_config, ref_vel = dialaSetup
    
    config = wlkr.getConfig()
    
    assert np.array_equal(config, ref_config), "The configuration from the walker did not match the reference."
    
def test_getVel(dialaSetup):
    """
    2) check that the velocities read in from the data file and returned by getVel match the reference.
    """ 
    wlkr, ref_config, ref_vel = dialaSetup
    
    vel = wlkr.getVel()
    
    assert np.array_equal(vel, ref_vel), "The velocity of the walker does not match the one gotten from getVel()"
    
    
def test_setConfig(dialaSetup):        
    """
    3) test the setConfig and setVel routines by making a change to the velocities and sending that data to the walker. We then retrieve those vectors to confirm that they are set correctly. 
    """    
    wlkr, ref_config, ref_vel = dialaSetup
    
    # set up walker 
    config = wlkr.getConfig()
    
    # now scale the positions and resend to walker. 
    config = config * 2.0
    
    wlkr.setConfig(config)
    
    config_pert = wlkr.getConfig()
    
    assert np.array_equal(config_pert, ref_config * 2.0), "The setConfig routine did not correctly set the configuration in the walker."

def test_setVel(dialaSetup): 
    """
    4) now check that reversing the velocities restores them back to the original. 
    """
    wlkr, ref_config, ref_vel = dialaSetup
    
    vel = wlkr.getVel()
    
    # reverse walker velocities. 
    vel = vel * -1.0
    
    wlkr.setVel(vel)
    
    assert np.array_equal(vel, ref_vel * -1.0), "the set Velocities command did not return correctly."

def test_reverseVel(dialaSetup):
    """
    5) now test whether the reverse velocities routine works as expected. 
    """
    wlkr, ref_config, ref_vel = dialaSetup
    
    wlkr.reverseVel()
    
    vel = wlkr.getVel()
    
    assert np.array_equal(vel, ref_vel * -1.0), "the reverse Velocities command did not correctly reverse the velocities."
    
def test_drawVel(dialaSetup):
    """
    6) test whether the draw velocity routine has redrawn the velocities. 
    """    
    wlkr, ref_config, ref_vel = dialaSetup
    
    ref_vel = np.loadtxt("gaussian_vel")
    
    wlkr.drawVel(seed=123456)
    
    vel = wlkr.getVel()
    
    assert np.array_equal(vel, ref_vel), "The draw velocities with seed 123456 did not return correctly."
    
    
def test_colvars():
    """
    This test checks the colvars routine using a set of test collective varaibles on the alanine dipeptide. 
    """

    
def test_equilibration():
    """
    This test checks the equilibration routine for the walker using a simple dynamics setup on a single particle system on a 2D potential. 
    """
    
    #import lammpsWalker
    
    #wlkr = lammpsWalker.lammpsWalker()
    
    #wlkr.close()
    
def test_minimize():
    """
    This test checks the minimization routine for the walker class. This routine uses a simple point particle on a 2D surface problem so that we have a good handle on how the minimization will work. 
    """
    # 