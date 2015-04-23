# -*- coding: utf-8 -*-
"""
This file implements the LAMMPS walker abstraction layer. The core of this idea is that 
"""

import random
#import sys
from walker_base import walker
import collectiveVariables
import outputClass
import numpy as np
import ctypes
import os.path
#from lammps import lammps


class lammpsWalker(walker):
    """
    This class implements the enhanced sampling walker API for the bindings to the LAMMPS package. 
    
    Some usage issues to note:
    
    1) Collective variables (CVs) are defined to the walker by constructing a list of
    CVs internally in the walker.colvars object. These CVs list takes the following
    format:

    
    ["coordinateType", atomids....]
    
    The coordinate type specifies which type of coordinate the CV is 
    (i.e. bond, angle, dihedral, etc.). The next items in the list are the atom 
    indecies involved in this specific instance of the CV.
    
    The walker will use this list to initialize them to the underlying LAMMPS objects. 
    """
    
    def __init__(self, inputFilename, logFilename, index = 0, debug = False):
        """
        Initializes a walker object. Takes the following input:
        
        * a filename that builds the model from a LAMMPS script.
        * a filename to pipe standard LAMMPS output into.
        * an index for the walker. (default=0)
        * debug flag to open writing verbose output. Importantly, this pipes LAMMPS output to standard output. (default=False)
        
        The __initLAMMPS__() routine will check for LAMMPS being an importable library from Python and raise an execption if it cannot be loaded. 
        
        """
        if inputFilename is not None:
            # check to see if there is actually a filename at that location.
            assert os.path.isfile(inputFilename), "The input file for the LAMMPS walker could not be found."
            
        # a filename for storing the data, note the appending of the index
        self.filename = inputFilename + "." + str(index)
        
        # walker index         
        self.index = index
        
        # start walker object from the passed lammps input file name
        self.lmp = self.__initLAMMPS__(inputFilename = inputFilename, logFilename = ".".join([logFilename, str(index)]), verbose = debug)
        
        # a list of the relevant collective variables
        self.colvars = []
        
        # we're constructing a list for collective output classes. 
        self.output = []
        
        self.dynamics = None
        
        # a list of commands used to equilibrate the system 
        self.equilCmds = []
        

    def close(self):
        """
        This function closes the LAMMPS object. 
        """
        self.lmp.close()
        return 0 
        
    def __initLAMMPS__(self, inputFilename=None, logFilename="log", verbose = False):
        """
        This function initializes a LAMMPS simulation given an input file. 
        
        Takes the following arguments:
        
        * input file name
        * 
        
        
        """
        try:
            from lammps import lammps
        except:
            print "ERROR: The LAMMPS python interface could not be found. Check your PYTHONPATH and compiled LAMMPS directory to make sure it is compiled and importable." 
        
        # initialize the lammps object
        if verbose == True:
            self.lmp = lammps()
        else:
            args = ["-sc", "none", "-echo", "log"]
            self.lmp = lammps("", args)
        
        # after the lammps object is created, initialize the lammps simulation. 
        # specify general log file parameters
        self.command("echo none")

        # if there is a specified filename, use it to set up the simulation.         
        if inputFilename != None:
            self.command("log " + logFilename)
            self.lmp.file(inputFilename)
        else:
            self.command("log " + logFilename)
        
        #self.__setupLAMMPS__(filename)

        return self.lmp    

    def addColvars(self, name, cvType, atomIDs):
        """
        Implements the addition of a collective variable to the list of collective variables held by this walker.
        
        The currently support collective variables are:
        * bonds
        * angles
        * dihedrals
        * x, y, z positions
        * x, y, z velocity components
        
        The collective variables are added as an instance of a collective variable object which contains the necessary fields defining that collective variable. 

        """
        __knownCVs__ = ['bond', 'angle', 'dihedral', 'x', 'y', 'z', 'vx', 'vy', 'vz']
        assert cvType in __knownCVs__, "cvType that was provided was not understood." 
        
        for cv in self.colvars:
            assert cv.name != name, "Collective variable names must be unique."
                
        # now append the collective variable to the walker list.      
        self.colvars.append(collectiveVariables.collectiveVariables(name, cvType, atomIDs))
        
        cv = self.colvars[-1]
        # first set the group for the colvar
        self.command("group " + cv.name + " id " + " ".join(map(str,cv.atomIDs)))
            
        # now set the appropriate colvar
        if cv.type == 'bond':
            self.command("compute " + cv.name + " " + cv.name + " bond/local dist" )
        elif cv.type == 'angle':
            self.command("compute " + cv.name + " " + cv.name + " angle/local theta")
        elif cv.type == 'dihedral':
            self.command("compute " + cv.name + " " + cv.name + " dihedral/local phi")        
        elif cv.type == 'x':
            self.command("compute " + cv.name + " " + cv.name + " property/atom x")
        elif cv.type == 'y':
            self.command("compute " + cv.name + " " + cv.name + " property/atom y")
        elif cv.type == 'z':
            self.command("compute " + cv.name + " " + cv.name + " property/atom z")
        elif cv.type == "vx":
            self.command("compute " + cv.name + " " + cv.name + " property/atom vx")
        elif cv.type == "vy":
            self.command("compute " + cv.name + " " + cv.name + " property/atom vy")
        elif cv.type == "vz":
            self.command("compute " + cv.name + " " + cv.name + " property/atom vz")  

        
        return 0
        
    def destroyColvars(self):
        """
        This function removes the colvars set by setColvars(). By default, it removes all of the collective variables in the list. It does not remove them from the collective variables list. 
        """        
        for cv in self.colvars:
            self.command("uncompute " + cv.name)
            
        self.colvars = []
            
        return 0    
                
    def equilibrate(self, center, restraint, numSteps):
        """
        This function prepares a LAMMPS image to be at the specified target position given by the vector 'center' passed and an arguments. 
        """
        #print "Equilibrating walker."
        
        assert len(center) == len(restraint), "The dimension of the center array and the restraint array do not match."
        assert len(center) == len(self.colvars), "The dimension of the center array and the number of collective variables do not match."
        
        
        # first enter the restraints based on computes data structure
        restCommand = "fix REST all restrain "
        
        # here we apply the constraints based on the collective variable definition
        # check to make sure we have collective variables defined. 
        if len(self.colvars) == 0:
            print "There are no collective variables defined in walker " + str(self.index)
            print "Aborting equilibration of walker " + str(self.index)
            return 0 
            
        # let's get a string representation of the colvars, it's friendly
        """
        cv = []
        for entry in self.colvars:
            cv.append(map(str, entry))
        """ 
        # now loop through each 
        for index, entry in enumerate(self.colvars):
            if entry.type == 'dihedral':
                #restCommand += " " + entry.type + " " + " ".join(map(str, entry.atomIDs)) + " " + " ".join(map(str, restraint[index])) + " " + str(center[index])
                restCommand += " " + entry.type + " " + " ".join(map(str,entry.atomIDs)) + " " + " ".join(map(str, restraint[index])) + " " + str(center[index] + 180.0)
            else:
                restCommand += " " + entry.type + " " + " ".join(map(str, entry.atomIDs)) + " " + " ".join(map(str, restraint[index])) + " " + str(center[index])
        
        # now issue restraint definition to the lammps object 
        self.command(restCommand)
        
        self.command("fix_modify REST energy yes")
                    
        # now run the equilibration dynamics     
        self.command("run " + str(numSteps) + " post no")
        
        """
        # apply SHAKE if used
        if self.shakeH:    
            self.lmp.command("fix 10 all shake 0.0001 500 0 m 1.008")
        """ 
        
        # now remove constraints for subsequent dynamics         
        self.command("unfix REST")
        
        # this resets the dynamics environment after the equilibration run
        self.command("run 0 post no")
        
        return 0
        
    def getConfig(self):
        """
        This function returns the current position of the LAMMPS simulation.
        """
        config = np.asarray(self.lmp.gather_atoms("x",1,3))
        
        return config
        
    def getVel(self):
        """
        This function returns the current velocities from the LAMMPS simulation.
        """
        vel = np.asarray(self.lmp.gather_atoms("v", 1, 3))
        
        return vel
        
    def setVel(self, vel):
        """
        This function sets the velocity to the lammps simulation. 
        """
        
        self.lmp.scatter_atoms("v", 1, 3, vel.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        
        return 0

    def getColvars(self):
        """
        This function returns the current position of the LAMMPS simulation in 
        colvars space.
        """
        # get an empty array with placeholders
        cvarray = []
        
        # now get cv's one by one from each compute defined
        for cv in self.colvars:
            if cv.type in ['x', 'y', 'z', 'vx', 'vy', 'vz']: 
                cvarray.append(self.lmp.extract_compute(cv.name, 1, 1)[0])
            else:
                #*** We REALLY need assurance here that what we are getting here is in fact not a NULL POINTER.
                cvarray.append(self.lmp.extract_compute(cv.name, 2, 1)[0])
                
        assert len(cvarray) == len(self.colvars), "Not all collective variables were added."

        return cvarray
    
    def setConfig(self, config):
        """
        This routine sets the internal configuration. 
        """                
        self.lmp.scatter_atoms("x", 1, 3, config.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        
        return 0 
    
    def drawVel(self, distType = 'gaussian', temperature = 310.0, seed = None):
        """
        This function redraws the velocities from a maxwell-boltzmann dist.
        """
        if seed is None:
            seed = random.randint(100000,999999)
        
        if distType == 'gaussian':
            self.command("velocity all create " + str(temperature) + " " + str(seed) + " dist gaussian")
        else:
            print "The drawVel() routine was passed a distribution Type that was not understood."
            
        return 0
    
    def reverseVel(self):
        """ 
        This function reverses the velocities of a given LAMMPS simulation
        """
        
        # set varibales for reversing velocities
        self.command("variable vx atom -vx")
        self.command("variable vy atom -vy")
        self.command("variable vz atom -vz")
        
        # set new velocities
        self.command("velocity all set v_vx v_vy v_vz")
        
        return 0     
        
    def propagate(self, numSteps, pre='no', post='no'):
        """
        This function issues a run command to the underlying dynamics to propagate
        the dynamics a given number of steps. 
        """
        
        self.command("run " + str(numSteps) + " pre " + str(pre) + " post " + str(post))
        
        return 0 
    
    def setDynamics(self, dynamics_instance):
        """
        This routine sets the dynamics for the walker. 
        """        
        __knownDynamics__ = ['langevin']
        
        assert dynamics_instance.type in __knownDynamics__, "Dynamics instance type was not recognized."        
        
        # first we should check to see if there is already a dynamics defined. 
        if self.dynamics is not None:
            # if there is a dynamics already present, remove it's fixes from LAMMPS
            for item in self.dynamics.fixes:
                self.command("unfix " + item)
        
        # now replace with new dynamics instance
        self.dynamics = dynamics_instance        
        
        self.dynamics.fixes = []
        
        # and set the required fixes
        if self.dynamics.type is 'langevin':                        
            # send the fixes to the underlying lammps object 
            self.command("fix 1 all nve")
            if self.dynamics.seed is None:
                self.command("fix 2 all langevin " + " ".join([str(self.dynamics.temperature), str(self.dynamics.temperature)]) + " " + str(self.dynamics.damping_coefficient) + " " + str(random.randint(100000, 999999)))
            else: 
                self.command("fix 2 all langevin " + " ".join([str(self.dynamics.temperature), str(self.dynamics.temperature)]) + " " + str(self.dynamics.damping_coefficient) + " " + str(self.dynamics.seed))
            
            # add the fixes to the internal list 
            self.dynamics.fixes.append("1")
            self.dynamics.fixes.append("2")
            
        # set shake
        if self.dynamics.shake is True:
            self.command("fix shk all shake 0.0001 20 0 m 1.0")
            self.dynamics.fixes.append("shk")
        
        self.propagate(0, pre='yes')
        
        return 0
    
    def command(self, command):
        """
        This function allows the user to issue a LAMMPS command directly to the
        LAMMPS object.
        """
        
        # issue the given command directly. 
        self.lmp.command(command)
        
        return 0 
        
    def minimize(self, args=None):
        """
        This function runs a minimization routine with the specified type.
        """
        # use default settings
        if args == None:
            self.command("minimize 1.0e-4 1.0e-6 100 1000")
        # else use the arguments provided. 
        else:
            self.command("minimize " + " ".join(args))
        
        return 0 
        
    def setTemperature(self, temp):
        """
        This function sets the temperature of the walker object. 
        """
        self.dynamics.temperature = temp
        
        self.setDynamics(self.dynamics)
        
        return 0
        
    def setTimestep(self, timestep):
        """
        This routine sets the dynamics time step.
        """
        self.command("timestep " + str(timestep))
        
        return 0
        
    def setOutput(self, name, outputType, filename, nSteps):
        """
        This routine sets up a mechanism for writing system information to file directly from the dynamics engine. This is equivalent to output constructed 
        """
        # assert outputType in __knownOutput__, "The output type specified was not recognized."        
        __knownOutputTypes__ = ['dcd', 'xyz', 'custom', 'xyz/mpiio', "atom"]
        assert outputType in __knownOutputTypes__, "outputType that was provided was not understood." 
        
        for out in self.output:
            assert out.name != name, "Output class names must be unique."
                
        # now append the collective variable to the walker list.      
        self.output.append(outputClass.outputClass(name, outputType, filename, nSteps))
        
        dump = self.output[-1]
        
        if dump.type == "custom":
            self.command("dump " + dump.name + " all " + " ".join([dump.type, str(dump.nSteps), dump.destination, "x y z"]))
        else:
            self.command("dump " + dump.name + " all " + " ".join([dump.type, str(dump.nSteps), dump.destination]))
        
        
    def removeOutput(self):
        """
        This routine removes the ouput pipes for information to be written to disk from the underlying dynamics engine. Right now this simply clears all existing output. 
        """
        for out in self.output:
            self.command("undump " + out.name)
            
        self.output = []
        
        return 0

        
walker.register(lammpsWalker)        

if __name__ == "__main__":
    print 'The lammpsWalker module is a subclass of walker:', issubclass(lammpsWalker, walker)
    
