# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 16:21:03 2014

This file implements the LAMMPS walker abstraction layer.


@author: Jeremy Tempkin
"""

import random
import sys
import walker

from lammps import lammps


class lammpsWalker(walker.velocityWalker):
    """
    This class implements the lammps bindings into the walker class.
    """
    def __init__(self, inputFilename, logFilename, index = 0):
        """
        Initializes the walker based on an input file for lammps.
        """
        # a general filename for storing the data, note the appending of the index
        self.filename = inputFilename + "." + str(index)
        
        # start walker object from the passed lammps input file name
        self.lmp = self.initLAMMPS(inputFilename=inputFilename, logFilename=".".join([logFilename, str(index)]))
        
        # a list of the relevant collective variables
        self.colvars = []
        
        # a list of commands used to equilibrate the system 
        self.equilCmds = []
        
        # by default, execute shake code. This flag lets the walker know
        # that shakeH will be used in the dynamics
        self.shakeH = True
        
        # here is the list of dynamics properties the walker should know:
        
        self.temperature = 310.0
        
    
        # walker index         
        self.index = index
        
    
    def close(self):
        """
        This function closes the lammps object. 
        """
        self.lmp.close()
        return 0 
        
    def initLAMMPS(self, inputFilename=None, logFilename="log"):
        """
        This function initializes a lammps simulation given an input file
        and the system parameters. Calls the setupLAMMPS internal routine. 
        
        """
        # initialize the lammps object
        #args = ["-sc", "none", "-echo", "log"]
        #self.lmp = lammps("", args)
        
        self.lmp = lammps()
        
        # after the lammps object is created, initialize the lammps simulation. 
        # specify general log file parameters
        self.lmp.command("echo none")

        # if there is a specified filename, use it to set up the simulation.         
        if inputFilename != None:
            self.lmp.command("log " + logFilename)
            self.lmp.file(inputFilename)
        else:
            self.lmp.command("log " + logFilename)
        
        #self.__setupLAMMPS(filename)

        return self.lmp    

    def __setupLAMMPS(self, filename=None):
        """
        The input file sets up and initializes the LAMMPS system from a setup
        file. The setup file should be a simple text file written in LAMMPS
        input file syntax.
        """
        # specify general log file parameters
        self.lmp.command("echo none")

        # if there is a specified filename, use it to set up the simulation.         
        if filename != None:
            self.lmp.command("log " + filename + ".log")
            self.lmp.file(filename)
        else:
            self.lmp.command("log lammps.log")
        
        """
        # set default values for keywords:
        self.dynaType = 'langevin'
        self.temperature = 310.0
        self.outputTrajFreq = 1000                  
        
        with open(sysParams['datadir'] + "/" + sysParams['setupFile'], "r") as ifile:
            while True:
                # get the next line
                line = ifile.readline()
                # check to see if the next line is the EOF, if so break the loop
                if line == '':
                    break
                # ignore emtpy lines
                if not line.strip():
                    continue
                # ignore comments 
                if line.strip()[0] == '#':
                    continue
                # if not skipped from above, parse the information 
                linesplit = line.split()
                
                # first parse for general parameter flags       
                # if you hit the equilibration section, append the lines
                if linesplit[0] == 'equilibrate':
                    while '}' not in line:
                        # get the following line
                        line = ifile.readline()
                        #check for comments
                        if line.split()[0] == "#":
                            continue
                        # get the next line and append it
                        if line.strip().strip('}'):
                            self.equilCmds.append(line.strip().strip('}'))
                elif linesplit[0] == 'colvars':
                    while '}' not in line:
                        # get the following line
                        line = ifile.readline()
                        #check for comments
                        if line.split()[0] == "#":
                            continue
                        # get the next line and append it
                        if line.strip().strip('}'):
                            self.colvars.append(line.strip().strip('}').split())
                # parse for specific keywords
                elif linesplit[0] == 'dynamics':
                    self.dynaType = linesplit[1]
                elif linesplit[0] == 'temperature':
                    self.temperature = float(linesplit[1])
                elif linesplit[0] == 'outputTrajFreq':
                    self.outputTrajFreq = int(linesplit[1])
                # note that if the read data command is found, read in the molecular data
                # this ordering was done to fix errors in the sequence of commands read into LAMMPS
                elif linesplit[0] == 'read_data':
                    self.lmp.command(line)
                elif linesplit[0] == 'shake':
                    if linesplit[1] == 'off':
                        self.shakeH = False
                else: 
                    # if it is a general lammps command, execute
                    self.lmp.command(line)
        
        # now set up dynamics fixes for the specified dynamics
        if self.dynaType == 'langevin':
            # langevin dynamics fixes
            self.lmp.command("fix 1 all nve")
            self.lmp.command("fix 2 all langevin " + str(self.temperature) + " " + str(self.temperature) + " 30.0 " + str(random.randint(100000,999999)))
            # initialize velocities
            self.lmp.command("velocity all create " + str(self.temperature) + " " + str(random.randint(100000,999999)) + " dist gaussian")
        else:
            print "Dynamics type not understood from", sysParams['setupFile']
            sys.exit(0)
        
        # set colvars and computes
        self.initColVars()
        # write out the configurations to a file. 
        self.lmp.command("dump 3 all dcd " + str(self.outputTrajFreq) + " " + filename + ".dcd")
        """
        
        
        
        return 0 
        
    def setColvars(self):
        """
        This function initializes the collective variable and trajectory output 
        for a LAMMPS simulation that is handed to this object. 
        """        
        cv = []
        for entry in self.colvars:
            cv.append(map(str, entry))
        
        # initialize the groups and the computes associated with them
        for index,entry in enumerate(cv):  
            if entry[0] == 'bond':
                self.lmp.command("group " + "b" + str(index) + " id " + " ".join(entry[1:]))
                self.lmp.command("compute " + str(index) + " b" + str(index) + " bond/local dist" )
            elif entry[0] == 'angle':
                self.lmp.command("group " + "a" + str(index) + " id " + " ".join(entry[1:]))
                self.lmp.command("compute " + str(index) + " a" + str(index) + " angle/local theta")
            elif entry[0] == 'dihedral':
                self.lmp.command("group " + "d" + str(index) + " id " + " ".join(entry[1:]))
                self.lmp.command("compute " + str(index) + " d" + str(index) + " dihedral/local phi")        
        
        return 0

        
    def destroyColvars(self):
        """
        This function removes the colvars set by initColVars().
        """        
        for cv in self.colvars:
            self.lmp.command("uncompute " + str(self.colvars.index(cv)))
            
        return 0    
                
    def equilibrate(self, center, restraint, numSteps):
        """
        This function prepares a LAMMPS image to be at the specified target 
        position given by the vector 'center' passed and an arguments. 
        """
        
        # first enter the restraints based on computes data structure
        restCommand = "fix REST all restrain "
        
        # here we apply the constraints based on the collective variable definition
        # check to make sure we have collective variables defined. 
        if len(self.colvars) == 0:
            print "There are no collective variables defined in walker " + str(self.index)
            print "Aborting equilibration of walker " + str(self.index)
            return 0 
            
        # let's get a string representation of the colvars, it's friendly
        cv = []
        for entry in self.colvars:
            cv.append(map(str, entry))
            
        # now loop through each 
        for index, entry in enumerate(cv):
            if entry[0] == 'dihedral':
                restCommand += " " + " ".join(entry) + " " + " ".join(map(str, restraint[index])) + " " + str(center[index] + 180.0)
            else:
                restCommand += " " + " ".join(entry) + " " + " ".join(map(str, restraint[index])) + " " + str(center[index])
        
        # now issue restraint definition to the lammps object 
        self.command(restCommand)
                    
        # now run the equilibration dynamics     
        self.command("run " + str(numSteps) + " post no")
        
        """
        # apply SHAKE if used
        if self.shakeH:    
            self.lmp.command("fix 10 all shake 0.0001 500 0 m 1.008")
        """ 
        
        # now remove constraints for subsequent dynamics         
        self.command("unfix REST")
        
        return 0
        
    def getConfig(self):
        """
        This function returns the current position of the LAMMPS simulation.
        """

        return self.lmp.gather_atoms("x",1,3)

    def getColvars(self):
        """
        This function returns the current position of the LAMMPS simulation in 
        colvars space.
        """
        # get an empty array with placeholders
        cvarray = [None]*len(self.colvars)
        
        # now get cv's one by one from each compute defined
        for i in range(len(self.colvars)):
            cvarray[i] = self.lmp.extract_compute(str(i), 2, 1)[0]

        return cvarray
    
    def setConfig(self, config):
        """
        This routine sets the internal configuration. 
        """        
        self.lmp.scatter_atoms("x", 1, 3, config)
        return 0 
    
    def drawVel(self, distType = 'gaussian', temperature = 310.0):
        """
        This function redraws the velocities from a maxwell-boltzmann dist.
        """
        if distType == 'gaussian':
            self.lmp.command("velocity all create " + str(temperature) + " " + str(random.randint(100000,999999)) + " dist gaussian")
        else:
            print "The drawVel routine was passed distType that was not understood."
            
        return 0
    
    def reverseVel(self):
        """ 
        This function reverses the velocities of a given LAMMPS simulation
        """
        
        # set varibales for reversing velocities
        self.lmp.command("variable vx atom -vx")
        self.lmp.command("variable vy atom -vy")
        self.lmp.command("variable vz atom -vz")
        
        # set new velocities
        self.lmp.command("velocity all set v_vx v_vy v_vz")
        
        return 0     
        
    def propagate(self, numSteps, pre='no', post='no'):
        """
        This function issues a run command to the underlying dynamics to propagate
        the dynamics a given number of steps. 
        """
        
        self.lmp.command("run " + str(numSteps) + " pre " + str(pre) + " post " + str(post))
        
        return 0 
    
    def setDynamics():
        """
        This routine sets the dynamics for the walker. 
        """
        
        return 0
    
    def command(self, command):
        """
        This function allows the user to issue a LAMMPS command directly to the
        LAMMPS object.
        """
        
        # issue the command directly. 
        self.lmp.command(command)
        
        return 0 
        
    def minimize(self, args=None):
        """
        This function runs a minimization routine with the specified type.
        """
        # use default settings
        if args == None:
            self.lmp.command("minimize 1.0e-4 1.0e-6 100 1000")
        else:
            self.lmp.command("minimize " + " ".join(args))
        
        return 0 
        
    def setTemperature(self, temp):
        """
        This function sets the temperature of the walker object.
        
        NOTE THAT THIS DOES NOT ALTER THE DYNAMICS THERMOSTAT. LAMMPS REQUIRES
        RESETING THIS THERMOSTAT. WE WILL LOOK INTO HOW TO DO THIS. 
        """
        self.temperature = temp
        
        return 0
        
    def setTimestep(self, timestep):
        """
        This routine sets the dynamics time step.
        """
        self.lmp.command("timestep " + str(timestep))
        
        return 0
        
        
        
    