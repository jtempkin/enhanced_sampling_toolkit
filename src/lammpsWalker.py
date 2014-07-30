# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 16:21:03 2014

@author: jtempkin
"""

import random
import sys
import walker
from lammps import lammps


class lammpsWalker(walker.velocityWalker):
    """
    This class implements the lammps bindings into the walker class.
    """
    def __init__(self, filename, sysParams , index = 0):
        """
        Initializes the walker based on an input file for lammps.
        """
        
        # a list of the relevant collective variables
        self.colvars = []
        # a list of commands used to equilibrate the system 
        self.equilCmds = []
        # a general filename for storing the data
        self.filename = filename
        # by default, execute shake code. This flag lets the walker know
        # that shakeH will be used in the dynamics
        self.shakeH = True
        # start walker object
        self.lmp = self.initLAMMPS(filename, sysParams)
        self.index = index
        
    
    def close(self):
        """
        This function closes the lammps object. 
        """
        self.lmp.close()
        return 0 
        
    def initLAMMPS(self, filename, sysParams):
        """
        This function initializes a lammps simulation given an input file
        and the system parameters. Calls the setupLAMMPS internal routine. 
        
        """
        # initialize the lammps object
        args = ["-sc", "none", "-echo", "log"]
        self.lmp = lammps("", args)
        # after the lammps object is created, initialize the lammps simulation. 
        self.setupLAMMPS(filename, sysParams)

        return self.lmp    
    
    def setupLAMMPS(self, filename, sysParams):
        """
        The input file sets up and initializes the LAMMPS system from a setup
        file. 
        """
        # specify general log file parameters
        self.lmp.command("echo none")
        self.lmp.command("log " + filename + ".log")
        
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
        
        return 0 
        
    def initColVars(self):
        """
        This function initializes the collective variable and trajectory output 
        for a LAMMPS simulation that is handed to this object. 
        """        
        
        # initialize the groups and the computes associated with them
        for cv in self.colvars:       
            if cv[0] == 'bond':
                self.lmp.command("group " + "b" + str(self.colvars.index(cv)) + " id " + " ".join(cv[1:]))
                self.lmp.command("compute " + str(self.colvars.index(cv)) + " b" + str(self.colvars.index(cv)) + " bond/local dist" )
            elif cv[0] == 'angle':
                self.lmp.command("group " + "a" + str(self.colvars.index(cv)) + " id " + " ".join(cv[1:]))
                self.lmp.command("compute " + str(self.colvars.index(cv)) + " a" + str(self.colvars.index(cv)) + " angle/local theta")
            elif cv[0] == 'dihedral':
                self.lmp.command("group " + "d" + str(self.colvars.index(cv)) + " id " + " ".join(cv[1:]))
                self.lmp.command("compute " + str(self.colvars.index(cv)) + " d" + str(self.colvars.index(cv)) + " dihedral/local phi")        
        
        return 0

        
    def destroyColvars(self):
        """
        This function removes the colvars set by initColVars().
        """        
        for cv in self.colvars:
            self.lmp.command("uncompute " + str(self.colvars.index(cv)))
            
        return 0
        
    def equilibrate(self, center):
        """
        This function prepares a LAMMPS image to be at the specified target 
        position given by umb.
        """
        
        # first enter the restraints based on computes data structure
        restCommand = "fix REST all restrain "
        for index, cv in enumerate(self.colvars):
            if cv[0] == 'dihedral':
                restCommand += " " + " ".join(cv) + " " + str(center[index] + 180.0)
            else:
                restCommand += " " + " ".join(cv) + " " + str(center[index])
        
        # now issue restraint definition to the lammps object 
        self.lmp.command(restCommand)
                    
        # now run the equilibration lines from the equilibrate file.     
        for line in self.equilCmds:
            if line[0] == '#':
                continue
            else:
                self.lmp.command(line)

        # apply SHAKE if used
        if self.shakeH:    
            self.lmp.command("fix 10 all shake 0.0001 500 0 m 1.008")
            
        # now remove constraints for subsequent dynamics         
        self.lmp.command("unfix REST")
        
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
    
    def drawVel(self, distType = 'gaussian'):
        """
        This function redraws the velocities from a maxwell-boltzmann dist.
        DOES THIS WORK FOR SIMULATIONS THAT DO NOT HAVE A VELOCITY?
        """
        if distType == 'gaussian':
            self.lmp.command("velocity all create " + str(self.temperature) + " " + str(random.randint(100000,999999)) + " dist gaussian")
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
        
    