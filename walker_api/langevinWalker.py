# -*- coding: utf-8 -*-
"""
This module implements the Walker API for the LAMMPS MD engine. See walker_base.py for a specification of the API. For details concerning the usage of the LAMMPS MD package see the excellent documentation at the LAMMPS webpage:

http://lammps.sandia.gov/


In particular, you may want to see how the Python wrapper to LAMMPS on which this implementation is based:

http://lammps.sandia.gov/doc/Section_python.html

Here we will outline basic usage guides for the walker API usage in LAMMPS.
"""

import random, copy
from walker_base import walker
import collectiveVariables
import outputClass
import numpy as np


class langevinWalker(walker):
    """
    """

    def __init__(self, logFilename, inputFilename=None, index = 0, debug = False):
        """
        Initializes a walker object.

        Arguments:
        -----------------------
        inputFilename -
            A string with a system path to an input file. This file is
            ignored now.

        LogFilename -
            A string to pipe standard LAMMPS output into.

        index (default: 0) -
            A optional index which is appended to the log files and output
            files to protect overwriting multiple copies

        debug (default: False) -
            A boolean flag for debugging purposes. If set to True, will direct
            LAMMPS output to stdout instead of to logFilename.

        """

        # a filename for storing the data, note the appending of the index
        self.filename = inputFilename + "." + str(index)

        self.logFilename = logFilename

        if inputFilename is not None:
            from zuckerman_potential import U
            self.U = U

        # walker index
        self.index = index

        # a list of the relevant collective variables
        self.colvars = []

        # we're constructing a list for collective output classes.
        self.output = []

        self.dynamics = None

        # a list of commands used to equilibrate the system
        self.equilCmds = []

        # set an initial array for coordinates and velocities
        self.coordinates = np.array([0.0])
        self.velocity = np.array([0.0])



    def close(self):
        """
        This function closes the object. Right now this does nothing.
        """

        return 0


    def addColvars(self, name, cvType, atomIDs):
        """
        Implements the addition of a collective variable to the list of collective variables held by this walker.

        Arguments:
        -----------------------
        name -
            A internal string used to reference this collective variable. Collective variable names must be unique.

        cvType -
            A string refering to the type of collective variable. Currently the following variables are supported:
                * 'bond'
                * 'angle'
                * 'dihedral'
                * 'x', 'y', 'z' positions
                * 'x', 'y', 'z' velocity components

        atomIDs -
            A list of the atom indecies involved in the collective variable. Should provide the right number of atom indecies for the collective variable:
                * 'bond' -> 2
                * 'angle'-> 3
                * 'dihedral' -> 4
                * position or velocity component -> 1

        """
        # make sure we know what the cv type is.
        __knownCVs__ = ['x', 'y', 'z', 'vx', 'vy', 'vz']
        assert cvType in __knownCVs__, "cvType that was provided was not understood."

        # check to make sure the name is unique.
        for cv in self.colvars:
            assert cv.name != name, "Collective variable names must be unique."

        # now append the collective variable to the walker list. Initialize a collective variable object.
        self.colvars.append(collectiveVariables.collectiveVariables(name, cvType, atomIDs))

        return 0

    def destroyColvars(self):
        """
        This function removes the colvars set by setColvars(). By default, it removes all of the collective variables in the list. It does not remove them from the collective variables list.
        """

        self.colvars = []

        return 0

    def equilibrate(self, center, restraint, numSteps):
        """
        This function prepares a LAMMPS image to be at the specified target position given by the vector 'center' passed and an arguments.

        Arguments:
        -----------------------

        """
        #print "Equilibrating walker."

        self.setConfig(center)

        return 0

    def getConfig(self):
        """
        This function returns the current position of the simulation.
        """
        config = copy.deepcopy(self.coordinates)

        return config

    def getVel(self):
        """
        This function returns the current velocities from the LAMMPS simulation.
        """
        vel = copy.deepcopy(self.velocity)

        return vel

    def setVel(self, vel):
        """
        This function sets the velocity to the lammps simulation.
        """

        if self.velocity.shape != vel.shape: raise Exception('The velocity provided is not the correct dimension.')

        self.velocity = copy.deepcopy(vel)

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
            cvarray.append(self.coordinates)

        assert len(cvarray) == len(self.colvars), "Not all collective variables were added."

        return copy.deepcopy(cvarray)

    def setConfig(self, config):
        """
        This routine sets the internal configuration.
        """

        if self.coordinates.shape != config.shape: raise Exception('The coordinates array provided is not the correct dimension.')

        self.coordinates = copy.deepcopy(config)

        return 0

    def drawVel(self, distType = 'gaussian', temperature = 310.0, seed = None):
        """
        This function redraws the velocities from a maxwell-boltzmann dist.
        """

        return 0

    def reverseVel(self):
        """
        This function reverses the velocities of a given LAMMPS simulation
        """

        return 0

    def propagate(self, numSteps, pre='yes', post='no'):
        """
        This function issues a run command to the underlying dynamics to propagate
        the dynamics a given number of steps.
        """

        for i in range(numSteps):
            # generate proposal dx in [-pi/2, pi/2]
            dx = (np.pi/2.0) * (random.random()*2 - 1)

            rand_val = random.random()

            if rand_val < min(1.0, np.exp(-self.U(self.coordinates+dx))/ np.exp(-self.U(self.coordinates))):
                self.coordinates += dx


        return 0

    def setDynamics(self, dynamics_instance):
        """
        This routine sets the dynamics for the walker.
        """

        return 0

    def minimize(self, args=None):
        """
        This function runs a minimization routine with the specified type.
        """

        return 0

    def setTemperature(self, temp):
        """
        This function sets the temperature of the walker object.
        """

        return 0

    def setTimestep(self, timestep):
        """
        This routine sets the dynamics time step.
        """

        return 0

    def setOutput(self, name, outputType, filename, nSteps):
        """
        This routine sets up a mechanism for writing system information to file directly from the dynamics engine. This is equivalent to output constructed
        """
        # assert outputType in __knownOutput__, "The output type specified was not recognized."
        __knownOutputTypes__ = ['text', 'h5py', 'binary']
        assert outputType in __knownOutputTypes__, "outputType that was provided was not understood."

        for out in self.output:
            assert out.name != name, "Output class names must be unique."

        # now append the collective variable to the walker list.
        self.output.append(outputClass.outputClass(name, outputType, filename, nSteps))


    def removeOutput(self):
        """
        This routine removes the ouput pipes for information to be written to disk from the underlying dynamics engine. Right now this simply clears all existing output.
        """

        self.output = []

        return 0


walker.register(langevinWalker)

if __name__ == "__main__":
    print 'The lanegevinWalker module is a subclass of walker:', issubclass(langevinWalker, walker)

