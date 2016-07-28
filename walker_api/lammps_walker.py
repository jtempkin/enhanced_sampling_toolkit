# -*- coding: utf-8 -*-
"""
This class implements the enhanced sampling walker API for the bindings to the
LAMMPS package.
"""

import random
from walker_base import walker
import collective_variables
import output_class
import numpy as np
import ctypes
import os


class Lammps(walker):
    """Implements the Walker API using the LAMMPS engine.
    """
    def __init__(self, inputFilename, logFilename, index=0, debug=False):
        """Instantiate a LAMMPS walker object.

        Parameters
        -----------------------
        inputFilename : string
            A string with a system path to a LAMMPS input file. This file is
            read and interpreted by LAMMPS's own interpreter so the syntax
            should be readable by LAMMPS. If None is provided, the input file
            is ignored.

        LogFilename : string
            A string to pipe standard LAMMPS output into.

        index : int (0)
            A optional index which is appended to the log files and output
            files to protect overwriting multiple copies

        debug : bool (False)
            A boolean flag for debugging purposes. If set to True, will direct
            LAMMPS output to stdout instead of to logFilename.

        """
        if inputFilename is not None:
            # check to see if there is actually a filename at that location.
            assert os.path.isfile(inputFilename), "The input file for the LAMMPS walker could not be found."

            # a filename for storing the data, note the appending of the index
            self.filename = inputFilename + "." + str(index)

        # walker index
        self.index = index

        if logFilename is not None:
            logfile = ".".join([logFilename, str(index)])
        else:
            logfile = None


        # start walker object from the passed lammps input file name
        self.lmp = self.__initLAMMPS__(inputFilename = inputFilename, logFilename = logfile, verbose = debug)

        # a list of the relevant collective variables
        self.colvars = []

        # we're constructing a list for collective output classes.
        self.output = []

        self.dynamics = None

        #self.Y_s = (self.get_position(), self.get_velocity())

        self.time = 0.0

        # a list of commands used to equilibrate the system
        self.equilCmds = []


    def close(self):
        """Close the underlying LAMMPS object.
        """
        self.lmp.close()
        return 0

    def __initLAMMPS__(self, inputFilename=None, logFilename=None, verbose = False):
        """Initialize a LAMMPS object using the provided input file.

        Parameters
        -----------------------
        inputFilename : string (None)
            Path to input file.

        logFilename : string (None)
            Path to log file.

        verbose : bool
            If true, write output to standard out.
        """

        try:
            from lammps import lammps
        except:
            print "ERROR: The LAMMPS python interface could not be found. Check your PYTHONPATH and compiled LAMMPS directory to make sure it is compiled and importable."


        # check to see we are handing an absolute pathname.
        if not os.path.isabs(inputFilename):
            inputFilename = os.path.abspath(inputFilename)

        # we're going to assume here that the user is taking care of the fact that LAMMPS needs to know where to find things but we will issue a worning in case they aren't.
        if os.path.dirname(inputFilename) != os.getcwd():
            print "WARNING: Getting lammps input file from directory other than current working directory. Be sure pathnames are correct in input files."

        # initialize the lammps object
        if verbose == True:
            self.lmp=lammps()
            #self.lmp = lammps("", ["-c", "on", "-sf", "cuda"])
        else:
            args = ["-sc", "none", "-echo", "none", "-log", "none"]
            self.lmp = lammps("", args)

        # after the lammps object is created, initialize the lammps simulation.
        # specify general log file parameters
        #self.command("echo none")

        # write out log file if needed
        if logFilename is not None:
            if not os.path.isabs(logFilename):
                logFilename = os.path.abspath(logFilename)

            self.command("log " + logFilename)
        else:
            self.command("log none")


        # if there is a specified filename, use it to set up the simulation.
        if inputFilename != None:
            #self.command("log " + logFilename)

            self.lmp.file(inputFilename)

        #self.__setupLAMMPS__(filename)

        return self.lmp

    def add_colvars(self, name, cvType, atomIDs):
        """Add a collective variable to the walker.

        The following variables are currently supported:

        * 'bond'
        * 'angle'
        * 'dihedral'
        * 'x', 'y', 'z' position components
        * 'x', 'y', 'z' velocity components

        The atomIDs must list the atoms indecies involved and should provide the right number of atom indecies for the collective variable:

        * 'bond' -> 2
        * 'angle'-> 3
        * 'dihedral' -> 4
        * position or velocity component -> 1

        Note that due to feature of LAMMPS's implementation, the bond, angle or dihedral must be know in the topology of the LAMMPS data file. I.e. LAMMPS must be able to build this compute without an error.

        Parameters
        -----------------------
        name : string
            A string used to reference this collective variable. Collective variable names must be unique.

        cvType : string
            A string refering to the type of collective variable.

        atomIDs : list
            A list of the atom indecies involved in the collective variable.
        """
        # make sure we know what the cv type is.
        __knownCVs__ = ['bond', 'angle', 'dihedral', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'fe']
        assert cvType in __knownCVs__, "cvType that was provided was not understood."

        # check to make sure the name is unique.
        for cv in self.colvars:
            assert cv.name != name, "Collective variable names must be unique."

        # now append the collective variable to the walker list. Initialize a collective variable object.
        self.colvars.append(collective_variables.CollectiveVariables(name, cvType, atomIDs))

        # grab a handle to the new object
        cv = self.colvars[-1]

        # first set the group for the colvar
        self.command("group " + cv.name + " id " + " ".join(map(str,cv.atomIDs)))

        # now set the appropriate colvar as a compute to LAMMPS
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

        self.propagate(0, pre='yes')

        self.ref_cv = self.get_colvars()

        return 0

    def destroy_colvars(self):
        """Remove all collective variables from the walker instance.
        """
        for cv in self.colvars:
            self.command("uncompute " + cv.name)

        self.colvars = []

        self.propagate(0, pre='yes')

        return 0

    def equilibrate(self, center, restraint, numSteps):
        """Equilibrate the walker to a value of the collective variabales.

        This equilibrates the system to the target value in the collective variable space. It applies a harmonic restraint at the specified coordinates with a given strength and runs a fixed number of dynamics steps. Note that sequences of the equilibrate commands can be used to implement dragging protocols.

        The length of the center and restraint parameters must match the length of the length of the collective variable array.

        Parameters
        --------------------
        center : list
            A list of the centers of the harmonic restraint in collective variables space.

        restraint : list
            A list of the strength of the harmonic.

        numSteps : int
            Number of dynamics timesteps to integrate.
        """

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

    def get_position(self):
        """Return the spatial coordinates of the walker.

        Returns
        -----------
        numpy.ndarray
            A one-dimensional array of the coordinates in [x1, y1, z1, x2, y2, z2,...] format.
        """
        config = np.asarray(self.lmp.gather_atoms("x",1,3), dtype='float64')

        return config

    def get_velocity(self):
        """Return the velocity of the walker.

        Returns
        -----------
        numpy.ndarray
            A one dimensional array of the velocities in [vx1, vy1, vz1, vx2, vy2, vz2,...] format.
        """
        vel = np.asarray(self.lmp.gather_atoms("v", 1, 3), dtype='float64')

        return vel

    def set_velocity(self, vel):
        """Set the velocity of the walker.

        Parameters
        --------------
        vel : numpy.ndarray
            A one-dimentional array of the velocities in [vx1, vy1, vz1, vx2, vy2, vz2,...] format.
        """

        self.lmp.scatter_atoms("v", 1, 3, vel.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

        self.propagate(0, pre='yes')

        return 0

    def get_colvars(self):
        """Return the values of the collective variables.

        Returns
        -------------
        numpy.ndarray
            A numpy array of the collective variables computed by this walker
        """
        # get an empty array with placeholders
        cvarray = []

        # now get cv's one by one from each compute defined
        for cv in self.colvars:
            if cv.type in ['x', 'y', 'z', 'vx', 'vy', 'vz']:
                cvarray.append(self.lmp.extract_compute(cv.name, 1, 1)[0])
            elif cv.type in ['fe']:
                cvarray.append(self.lmp.extract_compute(cv.name, 0, 1)[0])
            else:
                #*** We REALLY need assurance here that what we are getting here is in fact not a NULL POINTER.
                cvarray.append(self.lmp.extract_compute(cv.name, 2, 1)[0])

        assert len(cvarray) == len(self.colvars), "Not all collective variables were added."

        return np.array(cvarray)

    def set_position(self, config):
        """Set the spatial coordinates of the walker.

        Parameters
        -------------
        config : numpy.ndarray
            A one-dimentional numpy array of the position coordinates in [x1, y1, z1, x2, y2, z2,...] format.
            
        """
        self.lmp.scatter_atoms("x", 1, 3, config.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

        self.propagate(0, pre='yes')

        return 0

    def draw_velocity(self, distType = 'gaussian', temperature = 310.0, seed = None):
        """Redraw the velocity of the walker.

        Current supports the following distributions:

        * 'gaussian'
        """
        if seed is None:
            seed = random.randint(100000,999999)

        if distType == 'gaussian':
            self.command("velocity all create " + str(temperature) + " " + str(seed) + " dist gaussian")
        else:
            print "The draw_velocity() routine was passed a distribution Type that was not understood."

        return 0

    def reverse_velocity(self):
        """Reverse the velocity of the walker.

        This routine reverses the sign of the velocities.
        """

        # set varibales for reversing velocities
        self.command("variable vx atom -vx")
        self.command("variable vy atom -vy")
        self.command("variable vz atom -vz")

        # set new velocities
        self.command("velocity all set v_vx v_vy v_vz")

        self.propagate(0, pre='yes')

        return 0

    def propagate(self, numSteps, pre='no', post='no'):
        """Integrate the dynamics forward in time. 

        This routine integrates the dynamics of the underlying model. It expects an integer number specifying the number of time steps to take.

        Parameters
        -----------
        numSteps : int
            The number of model timesteps to take.

        pre : string ('no')
            Run the pre-dynamics setup in LAMMPS if set to 'yes'.

        post : string ('no')
            Run the post-dynamics timings in LAMMPS if set to 'yes'.
        """

        self.command("run " + str(numSteps) + " pre " + str(pre) + " post " + str(post))

        return 0

    def set_dynamics(self, dynamics_instance):
        """Apply a dynamics model to the walker.

        Parameters
        --------------
        dynamics_instance : dynamics
            Set up the dynamics model from the dynamics instance passed.
        """
        __knownDynamics__ = ['langevin', 'baoab']

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

        elif self.dynamics.type is 'baoab':
            if self.dynamics.seed is None:
                self.command("fix 1 all baoab " + " ".join([str(self.dynamics.damping_coefficient), str(self.dynamics.temperature), str(random.randint(100000, 999999))]))
            else:
                self.command("fix 1 all baoab " + " ".join([str(self.dynamics.damping_coefficient), str(self.dynamics.temperature), str(self.dynamics.seed)]))

        # set shake
        if self.dynamics.shake is True:
            self.command("fix shk all shake 0.0001 20 0 m 1.0")
            self.dynamics.fixes.append("shk")

        if self.dynamics.linear_momentum == False:
            self.command("fix mom all  momentum 1 linear 1 1 1")
            self.dynamics.fixes.append("mom")

        self.propagate(0, pre='yes')

        return 0

    def command(self, command):
        """Issue a LAMMPS command directly to the LAMMPS code.

        Parameters
        ------------
        command : string
            A command to pass directy to LAMMPS's command interpreter.
        """

        # issue the given command directly.
        self.lmp.command(command)

        return 0

    def minimize(self, args=None):
        """Perform a steepest descent minimization.
        """
        # use default settings
        if args is None:
            self.command("minimize 1.0e-4 1.0e-6 100 1000")
        # else use the arguments provided.
        else:
            self.command("minimize " + " ".join(args))

        return 0

    def set_temperature(self, temp):
        """Set the temperature.

        This command works only if the dynamics was set with the set_dynamics routine. Otherwise, dynamcis is set by input script provided at initialization.

        Parameters
        -------------
        temp : float
            The temperature in LAMMPS units.
        """
        self.dynamics.temperature = temp

        self.set_dynamics(self.dynamics)

        return 0

    def set_timestep(self, timestep):
        """Set the timestep of the integrator for the dynamics model.

        This command sets the timestep used in the model's integration routine.

        Parameters
        -------------
        timestep : float
            The integration timestep.
        """
        self.command("timestep " + str(timestep))

        return 0

    def set_output(self, name, outputType, filename, nSteps):
        """Set an output file for coordinates from the dynamics model.

        Parameters
        ------------
        name : string
            A string identifying this output type.

        outputType : string
            A string identifying the output file type.

        filename : string
            A string with the filename to write dynamics output.

        nSteps : int
            The number of model timesteps bewteen which to write.s
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


    def remove_output(self):
        """Remove all output for the model.
        """
        for out in self.output:
            self.command("undump " + out.name)

        self.output = []

        return 0

    def get_time(self):
        """Return the current time in number of timesteps.

        Returns
        ---------
        time : integer
            Returns the number of model time steps taken.
        """
        return self.time

    def set_time(self, t):
        """Set the time of the walker.

        This routine expects an integer time in the number of model timesteps.

        Parameters
        -----------
        t : integer
            An integer specifying the number of model timesteps taken.
        """

        self.time = t

        self.lmp.command("reset_timestep " + str(t))

        self.propagate(0, pre="yes")

        return None

walker.register(Lammps)

if __name__ == "__main__":
    print 'The lammpsWalker module is a subclass of walker:', issubclass(Lammps, walker)
