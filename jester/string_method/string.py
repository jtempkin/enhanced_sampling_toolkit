"""MS STRING
-------------------------------------
This is an implementation of the string object class used for the string method.

Some things to improve:

* need to review / refactor the structures of the classes here to ensure consistency and good python object behaviors.
* I should remove all extraneous methods from the object definition such that the objects are far simpler in construction.
* I should design how the objects should behave as Python objects.
* should write unit testing for the application package.

Jeremy Tempkin
11/26/14
"""

from scipy import interpolate
import math
import numpy as np
import lammpsWalker

class string:
    """
    A string object houses the properties and methods needed for executing a string method calculation.
    """

    def __init__(self, nimages, filename = None):
        """
        Initialize a string object.

        Input: 
        """
        self.nimages = nimages


    def __call__(self):
        """

        """

    def setImages(self, filename):
        """
        Here we initialize an array of walker objects for the execution of the
        string object.
        """
        # start an array of walkers to store the string configuration.
        self.images = []

        for i in range(self.nimages):
            wlkr = lammpsWalker.lammpsWalker(filename, index=i)
            self.images.append(wlkr)

        return 0


    def setColvars(self, cvs):
        """
        This function sets the colvars data structure needed for the string to operate. Also sets up a numpy array for storing the state of the string.
        """

        self.ncvs = len(cvs)

        self.state = np.zeros((self.nimages, self.ncvs))

        return 0


    def getState(self):
        """
        Returns the current state of the string in collective variable space.
        """

        return self.state

    def reparameterize(self):

        return 0

    def smoothString(self, kappa):
        """
        This string smooths the images of the string through averaging adjacent images.
        """
        nimages = len(colvars)
        ncvs = len(colvars[0])

        #print "colvars :"
        #for i in range(0, len(colvars[0]), 1):
            #print colvars[1][i][-1]

        for image in range(1,nimages-1,1):
            for cv in range(0, ncvs, 1):
                colvars[image][cv][-1] = colvars[image][cv][-1] + kappa * (colvars[image-1][cv][-1] + colvars[image+1][cv][-1] - 2.0 * colvars[image][cv][-1])

        #print "smooth colvars :"
        #for i in range(0, len(colvars[0]), 1):
            #print colvars[1][i][-1]

        return colvars
