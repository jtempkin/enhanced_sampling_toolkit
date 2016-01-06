# -*- coding: utf-8 -*-
"""
DEVELOPER NOTE: we may want to use abstract properties in the future. 
"""
import abc

class walker(object):
    """
    This object acts as an abstract base class that sets the core specification defining the Walker API. These routines 
    """
    __metaclass__=abc.ABCMeta
        
    @abc.abstractmethod
    def destroyColvars(self):
        """
        Removes all collective variables from the walker. 
        """
        return 
 
    @abc.abstractmethod
    def addColvars(self, cv):
        """
        Adds a new collective variable to the walker. Subsequent calls to routines that act on the defined collective variables will use the complete list of collective variables in the order they are added. To wipe this list of variables, see the destroyColvars() routine. 
        """
        return

    @abc.abstractmethod
    def getConfig(self):
        """
        This function should return the postition of the walker in configuration space as a one-dimensional numpy array. The coordinates are returned as [x1, y1, z1, x2, y2, z2, ...]. 
        """
        return 
        
    @abc.abstractmethod
    def setConfig(self, configuration):
        """
        Sets the configuration of the At minimum, it should take some sort of specification of the configuration.
        """
        return

    @abc.abstractmethod
    def getVel(self): 
        """
        This function returns the velocities of the system as a one-dimensional numpy array. The velocities are returned in a format that matches the getConfig() routine, specifically [v_x1, v_y1, v_z1, v_x2, v_y2, v_z2, ...]. 
        """
        return

    @abc.abstractmethod
    def setVel(self):
        """
        This routine sets the velocities of the system.
        """
        return

    @abc.abstractmethod
    def drawVel(self):
        """
        Draws a new value of the velocities for the walker. Can redraw velocities according to a uniform or Gaussian distributions. 
        """
        return
    
    @abc.abstractmethod        
    def reverseVel(self): 
        """
        This function reverses the velocity of the walker. The new velocities are given as 
        """
        return
    
    @abc.abstractmethod
    def equilibrate(self, center, restraint, numSteps):
        """
        """
        return
            
    @abc.abstractmethod
    def getColvars(self):
        """
        This function returns the location of the walker in the collective variable space. Values are returned as a one-dimensional numpy array 
        """
        return
        
    @abc.abstractmethod
    def propagate(self,numsteps):
        """
        Integrates the dynamics of the model forward in time. Takes the numSteps argument specifying the number of time steps to take. 
        """
        return
    
    @abc.abstractmethod
    def close(self):
        """
        Destroys the walker. 
        """
        return
        
    @abc.abstractmethod
    def setOutput(self):
        """
        This routine adds a source of output for the walker to write information to disk.
        """
        return 
        
    @abc.abstractmethod
    def removeOutput(self):
        """
        This routine adds a source of output for the walker to write information to disk.
        """
        return 
