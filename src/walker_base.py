# -*- coding: utf-8 -*-
"""
This section describes the base API used to define the walker objects. The base class is declared internally using Python's abstract base class module. The implementations of the bindings specific to the dynamics packages are registered to this base class definition through direct subclassing of the modules. A key feature of this decision is that it enforces the implementation of the dynamics bindings to implement each of the following methods in their class declarations. However, many of the methods can be overridden with empty methods if an implementation is incomplete or does not require all of the base class methods. 

DEVELOPER NOTE: we may want to use abstract properties in the future. 
"""
import abc

class walker(object):
    """
    Defines an abstract walker object.  It functions as an interface, providing methods that other
    walkers should implement.
    """
    __metaclass__=abc.ABCMeta
        
    @abc.abstractmethod
    def destroyColvars(self):
        """
        Removes colvars from the walker. 
        """
        return 
      
 
    @abc.abstractmethod
    def addColvars(self, cv):
        """
        Updates internal walker data structure containing a collective variable.
        """
        return

    @abc.abstractmethod
    def getConfig(self):
        """
        This function should return the configuration of the walker in configurtion space.
        """
        return 
        
    @abc.abstractmethod
    def setConfig(self, configuration):
        """
        This function should set the system at a specific place in configurtin space.
        At minimum, it should take some sort of specification of the configuration.
        """
        return
    
    @abc.abstractmethod
    def equilibrate(self, center, restraint, numSteps):
        """
        This function sets the walker at a specific configuration in the collective variable space.
        At minimum, any implementation will need to give it the configuration.
        """
        return
            
    @abc.abstractmethod
    def getColvars(self):
        """
        This function returns the location of the walker in the collective variable space.
        """
        return
        
    @abc.abstractmethod
    def propagate(self,numsteps):
        """
        This function propagates the simulation forward a given number of steps.
        It takes as an argument at least numsteps, the number of steps to progagate forward.
        """
        return
    
    @abc.abstractmethod
    def close(self):
        """
        Destroys the walker.
        """
        return

    @abc.abstractmethod
    def drawVel(self):
        """
        Draws a new value of the velocities for the walker.
        """
        return
    
    @abc.abstractmethod        
    def reverseVel(self): 
        """
        This function reverses the velocity of the walker.
        """
        return
        
    @abc.abstractmethod
    def getVel(self): 
        """
        This function returns the velocities of the system.
        """
        return
        
    @abc.abstractmethod
    def setVel(self):
        """
        This routine sets the velocities of the system.
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
