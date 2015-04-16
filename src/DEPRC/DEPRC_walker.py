# -*- coding: utf-8 -*-
"""
NOTE: we may want to use abstract properties in the future. 
"""
import abc

class walker(object):
    """
    Defines an abstract walker object.  It functions as an interface, providing methods that other
    walkers should implement.
    """
    __metaclass__=abc.ABCMeta
    
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
    def equilibrate(self, colvar):
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
        
class velocityWalker(walker):
    """
    This is an abstract class for dynamics which have a velocity component, such as 
    a protein under Langevin dynamics.  It extends the walker abstract class by
    adding additional abstract methods that are necessary for a walker that has
    a velocity.
    
    """
    
    __metaclass__=abc.ABCMeta

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
