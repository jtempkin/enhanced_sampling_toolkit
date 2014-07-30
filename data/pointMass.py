# -*- coding: utf-8 -*-
"""
Created on Thu May 22 17:24:19 2014

@author: erikthiede
"""
import numpy as np
import copy
import walker


class met_walker(walker.walker):
    """
    This gives a walker which propagates according to Metropolis Monte Carlo.
    Please note that the collective variable here is assumed to be synonomous with the spatial coordinate.
    """
    
    
    def __init__(self, U, kT=1, dmax=1):
        """
        This is the initialization method for the Metropolis Monte Carlo walker.
        Coords should be a numpy array, U should be a potential object.  kT and dx should both be scalars.
        """
        self.potential = U
        self.kT = kT
        self.config = None
        self.dmax = dmax
        self.currentweight = None 
        
    def getColvars(self):
        return np.copy(self.config)
        
    def equilibrate(self, colvar):
        """
        This sets the configuration of the walker to a configuration which fits
        under a certain collective variable.
        """
        self.config=np.copy(np.asarray(colvar))
        self.currentweight = np.exp(-1.0*self.potential.evaluate(self.config)/self.kT)
    
    def getConfig(self):
        return np.copy(self.config)
    
    def setConfig(self,coords):
        """
        We set the config to the old config.
        """
        self.config = np.copy(np.asarray(map(float,coords)))
        self.currentweight = np.exp(-1.0*self.potential.evaluate(self.config)/self.kT)
    
    def propagate(self,n):
        """
        Propagates the dynamics forward using the leapfrog algorithm and Langevin dynamics.
        """
        # We create a list of random gaussian numbers to draw from.
        for i in xrange(n):
            #  First, we calculate a random displacement
            step=(2.0*np.random.random_sample(self.config.shape)-1.0)*self.dmax
            # We calculate U at the proposed point.
            newweight=np.exp(-1.0*self.potential.evaluate(step+self.config)/self.kT)
            # We accept the move if we meet the metropolis test:  otherwise, we just stay at the same location.
            if(np.random.rand()<newweight/self.currentweight):
                self.config+=step
                self.currentweight=newweight
        return 0
        
        
    def close(self):
        
        return
        
class lan_walker(walker.velocityWalker):
    """
    Creates a walker that propagates forward according to Langevin dynamics on a surface.
    Please note that the collective variable here is assumed to be synonomous with the spatial coordinate.
    """
    
    def __init__(self, U,dt=1.0,m=1.0,gamma=1.0,kT=1.0):
        """
        The initialization method for a Langevin equation in Cartesian coords.
        Here, we expect coords to be an array or list, and U to be a potential object.
        m,v0,D, and kT should be numbers.  Coords and v0 should be passed as 
        num/scipy arrays.  Note that vLast here is the speed at V-dt/2
        
        In the absence of any direction, it defaults to natural units for dt, m, 
        gamma, and kT
        """

        self.config = None
        # We create the plethora of constants necessary for Langevin dynamics.
        self.m=m
        self.vlast = None
        self.potential = U
        self.gamma = gamma    
        self.kT = kT
        self.dt = dt

    def getColVar(self):
        return np.copy(self.config)
        
    def equilibrate(self, colvar):
        """
        This sets the configuration of the walker to a configuration which fits
        under a certain collective variable.
        """
        self.config=np.copy(np.asarray(colvar))
        self.drawVel()
        
    def drawVel(self):
        self.vlast=np.array([np.random.normal(0.0,np.sqrt(self.kT/self.m)) for i in self.vlast])
        
    def reverseVel(self):
        """
        Reverses the velocity the walker, by multiplying vlast by -1.  
        Note that used improperly, this could violate conservation of energy.
        The proper way to bounce this is:  save a configuration, and take a step.
        Set the configuration to the old configuration, and then reverse the velocity.
        """
        self.vlast*=-1
        
    def setConfig(self,coords):
        """
        We set the config to the old config, clear out the old lastconfig, and then randomize the velocity.
        """
        self.config=np.copy(np.asarray(coords))
        self.drawVel()
      
    def getConfig(self):
        return np.copy(self.config)
    def propagate(self,n):
        """
        Propagates the dynamics forward using the leapfrog algorithm and Langevin dynamics.
        """
        # We create a list of random gaussian numbers to draw from.
        for i in xrange(n):
            # We calculate the change in the velocity, dv.
            # We do this by adding the force from the potential, a frictional term (gamma*v/m), and
            # by adding a random number scaled by sqrt(2*gamma*kT/(m*dt)).
            random = np.random.normal(0,1,(1,len(self.config)))
            dv=((np.asarray(self.potential.force(self.config))-self.gamma*self.vlast)+np.sqrt(2*self.gamma*self.kT/(self.m*self.dt))*random)*self.dt/self.m
            # We update the velocity.
            self.vlast+=dv
            # We update the configuration.
            self.config+=self.vlast*self.dt
        return 0
        

        
