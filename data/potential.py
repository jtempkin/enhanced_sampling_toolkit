# -*- coding: utf-8 -*-
"""
Created on Fri May 23 10:22:09 2014

@author: erikthiede
"""
import scipy as sp


class potential:
    """
    This is a template file, which defiles the minimal methods that any potential object should implement.
    They should return the potential at a coordinate (evaluate), and calculate the force there (force)
    """
    def evaluate(self, coord):
        """
        This function should return the potential of the object at a specific coordinate.  
        As a placeholder, it returns 0.0: the value of a free particle.
        It takes as input coord, 
        """
        return 0.0
    def force(self,coord):
        """
        Returns the force on the particle from the potential. 
        As a default, it returns the value of a free particle in one dimension, 0.0
        """
        return 0.0

class rulepotential(potential):
    """
    This class implements a potential which is a function of the coordinate based on some coordinate.
    For instance, these could be defined using the lambda feature.  An example might be a particle in 
    a harmonic well, or on a Muller Brown surface.
    """
    def __init__(self, uRule,fRule):
        """
        The initialize statement takes in a rule for the potential, and for the force.
        """
        self.__u__=uRule
        self.__f__=fRule
    def evaluate(self, coord):
        """
        Returns the potential at a location, using the rule created at initialization.
        """
        return self.__u__(coord)
    def force(self,coord):
        """
        Returns the force on the particle at a coordinate, using the rule created at initialization.
        """
        return self.__f__(coord)

class MullerBrown(rulepotential):
    """
    This class is an implementation of the ruled potential for the Muller Brown surface.
    The Muller Brown surface is a sum of 4 gaussians that try to model the PMF surfaces of the kind
    that are often seen in chemical reactions.
    """
    A=sp.array([-200.0,-100.0,-170.0,15])
    a=sp.array([-1.,-1.,-6.5,0.7])
    b=sp.array([0,0,11,0.6])
    c=sp.array([-10.,-10,-6.5,0.7])
    x_0=sp.array([1.0,0.0,-0.5,-1])
    y_0=sp.array([0.0,0.5,1.5,1.0])
    
    def __init__(self):
        # We overwrite the constructor for rulepotential, since we already know the rule.
        return
    
    def __u__(self,coord):
        """
        Returns the potential energy of the Muller Brown surface.  
        Here, we have assumed that coord is a 2D iterable (typically, a scipy or numpy array object).
        """
        V=0.0
        # We rename the first and second coordinate to x and y, for readability.
        x=coord[0]
        y=coord[1]
        # The Muller Brown potential is calculated using a sum of 4 Gaussians.
        for i in xrange(4):
            V+=self.A[i]*sp.exp(self.a[i]*(x-self.x_0[i])**2+self.b[i]*(x-self.x_0[i])*(y-self.y_0[i])+self.c[i]*(y-self.y_0[i])**2)
        return V 
        
        
    def __f__(self,coord):
        """
        Returns the force at a location, from the Muller Brown Potential as a scipy array.
        We have assumed that coord is a 2D iterable (typically, a scipy or numpy array object).
        """
        fx=0.0
        fy=0.0
        # We rename the first and second coordinate to x and y, for readability.
        x=coord[0]
        y=coord[1]
        # The Muller Brown potential is defined as a sum of 4 Gaussians.
        # To calculate the force, we add up the negative partial deratives in both x and y.
        # This becomes the x and y components of our force, respectively.
        for i in xrange(4):
            factor=-1.0*self.A[i]* sp.exp(self.a[i]*(x-self.x_0[i])**2+self.b[i]*(x-self.x_0[i])*(y-self.y_0[i])+self.c[i]*(y-self.y_0[i])**2)
            fx+=(2.0*self.a[i]*(x-self.x_0[i])+self.b[i]*(y-self.y_0[i]))*factor
            fy+=(2.0*self.c[i]*(y-self.y_0[i])+self.b[i]*(y-self.y_0[i]))*factor
        return sp.array([fx,fy])
        