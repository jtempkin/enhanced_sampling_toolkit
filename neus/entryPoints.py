# -*- coding: utf-8 -*-
"""
This module acts as a wrapper for the entry point object.
"""

class entryPoints(object):
    """
    This class acts as a wrapper for an entry point object. Should have at least the structures to contain a phase space point.
    """
    def __init__(self, config, vel, time, Y_s = None):
        """
        Constructs an object containing at least a point in phase space. Can be initialized to contain other information like a reference point as well but we leave this to be handled dynamically. 
        
        We should note that the Y_s is a tuple of (X(s), V(s), s_index) where s_index specifies to which window this entry point is associated. 
        """
        self.config = config
        self.vel = vel
        self.time = time 
        self.Y_s = Y_s
        