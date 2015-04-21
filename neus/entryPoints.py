# -*- coding: utf-8 -*-
"""
This module acts as a wrapper for the entry point object.
"""

class entryPoints(object):
    """
    This class acts as a wrapper for an entry point object. Should have at least the structures to contain a phase space point.
    """
    def __init__(self, config, vel, time):
        """
        Constructs an object containing at least a point in phase space. Can be initialized to contain other information like a reference point as well but we leave this to be handled dynamically. 
        """
        self.config = config
        self.vel = vel
        self.time = time 
        