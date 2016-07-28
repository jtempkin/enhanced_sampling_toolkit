# -*- coding: utf-8 -*-
"""
A container class for collective variable definitions. We provide here a container class for
"""

class CollectiveVariables(object):
    def __init__(self, name, cvType, atomIDs):
        """
        Here we declare a name, type and list of atomID's as the parameters. These are all required parameters.
        """
        self.name = name
        self.type = cvType
        self.atomIDs = atomIDs
        
