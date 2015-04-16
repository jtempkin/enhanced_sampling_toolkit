# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 15:03:56 2015

@author: jeremytempkin
"""

class outputClass(object):
    """
    This class acts a container object for an output for a walker object.
    """
    def __init__(self, name, outputType, filename, nSteps=1000): 
        """
        Here we construct the base necessary fields. Note that nSteps has a default value. 
        """        
        self.name = name 
        self.outputType = outputType
        self.destination = filename
        self.nSteps = nSteps 
        