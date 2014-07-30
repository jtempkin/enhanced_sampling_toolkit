# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:26:00 2014

@author: jtempkin
"""

class DynamicsError(Exception):
    """
    A dynamics error exception class. 
    """
    def __init__(self, value):
        """
        init for dynamics error. 
        """
        print "Dynamics Error in umbrella", value, "occured."
        
def hangupHandler(signum, stack):
    """
    This function is designed to handle signal 1 errors sent from a subprogram.
    """
    print "Received: ", signum
    raise DynamicsError("An error occured in the Dynamics routine of ")