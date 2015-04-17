# -*- coding: utf-8 -*-
"""
This module contains definitions of the dynamics classes using to wrap information about driving dynamics inside the walkers. 
"""

class langevin(): 
    """
    This class specifies the langevin dynamcis parameters needed. 
    """
    
    def __init__(self, temperature, damping_coefficient, shake=False, seed=None):
        """
        This constructs a langevin object. Requires that you specify a temperature and a damping coefficient. 
        """
        self.type = "langevin"
        self.damping_coefficient = damping_coefficient
        self.temperature = temperature
        self.shake = shake
        self.seed = seed