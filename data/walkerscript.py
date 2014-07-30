# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 16:16:02 2014

@author: erikthiede
"""

import potential
import pointMass

kT = 40
dmax = 0.1

def getWalker():
    return pointMass.met_walker(getPotential(),kT,dmax)
    
def getPotential():
    return potential.MullerBrown()