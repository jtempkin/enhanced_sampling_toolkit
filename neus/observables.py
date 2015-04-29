# -*- coding: utf-8 -*-
"""
This module contains a list of observable routines. These functions serve to calculation and record the value of observables in the course of a calculation.

These routines are initialized as their own class and subsequently can be invoked during other calculations or called by other routines. We may consider updating these classes with a base class prototype similar to how the walker class was initialized. 
"""

import numpy as np
import copy
    
# ------------- OBSERVABLES FUNCTIONS ---------------------
class pmf:
    """
    This class represents a PMF observable.
    """
    def __init__(self, name, data):
        """
        Constructor for the time correlation function for the end to end distance.
        """
        self.name = name
        self.data = data
        self.nsamples = np.zeros(data.shape)

    def __call__(self, sample, colvars):

        temp_sample = copy.deepcopy(sample)
        # make an array of
        indx = np.zeros(len(temp_sample), dtype = np.int8)
        
        assert len(temp_sample) == self.data.ndim, "The sample received does not match the dimensionality of the data array for this observable."
        
        assert len(temp_sample) == len(colvars), "The sample dimentionality does not match the colvars dimensions"

        # we should add a check here to make sure that we are going to appropriately match the data type handed to us
        #assert

        for i,cv in enumerate(colvars):
            if cv.type == 'dihedral':
                # shift dihedral values up to ranges between [0.0, 360.0]
                temp_sample[i] += 180.0

                indx[i] = self.data.shape[i] - 1 - int(np.floor(temp_sample[i] / (360.0 / self.data.shape[i] ) ) )

            else:
                print "WARNING: accumulatePMF() does not support given collective variable."

        self.data[tuple(indx)] += 1.0
        self.nsamples[tuple(indx)] += 1.0

        return 0

class P1:
    """
    This routine returns the TCF of the end to end distance.
    """
    def __init__(self, name, s, stepLength, atomids, data):
        """
        Constructor for the time correlation function for the end to end distance.
        """
        self.s = s
        self.stepLength = stepLength
        self.atomids = atomids
        self.name = name
        self.data = data
        self.nsamples = np.zeros(data.shape)

    def __call__(self, wlkr):
        """
        This function takes the current position of the walker and updates the
        local autocorrelation function estimation.
        """
        time_indx = (wlkr.simulationTime - np.floor(wlkr.simulationTime / self.s) * self.s ) / self.stepLength

        # current configuration
        config = wlkr.getConfig()
        p1 = config[self.atomids[0]:self.atomids[0]+3]
        p2 = config[self.atomids[1]:self.atomids[1]+3]
        l1 = p1 - p2

        # reference configuration
        config = wlkr.Y_s[0]
        y1 = config[self.atomids[0]:self.atomids[0]+3]
        y2 = config[self.atomids[1]:self.atomids[1]+3]
        l2 = y1 - y2

        temp_val = np.sqrt(np.dot(l1, l2))

        self.data[time_indx] = (self.data[time_indx] * self.nsamples[time_indx] + temp_val) / (self.nsamples[time_indx] + 1.0)
        self.nsamples[time_indx] += 1.0

        return 0
