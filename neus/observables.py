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
    def __init__(self, name, data, data_width):
        """
        Constructor for the time correlation function for the end to end distance.
        """
        self.name = name
        self.data = data
        self.data_width = data_width
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

            elif cv.type == "bond":
                temp_sample[i] -= self.data_width[i][0]
            
                # get the index to accumulate 
                indx[i] = self.data.shape[i] - 1 - int(np.floor(temp_sample[i] / ((self.data_width[i][1] - self.data_width[i][0]) / self.data.shape[i] )))
            
            else: 
                print "WARNING: accumulatePMF() does not support given collective variable."

        self.data[tuple(indx)] += 1.0
        self.nsamples[tuple(indx)] += 1.0

        return 0

class P1:
    """
    This routine returns the TCF of the end to end distance.
    """
    def __init__(self, name, s, stepLength, atomids, data, cellDim):
        """
        Constructor for the time correlation function for the end to end distance.
        
        NOTE: Step Length is the time between samples taken, not necessarily the time the walker advances each timestep. 
        """
        self.s = s
        self.stepLength = stepLength
        self.atomids = atomids
        self.name = name
        self.data = data
        self.nsamples = np.zeros(data.shape)
        self.cellDim = cellDim

    def __call__(self, wlkr):
        """
        This function takes the current position of the walker and updates the
        local autocorrelation function estimation.
        """
        # check to see that we are accumulating a value at a time when it is appropriate. return otherwise 
        """
        if wlkr.simulationTime % self.stepLength:     
            return 0 
        """ 
        if not (wlkr.simulationTime - np.floor(wlkr.simulationTime / self.s) * self.s ) % (self.s / self.data.shape[0]) == 0.0:
            return 0 
        
        time_indx = (wlkr.simulationTime - np.floor(wlkr.simulationTime / self.s) * self.s ) / (self.s / self.data.shape[0])
        
        time_indx = int(time_indx)
        
        # current configuration
        config = wlkr.getConfig()
        config = np.reshape(config, (-1,3))
        p1 = config[self.atomids[0]-1]
        p2 = config[self.atomids[1]-1]
        l1 = p2 - p1
        
        # apply minimum image to l1 
        for dim in range(3):
            if l1[dim] > self.cellDim[dim] / 2.0: 
                l1[dim] -= self.cellDim[dim]
            elif l1[dim] < -self.cellDim[dim] / 2.0: 
                l1[dim] += self.cellDim[dim]

        # reference configuration
        config = wlkr.Y_s[0]
        config = np.reshape(config, (-1,3))
        y1 = config[self.atomids[0]-1]
        y2 = config[self.atomids[1]-1]
        l2 = y2 - y1
        
        # apply minimum image to l2
        for dim in range(3):
            if l2[dim] > self.cellDim[dim] / 2.0: 
                l2[dim] -= self.cellDim[dim]
            elif l2[dim] < -self.cellDim[dim] / 2.0: 
                l2[dim] += self.cellDim[dim]

        temp_val = np.dot(l1, l2)

        self.data[time_indx] = (self.data[time_indx] * self.nsamples[time_indx] + temp_val) / (self.nsamples[time_indx] + 1.0)
        self.nsamples[time_indx] += 1.0

        return 0
        

class dist_fluctuation_correlation:
    """
    This routine returns the TCF of the end to end distance.
    """
    def __init__(self, name, s, stepLength, atomids, data, cellDim, mean):
        """
        Constructor for the time correlation function for the end to end distance.
        
        NOTE: Step Length is the time between samples taken, not necessarily the time the walker advances each timestep. 
        """
        self.s = s
        self.stepLength = stepLength
        self.atomids = atomids
        self.name = name
        self.data = data
        self.nsamples = np.zeros(data.shape)
        self.cellDim = cellDim
        self.mean = mean

    def __call__(self, wlkr):
        """
        This function takes the current position of the walker and updates the
        local autocorrelation information using a fluctuation correlation function
        i.e. < (A(0) - <A>) * (A(t) - <A>) > 
        
        """
        # check to see that we are accumulating a value at a time when it is appropriate. return otherwise 
        """
        if wlkr.simulationTime % self.stepLength:     
            return 0 
        """
        if not (wlkr.simulationTime - np.floor(wlkr.simulationTime / self.s) * self.s ) % (self.s / self.data.shape[0]) == 0.0:
            return 0 
        
        time_indx = (wlkr.simulationTime - np.floor(wlkr.simulationTime / self.s) * self.s ) / (self.s / self.data.shape[0]) - 1
        
        time_indx = int(time_indx)
        
        # current configuration
        config = wlkr.getConfig()
        config = np.reshape(config, (-1,3))
        # get the atomic coordinates 
        p1 = config[self.atomids[0]-1]
        p2 = config[self.atomids[1]-1]
        # compute the norm of the displacement
        l1 = p2 - p1
        
        # apply minimum image to l1 
        for dim in range(3):
            if l1[dim] > self.cellDim[dim] / 2.0: 
                l1[dim] -= self.cellDim[dim]
            elif l1[dim] < -self.cellDim[dim] / 2.0: 
                l1[dim] += self.cellDim[dim]

        # reference configuration
        config = wlkr.Y_s[0]
        config = np.reshape(config, (-1,3))
        y1 = config[self.atomids[0]-1]
        y2 = config[self.atomids[1]-1]
        l2 = y2 - y1
        
        # apply minimum image to l2
        for dim in range(3):
            if l2[dim] > self.cellDim[dim] / 2.0: 
                l2[dim] -= self.cellDim[dim]
            elif l2[dim] < -self.cellDim[dim] / 2.0: 
                l2[dim] += self.cellDim[dim]

        # now get the norms of the displacements
        norm1 = np.linalg.norm(l1)
        norm2 = np.linalg.norm(l2)
        
        # now get the fluctuation correlation value 
        temp_val = np.dot(norm1 - self.mean, norm2 - self.mean)

        self.data[time_indx] = (self.data[time_indx] * self.nsamples[time_indx] + temp_val) / (self.nsamples[time_indx] + 1.0)
        self.nsamples[time_indx] += 1.0

        return 0
