# -*- coding: utf-8 -*-
"""
Created on Mon May  5 18:04:54 2014

@author: jtempkin
"""

"""
This file contains the class definition used to structure the basis functions.
Not surprisingly, the parent class "basisFunction" contains information that 
is pertinent to all basis functions regardless of size or shape.

The child classes then define characteristics and functions necessary to 
interact with the basis function (i.e. define functions that evaluate the 
value of the basis function at a coordinate).

The advantage here is that new basis functions can be added with easy, simply
add a new subclass defining that particular basis function.
"""

import numpy as np
from scipy import linalg as LA
import scipy as sp
import fileIO
import sys
import walker
#import acor
import errors
import random
import h5py

class basisFunction:
    """
    It is important to note that the configurations associated with 
    sampling of this box are stored at this level of the class structure.
    
    At this level, the characteristics defining the basisFunction class is
    a set of configurations inherent to the class. 
    """
    
    def flushDataToFile(self, filename):
        """
        This function flushes the internal data buffers configs and samples to
        files.
        
        This function uses the HDF5 python implementation to store data to disk.
        It should be noted that it creates two groups in the top-level called 
        'colvars' and 'timeSeries'. Because appending data is not as easy for 
        this file format, it stores each flush as a new dataset indexed in each
        group. 
        """
        if len(self.samples) == 0:
            print "Tried to flush an empty samples buffer."
            return 0 
        
        # get the last sample in the list
        lastSample = self.samples.pop()
        #lastBasisFunction = self.basisFnxTimeSeries.pop()
        #lastConfig = self.configs.pop() 
        
        # write out the remaining data structure
        with h5py.File(filename + ".hdf5", "a") as f_handle:
            # check to see if we've created a group for the colvars writes
            if 'colvars' in f_handle:
                # if so, add in a dataset for this flush
                dset = f_handle['colvars'].create_dataset("cv_" + str(len(f_handle['colvars'].keys())), np.asarray(self.samples).shape, dtype="f")
                # now write out the data to the data set
                dset = self.samples
            else:
                # we'll need to create a group colvars fisrt
                f_handle.create_group("colvars")
                # now flush the dataset
                dset = f_handle['colvars'].create_dataset("cv_" + str(len(f_handle['colvars'].keys())), np.asarray(self.samples).shape, dtype="f")
                dset = self.samples
         
       
        with open(filename + ".colvars", "a") as f_handle:
            np.save(f_handle, self.samples)
        """
        # write out the basis function time series as well            
        with h5py.File(filename + ".hdf5", "a") as f_handle:
            # check to see if we've created a group for the colvars writes
            if 'timeSeries' in f_handle:
                # if so, add in a dataset for this flush
                dset = f_handle['timeSeries'].create_dataset("ts_" + str(len(f_handle['timeSeries'].keys())), np.asarray(self.basisFnxTimeSeries).shape, dtype="f")
                dset = self.basisFnxTimeSeries
            else:
                # we'll need to create a group colvars fisrt
                f_handle.create_group("timeSeries")
                # now flush the dataset
                dset = f_handle['timeSeries'].create_dataset("ts_" + str(len(f_handle['timeSeries'].keys())), np.asarray(self.basisFnxTimeSeries).shape, dtype="f")
                dset = self.basisFnxTimeSeries
        # so now we will also cast text output files
        with open(filename + ".timeSeries", "a") as f_handle:
            np.savetxt(f_handle, self.basisFnxTimeSeries)

        """

        # delete current reference to the list
        del self.samples
        #del self.basisFnxTimeSeries
        #del self.configs
            
        #now replace the data structure with the endpoint
        self.samples = [lastSample]
        #self.basisFnxTimeSeries = [lastBasisFunction]
        #self.configs = [lastConfig]

        return 0 
    
     
class Box(basisFunction):
    """
    This class impements a rectangle in CV space. 
    """
    def __init__(self, center, width, periodicLength = None):
        """
        Initialize the box object.
        This takes as input:
            center            The center of the box in the collective variable space.  
                                This should ideally be a numpy/scipy array or list, but other iterables MIGHT work.
            width              This is the width out from the center in each collective coordinate.
                                The actual width of the box in each coordinate is twice this vector.
            periodicLength        The wrapping object takes a little explanation.  It is an optional array, 
                                where the i'th element corresponds to the i'th collective variable.  If the 
                                collective variable wraps around, the corresponding element is the range 
                                of the collective variable, e.g. 360 for an angle going from -180 to 180.  If the 
                                collective varable does not wrap around, the corresponding element is 
                                just 0.  If none of the collective variables 
                                wrap, leaving wrapping as None is perfectly fine.
        
        """
        # construct data structures for the basis function. 
        # center and width
        self.center = np.asarray(center)
        self.width = np.asarray(width)
        self.dimension = len(center)
      
        # set the radius to use for computing neighbors.
        self.radius = np.linalg.norm(self.width)**2
        self.neighborList = []
        
        
        # for storing CV's and configs
        self.samples = []
        self.numSamples = 0
        self.configs = []
        
        # if the cv is periodic, build a wrapping array 
        if periodicLength != None:
            self.wrapping = np.asarray(periodicLength)
        else:
            # set to the default value of None if not periodic
            self.wrapping = periodicLength
        
        #print "Box created at", center, "with widths", width
        
    def __call__(self, coord, umb):
        """
        Return the value of the basis function at this coordinate. 
        """
        # get sum of box indicators. 
        if self.indicator(coord) == 0.0:
            return 0.0
        else:
            norm = 0.0
            
            if len(self.neighborList) > 0:
                for i in range(len(self.neighborList)):
                    norm += umb[self.neighborList[i]].indicator(coord)
            else:
                for i in range(len(umb)):
                    norm += umb[i].indicator(coord)
                    
            return 1.0 / norm
            
    def indicator(self, coord):
        """
        Return the value of the indicator of the box at this point.  This will be 1.0 if coord is contained inside the box, and 0.0 if it is not.
        """
        # create a distance vector
        distancevec = sp.asarray(coord) - sp.asarray(self.center)
        # if any collective variable is periodic, construct dr, the adjuct for minimum image convetion for the periodic cv's
        if self.wrapping != None:
            # build dr 
            dr = np.zeros(distancevec.shape)
            # add values to dr if the CV wraps
            for i in xrange(len(self.wrapping)):
                if self.wrapping[i] != 0.0:
                    # This is an old trick from MD codes to find the minimum distance between two points.
                    dr[i] = self.wrapping[i] * np.rint(distancevec[i]/self.wrapping[i])
            # add min image vector
            distancevec -= dr
           # print distancevec
            
        # We return 1.0 if all the distances are smaller than the width of the box from the center, 0.0 otherwise.
        return float(np.prod(self.width > np.abs(distancevec)))
        
        
class Gaussian(basisFunction):
    """
    This class implements a Gaussian function. 
    """
    def __init__(self, mu, sig):
        """
        Initializes a Gaussian at a point with a given width in each dimension. 
        """
        # the center of the gaussian 
        self.center = np.array(mu)
        # defines the stdev in each dimension
        self.width = np.array(sig)
        # a list for samples of CV space
        self.samples = []
        # configurations buffer
        self.configs = []
        # stores the dimension of the CV space. 
        self.dimension = self.center.size
        
        #print "Gaussian created at ", mu, " with stdev ", sig
        
    def __call__(self, coord, umb):
        """
        Return the value of the basis function at this coordinate. 
        """
        # get sum of box indicators. 
        norm = 0.0
        for i in range(len(umb)):
            norm += umb[i].indicator(coord)
            
        return self.indicator(coord) / norm
            
    def indicator(self, coord):
        """
        Function that returns the value of the Gaussian at a point.
        """
        # first enforce nearest image convention
        dist = (coord-self.center)-360.0*np.rint((coord-self.center)/360.0)
        # compute the exponentials
        vals = np.exp(-(dist)**2/ (2.0*self.width**2))
        # return the product
        return np.prod(vals)


class Pyramid(basisFunction):
    """
    This class implements a pyramidal basis function
    - Erik 

    """
    def __init__(self, center, width, periodicLength = None):
        """
        Initializes the pyramid at a center, given a width in each dimension.
        """
        self.center = np.asarray(center)
        self.width = np.asarray(width)
        self.dimension = len(center)
        self.radius = np.sqrt(np.sum(self.width**2))
        
        # We calculate the slope of the pyramid.
        self.slopes = 1.0/self.width
        
        # We check if the box wraps around.
        if periodicLength != None:
            self.wrapping=np.asarray(periodicLength)
        else:
            self.wrapping = None
        #print "Gaussian created at ", mu, " with stdev ", sig
    
            
    def __call__(self, coord, umb):
        """
        Return the value of the basis function at this coordinate. 
        """
        # get sum of box indicators. 
        if self.indicator(coord) == 0.0:
            return 0.0
        else:
            norm = 0.0
            for i in range(0, len(umb), 1):
                norm += umb[i].indicator(coord)
            
            return 1.0 / norm
            
    def indicator(self, coord):
        """
        Return the value of the indicator of the box at this point.  This will be 1.0 if coord is contained inside the box, and 0.0 if it is not.
        """
        # create a distance vector
        distancevec = np.asarray(coord) - np.asarray(self.center)
        # if any collective variable is periodic, construct dr, the adjuct for minimum image convention for the periodic cv's
        if self.wrapping != None:
            # build dr 
            dr = np.zeros(distancevec.shape)
            # add values to dr if the CV wraps
            for i in xrange(len(self.wrapping)):
                if self.wrapping[i] != 0.0:
                    # This is an old trick from MD codes to find the minimum distance between two points.
                    dr[i] = self.wrapping[i] * np.rint(distancevec[i]/self.wrapping[i])
            # add min image vector
            distancevec -= dr

        # We calculate the value of 
        psiparts = 1.0-self.slopes*np.abs(distancevec)
         
        # We remove negative entries and return the minimum value.
        return min(psiparts.clip(min=0))

        
