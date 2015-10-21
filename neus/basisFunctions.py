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
import scipy as sp
import random
import copy

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

        if filename is not None:
            with open(filename + ".colvars", "a") as f_handle:
                np.savetxt(f_handle, self.samples)

        # delete current reference to the list
        del self.samples
            
        #now replace the data structure with the endpoint
        self.samples = [lastSample]

        return 0

    def getEntryPoint(self, key = None):
        """
        This routine returns an entry point from the entry point.
        """
        if not hasattr(self, 'entryPoints'): raise Exception('Window has not initialized entry point library.')

        if self.eptype == 'set':

            EP = random.sample(self.entryPoints, 1)

        elif self.eptype == 'dictionary':

            if key is None: raise Exception('key was not provided to draw entry point.')

            #EP = random.sample(self.entryPoints[key], 1)
            EP = random.sample(self.entryPointsList, 1)

        return EP[0]

    def updateEntryPointsFluxes(self, flux, key_map, flush=False):
        """
        This routine updates the sizes of the discretization of the entry point flux lists based on the maximum size of the entry point buffer.

        The flush keyword moves all stored points to the active buffer regardless of the fluxes
        """

        assert self.eptype == 'dictionary'

        self.entryPointsList.clear()

        # now build a collection of points whose composition is constructed based on the current estimate of the flux. 

        for i in xrange(flux.size):
            size = int(self.max_entrypoints * flux[i])
            if size == 0: 
                continue
            if size > self.max_entrypoints:
                size = self.max_entrypoints

            nentries = min(len(self.entryPoints[key_map[i]]), size)

            update_list = self.entryPoints[key_map[i]][0:nentries]

            self.entryPointsList.update(update_list)

        return 0 

    def flushAllToBuffer(self):
        """
        This routine adds all entry points to the selection buffer regardless of the fluxes.
        """
        assert self.eptype == 'dictionary'

        self.entryPointsList.clear()
    
        for key in self.entryPoints.keys():
            self.entryPointsList.update(self.entryPoints[key])

        return 0

    def addEntryPoint(self, ep, key = None):
        """
        This routine adds a new entry point to the current list
        """
        if not hasattr(self, 'entryPoints'): raise Exception('Window has not initialized entry point library.')

        if self.eptype is 'set':

            if len(self.entryPoints) > self.max_entrypoints:
                discard = self.entryPoints.pop()
                self.entryPoints.add(ep)

            else: 
                self.entryPoints.add(ep)

        elif self.eptype is 'dictionary':

            if key is None: raise Exception('key was not provided to add entry point.')

            #self.entryPoints[key].add(ep)

            #print key

            if len(self.entryPoints[key]) > self.max_entrypoints:
                # if we're too large, pop the first (oldest) entry and add to the end
                discard = self.entryPoints[key].pop(0)

                self.entryPoints[key].append(ep)

            else: 
                self.entryPoints[key].append(ep)

        return 0

    def addNewEntryPoint(self, ep, key_from):
        """
        This routine add an entry point to the list of proposed entry points.
        """
        if not hasattr(self, 'newEntryPoints'): raise Exception('Window has not initialized new entry point library.')

        self.newEntryPoints.add((ep, key_from))

        return 0

    def emptyNewEntryPoints(self):
        """
        This routine empties the new entry points data structure locally.
        """

        self.newEntryPoints.clear()

        return 0

    def initializeEntryPoints(self, eptype='set', keylist = None):
        """
        This routine initializes the entry point data structure in the umbrella.

        The structure of the entry Point library is designed to be flexible as to how one wants to track entries from
        neightbors. There are two schemes that are supported now. The first option is a single set containing all
        entry points to this window (we're actually using the Python set structure). The second type is a dictionary
        of sets. This allows one to group contributions to the entry point list for this window by assigning a key to
        each neighbor and add/draw entry points from these separate sets by passing a key value to the respective add/draw
        operations.

        To construct the second type of structure you specify the key structure with a list of keys that define the
        categories from which one can group entry points.
        """

        if eptype == 'set':

            self.eptype = eptype

            self.entryPoints = set()

        elif eptype == 'dictionary':

            self.eptype = eptype

            self.entryPoints = {}

            #if keylist is None: print "WARNING: no keys were initialized in the entry point structure in this window."

            if len(self.neighborList) > 0:
                for neighbor in self.neighborList:
                    self.entryPoints[neighbor] = []
            else:
                for key in keylist:
                    self.entryPoints[key] = []

            #for key in keylist:
                #self.entryPoints[key] = []

            self.keylist = keylist

            self.entryPointsList = set()

        else:
            raise Exception('eptype was not understood.')

        # now we need to create a deep copy called newEntryPoints
        self.newEntryPoints = set()

        return 0

    def addLocalObserbale(self, obs):
        """
        This routine adds a local observable to this window.
        """

        self.local_observables.append(obs)

        return 0

    def removeLocalObservales(self):
        """
        This routine removes all local observables contained by this window.
        """

        self.local_observables = []

        return 0

    def getLocalObservales(self):
        """
        This routine returns a list of the local observables.
        """

        return self.local_observables

    def getNumberOfEntryPoints(self, key=None):
        """
        This routine counts the size of the entrypoint library and returns it.
        """
        count = 0

        if self.eptype == 'set':
            count = len(self.entryPoints)
        elif self.eptype == 'dictionary':
            if key is None:
                for key in self.entryPoints.keys():
                    count += len(self.entryPoints[key])
            else:
                #count = len(self.entryPoints[key])
                count = len(self.entryPointsList)

        return count

class Box(basisFunction):
    """
    This class impements a rectangle in CV space. 
    """
    def __init__(self, center, width, periodicLength = None, max_entrypoints = 10000):
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
        self.radius = np.sqrt(np.sum(self.width**2)) * 2.0
        #self.radius = np.linalg.norm(self.width)**2
        self.neighborList = []
        self.max_entrypoints = max_entrypoints
        
        
        # for storing CV's and configs
        self.samples = []
        self.numSamples = 0
        self.configs = []
        
        # if the cv is periodic, build a wrapping array 
        if periodicLength is not None:
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
        if self.wrapping is not None:
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
    def __init__(self, mu, sig, max_entrypoints = 10000):
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
        self.max_entrypoints = max_entrypoints
        
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
    def __init__(self, center, width, ref_center=None, ref_width=None, time=None, periodicLength = None, max_entrypoints = 500):
        """
        Initializes the pyramid at a center, given a width in each dimension.
        """
        self.center = np.asarray(center)
        self.width = np.asarray(width)
        self.ref_center = ref_center
        self.ref_width = ref_width
        # if we get passed a time boundary, let's save the start and end points
        if time is not None:
            self.time_start = time[0]
            self.time_end = time[1]
        else:
            self.time_start = None
        self.dimension = len(center)
        self.radius = np.sqrt(np.sum(self.width**2))
        self.max_entrypoints = max_entrypoints
        self.walker_restart = None
        
        self.neighborList = []
        
        # for storing CV's and configs
        self.samples = []
        self.numSamples = 0
        self.configs = []
        
        # We calculate the slope of the pyramid.
        self.slopes = 1.0/self.width
        
        # We check if the box wraps around.
        if periodicLength is not None:
            self.wrapping=np.asarray(periodicLength)
        else:
            self.wrapping = None
        #print "Gaussian created at ", mu, " with stdev ", sig
    
            
    def __call__(self, wlkr, umb):
        """
        Return the value of the basis function at this coordinate. 
        """
        coord = wlkr.getColvars()
        ref_coord = wlkr.Y_s[2]

        # adjust for the current 
        coord = coord[0:self.dimension]
        ref_coord = coord[0:self.dimension]

        time = wlkr.simulationTime

        # get sum of box indicators for 
        if self.indicator(coord) == 0.0:
            return 0.0

        elif len(self.neighborList) > 0:
            norm = 0.0
            for i in self.neighborList:
                norm += umb[i].indicator(coord)

        else: 
            norm = 0.0
            for win in umb:
                norm += win.indicator(coord)
            
        assert norm != 0.0

        if self.ref_center is not None:
            norm *= self.refIndicator(ref_coord)
        if self.time_start is not None:
            norm *= self.timeIndicator(time)

        if norm > 0.0:
            return self.indicator(coord) / norm
        else: 
            return 0.0

    def indicator(self, coord):
        """
        Return the value of the indicator of the box at this point.  This will be 1.0 if coord is contained inside the box, and 0.0 if it is not.
        """
        # create a distance vector
        distancevec = np.asarray(coord) - np.asarray(self.center)

        # if any collective variable is periodic, construct dr, the adjuct for minimum image convention for the periodic cv's
        if self.wrapping is not None:

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

    def refIndicator(self, coord):
        """
        Return the value of the indicator for the reference coordinates if appropriate.

        NOTE: Currently we will simply implement the reference discretization as non-overlapping boxes. 
        """
        # create a distance vector
        distancevec = sp.asarray(coord) - sp.asarray(self.ref_center)

        # if any collective variable is periodic, construct dr, the adjuct for minimum image convetion for the periodic cv's
        if self.wrapping is not None:

            # build dr 
            dr = np.zeros(distancevec.shape)

            # add values to dr if the CV wraps
            for i in xrange(len(self.wrapping)):
                if self.wrapping[i] != 0.0:
                    # This is an old trick from MD codes to find the minimum distance between two points.
                    dr[i] = self.wrapping[i] * np.rint(distancevec[i]/self.wrapping[i])

            # add min image vector
            distancevec -= dr
            
        # We return 1.0 if all the distances are smaller than the width of the box from the center, 0.0 otherwise.
        return float(np.prod(self.ref_width > np.abs(distancevec)))


    def timeIndicator(self, coord):
        """
        Return the value of the indicator function for the time coordinate.

        This will be implemented currently as non-overlapping discretization in time. 
        """ 

        if self.time_start <= coord < self.time_end:
            return 1.0
        else:
            return 0.0
