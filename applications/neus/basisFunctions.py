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
import collections
from entryPoints import entry_point


class basisFunction:
    """
    It is important to note that the configurations associated with 
    sampling of this box are stored at this level of the class structure.
    
    At this level, the characteristics defining the basisFunction class is
    a set of configurations inherent to the class. 
    """

    def __init__(self, center, width, ref_center=None, ref_width=None, time=None, periodicLength = None, max_entrypoints = 500):
        """
        Defines the initialization routine for the windows super class. Gets called by the subclass types.

        Initializes the pyramid at a center, given a width in each dimension.

        The parameters passed define a spatial current point, reference point and a time boundary. 

        These parameters support sophisticated definitions of the process J(t). 
        """

        # define the parameters passed. 
        self.center = np.asarray(center)
        self.width = np.asarray(width)

        # add the reference centers and widths if needed
        if ref_center is not None:
            self.ref_center = np.asarray(ref_center)
        else:
            self.ref_center = None

        if ref_width is not None:
            self.ref_width = np.asarray(ref_width)
        else: 
            self.ref_width = None

        # if we get passed a time boundary, let's save the start and end points
        if time is not None:
            self.time_start = time[0]
            self.time_end = time[1]
        else:
            self.time_start = None
            self.time_end = None

        # record the dimension of the space
        self.dimension = len(self.center)

        # define a radius the usual way in this space
        self.radius = np.sqrt(np.sum(self.width**2))

        # this parameter sets the last phase space point at each iteration such that 
        self.walker_restart = None

        self.local_observables = set()

        # this sets the initial distribution at time zero. Used for finite time calculations where a
        # time zero injection distribution might be used. 
        self.initial_distribution_prob = 0.0

        self.neighborList = None
        
        # We check if the box wraps around a periodic variable
        if periodicLength is not None:
            self.wrapping=np.asarray(periodicLength)
        else:
            self.wrapping = None

        self.active = False

        # maximum number of entry points stored.
        self.max_entrypoints = max_entrypoints

        self.neighbor_prob = None

        return None

    def __len__(self):
        """
        Returns the number of entry points in the active buffer.
        """

        count = 0 

        for key in self.entryPoints: 
            count += len(self.entryPoints[key])

        return count
    
    def flush_data_to_file(self, filename):
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

    def reinject(self, walker):
        """
        This function initializes a simulation from the entry point list in the 
        current umbrella.
        """
        # let's get the current estimate of the flux
        #prob = self.z * self.F[:,i]

        # zero out flux from i to i
        #prob[i] = 0.0

        """
        # now let's zero out any neighbors with potentially nonzero flux but no stored entry points
        for indx in range(prob.size):
            # lets get the key to this index
            key = self.index_to_key[indx]
            # check size of this neighbor specifically and zero out probability if zero
            if self._umbrellas[i].getNumberOfEntryPoints(key=key) == 0:
                prob[indx] = 0.0
        """

        # normalize probability
        #assert prob.sum() > 0.0
        #prob /= prob.sum()

        # now choose neighbor proportional to prob 
        #I = np.random.choice(np.arange(prob.size), p=prob)

        # get the entry point from the umbrella window

        #assert self._umbrellas[i].getNumberOfEntryPoints(key=self.index_to_key[I])

        # we choose from the initial distribution with prob. stored in this window
        if random.random() < self.initial_distribution_prob:
            dist = self.initial_distribution
            EP = random.sample(dist, 1)[0]

        else:
            EP = self.get_entry_point()

        walker.setConfig(EP.q)
        walker.setVel(EP.p)

        # set the lag component of the walker state
        walker.Y_s = (EP.ref_q, EP.ref_p)

        walker.simulationTime = EP.time

        return 0


    def get_entry_point(self):
        """
        This routine returns an entry point from the buffers draw proportional to the entry fluxes.
        """

        assert hasattr(self, 'entryPoints'), 'Window has not initialized entry point library.'
    
        try:
            choice = np.random.choice(self.flux[self.neighborList].size, p=self.flux[self.neighborList])
        except ValueError:
            print self.flux[self.neighborList]

        # assertion checking that the neighbor has entry points stored. 
        assert len(self.entryPoints[choice]) > 0, "Error in entry point fluxes. Attempted to draw from flux zero neighbor."

        # draw the point uniformly from that buffer
        EP = random.choice(self.entryPoints[choice])

        return EP

    def update_entry_points_fluxes(self, flux):
        """
        This routine updates the sizes of the discretization of the entry point flux lists based on the maximum size of the entry point buffer.

        The flush keyword moves all stored points to the active buffer regardless of the fluxes
        """

        # set the flux
        self.flux = flux

        for i in self.neighborList:
            # if we don't find entries, set the flux to zero so we don't draw from this window
            if len(self.entryPoints[i]) == 0: 
                if self.flux[i] > 0.0: print "zeroing out", i
                self.flux[i] = 0.0


        # remove negative fluxes
        #self.flux[self.flux < 0.0] = 0.0
        
        # normalize if needed, since this will be treated as a probability 
        if self.flux.sum() != 0.0: 
            self.flux /= self.flux.sum()
            #print "flux:", self.flux.sum(), self.flux
            #print "neighbor Flux:", self.flux[self.neighborList].sum(), self.flux[self.neighborList]

        return 0 

    def add_entry_point(self, ep, key):
        """
        This routine adds a new entry point to the current list
        """
        if not hasattr(self, 'entryPoints'): raise Exception('Window has not initialized entry point library.')

        if key is None: raise Exception('key was not provided to add entry point.')

        self.entryPoints[key].append(ep)

        return 0

    def add_new_entry_point(self, ep, key_from):
        """
        This routine add an entry point to the list of proposed entry points.
        """
        if not hasattr(self, 'newEntryPoints'): raise Exception('Window has not initialized new entry point library.')

        self.newEntryPoints[key_from].append(ep)

        return 0

    def empty_new_entry_points(self):
        """
        This routine empties the new entry points data structure locally.
        """

        for key in self.newEntryPoints:
            self.newEntryPoints[key].clear()

        return 0

    def initialize_entry_points(self, maxEntryPoints = 500, neighborList = None, size = None, keylist = None):
        """
        This routine initializes the entry point data structure in the umbrella.

        The entry point library contains two structures. The first structure is a library of accepted points which is used to draw from when get_entry_point() is called. The second is a buffer library in which new proposed points are stored. The first structure gets update from this second structure when the update_entry_points_fluxes() routine is called. 

        The first structure acts like a double ended queue, in that it has a maximum lengh (set by the max_entrypoints attribute) and new entries are added in a FIFO manner. If max_entrypoints is set to None, this list grows to arbitrary length. 

        This routine does not take max_entrypoints as an argument and sets the structure of these libraries to that size. This is destructive to 
        """

        if hasattr(self, 'entryPoints'): print "WARNING: Calling initialize_entry_points() is destructive to current library."

        self.entryPoints = {}

        #if keylist is None: print "WARNING: no keys were initialized in the entry point structure in this window."

        # if there is no neighborlist, add all windows 
        if self.neighborList is None:
            for neighbor in range(size):
                self.entryPoints[neighbor] = collections.deque(maxlen=maxEntryPoints)
        else:
            for neighbor in self.neighborList:
                self.entryPoints[neighbor] = collections.deque(maxlen=maxEntryPoints)

        self.keylist = keylist

        # now we need to create a deep copy called newEntryPoints
        self.newEntryPoints = copy.deepcopy(self.entryPoints)

        return 0

    def set_initial_distribution(self, initial_distribution, initial_distribution_prob):
        """
        This routine sets the initial distribution and support probability for this window. Used for reinjection routine.
        """

        self.initial_distribution = initial_distribution
        self.initial_distribution_prob = initial_distribution_prob

        return 0 

    def add_local_obserbale(self, obs):
        """
        This routine adds a local observable to this window.
        """

        self.local_observables.add(obs)

        return 0

    def remove_local_observales(self):
        """
        This routine removes all local observables contained by this window.
        """

        self.local_observables.clear()

        return 0

    def get_local_observales(self):
        """
        This routine returns a list of the local observables.
        """

        return self.local_observables


class Box(basisFunction):
    """
    This class impements a rectangle in CV space. 
    """
    def __init__(self, center, width, periodicLength = None, max_entrypoints = 250):
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

    def __repr__(self):
        """
        Definition of repr special method to support identification of the basis function.
        """

        id = "basisFunctions.Gaussian(" + str(self.center) + ", " + str(self.width) + ")"

        return id

    def __init__(self, mu, sig, max_entrypoints = 250):
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
    This class implements a pyramidal basis function. 
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

    def __repr__(self):
        """
        Definition of repr special method to support identification of the basis function.
        """

        id = "basisFunctions.Pyramid(" + str(self.center) + ", " + str(self.width) + ")"

        return id


    def __init__(self, center, width, ref_center=None, ref_width=None, time=None, periodicLength = None, max_entrypoints = 500):
        """
        Initializes the pyramid at a center, given a width in each dimension.
        """

        # call parent constructor
        basisFunction.__init__(self, center, width, ref_center=ref_center, ref_width=ref_width, time=time, periodicLength=periodicLength, max_entrypoints=500)
        
        # We calculate the slope of the pyramid.
        self.slopes = 1.0/self.width

        return None
            
    def __call__(self, wlkr):
        """
        Return the value of the basis function at this coordinate. 
        """
        coord = wlkr.getColvars()
        ref_coord = wlkr.ref_cv
        time = wlkr.simulationTime

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

        val = min(psiparts.clip(min=0))

        if self.ref_center is not None:
            val *=self.ref_indicator(ref_coord)

        if self.time_start is not None:
            val *= self.time_indicator(time)

         
        # We remove negative entries and return the minimum value.
        return val

        """

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

        """


    def __bool__(self):
        """
        Returns True if the window is active and False if it is inactive.
        """

        return self.active


    def ref_indicator(self, coord):
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


    def time_indicator(self, coord):
        """
        Return the value of the indicator function for the time coordinate.

        This will be implemented currently as non-overlapping discretization in time. 
        """ 

        if self.time_start <= coord < self.time_end:
            return 1.0
        else:
            return 0.0
