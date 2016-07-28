# -*- coding: utf-8 -*-
"""
This module contains the definition of the partition object. The partition object represents the definition of a spatial decomposition defined in umbrella samping / stratification algorithms.

The basic features of this partition and how one uses this data structure to perform umbrella sampling type calculations is described as follows.

The core class of this module is the partition class.
"""

import numpy as np
import copy
import scipy as sp
from scipy import linalg as LA
from scipy import sparse
from scipy.sparse import linalg as LA_sparse
import random
import basisFunctions
import entryPoints
import observables

try:
    from mpi4py import MPI
except:
    print "assuming serial job."
    pass

class partition:
    """
    This class defines a set of windows that partition the sampling space.
    This class will contain an array of basisFunction objects.
    """
    def __init__(self, scratchdir=None, parallel=False):
        """
        Init routine. Initializes the following values:

        The obervables list is a bookkeeping of the observables of the self one
        wants to record during the simulaiton. This list needs to be set prior to sampling.

        The list shoud contain elements that are tuples of a function and a data array.

        Each element should take two arguments, one the walker object and one the umbrella index for which the sample is
        associated.
        """

        # given a list of center locations and widths, initialize a

        # create an umbrella list
        self._umbrellas = []

        # keep a list of any observable objects that should be estimated during sampling.
        self.observables = []

        return None

    def __call__(self, walker):
        """
        Returns an array of the support of the partition at the current CV state of the walker.
        """
        # build an array
        # ***** SHOULD LOOK TO SEE IF LIST COMPREHENSIONS ARE BETTER OR WORSE THAN CREATING AN ITERATOR OBJECT *******
        indicators = np.asarray([win(walker) for win in self])

        # normalize the values of the basis functions
        if np.sum(indicators) == 0.0:
            print walker.getColvars()
            print walker.simulationTime
            print walker.Y_s
        assert np.sum(indicators) != 0.0, str(walker.getColvars())

        return indicators / indicators.sum()


        coord = walker.getColvars()

        if umbrella_index is None:

            # go through the basis functions and construct the indicators array
            for i in xrange(len(self._umbrellas)):
                if self._umbrellas[i].indicator(coord) == 0.0:
                    continue
                else:
                    indicators[i] = self._umbrellas[i](walker, self._umbrellas)

        elif len(self._umbrellas[umbrella_index].neighborList) == 0:
            # go through the basis functions and construct the indicators array
            for i in xrange(len(self._umbrellas)):
                if self._umbrellas[i](walker, self._umbrellas) == 0.0:
                    continue
                else:
                    indicators[i] = self._umbrellas[i](walker, self._umbrellas)
        else:

            # make sure the partition has a neighborlist
            assert hasattr(self._umbrellas[umbrella_index], 'neighborList'), "There is no neighborlist defined."
            # now loop over neighbors and populate the indicator
            for i in self._umbrellas[umbrella_index].neighborList:

                indicators[i] = self._umbrellas[i](wlkr, self._umbrellas)

            # if we don't find any support, let's try this again and search the whole space.
            if np.sum(indicators) == 0.0:
                indicators = self.get_basis_function_values(wlkr, umbrella_index=None)

        # normalize the values of the basis functions
        if np.sum(indicators) == 0.0:
            print wlkr.getColvars()
            print wlkr.simulationTime
            print wlkr.Y_s
        assert np.sum(indicators) != 0.0, str(wlkr.getColvars())

        indicators = indicators / np.sum(indicators)

        return indicators

    def __getitem__(self, pos):
        """
        Definition of special method for returning window objects.
        """

        return self._umbrellas[pos]

    def __len__(self):
        """
        Definition of the size of the partition. Returns number of windows.
        """

        return len(self._umbrellas)

    def win_index(self, win):
        """
        Returns the index of teh given window in the partition.
        """
        return self._umbrellas.index(win)

    def set_umbrellas(self, umbrellas, neighborList=True,  s=None):
        """
        This routine replaces the current list of umbrellas with those provided
        and updates the matricies contained to match the new umbrella list. This
        behavior is destructive to the old matricies.
        """

        self._umbrellas = umbrellas

        # get the size of the new partition
        N = len(self)

        # build the neighborlist in the umbrellas by default.
        if neighborList:

            assert s is not None, "Stopping time needed to construct neighborlist."

            self.build_neighbor_list(umbrellas, s)
        # if there is no neighbor list, set variable to be every other window.
        else:

            for win in self:
                win.neighborList = np.arange(N)

        # initialize the matricies needed for the NEUS
        self.M = np.zeros((N,N))
        self.F = np.zeros((N,N))

        # initialize the F_list
        self.F_list = []

        [self.F_list.append([]) for i in range(N)]

        # set a and z
        self.a = np.zeros(N)
        self.z = np.zeros(N)

        # we should track how many samples are taken in the M matrix
        self.nsamples_M = np.zeros(N)

        # start a list of the windows associated with this partition.
        self.rank_window_index = None

        self.simulationTime = np.zeros(N)

        # here, the variable k represents how many times the transition matrix has been updated.
        self.k = np.zeros(N)

        # we'll also store a boolean of active windows
        self.active_windows = np.zeros(N, dtype=np.bool)

        return 0

    def get_basis_function_values(self, wlkr, umbrella_index = None):
        """
        This function takes a point in collective variable space and returns
        an array of the value of the basis functions at that point.

        If no umbrella index is passed, then search the whole space
        """
        # build an array
        indicators = np.zeros(len(self._umbrellas))

        coord = wlkr.getColvars()

        if umbrella_index is None:

            # go through the basis functions and construct the indicators array
            for i in xrange(len(self._umbrellas)):
                if self._umbrellas[i].indicator(coord) == 0.0:
                    continue
                else:
                    indicators[i] = self._umbrellas[i](wlkr, self._umbrellas)

        elif len(self._umbrellas[umbrella_index].neighborList) == 0:
            # go through the basis functions and construct the indicators array
            for i in xrange(len(self._umbrellas)):
                if self._umbrellas[i](wlkr, self._umbrellas) == 0.0:
                    continue
                else:
                    indicators[i] = self._umbrellas[i](wlkr, self._umbrellas)
        else:

            # make sure the partition has a neighborlist
            assert hasattr(self._umbrellas[umbrella_index], 'neighborList'), "There is no neighborlist defined."
            # now loop over neighbors and populate the indicator
            for i in self._umbrellas[umbrella_index].neighborList:

                indicators[i] = self._umbrellas[i](wlkr, self._umbrellas)

            # if we don't find any support, let's try this again and search the whole space.
            if np.sum(indicators) == 0.0:
                indicators = self.get_basis_function_values(wlkr, umbrella_index=None)

        # normalize the values of the basis functions
        if np.sum(indicators) == 0.0:
            print wlkr.getColvars()
            print wlkr.simulationTime
            print wlkr.Y_s
        assert np.sum(indicators) != 0.0, str(wlkr.getColvars())
        indicators = indicators / np.sum(indicators)

        return indicators



    def build_neighbor_list(self, umbrellas, s, debug=False):
        """
        This routine constructs a neighborlist based on the radius of the basis
        function.

        L is the vector specifying the periodic lengths in each dimension.
        """
        # we should make sure that ever basisFunction has a radius defined.
        for win in self:
            assert hasattr(win, 'radius')

        # now let's find all of the neighbors and populate a list of neighbors for each
        for i in range(len(self)):
            # add the current window to it's own neighborList
            umbrellas[i].neighborList.append(i)
            # now search the other windows

            for j in range(i+1, len(umbrellas)):
                # get the distance between the centers
                dr = umbrellas[i].center - umbrellas[j].center

                # apply minimum image convention if the dimension wraps
                for indx,dim in enumerate(umbrellas[i].wrapping):
                    if dim == -1.0:
                        continue
                    else:
                        # now apply the minimum image criteria
                        if abs(dr[indx]) > dim / 2.0:
                            if dr[indx] > 0.0:
                                dr[indx] -= dim
                            elif dr[indx] < 0.0:
                                dr [indx] += dim
                        else:
                            continue

                dist = np.linalg.norm(dr)
                # append i,j to each other's list
                #if dist <= 2.0 * (umbrellas[i].radius + umbrellas[j].radius):

                i_start = umbrellas[i].time_start % s
                i_end = umbrellas[i].time_end % s
                j_start = umbrellas[j].time_start % s
                j_end = umbrellas[j].time_end % s

                # we should also check here that we connect in time as well.
                if i_start == j_end:
                    if debug: print "added pair:", i, umbrellas[i].center, umbrellas[i].time_start, umbrellas[i].time_end, j, umbrellas[j].center, umbrellas[j].time_start, umbrellas[j].time_end
                    #umbrellas[j].neighborList.append(i)
                    umbrellas[i].neighborList.append(j)

                elif i_end == j_start:
                    if debug: print "added pair:", i, umbrellas[i].center, umbrellas[i].time_start, umbrellas[i].time_end, j, umbrellas[j].center, umbrellas[j].time_start, umbrellas[j].time_end
                    #umbrellas[i].neighborList.append(j)
                    umbrellas[j].neighborList.append(i)

                elif i_start == j_start:
                    if debug: print "added pair:", i, umbrellas[i].center, umbrellas[i].time_start, umbrellas[i].time_end, j, umbrellas[j].center, umbrellas[j].time_start, umbrellas[j].time_end
                    umbrellas[i].neighborList.append(j)
                    umbrellas[j].neighborList.append(i)


        return 0

    def build_keylist_to_index_map(self, keylist):
        """
        This routine takes the input keylist and constructs internal dictionaries that convert between indicies of F
        and elements of the keylist.
        """
        # lets put the keylist in the partition object
        self.keylist = keylist

        # these dictionaries convert between the index of F and the key
        self.key_to_index = {}
        self.index_to_key = {}

        # populate the dictionaries
        for i,key in enumerate(keylist):
            self.key_to_index[key] = i
            self.index_to_key[i] = key

        return 0
