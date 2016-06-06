"""
An implementation of the window object in the NEUS toolkit.
"""

import numpy as np
from entryPoints import entry_point
import collections
import copy
import random

class window:
    """
    The window object defines a single piece of the discretization defined by the NEUS J(t) process. An instance of the window object corresponds to a single value of the J. The object defines a data structure for storing the entry point flux lists defined by ~pi. It also defines a routine for returning an entry point based on the defined fluxes.

    This object contains the following items:
 
    * data structure for containing flux lists
    * data structure for the initial conditions distribution at J(0)
    * a list of physical fluxes used to draw entry points proportionally
    * routines for reinjection
    * routines for updating fluxes
    * routines for updating entry point lists

    """

    def __init__(self, center, width, ref_center=None, ref_width=None, time=None, periodic_length = None, max_list_size=100, initial_conditions=[], initial_conditions_probability=0.0, a=0.0):
        """
        Create an intsance of a window object.
        """
        # initialize the handed parameters
        self.center = np.asarray(center)
        self.width = np.asarray(width)

        if ref_center is not None:
            self.ref_center = np.asarray(ref_center)
        else:
            self.ref_center = ref_center

        if ref_width is not None:
            self.ref_width = np.asarray(ref_width)
        else:
            self.ref_width = ref_width

        # set time boundary. None represents no time discretization.
        if time is not None:
            self.time_start = time[0]
            self.time_end = time[1]

        # set periodicity of coordinates.
        if periodic_length is not None:
            self.wrapping=np.asarray(periodicLength)
        else:
            self.wrapping = None

        # set initial t=0 distribution and probability with which to draw from this distribution
        self.initial_conditions = initial_conditions
        self.initial_conditions_probability = 0.0

        # the variable "a" here represents the entry in the initial distribution vector provided
        self.a = a

        # store maximum size of flux list
        self.max_list_size = max_list_size

        # create a dictionary of neighboring flux lists. We will store sets of entry point objects here
        self.flux_lists = {}

        return None

    def __len__(self):
        """
        Return the total number of entry points stored in this window.
        """

        size = 0

        # iterate through keys in flux list dictionary and add up size
        for key in self.flux_lists:
            size += len(self.flux_lists[key])

        return size

    def __nonzero__(self):
        """
        Return true if this window has at least one entry point to draw from in the initial distribution or the flux lists.
        """

        if (len(self) > 0) or len(self.get_initial_conditions()) > 0:
            return True
        else:
            return False

    def __repr__(self):
        """
        Return the string identifying this window instance.
        """

        string = "window(" + str(self.center) + ")"

        return "window(" + str(self.center) + ")"

    def reinject(self):
        """
        Return an entry point drawn proportional to the neighbor fluxes.
        """

        # choose from initial distribution with probability stored
        if random.random() < self.initial_conditions_probability:
            ep = random.choice(self.initial_conditions)
        # draw from teh flux lists with remaining probability
        else:
            ep = self.get_entry_point()

        return ep

    def get_entry_point(self):
        """
        Return an entry point from the store flux lists proportional to the entry point fluxes.
        """

        assert self.fluxes is not None

        # create probability array
        p = np.zeros(self.fluxes.size)

        # get fluxes for the keys in the flux list
        p[self.flux_lists.keys()] = self.fluxes[self.flux_lists.keys()]

        # normalize
        p /= p.sum()

        # now choose neighbor with probability p
        neighbor = np.random.choice(p.size, p=p)

        # now draw an entry point from the flux list
        ep = random.choice(self.flux_lists[neighbor])

        return ep

    def update_fluxes(self, fluxes):
        """
        Set neighbor fluxes as numpy array.
        """

        self.fluxes = fluxes

        return None

    def add_entry_point(self, ep, key):
        """
        Add a new entry point to the specified flux list. If key is not already know, create a new set with label "key" and add it to the dictionary.
        """

        #assert key is isinstance(key, int), "Supplied key was not an integer."

        # check to see if dictionary key exists for provided key. if not create new dictionary key implemented as a double ended queue from collections
        if key in self.flux_lists:
            self.flux_lists[key].append(ep)
        else:
            self.flux_lists[key] = collections.deque(maxlen=self.max_list_size)
            self.flux_lists[key].append(ep)

        return None

    def set_initial_conditions(self, distribution):
        """
        Set the distribution of initial conditions to this window. This corresponds the distribution of the process at J(0). Distribution is expected to be an iterable of entry point objects.
        """

        self.initial_conditions = copy.deepcopy(distribution)

        return None

    def set_initial_conditions_probability(self, p):
        """
        Set the probability for drawing from the initial conditions.
        """

        self.initial_conditions_probability = p

        return None

    def get_initial_conditions(self):
        """
        Return the distribution of initial conditions.
        """

        return self.initial_conditions

    def clear_flux_list(self):
        """
        Clear the flux list of all entries.
        """

        self.flux_lists.clear()

        return None

    def get_flux_lists(self):
        """
        Return the dictionary of the flux lists.
        """
        return self.flux_lists
