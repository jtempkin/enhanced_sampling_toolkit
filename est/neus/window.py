"""The window object defines a single spatiotemporal restriction of the discretization defined by the NEUS :math:`J^{(t)}` process. An instance of the window object corresponds to a single value of the :math:`J^{(t)}`. The object defines a data structure for storing the entry point flux lists defined by :math:`\pi_{i}(t,xd)`. It also defines the routines require for constructing the flux lists and returning physically weighted elements. 

This object contains the following attributes:

* data structure for containing flux lists
* data structure for the initial conditions distribution at :math:`J^{(0)}`
* a list of physical fluxes used to draw entry points proportionally
* routines for reinjection
* routines for updating fluxes
* routines for updating entry point lists
"""

import numpy as np
from entry_points import entry_point
import collections
import copy
import random


class Window:
    def __init__(self, center, width, ref_center=None, ref_width=None, time=None, periodic_length = None, max_list_size=100, initial_conditions=[], initial_conditions_probability=0.0):
        """Create an intsance of a window object.

        Note that the width parameter refers to the maximum distance in each direction for which this window is defined to have nonzero support.

        Parameters
        ---------------
        center : numpy.ndarray, list
            The coordinates of the center of the window support.

        width : numpy.ndarray, list
            The width of the support of the window.

        ref_center : numpy.ndarray, list (None)
            The center of the support in the reference configuration at time :math:`t=0`

        ref_width : numpy.ndarray, list (None)
            The width of the support of the reference configuration at time :math:`t=0`

        time : numpy.ndarray, list (None)
            The time interval for which this window has nonzero support.

        periodic_length : numpy.ndarray, list (None)
            The length of the periodicity for each coordinate this window supports.

        max_list_size : interval (100)
            The maximum size of each :math:`\gamma_{ij}` distribution.

        initial_conditions : iterable
            The iterable of entry point objects at time :math:`t=0`.

        initial_conditions_probability : float
            A float between 0.0 and 1.0.
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
        else:
            self.time_start = None
            self.time_end = None

        # set periodicity of coordinates.
        if periodic_length is not None:
            self.wrapping=np.asarray(periodicLength)
        else:
            self.wrapping = None

        # set initial t=0 distribution and probability with which to draw from this distribution
        self.initial_conditions = copy.deepcopy(initial_conditions)
        self.initial_conditions_probability = 0.0

        # store maximum size of flux list
        self.max_list_size = max_list_size

        # create a dictionary of neighboring flux lists. We will store sets of entry point objects here
        self.flux_lists = {}
        self.flux_weights = {}

        return None

    def __len__(self):
        """Return the total number of entry points stored in this window.

        Returns
        -------------
        len : int
            The total number of entry points stored in the flux lists.
        """

        size = 0

        # iterate through keys in flux list dictionary and add up size
        for key in self.flux_lists:
            size += len(self.flux_lists[key])

        return size

    def __nonzero__(self):
        """Return boolean value of the window.

        Returns True if there is at least one entry point stored in the flux lists or if there is at least one entry point stored in the initial conditions library. Returns False otherwise.

        Returns
        ----------
        bool
            True if there is at least one entry point stored in the flux lists. False otherwise.
        """

        if (len(self) > 0) or len(self.get_initial_conditions()) > 0:
            return True
        else:
            return False

    def __repr__(self):
        """Return the string representation of this window instance.

        Returns
        -----------
        repr : string
            Returns the string represenation of the window.
        """

        string = "window(" + str(self.center) + ")"

        return "window(" + str(self.center) + ")"

    def reinject(self):
        """Return an entry point drawn proportional to the neighbor fluxes.

        This routine will attempt to return an entry point from the initial distribution attribute with probability equal to value stored in the initial_conditions_probability attribute. With remaining probability, it attempts to draw an entry point from the weighted flux lists.

        Returns
        ------------
        ep : entry_point
            Returns an entry point.
        """

        # choose from initial distribution with probability stored
        if random.random() < self.initial_conditions_probability:
            ep = random.choice(self.initial_conditions)
        # draw from teh flux lists with remaining probability
        else:
            ep = self.get_entry_point()

        return ep

    def get_entry_point(self):
        """Return an entry point from the :math:`\gamma_{ij}`.

        This routine selects a :math:`\gamma_{ij}` with probability proportional to the value stored in the fluxes array attribute. Once the :math:`\gamma_{ij}` is chosen, the entry point is selected from :math:`\gamma_{ij}` proportional to the weight of each entry point.

        Returns
        -----------
        entry_point
            An entry point object.
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

        # now draw an entry point from the flux list according to entry weight
        p = np.array(self.flux_weights[neighbor])
        p /= p.sum()
        # we have to draw this value from the deque index
        indx = np.random.choice(np.arange(p.size), p=p)
        # draw the ep value
        ep = self.flux_lists[neighbor][indx]

        return ep

    def update_fluxes(self, fluxes):
        """Set neighbor fluxes.

        Sets the array of fluxes observed from the neighboring windows. 

        Parameters
        --------------
        fluxes : ndarray
            A numpy array of the fluxes from the neighbors.
        """

        self.fluxes = fluxes

        return None

    def add_entry_point(self, ep, key, weight=1.0):
        """Add a new entry point to the specified flux list.

        If key is not already know, create a new set with label "key" and add it to the dictionary. Keys should be the index of the corresponding neighbor window.

        Parameters
        ------------
        ep : entry_point
            The entry point to add to :math:`\gamma_{ij}`

        key : int
            The :math:`i` index of the neighbor.

        weight : float (1.0)
            The weight of the entry_point object.
        """

        #assert key is isinstance(key, int), "Supplied key was not an integer."

        # check to see if dictionary key exists for provided key. if not create new dictionary key implemented as a double ended queue from collections
        if key in self.flux_lists:
            self.flux_lists[key].append(ep)
            self.flux_weights[key].append(weight)
        else:
            self.flux_lists[key] = collections.deque(maxlen=self.max_list_size)
            self.flux_lists[key].append(ep)
            # create and add weight
            self.flux_weights[key] = collections.deque(maxlen=self.max_list_size)
            self.flux_weights[key].append(weight)

        return None

    def set_initial_conditions(self, distribution):
        """Set the distribution of initial conditions to this window. 

        This set corresponds the distribution of the process at J(0). Distribution is expected to be an iterable of entry point objects.

        Parameters
        ---------------
        distribution : iterable
            An iterable of entry_point objects.
        """

        self.initial_conditions = copy.deepcopy(distribution)

        return None

    def set_initial_conditions_probability(self, p):
        """Set the probability for drawing from the initial conditions.

        Parameters
        -------------
        p : float
            A float between 0.0 and 1.0.
        """

        assert p <= 1.0, "Probability is not valid."
        assert p >= 0.0, "Probability is not valid."

        self.initial_conditions_probability = p

        return None

    def update_initial_conditions(self, items):
        """Adds elements of items to initial conditions array
        """
        for i in items:
            self.initial_conditions.append(items)

        return None

    def get_initial_conditions(self):
        """Return the distribution of initial conditions.

        Returns
        ----------
        iterable
            An iterable of the initial conditions.
        """

        return self.initial_conditions

    def clear_flux_list(self):
        """Clear the :math:`\gamma_{ij}` of all entries.
        """

        self.flux_lists.clear()
        self.flux_weights.clear()

        return None

    def get_flux_lists(self):
        """Return the dictionary of the flux lists.

        Returns
        -------------
        dictionary
            The dictionary of the :math:`\gamma_{ij}`
        """
        return self.flux_lists

    def get_flux_weights(self):
        """Return the dictionary of the entry_point weights.

        Returns
        -----------------
        dictionary
            A dictionary of the entry point weights.
        """
        return self.flux_weights
