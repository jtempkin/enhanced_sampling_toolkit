"""
Class definition for Pyramid shaped window. Inherets from Window.
"""

from neus.window import Window
import numpy as np
import collections


class Pyramid(Window):
    """Class definition of Pyramid which inherits from window. Implements the pyramid shaped basis function.
    """

    def __init__(self, center, width, ref_center=None, ref_width=None, time=None, periodic_length=None, max_list_size=100, initial_conditions=[], initial_conditions_probability=0.0):
        """Create an intsance of the Pyramid object.

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

        # call parent constructor
        Window.__init__(self, center, width, ref_center=ref_center, ref_width=ref_width, time=time, periodic_length=periodic_length, max_list_size=max_list_size, initial_conditions=initial_conditions, initial_conditions_probability=initial_conditions_probability)

        # We calculate the slope of the pyramid.
        self.slopes = 1.0/self.width

        return None

    def __repr__(self):
        """Return the string represetation of the Pyramid instance.

        Returns
        ----------
        string
            A string representation of the Pyramid instance.
        """
        id = "Pyramid(" + str(self.center) + ", " + str(self.width) + ")"

        return id

    def __call__(self, walker):
        """Return the value of the support for this Pyramid object.

        Parameters
        -------------
        walker : walker instance
            The walker instance for which to evaluate the support of the Pyramid.

        Returns
        ------------
        float
            The value of the support.
        """
        # DEVLEOPER Note: this is the key for how the NEUS application module relies on the walker object definition. We have to think carefully about how exactly we want the window object accept input in the call function. Should this explicitly make assumptions about the structure and callable functions of the walker class or should we try to generalize this to act through something like a numpy array?

        # second comment: This question arises separately really from how we wish this call function is implemented from a partition object. The partition object will

        # ok, here we are going to implement the window on top if the walker object definition. But we will enforce that the behavoir of Pyramid's call will depend on both the definition of the reference and the information it can get from walker.

        # check to see that we've received a walker object
        try:
            coord = walker.get_colvars()
        except AttributeError:
            coord = np.array(walker)

        assert coord.shape == self.center.shape, "walker collective variables array does not match the shape of this window instance definition."

        # create a distance vector
        distancevec = coord - self.center

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

            # check to see that
            if not hasattr(walker, "get_initial_colvars"):
                raise Exception("Walker object passed to Pyramid __call__ does not have support for getting refernce collective variable value.")

            # return initial state of the collective variable
            ref_coord = walker.get_initial_colvars()

            val *=self.ref_indicator(ref_coord)

        if self.time_start is not None:
            # check that the passed walker object has a time coordinate
            if not hasattr(walker, "get_time"):
                raise Exception("Walker object passed to Pyramid __call__ does not have support for getting time value.")

            # get time
            time = walker.get_time()

            # return indicator and multiply against support value
            val *= self.time_indicator(time)


        #  return the minimum value.
        return val

    def ref_indicator(self, coord):
        """Return the value of the support for the reference phase space point.

        Parameters
        -------------
        coord : numpy.ndarray
            The coordinates of the reference phase space point.

        Returns
        ------------
        float
            The value of the support on the reference coordinate.
        """

        assert coord.shape == self.ref_center.shape, "walker reference collective variables array does not match the shape of this window reference center."

        # create a distance vector
        distancevec = coord - self.ref_center

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

    def time_indicator(self, time):
        """Return the indicator function on the time interval.

        Takes the value 1.0 if the time provided is in the time interval for which this window has nonzero support.

        Parameters
        ------------
        time : int
            A time to evaluate.

        Returns
        ----------
        float
            The indicator value.
        """
        if self.time_start <= time < self.time_end:
            return 1.0
        else:
            return 0.0
