# -*- coding: utf-8 -*-
import abc


class walker(object):
    """Basic API specification of the walker object.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        """Instantiate a walker object.
        """
        return None

    @abc.abstractmethod
    def get_position(self):
        """Return the postiton of the walker.

        The return array is provided in the format [x1, y1, z1, x2, ...].

        Returns
        ---------
        numpy.ndarray
            A one-dimensional array of the current coordinates.
        """
        return None

    @abc.abstractmethod
    def set_position(self, config):
        """Set the configuration of the system.

        Takes an input array in the form of [x1, y1, z1, x2,..].

        Parameters
        --------------

        configuration : numpy.ndarray
            A one-dimentional numpy array of system coordinates.
        """
        return None

    @abc.abstractmethod
    def get_velocity(self):
        """Return the velocities of the system.

        The velocities are returned in a format that matches the getConfig() routine, specifically [v_x1, v_y1, v_z1, v_x2, v_y2, v_z2, ...] where v_x1 represents the x component of the velocity on the first particle, etc.

        Returns
        ------------
        numpy.ndarray
            A one dimensional numpy array of the current velocities.
        """
        return None

    @abc.abstractmethod
    def set_velocity(self, velocity):
        """Set the velocities of the system.

        Takes a one dimensional array of the velocities in the format [v_x1, v_y1, v_z1, v_x2, v_y2, v_z2, ...] where v_x1 represents the x component of the velocity on the first particle, etc. 

        Parameters
        -------------
        velocity : numpy.ndarray
            A one-dimentional numpy array of velocities.
        """
        return None

    @abc.abstractmethod
    def draw_velocity(self, distType='uniform'):
        """Draw a new value of the velocities.

        Redraws velocities from a specified distribtion.

        Parameters
        ------------
        distType : string, optional
            Specifies the type of distribution from which to draw the velocities. Currently supports 'uniform' and 'gaussian'.
        """
        return None

    @abc.abstractmethod
    def reverse_velocity(self, multFactor=-1.0):
        """Reverse the velocity of the walker.

        Sets the velocity to multFactor * vel.

        Parameters
        ------------
        multFactor : float
            Factor to scale the velocities. Takes -1.0 as default.
        """
        return None

    @abc.abstractmethod
    def equilibrate(self, center, restraint, numSteps):
        """
        """
        return None

    @abc.abstractmethod
    def get_colvars(self):
        """Return the location of the walker in the collective variable space.

        Returns
        ---------
        numpy.ndarray
            A one-dimensional numpy array of the current collective variables.
        """
        return None

    @abc.abstractmethod
    def add_colvars(self):
        """
        """
        return None

    @abc.abstractmethod
    def propagate(self, nSteps):
        """Integrate the dynamics of the model forward in time.

        Parameters
        -----------
        nSteps : int
            The number of time steps to integrate forward in time.
        """
        return None

    @abc.abstractmethod
    def close(self):
        """Destroy the walker.
        """
        return None

    @abc.abstractmethod
    def set_temperature(self, temp):
        """Set the temperature of the system.
        """
        return None

    @abc.abstractmethod
    def set_timestep(self, timestep):
        """Set the number of timesteps.
        """
        return None

    @abc.abstractmethod
    def get_time(self):
        """Return the time in number of model timesteps.
        """
        return None

    @abc.abstractmethod
    def set_time(self, time):
        """Set the time in number of model timesteps.
        """
        return None
