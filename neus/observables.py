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
        self.weights = np.zeros(data.shape)

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
        self.weights = np.zeros(data.shape)
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
        self.weights = np.zeros(data.shape)
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

class dihedral_fluctuation_correlation:
    """
    This routine returns the TCF of a dihedral angle. 
    """
    def __init__(self, name, s, stepLength, data, cvindex, mean):
        """
        Constructor for the time fluctuation correlation function of a dihedral.
        
        NOTE: Step Length is the time between samples taken, not necessarily the time the walker advances each timestep. 
        """
        self.s = s
        self.stepLength = stepLength
        self.name = name
        self.data = data
        self.nsamples = np.zeros(data.shape)
        self.weights = np.zeros(data.shape)
        self.mean = mean
        self.cvindex = cvindex

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
        
        time_indx = (wlkr.simulationTime - np.floor(wlkr.simulationTime / self.s) * self.s ) / (self.s / self.data.shape[0])
        
        time_indx = int(time_indx)
        
        # current colvars
        cv = wlkr.getColvars()
        d1 = cv[self.cvindex]

        config = wlkr.getConfig()

        # now let's get the reference dihedral
        wlkr.setConfig(wlkr.Y_s[0])
        wlkr.propagate(0, pre='yes')
        cv = wlkr.getColvars()

        d2 = cv[self.cvindex]

        wlkr.setConfig(config)

        wlkr.propagate(0, pre='yes')
        
        # now get the fluctuation correlation value 
        temp_val = np.dot(d1 - self.mean, d2 - self.mean)

        self.data[time_indx] = (self.data[time_indx] * self.nsamples[time_indx] + temp_val) / (self.nsamples[time_indx] + 1.0)
        self.nsamples[time_indx] += 1.0

        return 0

class dihedral_fluctuation_correlation_2:
    """
    This routine returns the TCF of a dihedral angle. 
    """
    def __init__(self, name, s, stepLength, data, atomids, mean):
        """
        Constructor for the time fluctuation correlation function of a dihedral.
        
        NOTE: Step Length is the time between samples taken, not necessarily the time the walker advances each timestep. 
        """
        self.s = s
        self.stepLength = stepLength
        self.name = name
        self.data = data
        self.nsamples = np.zeros(data.shape)
        self.weights = np.zeros(data.shape)
        self.mean = mean
        self.atomids = atomids

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
        
        time_indx = (wlkr.simulationTime - np.floor(wlkr.simulationTime / self.s) * self.s ) / (self.s / self.data.shape[0])
        
        time_indx = int(time_indx)
        
        # current colvars
        config = wlkr.getConfig()

        d1 = self.__get_dihedral__(config, self.atomids)

        ref_config = wlkr.Y_s[0]

        d2 = self.__get_dihedral__(ref_config, self.atomids)
        
        # now get the fluctuation correlation value 
        temp_val = np.dot(d1 - self.mean, d2 - self.mean)

        self.data[time_indx] = (self.data[time_indx] * self.nsamples[time_indx] + temp_val) / (self.nsamples[time_indx] + 1.0)
        self.nsamples[time_indx] += 1.0

        return 0

    def __get_dihedral__(self, config, atomids):
        """
        A routine to compute a dihedral angle from a molecular configuration.

        *** THIS IMPLEMENTATION IS NOT TESTED ***
        """

        xyz = np.reshape(config, (config.size/3, 3))

        p1 = xyz[atomids[0] - 1]
        p2 = xyz[atomids[1] - 1]
        p3 = xyz[atomids[2] - 1]
        p4 = xyz[atomids[3] - 1]

        b1 = p2 - p1
        b2 = p3 - p2
        b3 = p4 - p3

        n1 = np.cross(b1, b2)
        n1 /= np.linalg.norm(n1)

        n2 = np.cross(b2, b3)
        n2 /= np.linalg.norm(n2)

        m1 = np.cross(n1, b2 / np.linalg.norm(b2))

        x = np.dot(n1, n2)
        y = np.dot(m1, n2)

        dihed = np.arctan2(x, y) * -180.0 / np.pi

        return dihed


class electric_field:
    """
    This class constructs an observable that reports on the time-correlation 
    of the electric field at a point in space. 
    """
    
    def __init__(self, name, s, stepLength, atomids, data, cellDim, atom_exclusions, ref_atoms_ids = None):
        """
        Initializes the electric field observable. 
        
        Parameters:
        ------------------
        name - 
            name of the observable instance.
        
        s - 
            correlation time of the observable to compute up to.
            
        stepLength - 
            the step length of the dynamics propagated in the walker. 
            
        atomids - 
            the ids of the atoms on which to compute the electric field.
            
        data - 
            a refernce data structure used to store the values computed by the call
            
        cellDim - 
            dimensions of the periodic box
            
        atom_excusions - 
            a list of atom ids to exclude from the electric field calculation
        
        atom_charges - 
            a list of fixed atomic charges 
            
        ref_atoms_ids -
            a list of atom ids used to compute the unit vector along which the 
            electric field is computed. This is here to facilitate computation of 
            electric fields along a bond-axis. 
        
        """
        self.s = s
        self.stepLength = stepLength
        self.atomids = atomids
        self.name = name
        self.data = data
        self.nsamples = np.zeros(data.shape)
        self.weights = np.zeros(data.shape)
        self.cellDim = cellDim
        self.atom_exclusions = atom_exclusions
        self.ref_atoms_ids = ref_atoms_ids
        # note the units here are in atomic units 
        self.coulomb_constant = 1 / (4 * np.pi)
        
    def __call__(self, wlkr):
        """
        Updates internal data structure of the time-correlation function for the 
        electric field at a point.
        
        Reminder, this sums up the instantaneous force from atoms i in the system 
        by E = \sum_{i} k Q_i / r**2 where k = 1 / (4 * \pi * \eta_{0}) 
        """

        # btw, you can get the charges from the simulation from lmp.gather_atoms("q", 1, 1)[:]
        
        # first check we should update a point here. 
        
        if not (wlkr.simulationTime - np.floor(wlkr.simulationTime / self.s) * self.s ) % (self.s / self.data.shape[0]) == 0.0:
            return 0 
        
        time_indx = (wlkr.simulationTime - np.floor(wlkr.simulationTime / self.s) * self.s ) / (self.s / self.data.shape[0]) - 1
        
        time_indx = int(time_indx)

        self.atom_charges = wlkr.lmp.gather_atoms("q", 1, 1)[:]
        
        # current configuration
        config = wlkr.getConfig()
        config = np.reshape(config, (-1,3))
        
        efield = self.__get_efield__(config)
        
        # get reference configuration
        config = wlkr.Y_s[0]
        config = np.reshape(config, (-1, 3))
        
        ref_efield = self.__get_efield__(config)

        temp_val = efield * ref_efield
        
        self.data[time_indx] = (self.data[time_indx] * self.nsamples[time_indx] + temp_val) / (self.nsamples[time_indx] + 1.0)
        self.nsamples[time_indx] += 1.0
        
        # now get the exclusion
        
        
        return 0 
        
    def __min_image__(self, vec):
        """
        Applies minimum image to given vector.
        """        
        for dim in range(3):
            if vec[dim] > self.cellDim[dim] / 2.0: 
                vec[dim] -= self.cellDim[dim]
            elif vec[dim] < -self.cellDim[dim] / 2.0: 
                vec[dim] += self.cellDim[dim]
        
        return vec 
        
        
    def __get_efield__(self, config):
        """
        Returns the value of the electric field for the given configuration. 
        """
        point = config[self.atomids[0]-1]
        
        # now get the unit vector in the direction of the reference electric field
        unit_vec = config[self.ref_atoms_ids[1]-1] - config[self.ref_atoms_ids[0]-1]
        
        #apply minimum distance to the unit vector
        unit_vec = self.__min_image__(unit_vec)
        
        # normalize it 
        unit_vec /= np.linalg.norm(unit_vec)
        
        temp_val = 0.0        
        
        # now loop over the particles and add up the columbic force
        for atom_index, atom in enumerate(config):
            if (atom_index+1) in self.atom_exclusions: 
                continue
            
            # get the inter nuclear separation
            d = atom - point
            
            # apply minimum image to d 
            d = self.__min_image__(d)
                
            # now get the force vector multiplying the unit vector of d against the magnitude of the force 
            force_vec = (d / np.linalg.norm(d)) * self.coulomb_constant * self.atom_charges[atom_index] / np.linalg.norm(d)**2
            # now get the value along the bond axis
            temp_val += np.dot(force_vec, unit_vec)
            
        return temp_val
        
        
class cv_indicator_correlation:
    """
    This routine returns the value of an indicator function over a space in a given cv. 
    """
    def __init__(self, name, s, stepLength, data, cv_index, cv_range):
        """
        Constructor for the time fluctuation correlation function of a dihedral.
        
        NOTE: Step Length is the time between samples taken, not necessarily the time the walker advances each timestep. 
        """
        self.s = s
        self.stepLength = stepLength
        self.name = name
        self.data = data
        self.nsamples = np.zeros(data.shape)
        self.weights = np.zeros(data.shape)
        self.cv_index = cv_index
        self.cv_range = cv_range

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
        
        time_indx = (wlkr.simulationTime - np.floor(wlkr.simulationTime / self.s) * self.s ) / (self.s / self.data.shape[0])
        
        time_indx = int(time_indx)
        
        # current colvars
        cv = wlkr.getColvars()
        d1 = cv[self.cv_index]

        if self.cv_range[0] < d1 < self.cv_range[1]:
            d1 = 1.0
        else:
            d1 = 0.0

        config = wlkr.getConfig()

        # now let's get the reference colvar
        wlkr.setConfig(wlkr.Y_s[0])
        wlkr.propagate(0, pre='yes')
        cv = wlkr.getColvars()

        d2 = cv[self.cv_index]

        if self.cv_range[0] < d2 < self.cv_range[1]:
            d2 = 1.0
        else:
            d2 = 0.0

        wlkr.setConfig(config)

        wlkr.propagate(0, pre='yes')
        
        # now get the fluctuation correlation value 
        temp_val = d1 * d2

        self.data[time_indx] = (self.data[time_indx] * self.nsamples[time_indx] + temp_val) / (self.nsamples[time_indx] + 1.0)
        self.nsamples[time_indx] += 1.0

        return 0