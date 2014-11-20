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
"""
Potentially create a class here called "simulation" which will contain an array
of partitions (just starting with one) and an array of "walkers" (in serial
this will be just one). The code will then pair walkers (i.e. LAMMPS instances)
with windows contained in the partition. 
"""
class partition:
    """
    This class defines a set of windows that partition the sampling space.
    This class will contain an array of basisFunction objects. 
    """
    umbrellas = []
    
    def __init__(self, N):
        self.F = self.initializeFMat(N)
        self.F_error = np.zeros(self.F.shape)
        

    def computeFij(self, i, j):
        """
        This function computes and returns the value of the entry Fij according 
        to summing over the value of the jth basis function. (The overlap 
        version)
        """
        cv_sum = 0.0
        for k in range(0, len(self.umbrellas[i].samples), 1):
            cv_sum += self.umbrellas[j](self.umbrellas[i].samples[k], self.umbrellas)
            
        cv_sum /= len(self.umbrellas[i].samples)
        
        return cv_sum
    
    def accumulateFij(self, i, cvs):
        """
        This function will update the current value of Fij with the given cvs.
        Note that this function expects cvs to be a ctypes array that can be fed
        to umbrellas()
        """
        self.F[i,:] = (self.F[i,:] * (self.umbrellas[i].numSamples-1) + cvs) / (self.umbrellas[i].numSamples)

        """
        self.F[i][j] = (self.F[i][j] * self.umbrellas[i].numSamples + self.umbrellas[j](cvs, self.umbrellas)) / (self.umbrellas[i].numSamples + 1)
        """
        
        return 0
        
    def computeAcor(self, traj):
        """
        This function will compute the autocorrelation time and the std dev of the
        given data set. 
        """
        tau, mean, sigma = acor.acor(traj)
        
        return tau, mean, sigma

    def getBasisFunctionValues(self, coord):
        """
        This function takes a point in collective variable space and returns 
        an array of the value of the basis functions.
        """
        # build an array 
        indicators = np.zeros(len(self.umbrellas))
        # go through the basis functions and construct the indicators array
        for i in range(len(self.umbrellas)):
            if self.umbrellas[i].indicator == 0.0:
                continue
            else:
                indicators[i] = self.umbrellas[i].indicator(coord)
                
        # normalize the values of the basis functions
        indicators = indicators / np.sum(indicators)
        
        return indicators
    
    def initializeFMat(self, ncells):
        """
        This function constructs the matrix F.
        """
        F = np.zeros((ncells, ncells))
        
        return F   
    
    def sample(self, wlkr, numSteps, umbrellaIndex, walkerIndex, sysParams, rank=-1):
        """    
        This function takes a system lmp and propagates the dynamics to generate
        the required samples. 
        
        """
        # if we are passing in a value of the rank, we know it is MPI and report the rank
        if rank <= 0:
            print "Rank", rank            
            
        print "Sampling umbrella " + str(umbrellaIndex) + " walker " + str(walkerIndex) + " for " + str(numSteps) + " steps."
        # print "Umbrella: ", self.umbrellas[umbrellaIndex]
        
        # This block of code is for the implementation of the abstracted walker
        # assign an input filename for this walker. 
        inputFilename = sysParams['scratchdir'] + "/in_" + str(umbrellaIndex) + "_w" + str(walkerIndex) + ".diala"
        
        
        # assign an output file for the LAMMPS trajectory:
        wlkr.command("dump 1 all dcd 10 " + inputFilename + ".dcd")
        
        # if we are using a transition type, we cannot directly reconstruct F 
        # from post-processing only the cvs's trajectory so we create a data
        # structure here to store values of this transition
        self.umbrellas[umbrellaIndex].basisFnxTimeSeries = []
        
        
        """
        #print "Initializing walker."
        # now construct the type of walker needed 
        if sysParams['walkerType'] == "lammps":     
            # import the lammps walker module
            import lammpsWalker
            # build the walker 
            wlkr = lammpsWalker.lammpsWalker(inputFilename, sysParams)
        elif sysParams['walkerType'] == "fromPyFile":
            # to import a module from a data directory, we will use the importlib module
            import importlib
            # print sysParams['datadir']
            # enforce this directory at the head of the pythonpath 
            sys.path.insert(0,sysParams['datadir'])
            # ADD A CHECK TO SEE IF THE FILE HAS AN EXTENSION!
            # ALSO, CHECK IF THIS IS CONGRUENT WITH HOW THE DATA FILE IS USED FOR LAMMPS.
            # return an import module for the pyfile module. 
            walkermodule = importlib.import_module(sysParams['dataFile'])
            # create the walker from the custom directory
            wlkr = walkermodule.getWalker()
        """   
            
        # allocate an variable to compute the metropolis accenptance percentange
        self.umbrellas[umbrellaIndex].metropolis_acceptance = 0.0
        
        # get initial points 
        # note that configs[0] refers to Y_t, while configs[1] refers to Y_t+1 proposal 
        self.umbrellas[umbrellaIndex].configs.append(wlkr.getConfig())
        self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())
        # note to increment numSamples as a record keeping for the value of N
        self.umbrellas[umbrellaIndex].numSamples += 1
        # append the starting values of the basis functions 
        self.umbrellas[umbrellaIndex].basisFnxTimeSeries.append(self.getBasisFunctionValues(self.umbrellas[umbrellaIndex].samples[-1]))
        
        #print "Propagaing sampling."
        # now propagate the dynamics of the walker
        # attempt dynamics routine
        """
        try:         
            for i in range(0, numSteps, sysParams['stepLength']):
                # propagate the dynamics
                if i == 0:
                    wlkr.propagate(sysParams['stepLength'], pre='yes')
                else:
                    wlkr.propagate(sysParams['stepLength'])
                
                # shift values and get new configs
                self.umbrellas[umbrellaIndex].configs.append(wlkr.getConfig())
                self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())
                self.umbrellas[umbrellaIndex].numSamples += 1
                
                # now upddate the F matrix with the new sample
                if sysParams['Ftype'] == "transition":
                    self.umbrellas[umbrellaIndex].basisFnxTimeSeries.append(self.getBasisFunctionValues(self.umbrellas[umbrellaIndex].samples[-1]))
                    self.accumulateFij(umbrellaIndex, self.umbrellas[umbrellaIndex].basisFnxTimeSeries[-1])
                
                # enforce metropolisation
                self.metropolisMove(self.umbrellas[umbrellaIndex])
                
                if sysParams['Ftype'] == "overlap":
                    self.umbrellas[umbrellaIndex].basisFnxTimeSeries.append(self.getBasisFunctionValues(self.umbrellas[umbrellaIndex].samples[-1]))                
                    self.accumulateFij(umbrellaIndex, self.umbrellas[umbrellaIndex].basisFnxTimeSeries[-1])
                    
                # now check to see if the data buffer has become too large and flush buffer to file
                if not (i % (1000*sysParams['stepLength'])):
                    self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)
        
        except Exception:
            raise errors.DynamicsError(umbrellaIndex)
        """
        for i in range(0, numSteps, sysParams['stepLength']):
                # propagate the dynamics
            if i == 0:
                wlkr.propagate(sysParams['stepLength'], pre='yes')
                #wlkr.propagate(sysParams['stepLength'])
            else:
                wlkr.propagate(sysParams['stepLength'])
                
                # shift values and get new configs
            self.umbrellas[umbrellaIndex].configs.append(wlkr.getConfig())
            self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())
            self.umbrellas[umbrellaIndex].numSamples += 1
                
                # now upddate the F matrix with the new sample
            if sysParams['Ftype'] == "transition":
                self.umbrellas[umbrellaIndex].basisFnxTimeSeries.append(self.getBasisFunctionValues(self.umbrellas[umbrellaIndex].samples[-1]))
                self.accumulateFij(umbrellaIndex, self.umbrellas[umbrellaIndex].basisFnxTimeSeries[-1])
                
                # enforce metropolisation
            self.metropolisMove(self.umbrellas[umbrellaIndex], wlkr)
                
            if sysParams['Ftype'] == "overlap":
                self.umbrellas[umbrellaIndex].basisFnxTimeSeries.append(self.getBasisFunctionValues(self.umbrellas[umbrellaIndex].samples[-1]))                
                self.accumulateFij(umbrellaIndex, self.umbrellas[umbrellaIndex].basisFnxTimeSeries[-1])
                    
                # now check to see if the data buffer has become too large and flush buffer to file
            if not (i % (1000*sysParams['stepLength'])):
                self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)

        # flush the last data to file
        self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)
        
        # compute the metropolis acceptance if it was computed             
        if hasattr(self.umbrellas[umbrellaIndex], 'metropolis_acceptance'):
            self.umbrellas[umbrellaIndex].metropolis_acceptance = self.umbrellas[umbrellaIndex].metropolis_acceptance / self.umbrellas[umbrellaIndex].numSamples
            #print "Acceptance ratio: ", self.umbrellas[umbrellaIndex].metropolis_acceptance
        # now clean up this walker 
        wlkr.close()             
        
        return 0

    def metropolisMove(self, current_umb, wlkr): 
        """
        This function returns True if the proposed move is accepted and false if 
        rejected based on a metropolization rule for the current Index.
        """
        if current_umb.indicator(current_umb.samples[-2]) == 0.0:
            print "An error occured. umb(umb.samples[-2]) == 0.0."
            sys.exit(0)
        # probability of accepting new point.
        prob = min(1.0, current_umb(current_umb.samples[-1], self.umbrellas) / current_umb(current_umb.samples[-2], self.umbrellas) )
        # evaluate probability
        if random.random() <= prob:
            current_umb.metropolis_acceptance += 1.0
            return True
            
        else:
            
            # print "Simulation has exited, restarting with stored coordinates."
            # scatter old config to the walker
            wlkr.setConfig(current_umb.configs[-1])
            
            # set the "rejected" step to the previous config
            current_umb.configs[-1] = current_umb.configs[-2]
                    
            # now set the "rejected" step to the current CV point
            current_umb.samples[-1] = current_umb.samples[-2]   
                        
            # print "Simulation has exited, restarting with stored coordinates."
            # scatter old config to the walker
            wlkr.setConfig(current_umb.configs[-1])
                    
            # redraw the velocities
            if issubclass(type(wlkr),walker.velocityWalker):
                wlkr.drawVel()
        
            return False
        
        # if something went wrong, raise an error 
        raise Exception("Something went wrong in the Metropolis Move Section. We don't know what.")

    def getZ(self):
        """
        This function takes the matrix F and returns the eigenvector with 
        eigenvalue 1.
        
        """
        # compute via numpy interface to LAPACK the eigenvectors v and eigenvalues w
        # The linalg routine returns this as the first (i.e. largest) eigenvalue. 
        evals, evec = LA.eig(self.F, left=True, right=False)
        sort = np.argsort(evals)
        # normalize if needed. 
        self.z = evec[:,sort[-1]] / np.sum(evec[:,sort[-1]]) 
        
        return self.z 
        
    def getlogBound(self, transmatrix):
        """
        This is controller code that takes in a transition matrix, and returns the 
        logarithm of the bound 1/P_i[tj<ti].  It also checks for pathologically 
        small entries.
        - written by Erik, incorporated by Jeremy (date 5/16/2013)
        """
        pjimatrix = sp.zeros(transmatrix.shape)
        for i in range(transmatrix.shape[0]):
            for j in range(transmatrix.shape[1]):
                if i != j:
                    pjimatrix[i,j] = self.computeLogPji(transmatrix, j, i)
                    if pjimatrix[i,j] < -100000000:
                        #This checks for pathologically small entries 
                        print "Pathologically small entry suspected at i,j= %d, %d, with value %.9f" % (i, j, pjimatrix[i,j])
                            
        return -1.0 * pjimatrix
        
    def computeLogPji(self, transmatrix, j, i):
        """
        This code calculates the logarithm of P_i[tj<ti] (Pji for short)
        for a specific i and j, according to Brian's paper on Dropbox (as of 
        writing this comment, it is on page 5).  
        
        To do this, it first calculates the matrix (I-F_{j}^{i}) which we
        denote here as A, as well as the column vector F_{notj,j}.  The vector p_m 
        is the solution to the equation A*p_m=F_{notj,j}), and the i'th element of 
        this p_m is P_i[tj<ti].  From a linear algebra standpoint, the equation is 
        solved using Cramer's rule, combined with some log math.  
        
        ######Algorithm notes:######
        This routine uses logarithms to calculate Pji, in case of an 
        over/under/whateverflow error.  This is generally safer for large 
        senvitivities.
        
        There is no good reason to use Cramer's rule here, apart from the fact that
        Erik doesn't trust python's Linear Algebra solver.  This will hopefully be 
        replaced by a better method, once Erik gets off of his lazy ass and 
        rewrites all this so that it handles arbitrary precision.
        -- written by Erik, incorporated by Jeremy (date 5/16/2013)
        """
        
        numcolumns = len(transmatrix[0])
        
        # We form F_{j}^{i}, by taking the transition matrix, setting the i'th row 
        # equal to zero, and then taking the j'th principal minor.
        
        Fsupi = sp.copy(transmatrix)
        Fsupi[:,i] = 0
        Fnotjj = sp.delete(transmatrix[:,j],j)
        Fsupisubj = sp.delete(sp.delete(Fsupi,j,0),j,1)
        newi = i
        
        # Since we have deleted a column and row from our system, all the indices 
        # above j have been shifted by one.  Therefore, if i>j, we need change i 
        # correspondingly.
        if i > j:
            newi -= 1
            
        # We define the matrix A=(I-F_{j}^{i})
        Amatrix = sp.eye(numcolumns-1)-Fsupisubj
        
        ### We start calculating p_ji using Cramer's Rule. 
        # We calculate the top matrix in Cramer's Rule
        Aswapped = sp.copy(Amatrix)
        Aswapped[:,newi] = Fnotjj[:]
        
        # To take the determinants of the two matrices, we use LU decomposition.
        Pdenom, Ldenom, Udenom = sp.linalg.lu(Amatrix)
        Pnumer,Lnumer,Unumer = sp.linalg.lu(Aswapped)
        # The determinant is just the product of the elements on the diagonal.
        # Since Pji is guaranteed positive, we can also just take the absolute 
        # value of all the diagonal elements, and sum their logarithms.
        numerdiag = abs(sp.diagonal(Unumer))
        denomdiag = abs(sp.diagonal(Udenom))
        lognumdiag = sp.log(numerdiag)
        logdenomdiag = sp.log(denomdiag)
        logpji = sp.sum(lognumdiag) - sp.sum(logdenomdiag)
        
        # Note that this is returns in ln (log base e), not log base 10.
        return logpji
        
        
    def createUmbrellas(self, umbs):
        """
        This function creates a gridding of umbrellas based on a collective 
        variable specification passed as arguments. 
        
        OK, there are definitely better ways of doing this but I don't want to 
        implement them now. For now, I support just two dimentions and boxes. 
        """
        umbrellas = []
        
        # the stride is a list of the spacing of the centers of each box in each axis
        stride = [(umbs[0][1] - umbs[0][0]) / umbs[0][2], (umbs[1][1] - umbs[1][0]) / umbs[1][2]]
        
        for i in range(umbs[0][2]): 
            for j in range(umbs[1][2]):
                center = [umbs[0][0] + stride[0]*i + stride[0] / 2.0, umbs[1][0] + stride[1]*j + stride[1] / 2.0]
                width = [umbs[0][3], umbs[1][3]]
                umbrellas.append(Box(center, width))
        
        return umbrellas
     
    

# packages needed for the basisfunction objects
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
        # get the last sample in the list
        lastSample = self.samples.pop()
        lastBasisFunction = self.basisFnxTimeSeries.pop()
        lastConfig = self.configs.pop() 
        
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
            
        with open(filename + ".colvars", "a") as f_handle:
            np.savetxt(f_handle, self.samples)

        # delete current reference to the list
        del self.samples
        del self.basisFnxTimeSeries
        del self.configs
            
        #now replace the data structure with the endpoint
        self.samples = [lastSample]
        self.basisFnxTimeSeries = [lastBasisFunction]
        self.configs = [lastConfig]

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
            for i in range(0, len(umb), 1):
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
        
