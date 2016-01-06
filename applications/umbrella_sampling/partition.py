"""
This file contains the partition class definition that allows 


Potentially create a class here called "simulation" which will contain an array
of partitions (just starting with one) and an array of "walkers" (in serial
this will be just one). The code will then pair walkers (i.e. LAMMPS instances)
with windows contained in the partition. 
"""
from scipy import linalg as LA
import scipy as sp
import fileIO
import sys
import walker
#import acor
import errors
import random
import h5py
import numpy as np

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


    def sample_US(self, wlkr, numSteps, umbrellaIndex, walkerIndex, sysParams, debug=False):
        """
        This routine samples the given box via the equilibrium overlap method.
        """
        assert sysParams.has_key('scratchdir'), "Scratch directory was not specified in the sampling routine."

        assert sysParams['transitionMatrixType'] in ['transition','overlap']
        # assign an input filename for this walker.
        if debug:
            inputFilename = sysParams['scratchdir'] + "/" + str(umbrellaIndex) + "_w" + str(walkerIndex)
        else:
            inputFilename = None

        oldConfig = wlkr.getConfig()
        oldSample = wlkr.getColvars()

        f_handle = h5py.File(sysParams['scratchdir'] + "/ep." + str(umbrellaIndex) + ".h5py", "a")

        #print self.umbrellas[umbrellaIndex](oldSample, self.umbrellas), self.umbrellas[umbrellaIndex].indicator(oldSample)
        assert self.umbrellas[umbrellaIndex].indicator(oldSample) > 0.0, "The walker is not in the support of the current window."

        # get the sample from the initial state of the walker in CV space
        self.umbrellas[umbrellaIndex].samples.append(oldSample)

        # reset the local observables arrays for this round of sampling
        for obs in self.umbrellas[umbrellaIndex].local_observables:
            self.resetObservable(obs)

        assert sysParams.has_key('stepLength'), "StepLength was not specified in the sampling routine."

        self.accumulateObservables(wlkr, wlkr.getColvars(), wlkr.colvars, umbrellaIndex)

        # now we proceed with the sampling routine
        for i in range(0, numSteps, sysParams['stepLength']):

            # propagate the dynamics
            wlkr.propagate(sysParams['stepLength'])

            newConfig = wlkr.getConfig()
            newSample = wlkr.getColvars()

            if sysParams['transitionMatrixType'] == 'transition':
                # update the M matrix based on this sample
                self.M[umbrellaIndex,:] += self.getBasisFunctionValues(newSample, umbrellaIndex)
                # increment the number of samples
                self.nsamples_M[umbrellaIndex] += 1

            # get the new configuration
            if self.metropolisMove(self.umbrellas[umbrellaIndex], oldSample, newSample):
                # get the new sample position and append it to the samples list
                self.umbrellas[umbrellaIndex].samples.append(newSample)

                # set current position to "old position"
                oldConfig = newConfig
                oldSample = newSample

                # append the new sample to the observables
                self.accumulateObservables(wlkr, newSample, wlkr.colvars, umbrellaIndex)

            else:
                # if we reject the proposed position, append the old config
                self.umbrellas[umbrellaIndex].samples.append(oldSample)

                # set the walker configuration to the old state
                wlkr.setConfig(oldConfig)

                # redraw the velocities
                wlkr.drawVel(distType = 'gaussian', temperature = 310.0)

                # append the new sample to the observables
                self.accumulateObservables(wlkr, oldSample, wlkr.colvars, umbrellaIndex)

            # now check to see if the data buffer has become too large and flush buffer to file
            #if len(self.umbrellas[umbrellaIndex].samples) > 1000000: self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)

            if sysParams['transitionMatrixType'] == 'overlap':
                self.M[umbrellaIndex,:] += self.getBasisFunctionValues(oldSample, umbrellaIndex)
                # increment the number of samples
                self.nsamples_M[umbrellaIndex] += 1

            if i % 1000 == 0: f_handle.create_dataset(str(self.k) + "." + str(i) + ".config", data=wlkr.getConfig())


        # flush the last data to file after sampling has finished
        self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)

        f_handle.flush()
        f_handle.close()

        return 0
    
