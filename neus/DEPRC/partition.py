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
        
        # this matrix will store transition counts 
        self.trans = self.initializeFMat(N)
        
        # this will be a list that stores the total time spent sampling each
        # basis window
        self.samplingTimes = [0.0]*N
        

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
        
        # we will assign a time-index associated with this umbrella
        self.umbrellas[umbrellaIndex].time = 0.0
        
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

            # now we check to see if we've made an exit and reset if so 
            if self.umbrellas[umbrellaIndex].indicator(self.umbrellas[umbrellaIndex].samples[-1]) == 0.0:
                
                for ind in range(len(self.umbrellas)):
                    if self.umbrellas[ind].indicator(self.umbrellas[umbrellaIndex].samples[-1]) == 1.0:
                        #print "added entry from " + str(umbrellaIndex) + " to " + str(ind)
                        self.umbrellas[ind].entryPoints.append(np.asarray(self.umbrellas[umbrellaIndex].configs[-1][:]))
                        break
                        
                # print "Simulation has exited, restarting with stored coordinates."
                # scatter old config to the walker
                wlkr.setConfig(self.umbrellas[umbrellaIndex].configs[-1])
            
                # set the "rejected" step to the previous config
                self.umbrellas[umbrellaIndex].configs[-1] = self.umbrellas[umbrellaIndex].configs[-2]
                    
                # now set the "rejected" step to the current CV point
                self.umbrellas[umbrellaIndex].samples[-1] = self.umbrellas[umbrellaIndex].samples[-2]   
                        
                # print "Simulation has exited, restarting with stored coordinates."
                # scatter old config to the walker
                wlkr.setConfig(self.umbrellas[umbrellaIndex].configs[-1])
                    
                # redraw the velocities
                if issubclass(type(wlkr),walker.velocityWalker):
                    wlkr.drawVel()
                
                # now check to see if the data buffer has become too large and flush buffer to file
            if not (i % (1000*sysParams['stepLength'])):
                self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)

        # flush the last data to file
        self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)
        
        # compute the metropolis acceptance if it was computed             
        if hasattr(self.umbrellas[umbrellaIndex], 'metropolis_acceptance'):
            self.umbrellas[umbrellaIndex].metropolis_acceptance = self.umbrellas[umbrellaIndex].metropolis_acceptance / self.umbrellas[umbrellaIndex].numSamples
            #print "Acceptance ratio: ", self.umbrellas[umbrellaIndex].metropolis_acceptance
        
        return 0
        
    def sampleNeus(self, wlkr, numSteps, umbrellaIndex, walkerIndex, sysParams, rank=-1):
        """    
        This function takes a system lmp and propagates the dynamics to generate
        the required samples but storing and generating samples via NEUS. 
        
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
        
        # get initial points 
        # note that configs[0] refers to Y_t, while configs[1] refers to Y_t+1 proposal 
        #self.umbrellas[umbrellaIndex].configs.append(wlkr.getConfig())
        #self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())
        
        # note to increment numSamples as a record keeping for the value of N
        #self.umbrellas[umbrellaIndex].numSamples += 1
        
        for i in range(0, numSteps, sysParams['stepLength']):
                # propagate the dynamics
            if i == 0:
                wlkr.propagate(sysParams['stepLength'], pre='yes')
                #wlkr.propagate(sysParams['stepLength'])
            else:
                wlkr.propagate(sysParams['stepLength'])
                
                # shift values and get new configs
            #self.umbrellas[umbrellaIndex].configs.append(wlkr.getConfig())
            self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())
            self.umbrellas[umbrellaIndex].numSamples += 1
                
            # check for transition
            if self.umbrellas[umbrellaIndex].indicator(self.umbrellas[umbrellaIndex].samples[-1]) == 0.0:
                # find to which box the transition was made 
                for indx in range(len(self.umbrellas)):
                    
                    if self.umbrellas[indx].indicator(self.umbrellas[umbrellaIndex].samples[-1]) == 1.0:        
                        
                        # record a transition from the matrix
                        self.trans[umbrellaIndex][indx] += 1.0
                        
                        # append the entry point to the new window
                        self.umbrellas[indx].entryPoints.append(wlkr.getConfig())
                        
                        # drop the last point from the samples 
                        self.umbrellas[umbrellaIndex].samples.pop()
                        
                        # now reset the dynamics from an entry point
                        temp_indx = random.randint(0, len(self.umbrellas[umbrellaIndex].entryPoints)-1)
                        wlkr.setConfig(self.umbrellas[umbrellaIndex].entryPoints[temp_indx])
                        
                        # redraw velocities at this point
                        wlkr.drawVel(distType = 'gaussian', temperature = 310.0)
                        # we're done so let's stop the loop
                        break
                
                    
            # now check to see if the data buffer has become too large and flush buffer to file
            if not (i % (1000*sysParams['stepLength'])):
                self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)

        # flush the last data to file
        self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)
        
        # update the F matrix with the transitions 
        self.F[umbrellaIndex,:] = (self.F[umbrellaIndex,:] * self.samplingTimes[umbrellaIndex] + self.trans[umbrellaIndex,:]) / (self.samplingTimes[umbrellaIndex] + sysParams['timestep'] * sysParams['stepLength'] * self.umbrellas[umbrellaIndex].numSamples)
        
        # update the time spent sampling this window
        self.samplingTimes[umbrellaIndex] += sysParams['timestep'] * sysParams['stepLength'] * self.umbrellas[umbrellaIndex].numSamples    

        return 0

    def metropolisMove(self, current_umb, wlkr): 
        """
        This function returns True if the proposed move is accepted and false if 
        rejected based on a metropolization rule for the current Index.
        
        Here we add 
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
    
    def buildNeighborList(self):
        """
        This routine constructs a neighborlist based on the radius of the basis 
        function. 
        """
        
        for i in range(len(self.umbrellas)):
            for j in range(i+1, len(self.umbrellas), 1):
                dist = np.linalg.norm(self.umbrellas[i].center - self.umbrellas[j].center)
                if dist < (self.umbrellas[i].radius + self.umbrellas[j].radius):
                    self.umbrellas[i].neighborList.append(j)
                    self.umbrellas[j].neighborList.append(i)
    
            
        return 0
    
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
     
    