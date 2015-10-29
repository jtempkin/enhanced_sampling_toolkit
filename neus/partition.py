# -*- coding: utf-8 -*-
"""
This is an update partition module that contains the flexibility to perform both umbrella and NEUS sampling. 
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
    def __init__(self, N, scratchdir, parallel=False):
        """
        Init routine. Initializes the following values:

        The obervables list is a bookkeeping of the observables of the self one
        wants to record during the simulaiton. This list needs to be set prior to sampling.

        The list shoud contain elements that are tuples of a function and a data array.

        Each element should take two arguments, one the walker object and one the umbrella index for which the sample is
        associated.
        """

        # create an umbrella list
        self.umbrellas = []
        
        # initialize the matricies needed for the NEUS
        self.M = np.zeros((N,N))
        self.F = np.zeros((N,N))     
        self.a = np.zeros(N)
        # we should track how many samples are taken in the M matrix
        self.nsamples_M = np.zeros(N)
        self.z = np.zeros(N)

        # keep a list of any observable objects that should be estimated during sampling. 
        self.observables = []
        
        # start a list of the windows associated with this partition. 
        self.rank_window_index = None

        self.simulationTime = np.zeros(N)
        
        # here, the variable k represents how many times the transition matrix has been updated. 
        self.k = np.zeros(N)

        self.scratchdir = scratchdir

        self.active_windows = np.zeros(N, dtype=np.bool)
        
    def updateF(self, row, epsilon):
        """
        This routine update G according to:
            
            G_{ij}^{k+1} = (1 - \epsilon_{k}) * G_{ij}^{k} + \epsilon_{k} * M_{ij} / T_{i}

        for a finite time problem. 
        """

       
        temp_M = self.M[row]
        """
        # update elements
        temp_M[:] = self.M[row,:] / self.nsamples_M[row]
        temp_M[row] = 1 - temp_M.sum()
        """

        #self.F[row] = (self.k[row] * self.F[row] + temp_M) / (self.k[row] + 1)
        self.F[row] = ((1-epsilon) * self.F[row] + epsilon * temp_M)

        self.k[row] += 1

        return 0

    def addObservable(self, A, rank_index=None):
        """
        This routine adds an observable and initializes local copies of the observables in the basis windows.
        """
        assert hasattr(self, "observables"), "The partition observables lists were not properly initialized."

        # add the observable to both the partition and each window own by the ranks
        self.observables.append(A)
        
        if self.rank_window_index is None:
            for window in self.umbrellas:
                if not hasattr(window, "local_observables"): window.local_observables = []
                window.local_observables.append(copy.deepcopy(A))
        else:
            for window in self.rank_window_index:
                if not hasattr(self.umbrellas[window], "local_observables"): self.umbrellas[window].local_observables = []
                self.umbrellas[window].local_observables.append(copy.deepcopy(A))

        return 0

    def removeObservable(self):
        """
        This routine removes all observables from the list for the partition. 
        """
        # remove global version of observables
        self.observables = []
        
        # remove local version of observables from the windows as well. 
        for win in self.umbrellas:
            win.local_observables = []
            
        return 0 
        

    def computeObservables(self, rank=None):
        """
        This routine populates the partition observables with data averaged from the windows.
        """
        for o_indx,obs in enumerate(self.observables):
            temp = np.zeros(obs.data.shape)
            temp_weights = np.zeros(obs.nsamples.shape)
            
            if rank is None:
                for w_indx,win in enumerate(self.umbrellas):
                    temp += self.z[w_indx] * win.local_observables[o_indx].data
                    temp_weights += self.z[w_indx] * win.local_observables[o_indx].nsamples.astype(bool)
            else:
                for w_indx in self.rank_window_index[rank]:
                    temp += self.z[w_indx] * self.umbrellas[w_indx].local_observables[o_indx].data
                    temp_weights += self.z[w_indx] * self.umbrellas[w_indx].local_observables[o_indx].nsamples.astype(bool)


            obs.data[:] = temp[:]
            obs.weights[:] = temp_weights[:]

        return 0

    def computeZ(self, sparseSolve=True, finiteTime=False):
        """
        Solves for z vector given current G,a via solving the following linear 
        system:
            
            (I - G)^T z = a 
	
	    for a finite time process and
	
            zG = z 

	    for infinite time process.

        A = (np.identity(self.G.shape[0]) - self.G).transpose()
        
        self.z = np.linalg.solve(A, self.a)

        """
        # let's compute a z only for the list of active windows
        temp_F = self.F[self.active_windows, :][:, self.active_windows]

        if finiteTime:
            A = (np.identity(self.F.shape[0]) - self.F).transpose()
        
            self.z = np.linalg.solve(A, self.a)

            self.z /= self.z.sum()

            return 0 


        if sparseSolve:
            # compute via numpy interface to LAPACK the eigenvectors v and eigenvalues w
            # here will first convert F to sparse format before solving.
            # The linalg routine returns this as the first (i.e. largest) eigenvalue.
            
            F_sparse = sparse.coo_matrix(temp_F)
            evals, evec = LA_sparse.eigs(F_sparse.transpose())
            #evals, evec = LA.eig(self.F, left=True, right=False)
            sort = np.argsort(evals)
            # normalize if needed.

            self.z[self.active_windows] = evec[:,sort[-1]] / np.sum(evec[:,sort[-1]])
            self.z[np.logical_not(self.active_windows)] = 0.0

        else:
                
            evals, evec = LA.eig(self.F, left=True, right=False)
            sort = np.argsort(evals)
            # normalize if needed.

            #self.z[self.active_windows] = evec[:,sort[-1]] / np.sum(evec[:,sort[-1]])
            #self.z[np.logical_not(self.active_windows)] = 0.0

            self.z = evec[:,sort[-1]] / np.sum(evec[:,sort[-1]])
            


        return 0

    def accumulateObservables(self, wlkr, sample, colvars, indx):
        """
        This routine loops through the observables list and updates the samples in the corresponding windows.
        """

        for obs in self.umbrellas[indx].local_observables:
            # use the observable call routine to accumulate a sampling in the observables data structure
            if isinstance(obs, observables.pmf): 
                obs(sample, colvars)
            elif isinstance(obs, observables.P1):
                obs(wlkr)
            elif isinstance(obs, observables.dist_fluctuation_correlation):
                obs(wlkr)
            elif isinstance(obs, observables.electric_field):
                obs(wlkr)
            elif isinstance(obs, observables.dihedral_fluctuation_correlation):
                obs(wlkr)
            elif isinstance(obs, observables.dihedral_fluctuation_correlation_2):
                obs(wlkr)
            elif isinstance(obs, observables.cv_indicator_correlation):
                obs(wlkr)
            else:
                raise Exception('Called unknown observable type.')

        return 0

    def resetObservable(self, obs):
        """
        Set each element in the data of the passed observable to zero.
        """

        obs.data.fill(0.0)
        obs.nsamples.fill(0.0)

        return 0
        
    def reinject(self, wlkr, i):
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
            if self.umbrellas[i].getNumberOfEntryPoints(key=key) == 0:
                prob[indx] = 0.0
        """

        # normalize probability
        #assert prob.sum() > 0.0
        #prob /= prob.sum()

        # now choose neighbor proportional to prob 
        #I = np.random.choice(np.arange(prob.size), p=prob)

        # get the entry point from the umbrella window

        #assert self.umbrellas[i].getNumberOfEntryPoints(key=self.index_to_key[I])

        # we choose from the initial distribution with prob. stored in this window
        if random.random() < self.umbrellas[i].initial_distribution_prob:
            dist = self.umbrellas[i].entryPoints[i]
            EP = random.sample(dist, 1)[0]

        else:
            EP = self.umbrellas[i].getEntryPoint(self.index_to_key[i])

        # you should pass this argument as a ctypes array for now
        wlkr.setConfig(EP.config)
        wlkr.setVel(EP.vel)

        # set the lag component of the walker state
        wlkr.Y_s = EP.Y_s

        wlkr.simulationTime = EP.time

        return 0

    def getBasisFunctionValues(self, wlkr, umbrella_index = None):
        """
        This function takes a point in collective variable space and returns 
        an array of the value of the basis functions at that point. 
        
        If no umbrella index is passed, then search the whole space 
        """
        # build an array 
        indicators = np.zeros(len(self.umbrellas))

        coord = wlkr.getColvars()
        
        if umbrella_index is None:
    	
            # go through the basis functions and construct the indicators array
            for i in xrange(len(self.umbrellas)):
                if self.umbrellas[i].indicator(coord) == 0.0:
                    continue
                else:
                    indicators[i] = self.umbrellas[i](wlkr, self.umbrellas)
                    
        elif len(self.umbrellas[umbrella_index].neighborList) == 0: 
            # go through the basis functions and construct the indicators array
            for i in xrange(len(self.umbrellas)):
                if self.umbrellas[i](wlkr, self.umbrellas) == 0.0:
                    continue
                else:
                    indicators[i] = self.umbrellas[i](wlkr, self.umbrellas)
        else:
            
            # make sure the partition has a neighborlist
            assert hasattr(self.umbrellas[umbrella_index], 'neighborList'), "There is no neighborlist defined."
            # now loop over neighbors and populate the indicator
            for i in self.umbrellas[umbrella_index].neighborList:
                
                indicators[i] = self.umbrellas[i](wlkr, self.umbrellas)
                
            # if we don't find any support, let's try this again and search the whole space. 
            if np.sum(indicators) == 0.0:
                indicators = self.getBasisFunctionValues(wlkr, umbrella_index=None)
                
        # normalize the values of the basis functions
        if np.sum(indicators) == 0.0:
            print wlkr.getColvars()
            print wlkr.simulationTime
            print wlkr.Y_s
        assert np.sum(indicators) != 0.0, str(wlkr.getColvars())
        indicators = indicators / np.sum(indicators)
        
        return indicators
        
    def sample(self, wlkr, umbrellaIndex, numSteps, stepLength, walkerIndex, corrLength=None, debug=False):
        """    
        This function takes a system lmp and propagates the dynamics to generate
        the required samples but storing and generating samples via NEUS algorithm.
        
        Specifically, this performs an umbrella sampling routine using NEUS reinjection procedure for reinitializing the walker. 
        
        We should remove the need for the output to be specified internally here. 
        
        """
        # count the number of transitions made in this sampling routine if debug.
        if debug: ntransitions = 0
            
        assert hasattr(self, 'scratchdir'), "Scratch directory was not specified for this partition object."

        # assign an input filename for this walker.  
        if debug: 
            inputFilename = self.scratchdir + "/" + str(umbrellaIndex) + "_w" + str(walkerIndex)
        else:
            inputFilename = None

        # get the sample from the initial state of the walker in CV space
        self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())
        
        # reset the local observables arrays for this round of sampling
        for obs in self.umbrellas[umbrellaIndex].local_observables:
            self.resetObservable(obs)
        
        # now we proceed with the sampling routine            
        for i in range(0, numSteps, stepLength):
            
            # propagate the dynamics
            wlkr.propagate(stepLength)
            
            # update the record of the simulation time for the walker object. 
            wlkr.simulationTime += stepLength

            # if this is a finite time version, we should kill the walker if it hits the boundary
            if 30.0 < wlkr.getColvars()[0] < 90.0:
                #print "hit target", self.umbrellas[umbrellaIndex].center, wlkr.simulationTime, self.z[umbrellaIndex]
                self.reinject(wlkr, umbrellaIndex)

                self.umbrellas[umbrellaIndex].nhits += 1 

                # now update the statistics for hitting vs stopping. 

                continue
	    
            # now we check to see if we've passed the autocorrelation length
            # if we do, we reset the Y ( [t / s] * s) value to the current point
            if corrLength is not None:
                if (wlkr.simulationTime % corrLength) == 0.0:
                    #wlkr.Y_s = (wlkr.getConfig(), wlkr.getVel(), wlkr.getColvars())
                    #wlkr.simulationTime = 0.0
                    #print "stopping time hit"
                    # if we hit the stopping time, reset the walker
                    self.reinject(wlkr, umbrellaIndex)
                    continue


            
            # get the new sample position
            new_sample = wlkr.getColvars()
            self.umbrellas[umbrellaIndex].samples.append(new_sample)

            # check for a transition out of this index
            if self.umbrellas[umbrellaIndex](wlkr, self.umbrellas) == 0.0:
                if debug: ntransitions += 1 
                # choose the new j with probability {psi_0, ..., psi_N}

                indicators = self.getBasisFunctionValues(wlkr)
            
                # record a transition to the matrix
                self.M[umbrellaIndex,:] += indicators * stepLength
                
                # now we select a window to which to append this new entry point,use numpy to choose window index
                ep_targets = np.arange(indicators.size)[indicators.astype(bool)]
                        
                # create a new entry point and append the entry point to the new window
                newEP = entryPoints.entryPoints(wlkr.getConfig(), wlkr.getVel(), wlkr.simulationTime)
                newEP.Y_s = wlkr.Y_s

                for indx in ep_targets:
                    #if not self.active_windows[indx]: print "added entry point", indx
                    self.umbrellas[indx].addNewEntryPoint(newEP, self.index_to_key[umbrellaIndex])

                # drop the last point from the samples 
                self.umbrellas[umbrellaIndex].samples.pop()

                # reinject the walker into the current window 
                self.reinject(wlkr, umbrellaIndex)
            		
                # get the sample from the new starting point after reinjection
                self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())

            # if we do not detect a transition and handle that, we should add a count to M_ii
            else:
                self.M[umbrellaIndex, umbrellaIndex] += stepLength
            
            # let's accumulate a sample into the observables we are accumulating on this window
            self.accumulateObservables(wlkr, wlkr.getColvars(), wlkr.colvars, umbrellaIndex)

            # now check to see if the data buffer has become too large and flush buffer to file
            #if len(self.umbrellas[umbrellaIndex].samples) > 1000000: self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)
        
        # flush the last data to file after sampling has finished
        self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)

        # record the number of samples taken here
        self.nsamples_M[umbrellaIndex] = numSteps

        # now we compute the estimate of the flux from thsi iteration
        #self.M[umbrellaIndex,umbrellaIndex] = numSteps - self.M[umbrellaIndex, :].sum()

        self.M[umbrellaIndex, :] /= numSteps

        # here we will store the current position of the walker in an entry point structure
        newEP = entryPoints.entryPoints(wlkr.getConfig(), wlkr.getVel(), wlkr.simulationTime)
        newEP.Y_s = wlkr.Y_s
        self.umbrellas[umbrellaIndex].walker_restart = newEP

        if debug: print "row" ,umbrellaIndex, "of M:", self.M[umbrellaIndex,:]
    
        if debug: print "Recorded",ntransitions,"transitions"

        return 0

    def metropolisMove(self, current_umb, oldConfig, newConfig): 
        """
        This function returns True if the proposed move is accepted and false if 
        rejected based on a metropolization rule for the current Index.
        """
        assert current_umb(oldConfig, self.umbrellas) > 0.0, "The old walker configuration is not in the support of the window."
        
        # probability of accepting new point.
        prob = min(1.0, current_umb(newConfig, self.umbrellas) / current_umb(oldConfig, self.umbrellas) )
        
        # evaluate probability
        if random.random() <= prob:
        
            return True
            
        else:        
            
            return False
        
        # if something went wrong, raise an error 
        raise Exception("Something went wrong in the Metropolis Move Section. We don't know what.")
    
    def buildNeighborList(self, umbrellas, s, debug=False):
        """
        This routine constructs a neighborlist based on the radius of the basis 
        function.
        
        L is the vector specifying the periodic lengths in each dimension. 
        """      
        # we should make sure that ever basisFunction has a radius defined. 
        for win in umbrellas:
            assert hasattr(win, 'radius')
            
        # now let's find all of the neighbors and populate a list of neighbors for each 
        for i in range(len(umbrellas)):
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

    def createUmbrellas(self, colVarParams, wrapping, basisType="Box", neighborList=True, max_entryPoints = 500):
        """
        We create a grid of Umbrellas on the collective variable space.  Each collective variable is divided by evenly spaced umbrellas.
        It takes in as input:
            Params          A dictionary, which includes parameters which specify how to manipulate the collective parameters.
                                Currently, the following keywords are implemented:
                                    cvrangeN        Here, N is some integer, e.g. cvrange1, or cvrange2.
                                                        cvrange should have a value of a string of 4 scalars, separated by a comma, e.g. "-180,180,12,30"
                                                        These represent (respectively) the min and max values of the collective variable, the number of breaks in that axs, and 1/2 the width of the umbrella.
                                Theoretically, you can pass it the sysParams array in the main statement, and it should work.  However, this might be an unsafe practice.
                                It might be better to collect the collective variable params somewhere else, and just pass those as an argument.
        The function should return:
            um               A list of umbrella objects which cover the space of all the collective variables.
    
        """
        assert type(colVarParams) is np.ndarray, "The colvars array passed is not of the correct type."
        
        assert basisType in ['Box', 'Gaussian', 'Pyramid']
        
        # We check if the user provided any input to wrap the basis functions.
        if 'wrapping' is None:
            wrapping=np.zeros(len(colVarParams))
        
        # We make the following three lists:
        #       rangelist, a list where element i is a list of all the possible center values for 
        #                   collective variable i
        #       widthlist, a list where each element is the width of the Box in that dimension
        #       divisions, a list where element i is the number of divisions in collective variable i.
        centersList=[]
        widthlist=[]
        divisions=[]
        # We also create the boxwrap array.  This array will give the domain of each wrapping collective variable.
        # For example, it would contain 360 for an angle going from -180 to 180 degrees.
        boxwrap=[]
        for cvindex in xrange(len(colVarParams)):
            doeswrap = (wrapping[cvindex] != 0.0)
            v=colVarParams[cvindex]
            # We make an evenly spaced array containing all the points along a collective coordinate 
            # where we want to center a box.  We will have to do this differently idepending on whether
            # the coordinate wraps around:  the reason for this is that if the coordinate wraps around,
            # we will want overlap over the boundaries.
            if doeswrap:
                c1=(np.linspace(v[0],v[1],v[2]+1))
                centersList.append([(c1[i]+c1[i+1])/2 for i in xrange(len(c1)-1)])
                boxwrap.append(v[1]-v[0])
            else:
                isIncreasing=v[1]>v[0]
                if isIncreasing:
                    centersList.append(np.linspace(v[0]+v[3],v[1]-v[3],v[2]))
                else:
                    centersList.append(np.linspace(v[0]-v[3],v[1]+v[3],v[2]))
                boxwrap.append(0.0)
            widthlist.append(v[3])
            divisions.append(v[2])
        # We define some constants which will be useful:  the number of umbrellas numum, and the
        # cumulative product of the number of divisions for each successive collective variable, numthingie
        # (It's important for integer rounding)
        numum=np.product(divisions)
        numheirarchy=np.flipud(np.cumprod(np.flipud(divisions)))
        
        # We create the boxwrap array.  This array will give the domain of each wrapping collective variable.
        # For example, it would contain 360 for an angle going from -180 to 180 degrees.
        
        # We make um, the list of all our partition function objects.
        um=[]
        for i in xrange(int(numum)):
            umbrellaCoord=[]
            # We loop through the collective coordinates.
            for j in xrange(len(divisions)):
                # We figure out which center to use for the box, using integer rounding tricks and 
                # modular arithmetic.
                centerindex=int(np.floor(numheirarchy[j]*i/(numum))%divisions[j])
                umbrellaCoord.append(centersList[j][centerindex])
            # We make da box.
            wraparray=np.array(boxwrap)
            # We check if any elements of wraparray are nonzero.
            if wraparray.any():
   	        if basisType == "Box":
   	            um.append(basisFunctions.Box(umbrellaCoord,widthlist,boxwrap))
   	        elif basisType == "Gaussian":
   	            um.append(basisFunctions.Gaussian(umbrellaCoord,widthlist,boxwrap))
   	        elif basisType == "Pyramid":
   	            um.append(basisFunctions.Pyramid(umbrellaCoord,widthlist,boxwrap,max_entryPoints))
            else:
   	        if basisType =="Box":
  		    um.append(basisFunctions.Box(umbrellaCoord,widthlist))
   	        elif basisType == "Gaussian":
                    um.append(basisFunctions.Gaussian(umbrellaCoord,widthlist))
                elif basisType == "Pyramid":
                    um.append(basisFunctions.Pyramid(umbrellaCoord,widthlist,max_entryPoints))
                    
        # if we specify to build a neighborlist for the windows, let's build it here. 
        L = np.zeros(len(colVarParams))
        for i in range(len(L)):
            if wrapping[i] == 1:
                L[i] = colVarParams[i][1] - colVarParams[i][0]
            else:
                L[i] = -1.0
        
        if neighborList:
            self.buildNeighborList(L, um)
    
        return um

    def buildKeylistIndexMap(self, keylist):
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

    def communicateMPI(self, rank, comm, sparseSolve=False, finiteTime=False, debug=False):
        """
        This routine performs a round of MPI communication to synchronize information across all ranks.

        Let's note here that there is a sequence of communication / computation that takes place in this section.
        Therefore it makes the most sense to interleave the two components here.

        It should be noted that this being a initial version of the parallel neus code, an optimization of the
        MPI communication could be implemented at a later date.
        """

        """
        Step 1) Communication of M as an all reduce
        """

        # first we need to communicate the M matricies in which we've accumulated transition statistics

        # now reduce the M matrix at root, first making a send buffer,
        self.Mbuff = copy.deepcopy(self.M)
        comm.Allreduce([self.Mbuff, MPI.DOUBLE], [self.M, MPI.DOUBLE], op=MPI.SUM)

        comm.Barrier()

        #if rank == 0: print rank, "after", self.active_windows.sum()

        """
        Step 2) solution of the eigenvalue problem at each rank
        """
        # at rank, update F and compute z
        
        for row in range(self.F.shape[0]):
            if self.M[row].sum() > 0.0:
                self.updateF(row, self.epsilon)


        #if rank == 0: print self.z

        self.computeZ(sparseSolve=sparseSolve, finiteTime=finiteTime)

        #if rank == 0: print self.z



        """
        Step 3) estimation of observables on each rank
        """
        self.computeObservables()

        """
        Step 4) all reduction of averaged observables to each processor
        """

        for obs in self.observables:
            self.lbuff = copy.deepcopy(obs.data)
            comm.Allreduce([self.lbuff, MPI.DOUBLE], [obs.data, MPI.DOUBLE], op=MPI.SUM)

            self.lbuff = copy.deepcopy(obs.weights)
            comm.Allreduce([self.lbuff, MPI.DOUBLE], [obs.weights, MPI.DOUBLE], op=MPI.SUM)

            # normalize by the weights for nonzero weights only
            obs.data[obs.weights.astype(bool)] /= obs.weights[obs.weights.astype(bool)]

        if debug: print self.observables[0].data, rank


        return 0

