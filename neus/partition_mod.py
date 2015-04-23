# -*- coding: utf-8 -*-
"""
This is an update partition module that contains the flexibility to perform both umbrella and NEUS sampling. 
"""

import numpy as np
import copy
from scipy import linalg as LA
import scipy as sp
import random
import basisFunctions_neus_dipeptide as basisFunctions
import entryPoints
import h5py

class partition:
    """
    This class defines a set of windows that partition the sampling space.
    This class will contain an array of basisFunction objects. 
    """
    def __init__(self, N):
        """
        Init routine. Initializes the following values:

        The obervables list is a bookkeeping of the observables of the system one
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
        # we should track how many samples are taken in the M matrix
        self.nsamples_M = np.zeros(N)
        self.a = np.zeros(N)
        self.m = np.zeros(N)
        self.z = np.zeros(N)

        # keep a list of any observable objects that should be estimated during sampling. 
        self.observables = []

        self.stoppingTimes = np.zeros(N)
        self.simulationTime = np.zeros(N)
        
        # here, the variable k represents how many times the transition matrix has been updated. 
        self.k = 0
        
    def updateF(self, updateType=None):
        """
        This routine update G according to:
            
            G_{ij}^{k+1} = (1 - \epsilon_{k}) * G_{ij}^{k} + \epsilon_{k} * M_{ij} / T_{i}

        for a finite time problem. 
        """
        assert updateType in ['finiteTime', 'neus', 'overlap'], "Updatetype not understoond." 

        if updateType == 'finiteTime':
            
            assert np.any(self.stoppingTimes[:] != 0.0), "Stopping Times value was not properly initialized."
            
            temp_M = np.zeros(self.M.shape)
        
            for i in range(self.M.shape[0]):
                temp_M[i,:] = self.M[i,:] / self.stoppingTimes[i]
        
            self.F = (self.k * self.F + temp_M) / (self.k + 1)

        elif updateType == 'neus':
            
            temp_M = np.zeros(self.M.shape)

            for i in range(self.M.shape[0]):
                temp_M[i,:] = self.M[i,:] / self.stoppingTimes[i]

                # update diagonal elements of G in infinite time process
                temp_M[i,i] = 1 - temp_M[i,:].sum()

            self.F = (self.k * self.F + temp_M) / (self.k + 1)
     
        elif updateType == 'overlap':
            
            temp_M = np.zeros(self.M.shape)
            
            for i in range(self.M.shape[0]):
                temp_M[i,:] = self.M[i,:] / self.nsamples_M[i]
                
            self.F = (self.k * self.F + temp_M) / (self.k + 1)
                

        return 0
    
    def updateA(self):
        """
        This routine updates the a vector via:
            
            a_{i}^{k+1} = 1 - Sum_{l=1}^{n} G_{il}^{k+1}
            
        """
        for row in range(self.a.shape[0]):
            self.a[row] = 1 - self.G[row,:].sum()
        
        return 0

    def addObservable(self, A):
        """
        This routine adds an observable and initializes local copies of the observables in the basis windows.
        """
        assert hasattr(self, "observables"), "The partition observables lists were not properly initialized."

        # add the observable to both the partition and each window
        self.observables.append(A)
        for window in self.umbrellas:
            if not hasattr(window, "local_observables"): window.local_observables = []
            window.local_observables.append(copy.deepcopy(A))

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
        

    def computeObservables(self):
        """
        This routine populates the partition observables with data averaged from the windows.
        """
        for o_indx,obs in enumerate(self.observables):
            temp = np.zeros(obs.data.shape)

            for w_indx,win in enumerate(self.umbrellas):
                temp += self.z[w_indx] * win.local_observables[o_indx].data

            obs.data[:] = temp[:]

        return 0

    def computeZ(self, finiteTime=False):
        """
        Solves for z vector given current G,a via solving the following linear 
        system:
            
            (I - G)^T z = a 
	
	    for a finite time process and
	
            zG = z 

	    for infinite time process.
        """
        if finiteTime:
            A = (np.identity(self.F.shape[0]) - self.F).transpose()
        
            self.z = np.linalg.solve(A, self.a)

        else:
            # compute via numpy interface to LAPACK the eigenvectors v and eigenvalues w
            # The linalg routine returns this as the first (i.e. largest) eigenvalue.
            evals, evec = LA.eig(self.F, left=True, right=False)
            sort = np.argsort(evals)
            # normalize if needed.
            self.z = evec[:,sort[-1]] / np.sum(evec[:,sort[-1]])

        return 0

    def accumulateObservables(self, sample, colvars, indx):
        """
        This routine loops through the observables list and updates the samples in the corresponding windows.
        """

        for obs in self.umbrellas[indx].local_observables:
            # use the observable call routine to accumulate a sampling in the observables data structure
            obs(sample, colvars)

        return 0

    def resetObservable(self, obs):
        """
        Set each element in the data of the passed observable to zero.
        """
        for elem in np.nditer(obs.data, op_flags=['readwrite']):
            elem[...] = 0.0
            
        for elem in np.nditer(obs.nsamples, op_flags=['readwrite']):
            elem[...] = 0.0

        return 0
        
    def reinject(self, wlkr, i):
        """
        This function initializes a simulation from the entry point list in the 
        current umbrella.
        """
        assert len(self.umbrellas[i].entryPoints) != 0, "Reached reinjection routine with no entry points in the buffer."
    
        # now we initialize the starting coordinates from the entry points library
        temp_indx = random.randint(0, len(self.umbrellas[i].entryPoints)-1)
    
        # you should pass this argument as a ctypes array for now
        wlkr.setConfig(self.umbrellas[i].entryPoints[temp_indx].config)
        wlkr.setVel(self.umbrellas[i].entryPoints[temp_indx].vel)
        
        # set the other component of the walker state
        wlkr.Y_s = self.umbrellas[i].entryPoints[temp_indx].Y_s
              		
        wlkr.simulationTime = self.umbrellas[i].entryPoints[temp_indx].time

        
        return 0

    def getBasisFunctionValues(self, coord, umbrella_index = None):
        """
        This function takes a point in collective variable space and returns 
        an array of the value of the basis functions at that point. 
        
        If no umbrella index is passed, then search the whole space 
        """
        # build an array 
        indicators = np.zeros(len(self.umbrellas))
        
        if umbrella_index is None:
    	
            # go through the basis functions and construct the indicators array
            for i in range(len(self.umbrellas)):
                if self.umbrellas[i].indicator(coord) == 0.0:
                    continue
                else:
                    indicators[i] = self.umbrellas[i].indicator(coord)
                    
        elif len(self.umbrellas[umbrella_index].neighborList) == 0: 
            # go through the basis functions and construct the indicators array
            for i in range(len(self.umbrellas)):
                if self.umbrellas[i].indicator(coord) == 0.0:
                    continue
                else:
                    indicators[i] = self.umbrellas[i].indicator(coord)
        else:
            
            # make sure the partition has a neighborlist
            assert hasattr(self.umbrellas[umbrella_index], 'neighborList'), "There is no neighborlist defined."
            # now loop over neighbors and populate the indicator
            for i in self.umbrellas[umbrella_index].neighborList:
                
                indicators[i] = self.umbrellas[i].indicator(coord)
                
            # if we don't find any support, let's try this again and search the whole space. 
            if np.sum(indicators) == 0.0:
                indicators = self.getBasisFunctionValues(coord, umbrella_index=None)
                
        # normalize the values of the basis functions
        assert np.sum(indicators) != 0.0
        indicators = indicators / np.sum(indicators)
        
        return indicators
    
    def initializeMat(self, ncells):
        """
        This function constructs and returns an ncells by ncells matrix as a numpy array of zeros.
        """
        F = np.zeros((ncells, ncells))
        
        return F
        
    def sample_US(self, wlkr, numSteps, umbrellaIndex, walkerIndex, sysParams):
        """
        This routine samples the given box via the equilibrium overlap method. 
        """
        assert sysParams.has_key('scratchdir'), "Scratch directory was not specified in the sampling routine."
        
        assert sysParams['transitionMatrixType'] in ['transition','overlap']
        # assign an input filename for this walker. 
        #inputFilename = sysParams['scratchdir'] + "/" + str(umbrellaIndex) + "_w" + str(walkerIndex)
        inputFilename = None
    
        oldConfig = wlkr.getConfig()
        oldSample = wlkr.getColvars()
        
        f_handle = h5py.File(sysParams['scratchdir'] + "/ep." + str(umbrellaIndex) + ".h5py", "w")
        
        #print self.umbrellas[umbrellaIndex](oldSample, self.umbrellas), self.umbrellas[umbrellaIndex].indicator(oldSample)
        assert self.umbrellas[umbrellaIndex].indicator(oldSample) > 0.0, "The walker is not in the support of the current window."
        
        # get the sample from the initial state of the walker in CV space
        self.umbrellas[umbrellaIndex].samples.append(oldSample)
        
        # reset the local observables arrays for this round of sampling
        for obs in self.umbrellas[umbrellaIndex].local_observables:
            self.resetObservable(obs)

        assert sysParams.has_key('stepLength'), "StepLength was not specified in the sampling routine."
        
        # now we proceed with the sampling routine            
        for i in range(0, numSteps, sysParams['stepLength']):
            
            # propagate the dynamics
            wlkr.propagate(sysParams['stepLength'])
            
            newConfig = wlkr.getConfig()
            newSample = wlkr.getColvars()
            
            if sysParams['transitionMatrixType'] == 'transition':
                # update the M matrix based on this sample 
                self.M[umbrellaIndex,:] = self.getBasisFunctionValues(newSample, umbrellaIndex)
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
                self.accumulateObservables(newSample, wlkr.colvars, umbrellaIndex)
            
            else:
                # if we reject the proposed position, append the old config                
                self.umbrellas[umbrellaIndex].samples.append(oldSample)
                
                # set the walker configuration to the old state
                wlkr.setConfig(oldConfig)
                
                # redraw the velocities 
                wlkr.drawVel(distType = 'gaussian', temperature = 310.0)
                
                # append the new sample to the observables
                self.accumulateObservables(oldSample, wlkr.colvars, umbrellaIndex)

            # now check to see if the data buffer has become too large and flush buffer to file
            #if len(self.umbrellas[umbrellaIndex].samples) > 1000000: self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)
            
            if sysParams['transitionMatrixType'] == 'overlap':
                self.M[umbrellaIndex,:] = self.getBasisFunctionValues(oldSample, umbrellaIndex)
                # increment the number of samples 
                self.nsamples_M[umbrellaIndex] += 1
            
            if i % 10000 == 0: f_handle.create_dataset(str(i) + ".config", data=wlkr.getConfig())
        
        
        # flush the last data to file after sampling has finished
        self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)
        
        f_handle.flush()
        f_handle.close()
        
        return 0
        
    def sample_NEUS(self, wlkr, numSteps, umbrellaIndex, walkerIndex, sysParams):
        """    
        This function takes a system lmp and propagates the dynamics to generate
        the required samples but storing and generating samples via NEUS algorithm.
        
        Specifically, this performs an umbrella sampling routine using NEUS reinjection procedure for reinitializing the walker. 
        
        We should remove the need for the output to be specified internally here. 
        
        """      
        assert sysParams.has_key('scratchdir'), "Scratch directory was not specified in the sampling routine."
        # assign an input filename for this walker. 
        inputFilename = sysParams['scratchdir'] + "/" + str(umbrellaIndex) + "_w" + str(walkerIndex)
        
        # set the T_0 value, the value essentially acts as a "countdown" for how
        # much time evolution is left to be done in this box
        T_0 = self.stoppingTimes[umbrellaIndex]
    
        # reinject the walker to start in this box. 
        self.reinject(wlkr, umbrellaIndex)
        
        # get the sample from the initial state of the walker in CV space
        self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())
        
        # reset the local observables arrays for this round of sampling
        for obs in self.umbrellas[umbrellaIndex].local_observables:
            self.resetObservable(obs)

        assert sysParams.has_key('stepLength'), "StepLength was not specified in the sampling routine."
        
        # now we proceed with the sampling routine            
        for i in range(0, numSteps, sysParams['stepLength']):
            
            # propagate the dynamics
            wlkr.propagate(sysParams['stepLength'])
            
            # update the record of the simulation time for the walker object. 
            wlkr.simulationTime += sysParams['stepLength']

            # for infinite time processes let's simply count down
            T_0 -= sysParams['stepLength']
	    
            # now we check to see if we've passed the autocorrelation length
            # if we do, we reset the Y ( [t / s] * s) value to the current point
            if (wlkr.simulationTime % sysParams['corrLength']) == 0.0:
                wlkr.Y_s = (wlkr.getConfig(), wlkr.getVel())
            
            # get the new sample position
            self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())             

            # check for a transition out of this index
            if self.umbrellas[umbrellaIndex].indicator(self.umbrellas[umbrellaIndex].samples[-1]) == 0.0:
                
                # choose the new j with probability {psi_0, ..., psi_N}
                indicators = self.getBasisFunctionValues(self.umbrellas[umbrellaIndex].samples[-1])
            
                # record a transition to the matrix
                self.M[umbrellaIndex,:] += indicators
                
                # now we select a window to which to append this new entry point
                randVal = random.random()                

                for indx in range(len(indicators)):
                    if randVal < indicators[:indx+1].sum():
                        
                        # now we will need to reinject the walker
                        
                        # create a new entry point and append the entry point to the new window
                        newEP = entryPoints.entryPoints(wlkr.getConfig(), wlkr.getVel(), wlkr.simulationTime)
                        newEP.Y_s = wlkr.Y_s
                        self.umbrellas[indx].newEntryPoints.append([umbrellaIndex, newEP])
                        
                        # we're done so let's stop the loop
                        break
                
                
                # drop the last point from the samples 
                self.umbrellas[umbrellaIndex].samples.pop()

                # reinject the walker into the current window 
                self.reinject(wlkr, umbrellaIndex)
            		
                # get the sample from the new starting point after reinjection
                self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())

                    

            
            # let's accumulate a sample into the autocorrelation function
            self.accumulateObservables(wlkr.getColvars(), wlkr.colvars, umbrellaIndex)

            # now check to see if the data buffer has become too large and flush buffer to file
            #if len(self.umbrellas[umbrellaIndex].samples) > 1000000: self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)
        
        # flush the last data to file after sampling has finished
        self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)   

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
    
    def buildNeighborList(self, L, umbrellas):
        """
        This routine constructs a neighborlist based on the radius of the basis 
        function.
        
        L is the vector specifying the box lengths in each dimension. 
        """      
        # we should make sure that ever basisFunction has a radius defined. 
        for win in umbrellas:
            assert hasattr(win, 'radius')
            
        # now let's find all of the neighbors and populate a list of neighbors for each 
        for i in range(len(umbrellas)):
            for j in range(i+1, len(umbrellas)):
                # get the distance between the centers                    
                dr = umbrellas[i].center - umbrellas[j].center
                
                # apply minimum image convention if the dimension wraps 
                for indx,dim in enumerate(L):
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
                
                if dist <= (umbrellas[i].radius + umbrellas[j].radius):
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

     
    def createUmbrellas(self, colVarParams, wrapping, basisType="Box", neighborList=True):
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
   	            um.append(basisFunctions.Pyramid(umbrellaCoord,widthlist,boxwrap))
            else:
   	        if basisType =="Box":
  		    um.append(basisFunctions.Box(umbrellaCoord,widthlist))
   	        elif basisType == "Gaussian":
                    um.append(basisFunctions.Gaussian(umbrellaCoord,widthlist))
                elif basisType == "Pyramid":
                    um.append(basisFunctions.Pyramid(umbrellaCoord,widthlist))
                    
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