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
import basisFunctions_neus as basisFunctions
import acor
import errors
import random
import h5py
import numpy as np

class partition:
    """
    This class defines a set of windows that partition the sampling space.
    This class will contain an array of basisFunction objects. 
    """
    
    def __init__(self, N):
        
        # create an umbrella list
        self.umbrellas = []
        
        # initialize the matricies needed for the NEUS
        self.M = np.zeros((N,N))
        self.G = np.zeros((N,N))     
        self.a = np.zeros(N)
        self.m = np.zeros(N)
        self.z = np.zeros(N)

        self.stoppingTimes = np.zeros(N)
        self.simulationTime = np.zeros(N)
        
        self.k = 0
        
    def updateG(self, finiteTime=False):
        """
        This routine update G according to:
            
            G_{ij}^{k+1} = (1 - \epsilon_{k}) * G_{ij}^{k} + \epsilon_{k} * M_{ij} / T_{i}

        for a finite time problem. 
        """
	if finiteTime:
            temp_M = np.zeros(self.M.shape)
        
            for i in range(self.M.shape[0]):
                temp_M[i,:] = self.M[i,:] / self.stoppingTimes[i]
        
            self.G = (self.k * self.G + temp_M) / (self.k + 1)

	else:
	    temp_M = np.zeros(self.M.shape)

	    for i in range(self.M.shape[0]):
		temp_M[i,:] = self.M[i,:] / self.stoppingTimes[i]

		# update diagonal elements of G in infinite time process
		temp_M[i,i] = 1 - temp_M[i,:].sum()

	    self.G = (self.k * self.G + temp_M) / (self.k + 1)

        return 0
    
    def updateA(self):
        """
        This routine updates the a vector via:
            
            a_{i}^{k+1} = 1 - Sum_{l=1}^{n} G_{il}^{k+1}
            
        """
        for row in range(self.a.shape[0]):
            self.a[row] = 1 - self.G[row,:].sum()
        
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
            A = (np.identity(self.G.shape[0]) - self.G).transpose() 
        
            self.z = np.linalg.solve(A, self.a)

	else:
	    # compute via numpy interface to LAPACK the eigenvectors v and eigenvalues w
            # The linalg routine returns this as the first (i.e. largest) eigenvalue.
            evals, evec = LA.eig(self.G, left=True, right=False)
            sort = np.argsort(evals)
            # normalize if needed.
            self.z = evec[:,sort[-1]] / np.sum(evec[:,sort[-1]])

        return 0
        
    def computeAcor(self):
        """
        This function will accumulate the autocorrelation function from the weighted
        local estimations generated in each window. 
        """
        # first gather the local estimations of the acor function
        temp = np.zeros(self.acorr.shape)
        for indx, window in enumerate(self.umbrellas):
            temp += self.z[indx] * window.localCorrFunc
            
        # now update the global estimation
        self.acorr = (self.k * self.acorr + temp) / (self.k + 1)
        
        
        return 0
        
    def accumulateAcorr(self, wlkr, i, s, stepLength):
        """
        This function takes the current position of the walker and updates the 
        local autocorrelation function estimation. 
        """
        time_indx = (self.simulationTime[i] - np.floor(self.simulationTime[i] / s) * s ) / stepLength
        # we will need to accumulate a single sample for the given index
        temp_val = np.dot(wlkr.getVel(), wlkr.Y_s[1])
        self.umbrellas[i].localCorrFunc[time_indx] = (self.umbrellas[i].localCorrFunc[time_indx] * self.umbrellas[i].localCorrFunc_nsamples[time_indx] + temp_val) / (self.umbrellas[i].localCorrFunc_nsamples[time_indx] + 1.0)
        self.umbrellas[i].localCorrFunc_nsamples[time_indx] += 1.0
        
        return 0
        
    def reinject(self, wlkr, i, uniform=False):
        """
        This function initializes a simulation from the entry point list in the 
        current umbrella.
        """
        if len(self.umbrellas[i].entryPoints) == 0:
            # if there are no entry points recorded, set time to zero and 
            # reinject at the center of the box
            wlkr.setConfig(np.array([self.umbrellas[i].center[0], self.umbrellas[i].center[1], 0.0]))
	    wlkr.propagate(1) 
	    #vel = np.random.normal(0.0, 0.728, (3,0))
	    #wlkr.setVel(vel)
            # set the reference to the current location 
            wlkr.Y_s = (wlkr.getConfig(), wlkr.getVel())
            
            # don't redraw vel. 
            self.simulationTime[i] = 0.0

            
        else:    
            # now we initialize the starting coordinates from the entry points library
            temp_indx = random.randint(0, len(self.umbrellas[i].entryPoints)-1)
    
            # you should pass this argument as a ctypes array for now
            wlkr.setConfig(self.umbrellas[i].entryPoints[temp_indx][1])
            wlkr.setVel(self.umbrellas[i].entryPoints[temp_indx][2])
            # set the other component of the walker state
            wlkr.Y_s = (self.umbrellas[i].entryPoints[temp_indx][3], self.umbrellas[i].entryPoints[temp_indx][4])
              		
            self.simulationTime[i] = self.umbrellas[i].entryPoints[temp_indx][0]

        
        return 0

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
    
    def initializeMat(self, ncells):
        """
        This function constructs the matrix F.
        """
        F = np.zeros((ncells, ncells))
        
        return F   
        
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
        inputFilename = sysParams['scratchdir'] + "/" + str(umbrellaIndex) + "_w" + str(walkerIndex)
        
        # assign an output file for the LAMMPS trajectory:
        #wlkr.command("dump 1 all xyz " + str(sysParams['stepLength']) + " " + inputFilename + ".xyz")
        
        # set the T_0 value, the value essentially acts as a "countdown" for how
        # much time evolution is left to be done in this box
        T_0 = self.stoppingTimes[umbrellaIndex]
        
        # inject the walker 
        self.reinject(wlkr, umbrellaIndex)
        
        # get the sample from the reinjection point
        self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())
        
        # also record the initial simulation time in a separate variable
        t = self.simulationTime[umbrellaIndex]

	# reset local correlation function variables
	self.umbrellas[umbrellaIndex].localCorrFunc[:] = 0.0
	self.umbrellas[umbrellaIndex].localCorrFunc_nsamples[:] = 0.0
        
        #for i in range(0, numSteps, sysParams['stepLength']):
        while T_0 > 0.0:
            # propagate the dynamics
            wlkr.propagate(sysParams['stepLength'])
            self.simulationTime[umbrellaIndex] += sysParams['stepLength']

            # for infinite time processes let's simply count down 
	    T_0 -= sysParams['stepLength']
            
            # now we check to see if we've passed the autocorrelation length
            # if we do, we reset the Y ( [t / s] * s) value to the current point
            if (self.simulationTime[umbrellaIndex] % sysParams['corrLength']) == 0.0:
                wlkr.Y_s = (wlkr.getConfig(), wlkr.getVel())
            
            # get the new sample position
            self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())             

            """    
            # check to see if we have hit the stopping time:
            if self.simulationTime[umbrellaIndex] >= self.stoppingTimes[umbrellaIndex]:
                print "Stopping time reached."
                # update diagonal element of M
                self.M[umbrellaIndex,umbrellaIndex] += self.simulationTime[umbrellaIndex] - t - sysParams['stepLength']
                
                # now we update the T_0 value 
                T_0 = T_0 - self.simulationTime[umbrellaIndex] + t
                
                self.reinject(wlkr, umbrellaIndex)
                t = self.simulationTime[umbrellaIndex]
                
            # also check to see if we have run T_i length of simulation
            elif self.simulationTime[umbrellaIndex] >= (T_0 + t): 
                print "time greater than T_0 + t."
                # likewise, update the diagonal element of M
                self.M[umbrellaIndex,umbrellaIndex] += self.simulationTime[umbrellaIndex] - t - sysParams['stepLength']
                
                # now we update the T_0 value 
                T_0 = T_0 - self.simulationTime[umbrellaIndex] + t
                
                self.reinject(wlkr, umbrellaIndex)
                t = self.simulationTime[umbrellaIndex]
            """                                                                
            # check for a transition out of this index
            if self.umbrellas[umbrellaIndex].indicator(self.umbrellas[umbrellaIndex].samples[-1]) == 0.0:
                #print "recorded transition."
                
                # choose the new j with probability {psi_0, ..., psi_N}
                indicators = self.getBasisFunctionValues(self.umbrellas[umbrellaIndex].samples[-1])
            
                randVal = random.random()

                # now we update the T_0 value for making this transition
		#T_0 = T_0 - self.simulationTime[umbrellaIndex] + t 

                for indx in range(len(indicators)):
                    if randVal < indicators[:indx+1].sum():
                    #if self.umbrellas[indx].indicator(self.umbrellas[umbrellaIndex].samples[-1]) == 1.0:        
                        
                        # record a transition to the matrix
                        self.M[umbrellaIndex,:] += indicators
                        #self.M[umbrellaIndex][indx] += 1.0
                        
                        # now we will need to reinject the walker
                        
                        # append the entry point to the new window
                        self.umbrellas[indx].newEntryPoints.append((self.simulationTime[umbrellaIndex], wlkr.getConfig(), wlkr.getVel(), wlkr.Y_s[0], wlkr.Y_s[1]))
                        
                        # drop the last point from the samples 
                        self.umbrellas[umbrellaIndex].samples.pop()
                        
                        # reinject into this window
            		self.reinject(wlkr, umbrellaIndex)
            		
            		# get the sample from the new starting point
            		self.umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())

            		t = self.simulationTime[umbrellaIndex]

                        # we're done so let's stop the loop
                        break
                    

            
            # let's accumulate a sample into the autocorrelation function
            self.accumulateAcorr(wlkr, umbrellaIndex, sysParams['corrLength'], sysParams['stepLength'])

            # now check to see if the data buffer has become too large and flush buffer to file
            if len(self.umbrellas[umbrellaIndex].samples) > 1000000:
                self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)
        
        print "Transitions recorded: ", self.M[umbrellaIndex,:].sum() - self.M[umbrellaIndex,umbrellaIndex]
        # flush the last data to file
        self.umbrellas[umbrellaIndex].flushDataToFile(inputFilename)   

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

     
    def createUmbrellas(self, Params, basisType="Box"):
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
        #First, we pick up the parameters we care about from the systemParams.
        colVarParams = Params["cvrange"]
        # We check if the user provided any input to wrap the basis functions.
        if 'wrapping' in Params:
            wrapping=np.array(Params['wrapping'])
        else:
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
   	        elif basisType == "Cone":
   	            um.append(basisFunctions.Cone(umbrellaCoord,widthlist,boxwrap))
            else:
   	        if basisType =="Box":
  		    um.append(basisFunctions.Box(umbrellaCoord,widthlist))
   	        elif basisType == "Cone":
                    um.append(basisFunctions.Cone(umbrellaCoord,widthlist))
    
        return um