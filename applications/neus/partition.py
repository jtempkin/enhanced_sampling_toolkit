# -*- coding: utf-8 -*-
"""
This module contains the definition of the partition object. The partition object represents the definition of a 
spatial decomposition defined in umbrella samping / stratification algorithms. 

The basic features of this partition and how one uses this data structure to perform umbrella sampling type calculations
is described as follows. 
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
    def __init__(self, scratchdir=None, parallel=False):
        """
        Init routine. Initializes the following values:

        The obervables list is a bookkeeping of the observables of the self one
        wants to record during the simulaiton. This list needs to be set prior to sampling.

        The list shoud contain elements that are tuples of a function and a data array.

        Each element should take two arguments, one the walker object and one the umbrella index for which the sample is
        associated.
        """

        # given a list of center locations and widths, initialize a 

        # create an umbrella list
        self._umbrellas = []

        # keep a list of any observable objects that should be estimated during sampling. 
        self.observables = []

        return None

    def __call__(self, walker):
        """
        Returns an array of the support of the partition at the current CV state of the walker.
        """
        # build an array 
        indicators = np.asarray([win(walker) for win in self])

        # normalize the values of the basis functions
        if np.sum(indicators) == 0.0:
            print wlkr.getColvars()
            print wlkr.simulationTime
            print wlkr.Y_s
        assert np.sum(indicators) != 0.0, str(wlkr.getColvars())

        return indicators / indicators.sum()


        coord = wlkr.getColvars()
        
        if umbrella_index is None:
        
            # go through the basis functions and construct the indicators array
            for i in xrange(len(self._umbrellas)):
                if self._umbrellas[i].indicator(coord) == 0.0:
                    continue
                else:
                    indicators[i] = self._umbrellas[i](wlkr, self._umbrellas)
                    
        elif len(self._umbrellas[umbrella_index].neighborList) == 0: 
            # go through the basis functions and construct the indicators array
            for i in xrange(len(self._umbrellas)):
                if self._umbrellas[i](wlkr, self._umbrellas) == 0.0:
                    continue
                else:
                    indicators[i] = self._umbrellas[i](wlkr, self._umbrellas)
        else:
            
            # make sure the partition has a neighborlist
            assert hasattr(self._umbrellas[umbrella_index], 'neighborList'), "There is no neighborlist defined."
            # now loop over neighbors and populate the indicator
            for i in self._umbrellas[umbrella_index].neighborList:
                
                indicators[i] = self._umbrellas[i](wlkr, self._umbrellas)
                
            # if we don't find any support, let's try this again and search the whole space. 
            if np.sum(indicators) == 0.0:
                indicators = self.get_basis_function_values(wlkr, umbrella_index=None)
                
        # normalize the values of the basis functions
        if np.sum(indicators) == 0.0:
            print wlkr.getColvars()
            print wlkr.simulationTime
            print wlkr.Y_s
        assert np.sum(indicators) != 0.0, str(wlkr.getColvars())

        indicators = indicators / np.sum(indicators)
        
        return indicators

    def __getitem__(self, pos):
        """
        Definition of special method for returning window objects. 
        """

        return self._umbrellas[pos]

    def __len__(self):
        """
        Definition of the size of the partition. Returns number of windows.
        """

        return len(self._umbrellas)

    def win_index(self, win):
        """
        Returns the index of teh given window in the partition.
        """
        return self._umbrellas.index(win)

    def set_umbrellas(self, umbrellas, neighborList=True,  s=None):
        """
        This routine replaces the current list of umbrellas with those provided
        and updates the matricies contained to match the new umbrella list. This 
        behavior is destructive to the old matricies. 
        """ 

        self._umbrellas = umbrellas

        # get the size of the new partition 
        N = len(self)

        # build the neighborlist in the umbrellas by default. 
        if neighborList:

            assert s is not None, "Stopping time needed to construct neighborlist."

            self.build_neighbor_list(umbrellas, s)
        # if there is no neighbor list, set variable to be every other window. 
        else:

            for win in self:
                win.neighborList = np.arange(N)

        # initialize the matricies needed for the NEUS
        self.M = np.zeros((N,N))
        self.F = np.zeros((N,N))     

        # initialize the F_list 
        self.F_list = []

        [self.F_list.append([]) for i in range(N)]

        # set a and z
        self.a = np.zeros(N)
        self.z = np.zeros(N)

        # we should track how many samples are taken in the M matrix
        self.nsamples_M = np.zeros(N)
        
        # start a list of the windows associated with this partition. 
        self.rank_window_index = None

        self.simulationTime = np.zeros(N)
        
        # here, the variable k represents how many times the transition matrix has been updated. 
        self.k = np.zeros(N)

        # we'll also store a boolean of active windows 
        self.active_windows = np.zeros(N, dtype=np.bool)

        return 0 
        
    def update_F(self, row, epsilon):
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

        # update an total average
        self.F[row] = (self.k[row] * self.F[row] + temp_M) / (self.k[row] + 1)

        # keep a list of the most recent F 
        """
        if len(self.F_list[row]) >= epsilon:
            self.F_list[row].pop(0)

        self.F_list[row].append(temp_M)
        self.F[row].fill(0.0)
        for m in self.F_list[row]:
            self.F[row] += m
        
        self.F[row] /= len(self.F_list[row])
        """

        # the one with the decay parameter
        #self.F[row] = ((1-epsilon) * self.F[row] + epsilon * temp_M)

        if self.k[row] / (self.k[row] + 1) < epsilon:
            self.k[row] += 1

        return 0

    def add_observable(self, A, rank_index=None):
        """
        This routine adds an observable and initializes local copies of the observables in the basis windows.
        """
        assert hasattr(self, "observables"), "The partition observables lists were not properly initialized."

        # add the observable to both the partition and each window own by the ranks
        self.observables.append(A)
        
        if self.rank_window_index is None:
            for window in self._umbrellas:
                if not hasattr(window, "local_observables"): window.local_observables = []
                window.local_observables.append(copy.deepcopy(A))
        else:
            for window in self.rank_window_index:
                if not hasattr(self._umbrellas[window], "local_observables"): self._umbrellas[window].local_observables = []
                self._umbrellas[window].local_observables.append(copy.deepcopy(A))

        return 0

    def remove_observable(self):
        """
        This routine removes all observables from the list for the partition. 
        """
        # remove global version of observables
        self.observables = []
        
        # remove local version of observables from the windows as well. 
        for win in self._umbrellas:
            win.local_observables = []
            
        return 0 
        

    def compute_observables(self, rank=None):
        """
        This routine populates the partition observables with data averaged from the windows.
        """
        for o_indx,obs in enumerate(self.observables):
            temp = np.zeros(obs.data.shape)
            temp_weights = np.zeros(obs.nsamples.shape)
            
            if rank is None:
                for w_indx,win in enumerate(self._umbrellas):
                    temp += self.z[w_indx] * win.local_observables[o_indx].data
                    temp_weights += self.z[w_indx] * win.local_observables[o_indx].nsamples.astype(bool)
            else:
                for w_indx in self.rank_window_index[rank]:
                    temp += self.z[w_indx] * self._umbrellas[w_indx].local_observables[o_indx].data
                    temp_weights += self.z[w_indx] * self._umbrellas[w_indx].local_observables[o_indx].nsamples.astype(bool)


            obs.data[:] = temp[:]
            obs.weights[:] = temp_weights[:]

        return 0

    def compute_z(self, sparseSolve=True, finiteTime=False):
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

            """
            Construct G_bar as

            G_bar = [[ G     1-G1

                        a     0   ]]

            Then solve the eigenvalue problem choosing "pivot" off maximum z excluding the final row.

            We'll choose a max_z

            and solve the equation z_bar^T = z_bar^t G_bar

            using the following strategy:
            """

            G_bar = np.zeros((self.F.shape[0]+1, self.F.shape[0]+1))

            G_bar[:-1, :-1] = self.F
            G_bar[-1, :-1] = self.a
            G_bar[:-1, -1] = np.ones(self.a.shape) - np.dot(self.F, np.ones(self.a.shape))

            #evals, evec = LA.eig(G_bar, left=True, right=False)
            #sort = np.argsort(evals)

            #temp_z = evec[:,sort[-1]] / np.sum(evec[:,sort[-1]])

            #self.z = temp_z[:-1].real

            #return 0

            # get the maximum z value 
            max_z = np.argmax(self.z)

            # slice out the row corresponding to maximum z
            temp_F = np.delete(np.delete(G_bar, max_z, 0), max_z, 1)

            # now construct A and run through solver.
            A = (np.identity(temp_F.shape[0]) - temp_F).transpose()
            temp_z = np.linalg.solve(A, np.delete(G_bar[max_z], max_z).transpose())

            z_bar = np.zeros(self.z.size+1)

            z_bar[max_z] = 1.0

            z_bar[0:max_z] = temp_z[0:max_z]
            z_bar[max_z+1:] = temp_z[max_z:]

            # now we renormalize the z_bar so that the last entry is one and take the rest
            z_bar /= z_bar[-1]

            # set remainig weights and return
            self.z = z_bar[:-1]

            #self.z /= self.z.sum()   

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

            # here we will do a different eigenvalue solution procedure to solve for the weights
                
            evals, evec = LA.eig(self.F, left=True, right=False)
            sort = np.argsort(evals)
            # normalize if needed.

            #self.z[self.active_windows] = evec[:,sort[-1]] / np.sum(evec[:,sort[-1]])
            #self.z[np.logical_not(self.active_windows)] = 0.0

            self.z = evec[:,sort[-1]] / np.sum(evec[:,sort[-1]])
            


        return 0

    def accumulate_observables(self, wlkr, sample, colvars, indx):
        """
        This routine loops through the observables list and updates the samples in the corresponding windows.
        """

        for obs in self._umbrellas[indx].local_observables:
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

    def reset_observable(self, obs):
        """
        Set each element in the data of the passed observable to zero.
        """

        obs.data.fill(0.0)
        obs.nsamples.fill(0.0)

        return 0
        
    def reinject(self, wlkr, win):
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
            if self._umbrellas[i].getNumberOfEntryPoints(key=key) == 0:
                prob[indx] = 0.0
        """

        # normalize probability
        #assert prob.sum() > 0.0
        #prob /= prob.sum()

        # now choose neighbor proportional to prob 
        #I = np.random.choice(np.arange(prob.size), p=prob)

        # get the entry point from the umbrella window

        #assert self._umbrellas[i].getNumberOfEntryPoints(key=self.index_to_key[I])

        # we choose from the initial distribution with prob. stored in this window
        if random.random() < win.initial_distribution_prob:
            dist = win.initial_distribution
            EP = random.sample(dist, 1)[0]

        else:
            EP = win.get_entry_point()

        wlkr.setConfig(EP.config)
        wlkr.setVel(EP.vel)

        # set the lag component of the walker state
        wlkr.Y_s = EP.Y_s

        wlkr.simulationTime = EP.time

        return 0

    def get_basis_function_values(self, wlkr, umbrella_index = None):
        """
        This function takes a point in collective variable space and returns 
        an array of the value of the basis functions at that point. 
        
        If no umbrella index is passed, then search the whole space 
        """
        # build an array 
        indicators = np.zeros(len(self._umbrellas))

        coord = wlkr.getColvars()
        
        if umbrella_index is None:
    	
            # go through the basis functions and construct the indicators array
            for i in xrange(len(self._umbrellas)):
                if self._umbrellas[i].indicator(coord) == 0.0:
                    continue
                else:
                    indicators[i] = self._umbrellas[i](wlkr, self._umbrellas)
                    
        elif len(self._umbrellas[umbrella_index].neighborList) == 0: 
            # go through the basis functions and construct the indicators array
            for i in xrange(len(self._umbrellas)):
                if self._umbrellas[i](wlkr, self._umbrellas) == 0.0:
                    continue
                else:
                    indicators[i] = self._umbrellas[i](wlkr, self._umbrellas)
        else:
            
            # make sure the partition has a neighborlist
            assert hasattr(self._umbrellas[umbrella_index], 'neighborList'), "There is no neighborlist defined."
            # now loop over neighbors and populate the indicator
            for i in self._umbrellas[umbrella_index].neighborList:
                
                indicators[i] = self._umbrellas[i](wlkr, self._umbrellas)
                
            # if we don't find any support, let's try this again and search the whole space. 
            if np.sum(indicators) == 0.0:
                indicators = self.get_basis_function_values(wlkr, umbrella_index=None)
                
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
        self._umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())
        
        # reset the local observables arrays for this round of sampling
        for obs in self._umbrellas[umbrellaIndex].local_observables:
            self.resetObservable(obs)
        
        # now we proceed with the sampling routine            
        for i in range(0, numSteps, stepLength):
            
            # propagate the dynamics
            wlkr.propagate(stepLength)
            
            # update the record of the simulation time for the walker object. 
            wlkr.simulationTime += stepLength

            # if this is a finite time version, we should kill the walker if it hits the boundary
            if 30.0 < wlkr.getColvars()[0] < 90.0:
                #print "hit target", self._umbrellas[umbrellaIndex].center, wlkr.simulationTime, self.z[umbrellaIndex]
                self.reinject(wlkr, umbrellaIndex)

                self._umbrellas[umbrellaIndex].nhits += 1.0 

                # now update the statistics for hitting vs stopping. 

                continue
	    
            # now we check to see if we've passed the autocorrelation length
            # if we do, we reset the Y ( [t / s] * s) value to the current point
            if corrLength is not None:
                if (wlkr.simulationTime % corrLength) == 0.0:
                    #wlkr.Y_s = (wlkr.getConfig(), wlkr.getVel(), wlkr.getColvars())
                    #wlkr.simulationTime = 0.0
                    #print "stopping time hit"

                    self._umbrellas[umbrellaIndex].nstops += 1.0

                    # if we hit the stopping time, reset the walker
                    self.reinject(wlkr, umbrellaIndex)
                    continue
            
            # get the new sample position
            new_sample = wlkr.getColvars()
            self._umbrellas[umbrellaIndex].samples.append(new_sample)

            # check for a transition out of this index
            if self._umbrellas[umbrellaIndex](wlkr, self._umbrellas) == 0.0:
                if debug: ntransitions += 1 
                # choose the new j with probability {psi_0, ..., psi_N}

                indicators = self.get_basis_function_values(wlkr)
            
                # record a transition to the matrix
                self.M[umbrellaIndex,:] += indicators
                
                # now we select a window to which to append this new entry point,use numpy to choose window index
                ep_targets = np.arange(indicators.size)[indicators.astype(bool)]
                        
                # create a new entry point and append the entry point to the new window
                newEP = entryPoints.entryPoints(wlkr.getConfig(), wlkr.getVel(), wlkr.simulationTime)
                newEP.Y_s = wlkr.Y_s

                for indx in ep_targets:
                    #if not self.active_windows[indx]: print "added entry point", indx
                    self._umbrellas[indx].addNewEntryPoint(newEP, umbrellaIndex)

                # drop the last point from the samples 
                self._umbrellas[umbrellaIndex].samples.pop()

                # reinject the walker into the current window 
                self.reinject(wlkr, umbrellaIndex)
            		
                # get the sample from the new starting point after reinjection
                self._umbrellas[umbrellaIndex].samples.append(wlkr.getColvars())

            # if we do not detect a transition and handle that, we should add a count to M_ii
            else:
                self.M[umbrellaIndex, umbrellaIndex] += 1.0

            #print self.M[umbrellaIndex, :]
            
            # let's accumulate a sample into the observables we are accumulating on this window
            self.accumulateObservables(wlkr, wlkr.getColvars(), wlkr.colvars, umbrellaIndex)

            # now check to see if the data buffer has become too large and flush buffer to file
            #if len(self._umbrellas[umbrellaIndex].samples) > 1000000: self._umbrellas[umbrellaIndex].flushDataToFile(inputFilename)
        
        # flush the last data to file after sampling has finished
        self._umbrellas[umbrellaIndex].flushDataToFile(inputFilename)

        # record the number of samples taken here
        self.nsamples_M[umbrellaIndex] = numSteps / stepLength

        # now we compute the estimate of the flux from thsi iteration
        #self.M[umbrellaIndex,umbrellaIndex] = numSteps - self.M[umbrellaIndex, :].sum()

        self.M[umbrellaIndex, :] /= self.nsamples_M[umbrellaIndex]

        # here we will store the current position of the walker in an entry point structure
        newEP = entryPoints.entryPoints(wlkr.getConfig(), wlkr.getVel(), wlkr.simulationTime)
        newEP.Y_s = wlkr.Y_s
        self._umbrellas[umbrellaIndex].walker_restart = newEP

        if debug: print "row" ,umbrellaIndex, "of M:", self.M[umbrellaIndex,:]
    
        if debug: print "Recorded",ntransitions,"transitions"

        return 0

    
    def build_neighbor_list(self, umbrellas, s, debug=False):
        """
        This routine constructs a neighborlist based on the radius of the basis 
        function.
        
        L is the vector specifying the periodic lengths in each dimension. 
        """      
        # we should make sure that ever basisFunction has a radius defined. 
        for win in self:
            assert hasattr(win, 'radius')
            
        # now let's find all of the neighbors and populate a list of neighbors for each 
        for i in range(len(self)):
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

    def build_keylist_to_index_map(self, keylist):
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
                self.update_F(row, 0.9) #self.epsilon)


        #if rank == 0: print self.z

        self.compute_z(sparseSolve=sparseSolve, finiteTime=finiteTime)

        #if rank == 0: print self.z



        """
        Step 3) estimation of observables on each rank
        """
        self.compute_observables()

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

