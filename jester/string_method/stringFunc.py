# this file contains functions related to performing operations on the string images. 
# 
# Jeremy Tempkin, 11/21/13

import numpy as np
from scipy import interpolate
import math

def reparString(colvars_input):
    "This function reparameterizes the string images."
    
    
    # first unwrap the dihedral data to ensure correct repar spaces
    colvars_input = unwrapDihed(colvars_input)
    # now get a numpy array of the colvar data.
    colvars = ColvarToBuffer(colvars_input)
    # get the colvars data shape
    nimages = colvars.shape[0]
    ncvs = colvars.shape[1]
    # get the arc length of the string.
    arc = np.zeros(nimages)
    for i in range(1, nimages, 1):
        arc[i] = arc[i-1]
        seg = 0.0
        for j in range(0, ncvs, 1):
            seg += (colvars[i][j] - colvars[i-1][j])**2
        arc[i] += np.sqrt(seg)
    arc = (arc/arc[nimages-1])
    #print "arc ", arc
    # now set up cubic spline interpolation.
    unif = np.linspace(0.0, 1.0, nimages)
    colvars_new = np.zeros_like(colvars)
    for cv in range(0, ncvs, 1):
        # linear spline:
        #f_linear = interpolate.interp1d(arc, colvars[:,cv])
        #colvars_new[:,cv] = f_linear(unif)  
        # cubic spline:
        f_spline = interpolate.splrep(arc, colvars[:,cv])
        colvars_new[:, cv] = interpolate.splev(unif, f_spline)
        
    # now convert the data structure back to the one used by main
    colvars_input = BufferToColvar(colvars_input, colvars_new)
    # wrap the dihedrals for the correct data structure syntax
    colvars_input = wrapDihed(colvars_input)
    return colvars_input
    
def smoothString(kappa, colvars):
    "This string smooths the images of the string through averaging adjacent images."
    nimages = len(colvars)
    ncvs = len(colvars[0])
    
    #print "colvars :"
    #for i in range(0, len(colvars[0]), 1):
        #print colvars[1][i][-1]
    
    for image in range(1,nimages-1,1):
        for cv in range(0, ncvs, 1):
            colvars[image][cv][-1] = colvars[image][cv][-1] + kappa * (colvars[image-1][cv][-1] + colvars[image+1][cv][-1] - 2.0 * colvars[image][cv][-1])
    
    #print "smooth colvars :"
    #for i in range(0, len(colvars[0]), 1):
        #print colvars[1][i][-1]
        
    return colvars
    
def ColvarToBuffer(colvars):
    "this function takes the colvar data structure used by main and returns a numpy array with the colvar data for the string operations."
    nimages = len(colvars)
    ncvs = len(colvars[0])
    
    colvars_buf = np.zeros((nimages, ncvs))
    
    for i in range(0, nimages, 1):
        for j in range(0, ncvs, 1):
            colvars_buf[i][j] = colvars[i][j][-1]
    
    return colvars_buf
    
def BufferToColvar(colvars, colvars_buffer):
    "this function takes a colvar buffer array and a colvar data structure and populates the second with values from the first."
    nimages = colvars_buffer.shape[0]
    ncvs = colvars_buffer.shape[1]
    
    for i in range(0, nimages, 1):
        for j in range(0, ncvs, 1):
            colvars[i][j][-1] = colvars_buffer[i][j]
    
    return colvars

def zeroColvars(colvars_diff):
    "This function takes the colvars data structure and zeros out the colvars. Useful for difference vector."
    for i in range(0, len(colvars_diff), 1):
        for j in range(0, len(colvars_diff[i]), 1):
            colvars_diff[i][j][-1] = 0.0
    
    return colvars_diff

def getDiff(colvars, colvars_CG, colvars_diff, Delta):
    "This function takes the CG and FG colvars and determines the difference vector."
    for i in range(0, len(colvars), 1):
        for j in range(0, len(colvars[i]), 1):
            colvars_diff[i][j][-1] = colvars[i][j][-1] - Delta * colvars_CG[i][j][-1]
    
    return colvars_diff

def geomAvg(colvars_diff, stringParams, w_sum, w_norm):
    "This function applies the geometric averaging to the difference term."
    # increment the geometric normalization. 
    w_norm = w_norm * float(stringParams['w_decay']) + 1.0 
    
    # now apply w_decay to history term
    for image in range(0, len(w_sum), 1):
        for cv in range(0, len(w_sum[image]), 1):
            w_sum[image][cv][-1] = w_sum[image][cv][-1] * float(stringParams['w_decay'])
    
    # now add current difference term and apply normalization.
    for image in range(0, len(w_sum), 1):
        for cv in range(0, len(w_sum[image]), 1):
            colvars_diff[image][cv][-1] = (colvars_diff[image][cv][-1] + w_sum[image][cv][-1]) / w_norm
    
    return colvars_diff, w_sum, w_norm

def copyColvars(colvars_from, colvars_to):
    "This function takes the values from one data set and sets them as the values for the other."
    for i in range(0, len(colvars_from), 1):
        for j in range(0, len(colvars_from[i]), 1):
            colvars_to[i][j][-1] = colvars_from[i][j][-1]
            
    return colvars_to

def addDiff(colvars, colvars_diff, Delta):
    "This function adds the colvars_diff to the given colvars."
    for i in range(0, len(colvars), 1):
        for j in range(0, len(colvars[i]), 1):
            colvars[i][j][-1] = Delta * colvars[i][j][-1] + colvars_diff[i][j][-1]
            
    return colvars
    

def printRMSD(colvars):
    "This function prints the image-wise RMSD in CV space to STDOUT. Mainly for debugging purposes."
    dist = [0.0]*len(colvars)
    for i in range(1, len(colvars), 1):
        for j in range(0, len(colvars[i]), 1):
            dist[i] += (colvars[i][j][-1] - colvars[i-1][j][-1])**2
            dist[i] = math.sqrt(dist[i])
        #print "%d  %lf" % (i, dist)
    
    # get the total length of the string
    total_dist = 0.0
    for i in range(0, len(dist), 1):
        total_dist += dist[i]
    
    # now normalize the dist array by the total distance of the string.
    for i in range(0, len(dist), 1):
        dist[i] /= total_dist
        print "%d  %lf" % (i, dist[i])
    
    return dist

def unwrapDihed(colvars):
    "This is a function that checks to see if the dihedral values wrap around the modulus and adjusts the ranges to remove jumps across the boundary."
    diff = 0.0
    threshold = 100.0
    for cv in range(0, len(colvars[0]), 1):
        # check that the cv is actually a dihedral.
        if colvars[0][cv][0] == 'dihedral':
            for image in range(0, len(colvars)-1, 1):
                diff = colvars[image][cv][-1] - colvars[image+1][cv][-1]
                # check for neg to pos transition
                if diff > threshold:
                    colvars[image+1][cv][-1] += 360.0
                # check for pos to neg transition
                if diff < -threshold:
                    colvars[image+1][cv][-1] -= 360.0
            
    return colvars

def wrapDihed(colvars):
    "This function takes in a colvar data set and checks to see if the colvars are properly wrapped around the -180 to 180 ratio."
    for image in range(0, len(colvars), 1):
        for cv in range(0, len(colvars[image]), 1):
            # check to see if cv is a dihedral
            if colvars[image][cv][0] == 'dihedral':
                if colvars[image][cv][-1] > 180.0:
                    colvars[image][cv][-1] -= 360.0
                elif colvars[image][cv][-1] < -180.0:
                    colvars[image][cv][-1] += 360.0
            
    return colvars

def genInterp(colvars, colvars_old, nSteps):
    "This function takes two colvars sets and generates an interpolation between them."
    # building the needed data structure.
    colvars_interp = []
    for i in range(0, nSteps, 1):
        colvars_interp.append(colvars_old)
    # now generate each interpolation and attach it to the data structure. 
    for i in range(0, len(colvars), 1):
        for j in range(0, len(colvars[i]), 1):
            interp = np.linspace(colvars_old[i][j][-1], colvars[i][j][-1], nSteps)
            for k in range(0, nSteps, 1):
                colvars_interp[k][i][j][-1] = interp[k]
                
    return colvars_interp