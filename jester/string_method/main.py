# MS STRING - NAMD
# -------------------------------------
# The main file for the python script that runs MS string calculations using NAMD 
# as the engine.
# Jeremy Tempkin
# 11/6/13

#----------IMPORT------------
import sys 
import time
import copy
import numpy as np
from subprocess import *
from fileIO import * # import string configuration parameters
from initSimul import * # input initialization functions
from callNAMD import * # import NAMD subroutines
from stringFunc import * # import functions related to controlling the string
from callCHARMM import * # import CHARMM subroutines


#----------SET SYSTEM PARAMETERS------------

#----------SET UP STRING--------------------

#----------MAIN LOOP------------------------


#----------FILE INPUT------------

# the system parameters are stored in a dictionary, sysParam
stringParams = {}
# read in the string parameters. 
err = readStringParams(sys.argv[1], stringParams) 

# check for errors in reading the string parameters. 
if err == 0: 
    print "There was an error reading in the system parameters. Exiting now."
    sys.exit(0)

# report successful parameter read.
print "The following parameters were read in:"
print stringParams

# check to make sure the stringType was understood. 
if stringParams['stringType'] == 'FG' or 'CG' or 'MS':
    print "The stringType selected is: " + stringParams['stringType']
else:
    print "stringType was not understood. Exiting now."
    sys.exit(0)

# now read in the collective varibale definiions.
print "Now reading in collective variable definitions for selected stringType: " + stringParams['stringType']

colvars = readColVarsDef(stringParams, stringParams['scratchdir'] + stringParams['CVFile'])
#print colvars
    
# get length of colvars.
stringParams['ncvs'] = len(colvars)
print "Collective variables read complete."
print str(stringParams['ncvs']) + " colvars were read."

# read in initial string images from structure input files.
#structure = readPDBstructure("./ionized.pdb")
#print "Initial Structures read in:"
#print structure

#-----------SIMULATION SETUP-------------

# establish dictionary for containing the simulation parameters that are used to set up NAMD runs
SimulParam = {}
initSimulParam(SimulParam)
print "Simulation Parameters initialized."

# set working directory for current image to the scratchdir
stringParams['imagedir'] = stringParams['scratchdir'] + "image"

# now build the batch submission files.
writeBatchFiles(stringParams)
# now write out a collective variable config file.
# this file defines the colvars to NAMD or CHARMM when they are recorded in the release phase.
writeColVarsConfig(colvars, stringParams)
# now, build the colvar data structure and populate for each image.
colvars = [0]*int(stringParams['nimages'])
if stringParams['stringType'] == 'FG' or 'MS':
    for i in range(0, int(stringParams['nimages']), 1):
        colvars[i] = readColVarsDef(stringParams, stringParams['imagedir'] + str(i) + "/colvars_img" + str(i) + ".restraint")
        #print colvars[i]

if stringParams['stringType'] == 'CG':
    for i in range(0, int(stringParams['nimages']), 1):
        colvars[i] = readColVarsDef(stringParams, stringParams['imagedir'] + str(i) + "/charmm_colvars_img" + str(i) + ".restraint")
        #print colvars[i]

        
        
# if running a MS string, initialize the needed colvars data structures
if stringParams['stringType'] == 'MS':
    colvars_CG = copy.deepcopy(colvars)
    colvars_old_CG = copy.deepcopy(colvars)
    colvars_diff = copy.deepcopy(colvars)
    colvars_new = copy.deepcopy(colvars)
    colvars_CG_new = copy.deepcopy(colvars)
    
colvars_old = copy.deepcopy(colvars)

# initialize geometric averaging terms.
w_norm = 0.0
w_sum = copy.deepcopy(colvars)
w_sum = zeroColvars(w_sum)


# open up a system log file for the colvars.
if int(stringParams['firstStep']) == 0:
    ofile_colvars = open(stringParams['scratchdir'] + "colvars.out", "w")
    writeColVarsToFile(colvars, ofile_colvars, 0)
else:
    ofile_colvars = open(stringParams['scratchdir'] + "colvars.out", "a")
    print "Restart run was requested. Starting at step: " + stringParams['firstStep'] + "."

# read in the previous colvars step locations. 
if int(stringParams['firstStep']) > 0:
    readRestartColvars(colvars, colvars_old, stringParams['scratchdir'] + "colvars.out")


    
#---------MAIN-----------

start_time = time.time()

for itr in range(int(stringParams['firstStep']), int(stringParams['niter']), 1):
    
    print "Iteration " + str(itr)
    
    #print "Starting Colvars :", itr
    #for i in range(0, len(colvars[0]), 1):
        #print  colvars[0][i][-1]
    #for i in range(0, len(colvars[0]), 1):
        #print colvars[0][i][3]
    
    if itr != int(stringParams['firstStep']):
        writeColVarsToFile(colvars, ofile_colvars, itr)
    
    if stringParams['stringType'] == 'FG':
        #----------FG ITERATION-------------    
        print "Iter", itr, ": Starting FG iteration."
        print "Iter", itr, ": Writing restraint file."
        # equilibrate to the target image in steps. 
        
        # generate the interpolation. 
        colvars_interp = genInterp(colvars, colvars_old, int(stringParams['nDragSteps_FG']))
        # make a copy of the current colvars location
        colvars_old = copyColvars(colvars, colvars_old)
        
        for step in range(0, int(stringParams['nDragSteps_FG']), 1):
            for i in range(0, int(stringParams['nimages']), 1):
                # setting input files for constrained FG run
                writeColVarsRestraint(colvars_interp[step][i], stringParams['restraint'], stringParams['imagedir'] + str(i) + "/colvars_img" + str(i) + ".restraint")
            
            print "Iter", itr, ": Submitting restraint batch step" + str(step) + "."
            # call constrained NAMD run
            callNAMDConstraint(SimulParam, stringParams)
        
        
        # now release the image
        for i in range(0, int(stringParams['nimages']), 1):
            # setting input files for release FG run
            writeColVarsRestraint(colvars[i], stringParams['weakRestraint'], stringParams['imagedir'] + str(i) + "/colvars_img" + str(i) + "_weak.restraint")
        
        print "Iter", itr, ": Submitting release batch."
        # call unconstrained NAMD run
        callNAMD(SimulParam, stringParams)
        
        # now collect the colvar data from the release phase.
        for i in range(0, int(stringParams['nimages']), 1):
            print "Iter", itr, ": Reading colvar traj " + str(i)
            #print "Before Read " + str(i) + ": \n", colvars[i]
            colvars[i] = readColVarsTraj(stringParams, colvars[i], stringParams['imagedir'] + str(i) + "/img" + str(i) + ".colvars.traj")
            #print " \nAfter read " + str(i) + ": \n", colvars[i]
    
        #----------REPARAMETERIZE STRING VARIABLES---------------
        for i in range(0, int(stringParams['reparItr']), 1):
            print "Repar Iteration: " + str(i)
            smoothString(float(stringParams['kappa']), colvars)
            colvars = reparString(colvars)
            printRMSD(colvars)
        
        #for i in range(0, len(colvars), 1):
            #print "Image ", i
            #for j in range(0, len(colvars[i]), 1):
                #print colvars_old[i][j][-1], colvars[i][j][-1]
        #print "After Repar."
        #for i in range(0, len(colvars[0]), 1):
            #print colvars[0][i][-1]
        
        #-----------------WRITE COLVARS-------------------------
        # save out trajectory files for image steps. 
        if itr % int(stringParams['backupFreq']) == 0:
            saveTrajFiles(stringParams, itr)
        
       
    
    elif stringParams['stringType'] == 'CG': 
        #------------CG ITERATION-----------------
        print "Iter", itr, ": Starting CG iteration."
        print "Iter", itr, ": Writing restraint file."
        
        # generate the interpolation. 
        colvars_interp = genInterp(colvars, colvars_old, int(stringParams['nDragSteps_CG']))
        # make a copy of the current colvars location
        colvars_old = copyColvars(colvars, colvars_old)
        
        for step in range(0, int(stringParams['nDragSteps_CG']), 1):
            for i in range(0, int(stringParams['nimages']), 1):
                # setting input files for constrained CG run
                CHARMM_writeColVarsRestraint(colvars_interp[step][i], stringParams['restraint_CG'], stringParams['imagedir'] + str(i) + "/charmm_colvars_img" + str(i) + ".restraint")

            print "Iter", itr, ": Submitting restraint batch step " + str(step) + "."
            # call constrained CHARMM run
            callCHARMMConstraint(SimulParam, stringParams)
        
        for i in range(0, int(stringParams['nimages']), 1):
            # setting input files for release CG run
            CHARMM_writeColVarsRestraint(colvars[i], stringParams['weakRestraint_CG'], stringParams['imagedir'] + str(i) + "/charmm_colvars_img" + str(i) + "_weak.restraint")
        
        print "Iter", itr, ": Submitting release batch."
        # call unconstrained CHARMM run
        callCHARMM(SimulParam, stringParams)
    
        # now collect the colvar data from the release phase.
        for i in range(0, int(stringParams['nimages']), 1):
            #print "Before Read " + str(i) + ": \n", colvars[i]
            print "Iter", itr, ": Reading colvar traj " + str(i)
            colvars[i] = CHARMM_readColVarsTraj(stringParams, colvars[i], stringParams['imagedir'] + str(i) + "/charmm_img" + str(i) + ".colvars.traj") 
            #print " \nAfter read " + str(i) + ": \n", colvars[i]
            
        #----------REPARAMETERIZE STRING VARIABLES---------------
        for i in range(0, int(stringParams['reparItr']), 1):
            print "Repar Iteration: " + str(i)
            smoothString(float(stringParams['kappa']), colvars)
            colvars = reparString(colvars)
            printRMSD(colvars)
        
        #for i in range(0, len(colvars), 1):
            #print "Image ", i
            #for j in range(0, len(colvars[i]), 1):
                #print colvars_old[i][j][-1], colvars[i][j][-1]
        #print "After Repar."
        #for i in range(0, len(colvars[0]), 1):
            #print colvars[0][i][-1]
        
        #-----------------WRITE COLVARS-------------------------
        # save out trajectory files for image steps. 
        if itr % int(stringParams['backupFreq']) == 0:
            saveTrajFiles(stringParams, itr)
        
            
    elif stringParams['stringType'] == 'MS': 
        # zero out the difference vector
        colvars_diff = zeroColvars(colvars_diff)
    
        #----------FG ITERATION-------------    
        print "Iter", itr, ": Starting FG iteration."
        print "Iter", itr, ": Writing restraint file."
        
        # generate the interpolation. 
        colvars_interp = genInterp(colvars, colvars_old, int(stringParams['nDragSteps_FG']))
        # make a copy of the current colvars location
        colvars_old = copyColvars(colvars, colvars_old)
        colvars_new = copyColvars(colvars, colvars_new)
        
        for step in range(0, int(stringParams['nDragSteps_FG']), 1):
            for i in range(0, int(stringParams['nimages']), 1):
                # setting input files for constrained FG run
                writeColVarsRestraint(colvars_interp[step][i], stringParams['restraint'], stringParams['imagedir'] + str(i) + "/colvars_img" + str(i) + ".restraint")
            
            print "Iter", itr, ": Submitting restraint batch step" + str(step) + "."
            # call constrained NAMD run
            callNAMDConstraint(SimulParam, stringParams)
        
        for i in range(0, int(stringParams['nimages']), 1):
            # setting input files for release FG run
            writeColVarsRestraint(colvars[i], stringParams['weakRestraint'], stringParams['imagedir'] + str(i) + "/colvars_img" + str(i) + "_weak.restraint")
        
        print "Iter", itr, ": Submitting release batch."
        # call unconstrained NAMD run
        callNAMD(SimulParam, stringParams)
        
        # now collect the colvar data from the release phase.
        for i in range(0, int(stringParams['nimages']), 1):
            #print "Before Read " + str(i) + ": \n", colvars[i]
            print "Iter", itr, ": Reading colvar traj " + str(i)
            colvars_new[i] = readColVarsTraj(stringParams, colvars_new[i], stringParams['imagedir'] + str(i) + "/img" + str(i) + ".colvars.traj")
            #print " \nAfter read " + str(i) + ": \n", colvars[i]
            
        #----------REPARAMETERIZE FG STRING VARIABLES---------------
        for i in range(0, int(stringParams['reparItr']), 1):
            print "Repar Iteration: " + str(i)
            smoothString(float(stringParams['kappa']), colvars_new)
            colvars_new = reparString(colvars_new)
            printRMSD(colvars_new)

        ofile_colvars_FG = open(stringParams['scratchdir'] + "FG_colvars.out", "w")
        writeColVarsToFile(colvars_new, ofile_colvars_FG, itr)
        ofile_colvars_FG.close()
        #------------CG ITERATION-----------------
        print "Iter", itr, ": Starting CG iteration."
        print "Iter", itr, ": Writing restraint file."
        
        # generate the interpolation
        colvars_interp_CG = genInterp(colvars_CG, colvars_old_CG, int(stringParams['nDragSteps_CG']))
        # make a copy of the current colvars location
        colvars_old_CG = copyColvars(colvars_CG, colvars_old_CG)
        colvars_CG_new = copyColvars(colvars_CG, colvars_CG_new)
        
        for step in range(0, int(stringParams['nDragSteps_CG']), 1):
            for i in range(0, int(stringParams['nimages']), 1):
                # setting input files for constrained CG run
                CHARMM_writeColVarsRestraint(colvars_interp_CG[step][i], stringParams['restraint_CG'], stringParams['imagedir'] + str(i) + "/charmm_colvars_img" + str(i) + ".restraint")
            
            print "Iter", itr, ": Submitting restraint batch step" + str(step) + "."
            # call constrained CHARMM run
            callCHARMMConstraint(SimulParam, stringParams)
        
        for i in range(0, int(stringParams['nimages']), 1):
            # setting input files for release CG run
            CHARMM_writeColVarsRestraint(colvars_CG[i], stringParams['weakRestraint_CG'], stringParams['imagedir'] + str(i) + "/charmm_colvars_img" + str(i) + "_weak.restraint")
        
        print "Iter", itr, ": Submitting release batch."
        # call unconstrained CHARMM run
        callCHARMM(SimulParam, stringParams)
    
        # now collect the colvar data from the release phase.
        for i in range(0, int(stringParams['nimages']), 1):
            #print "Before Read " + str(i) + ": \n", colvars[i]
            print "Iter", itr, ": Reading colvar traj " + str(i)
            colvars_CG_new[i] = CHARMM_readColVarsTraj(stringParams, colvars_CG_new[i], stringParams['imagedir'] + str(i) + "/charmm_img" + str(i) + ".colvars.traj")
            #print " \nAfter read " + str(i) + ": \n", colvars[i]
        
        #----------REPARAMETERIZE CG STRING VARIABLES---------------
        for i in range(0, int(stringParams['reparItr']), 1):
            print "Repar Iteration: " + str(i)
            smoothString(float(stringParams['kappa']), colvars_CG_new)
            colvars_CG_new = reparString(colvars_CG_new)
            printRMSD(colvars_CG_new)
        
        
        ofile_colvars_CG = open(stringParams['scratchdir'] + "CG_colvars.out", "w")
        writeColVarsToFile(colvars_CG_new, ofile_colvars_CG, itr)
        ofile_colvars_CG.close()
        
        #------------DETERMINE CORRECTION TERM---------------------
        colvars_diff = getDiff(colvars_new, colvars_CG_new, colvars_diff, float(stringParams['Delta']))
        # apply geometric averaging.
        #for image in range(0, len(colvars_diff), 1):
            #for cv in range(0, len(colvars_diff[image]), 1):
                #print colvars_diff[image][cv][-1]
                
        colvars_diff, w_sum, w_norm = geomAvg(colvars_diff, stringParams, w_sum, w_norm)
        
        #ofile_colvars_diff = open(stringParams['scratchdir'] + "diff_colvars" + str(itr) + ".out", "w")
        #writeColVarsToFile(colvars_diff, ofile_colvars_diff, itr)
        #ofile_colvars_diff.close()
        
        ofile_colvars_CG_inner = open(stringParams['scratchdir'] + "CG_colvars_inner.out", "w")
        #-------------PROPOGATE INNER LOOP-------------------------
        for itr_inner in range(0, int(stringParams['niter_inner']), 1):
            
            writeColVarsToFile(colvars_CG, ofile_colvars_CG_inner, itr_inner)
            
            #------------CG ITERATION-----------------
            print "Iter", itr, ": Starting inner loop iteration " + str(itr_inner)
            print "Iter", itr, ": Writing restraint file."
            
            # generate the interpolation
            colvars_interp_CG = genInterp(colvars_CG, colvars_old_CG, int(stringParams['nDragSteps_CG']))
            # make a copy of the current colvars location
            colvars_old_CG = copyColvars(colvars_CG, colvars_old_CG)
            
            for step in range(0, int(stringParams['nDragSteps_CG']), 1):    
                for i in range(0, int(stringParams['nimages']), 1):
                    # setting input files for constrained CG run
                    CHARMM_writeColVarsRestraint(colvars_interp_CG[step][i], stringParams['restraint_CG'], stringParams['imagedir'] + str(i) + "/charmm_colvars_img" + str(i) + ".restraint")
                
                print "Iter", itr, ": Submitting restraint batch step " + str(step) + "."
                # call constrained CHARMM run
                callCHARMMConstraint(SimulParam, stringParams)
            
            for i in range(0, int(stringParams['nimages']), 1):
                # setting input files for release CG run
                CHARMM_writeColVarsRestraint(colvars_CG[i], stringParams['weakRestraint_CG'], stringParams['imagedir'] + str(i) + "/charmm_colvars_img" + str(i) + "_weak.restraint")
            
            print "Iter", itr, ": Submitting release batch."
            # call unconstrained CHARMM run
            callCHARMM(SimulParam, stringParams)
        
            # now collect the colvar data from the release phase.
            for i in range(0, int(stringParams['nimages']), 1):
                #print "Before Read " + str(i) + ": \n", colvars[i]
                print "Iter", itr, ": Reading colvar traj " + str(i)
                colvars_CG[i] = CHARMM_readColVarsTraj(stringParams, colvars_CG[i], stringParams['imagedir'] + str(i) + "/charmm_img" + str(i) + ".colvars.traj")
                #print " \nAfter read " + str(i) + ": \n", colvars[i]
            
            #----------REPARAMETERIZE CG STRING VARIABLES---------------
            
            for i in range(0, int(stringParams['reparItr']), 1):
                print "Repar Iteration: " + str(i)
                smoothString(float(stringParams['kappa']), colvars_CG)
                colvars_CG = reparString(colvars_CG)
                printRMSD(colvars_CG) 
                
            #-----------ADD DIFFERENCE TERM----------------------
            colvars_CG = addDiff(colvars_CG, colvars_diff, float(stringParams['Delta']))
        
            #-------------END OF INNER LOOP----------------------
        
        ofile_colvars_CG_inner.close()
        
        #--------------SET FG TO CG POSITION-----------------------
        colvars = copyColvars(colvars_CG, colvars)
        
        #-----------------WRITE COLVARS-------------------------
        # save out trajectory files for image steps. 
        if itr % int(stringParams['backupFreq']) == 0:
            saveTrajFiles(stringParams, itr)
        


writeColVarsToFile(colvars, ofile_colvars, int(stringParams['niter']))
# close output files.    
ofile_colvars.close()

end_time = time.time()
total_time = end_time - start_time
print "The total elapsed time is: " + str(total_time) + " sec."

print "Done."
