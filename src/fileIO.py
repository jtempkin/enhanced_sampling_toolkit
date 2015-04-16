# -*- coding: utf-8 -*-
"""
This file contains a number of file input/output routines. 

The following routines provide a number of standalone tools for processing common input/output tasks one may need in scripting enhanced sampling algorithms. 
"""

import os.path
import shutil
import subprocess
import numpy as np
import re

def readSysInput(filename, sysParams, trustMe = False):
    """
    This function reads in  system input file that specifies global 
    parameters for the umbrella sampling code.  It takes as input: 
        filename:   the name of the file to be read (a string)
        sysParams:  the dictionary which the parameters will be added to
        trustMe:    a boolean which specifies if the method should trust that
                    the user provided sensible input.

    
    Please note that this is a destructive method, since it modifies sysParams directly.
    """
    
    if not os.path.exists(filename):
        raise IOError("The input file was not found.")
        
    with open(filename, "r") as ifile:
        # We define a capture string, which looks for parts of the input file with
        # the form "parameter:  value".
        pattern = re.compile(r'(\w+)\s*\:\s*(.*)')
        # We find all instances of the capture string, and make them into a dict.
        tempDict = dict(pattern.findall(ifile.read()))
        

    # We parse the input files. This converts values to the correct type
    parseInput(tempDict)
    
    # add in the newly parsed dictionary
    sysParams.update(tempDict)
    
    # only call the check function if we don't trust the user. silly user. 
    if trustMe == False:
        checkInput(sysParams, trustMe)
    else:
        print "Input parser is naively trusting the input file."
    
    return 0
    
def parseInput(dictionary):
    """
    Parses the input data into the correct filetypes.
    It takes as input a dictionary: inputDict.
    """
    # The following lists are lists of keywords to be interpreted as various data types.
    
    # Keywords to be parsed as integers
    intKeywords  = ["nwalkers", "walkerSteps", "stepLength"]
    
    # Keywords to be parsed as paths to directories
    filepathKeywords = ["wkdir", "datadir"]
    
    # for floats 
    # floatKeywords = []
    
    # now convert integers
    for key in intKeywords:
        if key in dictionary:
            dictionary[key] = int(dictionary[key])
    
    # convert filepaths
    for key in filepathKeywords:
        if key in dictionary:
            dictionary[key] = dictionary[key].rstrip("/")
            
    # We handle more complex entries separately.
    # parse the list of cv's that are periodic
    if "wrapping" in dictionary:
        dictionary["wrapping"] = np.array(map(float, dictionary["wrapping"].split(',')))    
    
    # parse ranges arrays for periodic cv's 
    if "cvrange" in dictionary:
        dictionary["cvrange"] = np.array([map(float,entry.split(",")) for entry in dictionary["cvrange"].split(";")])
        dictionary["ncells"] = int(np.product(dictionary["cvrange"][:,2]))
 
       
def checkInput(inputDict, trustMe):
    """
    Parses the input strings into information usable by the program.
    For instance, it will convert pure numbers to floats, or comma separated values
    to lists.  It also checks if the input provided is sensible, if trustMe is false.  
    It takes two inputs:
        inputDict       A dictionary of input keywords
        
    This code is currently not finished.
    """
    
    try: 
        wktyp = inputDict["walkerType"]
    except KeyError: 
        if trustMe == False:
            raise IOError("No walker type specified.")
    if wktyp == "lammps":
        if "datadir" not in inputDict and trustMe == False:
            raise IOError("You have not provided a directory with the LAMMPS data")
    else:
        if trustMe == False:
            raise IOError("Walker type not recognized.  Please include a line that specifies the type of walker, e.g. walkerType: lammps")
        
    # THIS IS NOT DONE!!!!
    raise NotImplementedError("Dude, the input checker has not been fully implemnted!")
    return 0

    
def makeWkdir(dirname):
    """
    This function creates the scratch directory for the intermediate MD files. The default behavior is to check to see if the specified directory exists and if so, move it to a backupfile using the name plus ".backup" extension. *** THIS BEHAVIOR SHOULD CHANGE SINCE THERE IS STILL THE POSSIBILITY OF DATA LOSS HERE***
    """
    if os.path.exists(dirname):
        # print some warnings that the target dictory exists
        print "Working Directory already exists."        
        print "Making backup of old working directory."
        if os.path.exists(dirname + ".backup"):
            shutil.rmtree(dirname + ".backup")
        shutil.move(dirname, dirname + ".backup")
        subprocess.call(['mkdir', dirname])
    else:
        subprocess.call(['mkdir', dirname])
        
    return 0

def writeMat(mat, filename, binary=False, log=False):
    """
    Writes the given 2D matrix to the filename. 
    
    If binary arg is set to True, writes it as a numpy binary array at the chosen filename. 
    
    Acts effectively as a thin wrapper to numpy.save / numpy.savetxt routines. 
    
    We'll add some HDF5 support in future releases. 
    """
    if binary:
        if log:
            np.save(filename, np.log(mat.real))
        else:
            np.save(filename, mat.real)
    else:
        if log:
            np.savetxt(filename, np.log(mat.real))
        else:
            np.savetxt(filename, mat.real)
    
    return 0

def writeColvars(colvars, filename, binary=False):
    """ 
    This function writes a data file to filename consisting of the provided 
    list of lists given in colvars. 
    """
    if binary:
        np.save(filename, colvars)
    else:
        np.savetxt(filename, colvars)
    """
    ofile = open(filename, "w")
    # loop over the cvs
    for cv in range(0, len(colvars), 1):
        ofile.write("%lf  %lf \n" % tuple(colvars[cv]))
    
    ofile.close()
    """
    return 0
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
