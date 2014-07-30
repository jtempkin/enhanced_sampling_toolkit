# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 15:34:23 2014

This file contains a number of file input/output routines. 

@author: jtempkin
"""

import os.path
import basisFunctions
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

def createUmbrellas(Params):
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
    colVarParams=Params["cvrange"]
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
            um.append(basisFunctions.Box(umbrellaCoord,widthlist,boxwrap))
        else:
            um.append(basisFunctions.Box(umbrellaCoord,widthlist))

    return um
    
def makeWkdir(sysParams):
    """
    This function creates the wkdir for the intermediate MD files.
    """
    if os.path.exists(sysParams['wkdir']):
        print "Working Directory already exists."        
        print "Making backup of old working directory."
        if os.path.exists(sysParams['wkdir'] + ".backup"):
            shutil.rmtree(sysParams['wkdir'] + ".backup")
        shutil.move(sysParams['wkdir'], sysParams['wkdir'] + ".backup")
        #subprocess.call(['cp', '-r', sysParams['wkdir'], sysParams['wkdir'] + ".backup"])    
        #subprocess.call(['rm', '-r', sysParams['wkdir']])
        subprocess.call(['mkdir', sysParams['wkdir']])
    else:
        subprocess.call(['mkdir', sysParams['wkdir']])
        
    return 0

def writeMat(mat, filename, binary=False, log=False):
    """
    Writes the given 2D matrix to the filename. 
    
    If binary arg is set to True, writes it as a numpy binary array at the chosen filename. 
    
    >>> import numpy as np
    >>> z = np.zeros((1,4))
    >>> writeMat(z,"./debug/z.test.out")
    
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
