# -*- coding: utf-8 -*-
"""
Pytest module for unit testing the window module.

Specifically, this module tests the following functionality: 
 
* tests initialization of window objects and yields the correct sequence of parameters 
 
* tests the initilization of the initial conditions buffers and the initial condition probability parameters  
 
* tests the routines implementing the len() function 
 
* tests the representation function 
 
* tests routines that yield the initial conditions buffers 
 
* tests the routines that set and return the flux array  
 
* tests the routines that set and return the flux lists  
 
* tests the correct initialization of the flux list data structure 
 
* tests the routines that return entry point objects from the flux list data structures
"""

import pytest
import numpy as np
from est import neus
from est import entry_points


@pytest.fixture
def setup():
    """Sets up some structures and lists useful for the remaining test.
    """
    
    # let's create a window instance
	window = neus.window.Window([1.0, 1.0], [0.1, 0.1])
	
	# we're going to add an empty window
	empty_window = neus.window.Window([1.0, 1.0], [0.1, 0.1])
	
	# create an entry point for testing purposes
	entry = entry_points.entry_point( q=  [1.0, 1.0], 
                                      p= [1.0, 1.0],
                                      q_ref= [1.0, 1.0],
                                      p_ref= [1.0, 1.0],
                                      time= 0,
                                      cv= [1.0])
   
   # add an entry point to window                    
   window.add_entry_point()
                                      
                                      

def test_initialization():
    """This routine tests that the initialization routine returns a BF object with the correct parameters.s
    """

    #window = neus.window.Window([1.0, 1.0], [0.1, 0.1])

    assert isinstance(window, neus.window.Window), "Object returned did not match expected type."

    assert np.array_equal(window.center, [1.0, 1.0]), "Center was not set correctly."

    assert np.array_equal(window.width, [0.1, 0.1]), "Widths was not set correctly."

    # assert np.array_equal(), "Time not set correctly."


def test_fluxlists():
    """Test the behavior of the routines associated with the flux lists.
    """

	# test adding entry points to the flux lists
	

def test_window_len():
	"""Test the behavior of the window length function.
	"""
	

def test_window_bool():
	"""Test the truth value of the window.
	"""

def test_window_print():
	"""Test the string __repr__ functionality of the window object
	"""
