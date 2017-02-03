# -*- coding: utf-8 -*-
"""
This module unit tests the partition class. 
 
Specifically, this test suite makes sure that the partition acts like a callable list in Python. The test suite will cover the following behaviors: 
 
* partition object satisfies the basic behavior of a list 
 
* partition correctly enforces that elements of the partition are callable and return a float object 
 
* partition's call routine correcty returns a numpy array of floats normalized over the return value of it's elements call routines  
 
"""

import pytest 
from est import neus 
 
@pytest.fixture 
def setup(): 
    """Sets up an empty partition object used for further testing of the partition module. 
    """ 
 
    # construct some callable test objects that are not window objects.  
 
def test_adding_windows(): 
    """Test addition of window objects and safety assertions that windows are callable. 
    """ 
 
def test_partition_call(): 
    """Test call routine of partition by passing a callable object to partition.  
    """ 
 
    # add some callable test objects. 
 
def test_list(): 
    """Testing core routines that implement list-like behavior in partition. 
    """ 
 
    # test the pop and sequence behaviors 
 
    # test the index behaviors 
