# -*- coding: utf-8 -*-
"""
pytest module for testing basis function module. 
"""

import pytest
import numpy as np

@pytest.fixture
def setup():
    """
    Sets up some structures and lists useful for the remaining test. 
    """


def test_initialization():
    """
    This routine tests that the initialization routine returns a BF object with the correct parameters. 
    """
    import basisFunctions_neus_dipeptide
    window = basisFunctions_neus_dipeptide.Box([1.0,1.0],[0.1,0.1])
    
    assert isinstance(window, basisFunctions_neus_dipeptide.Box), "Object returned did not match expected type."
    
    assert np.array_equal(window.center,[1.0, 1.0]), "Center was not set correctly."
    
    assert np.array_equal(window.width,[0.1, 0.1]), "Widths was not set correctly."
    
    
def test_indicator():
    """
    This tests the indicators of the various shapes. 
    """
    