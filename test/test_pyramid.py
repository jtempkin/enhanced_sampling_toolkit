# -*- coding: utf-8 -*-
"""
Pytest module for unit testing the Pyramid subclass of window. 

Specifically, this module test the call function implementation of the Pyramid subclass.

* test that the call function correctly implements a pyramid of a given width and height in spatical coordinates
* tests the call function generalizes to multiple coordinates and correctly returns the call function output
* tests the function of the indicators on the time-coordinate
* tests the function of the indicators on a reference time-coordinate
"""

import pytest
import numpy as np
from est import neus


@pytest.fixture
def setup():
    """Sets up some structures and lists useful for the remaining test.
    """


def test_initialization():
    """This routine tests that the initialization routine returns a BF object with the correct parameters.s
    """

    window = neus.window.Window([1.0, 1.0], [0.1, 0.1])

    assert isinstance(window, neus.window.Window), "Object returned did not match expected type."

    assert np.array_equal(window.center, [1.0, 1.0]), "Center was not set correctly."

    assert np.array_equal(window.width, [0.1, 0.1]), "Widths was not set correctly."

    # assert np.array_equal(), "Time not set correctly."


def test_fluxlists():
    """Test the behavior of the routines associated with the flux lists.
    """
