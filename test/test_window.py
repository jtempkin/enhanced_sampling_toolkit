# -*- coding: utf-8 -*-
"""
Pytest module for unit testing the window module.
"""

import pytest
import numpy as np
import neus


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
