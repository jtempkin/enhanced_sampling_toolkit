.. Documentation for applications folder. Will cover the 

Applications
=================

Nonequilibrium Umbrella Sampling
--------------------------------

The Nonequilibrium Umbrella Sampling (NEUS) application module provides a flexible set of tools built on top of the Walker API that implements the found in **Tempkin et at, in prep**. 

More detailed examples of the specific implementations of the use of the can be found in the examples subdirectory. 

Partition Module
<<<<<<<<<<<<<<<<<<<

The partition module provides a simple structure for handling and manipulating a list of windows. The partition module is effectively a callable Python list. As such, it supports. 

To initialize a partition object:

>>> from est.neus import partition
>>> import numpy as np
>>> sys = partition.partition()

One can add elements to the partition one by one:

>>> 

Or hand partition a list at initialization:

>>> from est.neus import Box
>>> windows = [Box.box(i,j) for i in range(5) for j in range(10)]
>>> sys = partition.partition(windows)
>>> print sys

Note that partition requires the elements it contains to be callable items. Partition will enforce this behavior:

>>> sys.append("not callable")

See below for a more detailed specification.

.. automodule:: partition
    :members: 

Window Module
<<<<<<<<<<<<<<<<<<<<<<<

The window object defines a single spatiotemporal restriction of the discretization defined by the NEUS J(t) process. An instance of the window object corresponds to a single value of the J. The object defines a data structure for storing the entry point flux lists defined by ~pi. It also defintes the routines require for constructing the flux lists and returning physically weighted elements. 

To import and initialize a window object:

>>> from est.neus import Window
>>> win = Window.window(center=[], width=[])
>>> print "Hello"

To 

.. automodule:: window
		:members:

The Pyramid module implements a support function for the Window class of the form:

**equation**.

.. automodule:: pyramid
    :members:

The Box module implements a support function for the Window class of the form:

**Box indicator function**.

Entry Point Module
<<<<<<<<<<<<<<<<<<<<

.. automodule:: entryPoints
		:members:
