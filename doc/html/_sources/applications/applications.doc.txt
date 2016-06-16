.. Documentation for applications folder. Will cover the 

Applications
=================

Nonequilibrium Umbrella Sampling
--------------------------------

The Nonequilibrium Umbrella Sampling (NEUS) application module provides a flexible set of tools built on top of the Walker API that implements the found in **Tempkin et at, in prep**. 

The NEUS application tools developed here include three things

* an module called Window that implements the data structures and routines that represent a single spatiotemporal discretization in the NEUS scheme. In the paper, this would correspond to a single value of the J(t) process.
* a module called partition that represents a collection of these windows.
* a set of routines for solving the affine eigenvalue problem.

DEV: Put here a basic overview of the types of calculations one can do using this code base with a reference to the testable code in the examples/ subdirectory.

* A summary of some of the results or a summary of a basic implementation could be really useful here. 
* Include a LAMMPS walker implementation of NEUS on the dipeptide example (should illustration the usage model in a real application)

Partition Module Basic Usage
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

.. include:: neus/partition-basic-usage.rst

.. include:: neus/Untitled.rst

Partition Class
##################

.. autoclass:: partition.partition
    :members: __init__, __call__, index, append

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

Pyramid Class
################

.. automodule:: pyramid
    :members:

The Box module implements a support function for the Window class of the form:

**Box indicator function**.

Entry Point Module
<<<<<<<<<<<<<<<<<<<<

.. automodule:: entryPoints
		:members:
