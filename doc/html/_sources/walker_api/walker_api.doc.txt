.. Documentation for src folder. Leverages the autodoc functionallity
   for formatting the lammpsWalker docstrings. Uses sphinx to compile the
   documentation. 

Walker API
========================

Introduction
-------------

The central design principle of the Walker API is to provide an abstracted interface between the "algorithm" level code and the code that implements the underlying dynamical model. This section describes the base API used to define the walker object and the API routines. The implementations of the bindings specific to the dynamics packages are registered to this base class definition through direct subclassing of the modules. A key feature of this decision is that it enforces the implementation of the dynamics bindings to implement each of the following methods in their class declarations. However, many of the methods can be overridden with empty methods if an implementation is incomplete or does not require all of the base class methods. 


.. automodule:: walker_base
        :members: 

LAMMPS Walker Module
-------------------------

This module implements the Walker API for the LAMMPS MD engine. See walker_base.py for a specification of the API. For details concerning the usage of the LAMMPS MD package see the excellent documentation at the LAMMPS webpage:

http://lammps.sandia.gov/

In particular, you may want to see how the Python wrapper to LAMMPS on which this implementation is based:

http://lammps.sandia.gov/doc/Section_python.html

Here we will outline basic usage guides for the walker API usage in LAMMPS.

OpenMM Walker Module
-------------------------

We plan to implement an OpenMM walker API module in a future release. 