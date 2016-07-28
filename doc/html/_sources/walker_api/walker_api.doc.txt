.. Documentation for src folder. Leverages the autodoc functionallity
   for formatting the lammpsWalker docstrings. Uses sphinx to compile the
   documentation. 

Walker API
========================

The Walker API is built to provide an abstract interface between the "algorithm" level code and the code that implements the underlying dynamical model. 

This section describes some of the basic usage of the implementation of the walker API for the LAMMPS package. The example pages are written from Jupyter notebooks included in the 'data/' directory for interactive use. The full API is documented in the last section.

LAMMPS Walker Module
----------------------

This module implements the Walker API for the LAMMPS MD engine. See walker_base.py for a specification of the API. For details concerning the usage of the LAMMPS MD package see the excellent documentation at the LAMMPS webpage:

http://lammps.sandia.gov/

In particular, you may want to see how the Python wrapper to LAMMPS on which this implementation is based:

http://lammps.sandia.gov/doc/Section_python.html

To use the LAMMPS walker module provided in the toolkit, you will need to have the LAMMPS MD engine available as an importable Python module. Please follow in the instructions provided by the above link to download, compile and install LAMMPS for Python. To check to see if LAMMPS is installed correctly, try importing it interactively as

.. code:: Python
    
    from lammps import lammps

If this import returns no errors, the LAMMPS Walker module should be able to import and use the compiled LAMMPS build on your machine. 

.. include:: LAMMPS-Walker-Basic-Usage.rst

LAMMPS Walker API
--------------------

.. automodule:: walker_api.lammps_walker
    :members:

