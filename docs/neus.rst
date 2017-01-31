.. Documentation for NEUS modules 

Nonequilibrium Umbrella Sampling (NEUS)
==========================================

.. automodule:: jester.neus
    :members:

The NEUS toolkit provided in this package contains three modules:

* A module called Window that implements the data structures and routines that represent a single spatiotemporal discretization in the NEUS scheme. An instance of the Window object corresponds to a single value of the :math:`J^{(t)}` process.
* A module called partition acts as a collection of these windows objects that expresses the full :math:`J^{(t)}` index space.
* A module called entry points that provides a named tuple for storing phase space points as elements in the :math:`\tilde \gamma_{ij}`

Below we describe the basic usage of the NEUS toolkit. Please see the Jupyter notebooks provided in the data folder for an interactive version.

NEUS Module Basic Usage
--------------------------------

.. include:: NEUS-basic-usage.rst

Partition Module API
-----------------------

.. autoclass:: jester.neus.partition.Partition
    :members: 
    :special-members:

Window Module API
--------------------

.. automodule:: jester.neus.window
    :members: 

.. autoclass:: jester.neus.window.Window
    :members:
    :special-members:


Pyramid Module API
--------------------
.. automodule:: jester.neus.pyramid
 
.. autoclass:: jester.neus.pyramid.Pyramid
    :members: __init__, __call__, ref_indicator, time_indicator


Entry Point Module
---------------------

.. automodule:: jester.neus.entry_points


