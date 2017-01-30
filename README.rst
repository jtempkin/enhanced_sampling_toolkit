README FILE FOR THE ENHANCED SAMPLING TOOLKIT
---------------------------------------------------
Primary Author: Jeremy Tempkin
Created: 5/8/2014

Contributing Authors:
Seyit Kale
Erik Thiede

Description:

This codebase is developed in order to provide a flexible and extensible toolkit for rapidly prototyping enhanced sampling algorithms for use in molecular simulations. The codebase is written entirely in Python and acts as a wrapper to various well-established molecular dynamics engines. The design of this toolkit prioritzes facilitating rapid code development at the algorithm level. This priority is what motivates the decision to use Python as the primary language. This decision also facilitates simple connections between this code and well-developed Python tools already established. 

This code consists of two components. The first is what we call the Walker API. This part consists of an API specification that abstracts the interaction between an enhanced sampling algorithm and the underlying MD engine that performs the integration. The Walker API facilitates the rapid construction of higher-level algorithm level code built on top of these basic interactions. Implementations built on this API have the benefit of rapid portability between various MD codes that implement the models and dynamics a user may want by a single algorithm code. Ideally, since the API consists of fairly high-level interactions, this allows the programmer to focus solely on the algorithm design in a manner agnostic to the specifics of the MD engine. 

System Prerequisits:

1) Python v2.7.x (v2.7.9 or higher is prefered)

2) LAMMPS distribution built as a shared library. To install LAMMPS see: http://lammps.sandia.gov/doc/Section_python.html

3) mpi4py distribution for MPI parallelization.

Folders:

walker_api - The main source code for this project. This contains the core API of the toolkit in the *Walker.py files.

doc - documentation for the project. This folder serves to hold the build targets for the Sphinx builds. 

applications - container for application modules built using the Walker API toolkit.

test - unit testing for the toolkit. 