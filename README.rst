Enhanced Sampling Toolkit
---------------------------------------------------
Primary Author: Jeremy O. B. Tempkin
Created: 5/8/2014

Introduction 
*************************** 
 
The Enhanced Sampling Toolkit is a Python package designed to provide a flexible and extensible toolkit for rapidly prototyping enhanced sampling algorithms for use in molecular simulations. The codebase is written entirely in Python and acts as a wrapper to various well-established molecular dynamics engines.

This codebase is developed in order to provide a flexible and extensible toolkit for rapidly prototyping enhanced sampling algorithms for use in molecular simulations. The codebase is written entirely in Python and acts as a wrapper to various well-established molecular dynamics engines. The design of this toolkit prioritizes facilitating rapid code development at the algorithm level. This priority is what motivates the decision to use Python as the primary language. This decision also facilitates simple connections between this code and well-developed Python tools already established. 

This code consists of two components. The first is what we call the Walker API. This part consists of an API specification that abstracts the interaction between an enhanced sampling algorithm and the underlying MD engine that performs the integration. The Walker API facilitates the rapid construction of higher-level algorithm level code built on top of these basic interactions. Implementations built on this API have the benefit of rapid portability between various MD codes that implement the models and dynamics a user may want by a single algorithm code. Ideally, since the API consists of fairly high-level interactions, this allows the programmer to focus solely on the algorithm design in a manner agnostic to the specifics of the MD engine. 

Structure of the toolkit 
******************************* 

Installation 
**************************** 
 
To install the toolkit. 
 
To interact with  
 
Gitpages 
************************** 
 
 
System Requirements 
*************************** 

1) Python v2.7.x (v2.7.9 or higher is required.)

2) LAMMPS distribution built as a shared library that is importable by the Python interpreter. To install LAMMPS see: http://lammps.sandia.gov/doc/Section_python.html

3) mpi4py distribution for MPI parallelization.

4) h5py, numpy, scipy are required but are easily available in most standard Python distributions. 

Package contents 
*************************** 

est - The main source code for this project. This contains the core API of the toolkit in the walker submodule. It also contains various libraries that provide high-level tools for algorithm development. 

docs - documentation for the project. This folder serves to hold the build targets for the Sphinx builds. 

profiling - source for some simple profiling of the application libraries. 

test - unit testing suit for the toolkit. 