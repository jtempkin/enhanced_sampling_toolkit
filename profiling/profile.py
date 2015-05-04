# -*- coding: utf-8 -*-
"""
This is a root Python module for conducting profiling of the enhanced sampling toolkit. The profiling tests the set of primary components of the walker classes and the units of the application sets by running simple test systems. We'll set up a series of different systems of varying size and complexity to measure how the code is scaling with system size and complexity. The profiling information is generated with Python's own cProfile module. 
"""

import sys
import cProfile
import pstats

if sys.argv[1] == 'US':
    cProfile.run("profile_US.py", "profile_US")
    p = pstats.Stats("profile_US")

p.strip_dirs().sort_stats("time").print_stats()
    
