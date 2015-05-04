# -*- coding: utf-8 -*-
"""
This is a root Python module for conducting profiling of the enhanced sampling toolkit. The profiling tests the set of primary components of the walker classes and the units of the application sets by running simple test systems. We'll set up a series of different systems of varying size and complexity to measure how the code is scaling with system size and complexity. The profiling information is generated with Python's own cProfile module. 
"""

import sys
import cProfile
import pstats

with open(sys.argv[1], "r") as f_handle:
    cProfile.run(f_handle, sys.argv[1] + ".pstats")
    
p = pstats.Stats(sys.argv[1] + ".pstats")
p.strip_dirs().sort_stats("time").print_stats()
    
