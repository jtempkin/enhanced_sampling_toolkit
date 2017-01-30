# -*- coding: utf-8 -*-
"""
A profile for comparing the direct dynamics call. 
"""

import lammpsWalker
import dynamics 
import os

nsteps = 1000000
stepLength = 5
temperature = 310.0
damping_coefficient = 30.0
scratchdir = os.path.abspath("data")
inputFilename = os.path.abspath("data/input.enk")
logFilename = scratchdir + "/log.lammps"

os.chdir(scratchdir)

# lets instantiate a walker object to sample this window.

wlkr = lammpsWalker.lammpsWalker(inputFilename, logFilename, index=0, debug=False)

# set langevin dynamics object and pass it to the walker.
dyn = dynamics.langevin(temperature, damping_coefficient, shake = True)            
wlkr.setDynamics(dyn)
wlkr.setTimestep(1.0)
wlkr.drawVel()

wlkr.propagate(nsteps)
