# -*- coding: utf-8 -*-
"""
This module acts as a template for defining named tuples that store the state of a system at a particular point.

The named tuple provide the following fields:

* q - the position vector of the system at time t
* p - the velocity vector of the system at time t
* ref_q - the position vector of the system at previous time t'
* ref_p - the velocity vector of the system at previous time t'
* time - the time of the system
* cv - the value of the collective variable at time t

The NEUS module provides the entry_point definition as an independent module so that analysis and manipulation of these objects is possible.
"""

import collections

entry_point = collections.namedtuple('entry_point', ['q', 'p', 'ref_q', 'ref_p', 'time', 'cv'])
