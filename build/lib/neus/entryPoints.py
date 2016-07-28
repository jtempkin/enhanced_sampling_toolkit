# -*- coding: utf-8 -*-
"""
This module acts as a template for defining named tuples that store entry points.

It provides a consistent entry point contruction as a named tuple.

It's provided as an independent file so that analysis and manipulation of these objects can be possible.
"""

import collections

entry_point = collections.namedtuple('entry_point', ['q', 'p', 'ref_q', 'ref_p', 'time', 'cv'])
