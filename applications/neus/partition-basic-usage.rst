
The partition module provides a simple structure for handling and
manipulating a list of windows. The partition module is effectively a
callable Python list.

To initialize a partition object:

.. code:: python

    import partition
    import numpy as np
    sys = partition.partition()

One can add elements to the partition one by one:

.. code:: python

    # import a window object with a pyramid support
    import pyramid
    win = pyramid.Pyramid([0.0], [1.0])
    
    sys.append(win)
    
    print sys


.. code::

    partition([Pyramid([ 0.], [ 1.]), Pyramid([ 0.], [ 1.]), Pyramid([ 0.], [ 1.])])


Or hand partition a list at initialization:

.. code:: python

    windows = [pyramid.Pyramid(i,[1.0]) for i in np.linspace(0.0, 10.0, 11)]
    sys = partition.partition(windows)
    print sys


.. parsed-literal::

    partition([Pyramid(0.0, [ 1.]), Pyramid(1.0, [ 1.]), Pyramid(2.0, [ 1.]), Pyramid(3.0, [ 1.]), Pyramid(4.0, [ 1.]), Pyramid(5.0, [ 1.]), Pyramid(6.0, [ 1.]), Pyramid(7.0, [ 1.]), Pyramid(8.0, [ 1.]), Pyramid(9.0, [ 1.]), Pyramid(10.0, [ 1.])])


Further, partition supports the addition operation for merging two
partition objects:

.. code:: python

    #union = sys1 + sys2
    #print union

