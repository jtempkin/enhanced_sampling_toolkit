
The NEUS application toolkit provides three modules:

-  A Window object
-  A partition object
-  An Entry Point collection

The Window object represents the data strutures and routines that
represent a single spatiotemporal restriction, or a single "window" in
your NEUS scheme. In the language of NEUS, this corresponds to a single
restricted value of the :math:`J^{(t)}` process.

The partition object represents a callable list of Window instances.
Again, in the language of NEUS, it represents the full index space of
the :math:`J^{(t)}` process.

The entry point module provides a named tuple structure for storing the
state of the system at a particular point in phase space.

In the following, we will illustrate some basic usage of these objects
with the idea of familiarizing the user with their tools and syntax.

.. code:: python

    %matplotlib inline
    from matplotlib import pyplot as plt
    import numpy as np
    from neus import pyramid
    from neus import partition 

Window objects
^^^^^^^^^^^^^^

The NEUS application in the toolkit provides a parent module called
"window". The shape of the support of each window can be customized by
writing a class that inherents from the window parent class.

For example, the NEUS toolkit provides a pyramidal shaped module named
neus.pyramid.

.. code:: python

    # instantiate a window object
    center = [1.0]
    width = [0.5]
    win = pyramid.Pyramid(center, width)
    
    # plot the support of the pyramid object
    x = np.linspace(0.0, 2.0, 100)
    out = [win([i]) for i in x]
    plt.plot(x, out)




.. parsed-literal::

    [<matplotlib.lines.Line2D at 0x10f399c10>]




.. image:: output_3_1.png


Partition objects
^^^^^^^^^^^^^^^^^

The NEUS toolkit provides the partition class that represents a
collection of window objects. It's implementation is effectively a
callable list, whose call function returns a normlized vector of the
supports of it's elements. Note that it enforces that each element of a
partition must itself be a callable object.

.. code:: python

    # create an empty partition object
    sys = partition.Partition()
    
    # create a list of winodw objects
    width = 0.75
    
    centers = [x for x in np.arange(-3, 3)]
    
    # now add the windows to partition
    for center in centers:
        win = pyramid.Pyramid([center], [width])
        sys.append(win)
    
        
    print sys


.. parsed-literal::

    partition([Pyramid([-3], [ 0.75]), Pyramid([-2], [ 0.75]), Pyramid([-1], [ 0.75]), Pyramid([0], [ 0.75]), Pyramid([1], [ 0.75]), Pyramid([2], [ 0.75])])


Partition objects act like callable lists so some of the expected list
operations also work with partitons:

.. code:: python

    # element access
    print sys[0]
    
    # slicing
    print sys[2:4]
    
    # list concatenation
    print sys[2:4] + [sys[0]]


.. parsed-literal::

    Pyramid([-3], [ 0.75])
    [Pyramid([-1], [ 0.75]), Pyramid([0], [ 0.75])]
    [Pyramid([-1], [ 0.75]), Pyramid([0], [ 0.75]), Pyramid([-3], [ 0.75])]


The partition's call routine returns the normalized vector of supports:

.. code:: python

    # choose a point with nonzero support
    cv = [-2.5]
    
    # return normalized vector of the support
    print sys(cv)




.. parsed-literal::

    array([ 0.5,  0.5,  0. ,  0. ,  0. ,  0. ])


