"""

Or hand partition a list at initialization:

>>> from est.neus import Box
>>> windows = [Box.box(i,j) for i in range(5) for j in range(10)]
>>> sys = partition.partition(windows)
>>> print sys

Note that partition requires the elements it contains to be callable items. Partition will enforce this behavior:

>>> sys.append("not callable")

See below for a more detailed specification.

Here we will show some basic usage of the partition class.

.. code:: python

    import numpy as np
    import window

.. code:: python

    wins = [i for i in range(10)]
    print wins


.. parsed-literal::

    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
"""

import numpy as np
import errors


class Partition:
    def __init__(self, umbrellas=[]):
        """Construct an instance of the partition.
        """

        self._umbrellas = umbrellas

        # let's check to see if every object in the windows is callable.
        for win in self._umbrellas:
            if not hasattr(win, "__call__"):
                raise Exception("Items in partition must be callable.")

        return None

    def __getitem__(self, key):
        """Return an item in the list.
        """

        return self._umbrellas[key]

    def __setitem__(self, key, value):
        """Set item to value.
        """

        if not hasattr(value, '__call__'):
            raise Exception("Items in partition must be callable.")

        self._umbrellas[key] = value

        return None

    def index(self, item):
        """Return the first index in which item appears in partition.

        Extended description.

        Parameters
        ------------
        item : obj
            The window object to check for membership.

        Returns
        ------------
        index : int
            The index of item in the list. Returns ValueError if item is not
            contained in the partition.
        """
        return self._umbrellas.index(item)

    def __delitem__(self, key):
        """Delete window.
        """
        del self._umbrellas[key]

        # note that we really should be careful about how to ensure consistency in the window data structures here. If we delete a window in the middle of a list, we may offset the indicies of the keys for the fluxlists in other remaining windows.
        # We may consider that to alter the data structure of the windows at this level will require some checks that the objects in partition can support those operations.

        return None

    def __len__(self):
        """Return the number of windows in the partition.
        """

        return len(self._umbrellas)

    def __contains__(self, item):
        """Return boolean for if the partition contains the specified item.
        """

        return item in self._umbrellas

    def __iter__(self):
        """Return an iterator for the list of basis functions.
        """

        return iter(self._umbrellas)

    def __repr__(self):
        """Return string representation of the partition.
        """

        return "partition(" + repr(self._umbrellas) + ")"

    def append(self, value):
        """Append a new value to the partition.

        Note that this enforces that the input object is callable. Will throw an exception if value is not callable.

        Parameters
        ------------
        value : obj
            A callable Python object.

        Returns
        -------------
        None
        """

        if not hasattr(value, '__call__'):
            raise Exception("Items in partition must be callable.")

        self._umbrellas.append(value)

        return None

    def __call__(self, pos):
        """Return a numpy array of the normalized supports of the windows at the given.
        """
        # build an array
        indicators = np.asarray([win(pos) for win in self])

        # raise a support error if the coordinates provided have zero support on the partition.
        if np.sum(indicators) == 0.0:
            raise errors.SupportError("The coordinates provided to partition instace had zero support on the partition.")

        # return the normalized values of the windows as an array.
        return indicators / indicators.sum()
