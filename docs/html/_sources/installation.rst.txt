Installation
=================

To install the Enhanced Sampling Toolkit on your local machine, download the package from the `Github page <https://github.com/jtempkin/enhanced_sampling_toolkit>`_.

There is a simple Makefile implementation for installing the package with pip. To install the package using the Makefile, run:

.. code:: 

    make install

from the root of the toolkit package. If you plan on developing the files in the toolkit and want to install a soft-linked version that points back to the local source, run

.. code:: 

    make install-dev

Note that this will allow changes made locally to the package files to visible in your active Python environment. To uninstall the package:

.. code::

   make uninstall

If you'd like to reinstall the package for whatever reason, run:

.. code::

   make reinstall

which will run the uninstall and install lines. 
    