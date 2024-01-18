============
Installation
============


To begin using ``LARRY``, install from PYPI


Install via pip (recommended)
"""""""""""""""""""""""""""""
.. code-block:: python

    pip install LARRY-dataset


GitHub (Developer version)
""""""""""""""""""""""""""

To access the latest version of ``LARRY`` from GitHub, clone the 
repository, ``cd`` into the project's root directory and install the
editable version.

.. code-block:: python

    # Install the developer version via GitHub
    
    git clone https://github.com/scDiffEq/LARRY-dataset.git; cd ./LARRY-dataset;
    pip install -e .


Troubleshooting
"""""""""""""""

The ``pykeops`` library creates a cache. Sometimes, when you switch devices
though retain the same disc (common when using a VM, for example), this cache
will no longer be compatible with the installed drivers for that device. To
clear and rewrite this cached, we can perform the following:

.. code-block:: python

    import pykeops

    pykeops.clean_pykeops()