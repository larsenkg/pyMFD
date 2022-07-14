Usage
=====

Installation
------------
pyMFD was tested with Python 3.7.11 and has the following requirements::

    numpy      >= 1.21.2
    scipy      >= 1.7.1
    matplotlib >= 3.5.0
    bottleneck >= 1.3.2

The following are required to build the documentation::

    pytest   >= 7.1.2
    sphinx   >= 4.4.0
    numpydoc >= 1.4.0

    (optional)
    pydata-sphinx-theme >= 0.9.0  

The most up-to-date version of pyMFD is available on Github:

.. code-block:: console

    $ git clone https://github.com/larsenkg/pyMFD.git

pyMFD can also be installed using pip:

.. code-block:: console

    $ pip install pyMFD

Usage
-----
pyMFD consists of the following modules:

 - :mod:`pyMFD.nanoscope` for loading NanoScope AFM scan files
 - :mod:`pyMFD.scan_params` for loading JSON scan parameter files
 - :mod:`pyMFD.FV` for working with force-volume data
 - :mod:`pyMFD.summarize` for summarizing force ramps into e.g. compliance maps
 - :mod:`pyMFD.cantilever` for analyzing (micro)cantiliver compliance

Example
-------
See :mod:`example`.