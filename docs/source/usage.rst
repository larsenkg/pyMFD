Usage
=====

Installation
------------
pyMFD was tested with Python 3.7.11 and has the following requirements::

    numpy      >= 1.21.2
    scipy      >= 1.7.1
    matplotlib >= 3.5.0
    bottleneck >= 1.3.2

The following is required to run the tests::

    pytest   >= 7.1.2

The following are required to build the documentation::

    sphinx   >= 4.4.0
    numpydoc >= 1.4.0

    (optional)
    pydata-sphinx-theme >= 0.9.0  

The most up-to-date version of pyMFD is available on Github. If you would 
like to run the tests, please install this version:

.. code-block:: console

    $ git clone https://github.com/larsenkg/pyMFD.git
    $ cd pyMFD
    $ python -m venv env
    $ source env/bin/activate
    (env) $ python -m pip install numpy==1.21.2 scipy==1.7.1 matplotlib==3.5.0 \ 
            bottleneck==1.3.2 pytest==7.1.2
    (env) $ python -m pytest pyMFD/tests/

pyMFD can also be installed using pip, which will also install the required 
dependencies:

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