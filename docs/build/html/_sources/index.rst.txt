.. pyMFD documentation master file, created by
   sphinx-quickstart on Wed Jul 13 13:32:12 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyMFD Documentation
===================

pyMFD makes analyzing force-volume data for multipoint force-deflection (MFD) 
simple. MFD is a new method for characterizing suspended thin films. Small 
(e.g. 1 µm x 2 µm) cantilevers are cut (or etched) into a suspended thin 
film and then an atomic force microscope (AFM) is used to record many 
force-deflection ramps on and around the cantilever (using a force-volume
scan). pyMFD provides the capability to load the force-volume data, summarize 
each force-deflection ramp into a mechanical compliance, and fit the 
compliance data on the cantilever to estimate Young's modulus.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage
   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
