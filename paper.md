---
title: 'pyMFD: A package for multipoint force-deflection analysis'
tags:
  - Python
  - physics
  - materials science
  - material characterization
authors:
  - name: Kyle G. Larsen
    orcid: 0000-0001-9968-2133
    affiliation: 1 
  - name: Robert C. Davis
    orcid: 0000-0002-6165-4396
    affiliation: 1
affiliations:
 - name: Department of Physics and Astronomy, Brigham Young University, Provo, UT, USA
   index: 1
date: 14 July 2022
bibliography: paper.bib
---

# Summary

Multipoint force-deflection (MFD) is a technique for measuring Young’s modulus in suspended thin films. It involves using an atomic force microscope in force-volume mode to measure a 2-D array of force-deflection ramps (which are 1-D data) on and around a microcantilever in the thin film. `pyMFD` can load this 3-D force-volume data (from a NanoScope file), analyze it, and calculate the Young’s modulus of the film. The analysis consists of summarizing each force-deflection ramp into a single compliance value (where compliance is the inverse of stiffness). The compliance along the cantilever is then fit to a fixed-free Euler beam [@howell_compliant_2001] and Young’s modulus is extracted.

# Statement of need

Multipoint force-deflection of microcantilevers is a new technique with no previously available software for facilitating analysis. Gwyddion [@gwyddion] is an open source program for visualizing and analyzing atomic force microscopy data. Support was recently added for loading force-volume data, but no tool yet exists for extracting compliance and performing MFD analysis. `pyMFD` was created to allow researchers in materials science, physics, and engineering to easily analyze MFD measurements. MFD deals with relatively large amounts of data compared to standard beam bending [@weihs_mechanical_1988] where only one force-deflection ramp is acquired. `pyMFD` can load NanoScope force-volume files, which can contain 4096 force-deflection ramps for a scan with 64x64 “pixels.” It can summarize the data into a 2-D image (a compliance map), from which it extracts the compliance along the cantilever, and calculates Young’s modulus.

# References