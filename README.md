# Global Magnetic Field Models

This set of fortran codes uses synoptic maps located on the solar surface and calculates a global magnetic field. It calculates both a PFSS model (Altschuler and Newkirk, 1969; Schatten, Wilcox and Ness, 1969) and a MHS model (Neukirch, 1995). From the MHS model a plasma pressure, density etc. may be calculated. The details of these codes and the mathematical models behind them can be found in Williams (2018). Williams (2018) describes a PFSS code and an MHS code. This repository is an updated version of both of these models (those with finite boundary conditions) codes for ease of use. The original MHS code listed in the thesis is still available where the infinite boundary condition version can be found.

The radial synoptic maps must first be extracted into a raw binary file. This is possible using `synoptic_map.py`. Once the synoptic map data is in the appropriate format, running `make` should create the binaries for `pfss` and `mhs_finite`.

Both have similar calling structures:
```sh
pfss -i {filename of synmap} -l {number of harmonics}
mhs_finite -i {filename of synmap} -l {number of harmonics} -a {value of alpha} -d {value of d}
```
These will output a magnetic field into a directory which by default is called `data`.

There are several other options:
```
-o: specify the output directory
-f: value to use for the gaussian filter
-r or --rmax: change the outer radius of the model (R_max in Williams, 2018)
-e: add string onto the end of the filename, before .dat
--np: change number of processors to use by OpenMP
```

OpenMP can be turned off using `make openmp=off`. By default it is on.
