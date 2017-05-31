# Finite Volume Seismic Model Repository
Code used to generate results for the manuscript

[A High-Resolution Finite Volume Seismic Model to Generate Seafloor Deformation for Tsunami Modeling](https://link.springer.com/article/10.1007/s10915-017-0459-y) ([preprint](https://arxiv.org/abs/1701.01430))

## Dependencies

This code requires installation of the Clawpack software ([www.clawpack.org](https://www.clawpack.org)).  It was last tested with v5.4.0.

## Make Commands

The follow commands are available in the included Makefiles:

* figures: remove all previous output and plots, recompile code, make all figures
* .figures: check dependencies and remake only figures that need it

Note the 2d and 3d code exist in separate directories.
