# IRT Toolkit

Infrasound Ray Tracer 3D Toolkit, an in-house developed ray tracing algorithm (cast in spherical coordinates) that takes into account the full effect of the 3-D inhomogeneous wind and temperature fields.


## Installing pre-requisites using conda

Create and activate the irt environment:

```
conda create -y -n irt make compilers netcdf4 netcdf-fortran eccodes
conda activate irt
```

### macOS and 10.9 SDK

#### 10.9 SDK
Conda provides the clang compilers for macOS. But the macOS SDK is still required. The SDK license prevents it from being bundled in the conda package. The SDK has to be installed manually. For compatibility issue, conda packages are built with the 10.9 SDK.

The 10.9 SDK can be downloaded from:

https://github.com/devernay/xcodelegacy
https://github.com/phracker/MacOSX-SDKs

Download MacOSX10.9.sdk.tar.xz and untar it under `/opt/MacOSX10.9.sdk`.

#### Set environment variables

Before to be able to compile, three variables have to be set on macOS: `MACOSX_DEPLOYMENT_TARGET`, `CONDA_BUILD_SYSROOT`, and `SDKROOT`.


1. Locate the directory for the conda environment in your terminal window by running in the terminal echo `$CONDA_PREFIX`.
2. Enter that directory and create these subdirectories and files:

```
cd $CONDA_PREFIX
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
touch ./etc/conda/deactivate.d/env_vars.sh
```

3. Edit `./etc/conda/activate.d/env_vars.sh` as follows:

```
#!/bin/sh

export CONDA_BUILD_SYSROOT='/opt/MacOSX10.9.sdk'
export SDKROOT='/opt/MacOSX10.9.sdk'
export MACOSX_DEPLOYMENT_TARGET='10.9'
```

4. Edit `./etc/conda/deactivate.d/env_vars.sh` as follows:

```
#!/bin/sh

unset CONDA_BUILD_SYSROOT
unset SDKROOT
unset MACOSX_DEPLOYMENT_TARGET
```

5. Load environment variables:
```
conda activate irt
```


## Install source code

### Configure

The fortran compiler via conda is limited to GNU gfortran.
Make an alias to avoid compiler check issues in `configure` (only `gfortran` 
and not the full name `x86_64-apple-darwin13.4.0-gfortran` is currently allowed).

```
alias gfortran=`$CONDA_PREFIX/bin/nc-config --fc`

FCFLAGS="-O3 -m64 -funroll-all-loops -fpic -I${CONDA_PREFIX}/include"
FCLIBS="-L${CONDA_PREFIX}/lib -lnetcdf -lnetcdff -leccodes_f90 -leccodes"

FC=gfortran FCFLAGS=$FCFLAGS FCLIBS=$FCLIBS ./configure --prefix=$CONDA_PREFIX
```

### Compile and install

```
make
make install
```

Make check will give errors due to implementation errors. These can be ignored.


## Reference

> Smets, P. S. M., 2018. Infrasound and the Dynamical Stratosphere: A new application for operational weather and climate prediction, Dissertation, Delft University of Technology. Doi: [10.4233/uuid:517f8597-9c24-4d01-83ed-0f430353e905](https://doi.org/10.4233/uuid:517f8597-9c24-4d01-83ed-0f430353e905)

For additional background information see Section 2.2 and Appendix A.1.
