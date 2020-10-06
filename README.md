# IRT Toolkit

Infrasound Ray Tracer 3D Toolkit

## Install dependencies using conda

Create a new conda environment
```
conda create --name=irt
```

Activate `irt` environment
```
conda activate irt
```

Add dependencies
```
conda install -c conda-forge compilers netcfd4 netcdf-fortran eccodes
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
