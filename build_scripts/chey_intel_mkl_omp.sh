#!/bin/bash

cd ../fortran/build_help
ifort -o sizes -qopenmp omp_sizes.f90
python sub_sizes.py

cd ..
ifort ompgen.F90 -qopenmp -fpp -save-temps 2>/dev/null
mv ompgen.i90 omp.f90
f2py *.f90 -m _wrffortran -h wrffortran.pyf --overwrite-signature
cd ..

python setup.py clean --all
python setup.py config_fc --f90flags="-O3 -xHost -mkl -qopenmp" build --compiler=intelem --fcompiler=intelem
pip install .

