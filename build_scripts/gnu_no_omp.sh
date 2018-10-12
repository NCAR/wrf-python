#!/bin/bash

cd ../fortran
gfortran -E ompgen.F90 -cpp -o omp.f90
#f2py *.f90 -m _wrffortran -h wrffortran.pyf --overwrite-signature
cd ..

python setup.py clean --all
python setup.py config_fc --f90flags="-mtune=generic" build_ext build
pip install .
