#!/bin/bash

unset LDFLAGS

cd ../fortran
$FC -E ompgen.F90 -cpp -o omp.f90
cd ..

python setup.py clean --all
python setup.py config_fc --f90flags="-mtune=generic" build_ext build
pip install .
