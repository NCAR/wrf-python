#!/bin/bash

unset LDFLAGS

cd ../fortran/build_help
$FC -o sizes -fopenmp omp_sizes.f90
python sub_sizes.py

cd ..
$FC -E ompgen.F90 -fopenmp -cpp -o omp.f90
cd ..

python setup.py clean --all
python setup.py config_fc --f90flags="-mtune=generic -fopenmp" build_ext --libraries="gomp" build
pip install .

