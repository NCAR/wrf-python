#!/bin/bash

cd ../fortran
gfortran -E ompgen.F90 -cpp -o omp.f90
f2py *.f90 -m _wrffortran -h wrffortran.pyf --overwrite-signature
cd ..


if [$CONDA_BUILD == 1]
	if [ `uname` == Darwin ]; then
    	LDFLAGS="$LDFLAGS -undefined dynamic_lookup -bundle"
	fi
	$PYTHON setup.py clean --all
	$PYTHON setup.py config_fc --f90flags="-mtune=generic" build_ext build
	$PYTHON setup.py install --single-version-externally-managed --record=record.txt
else
	python setup.py clean --all
	python setup.py config_fc --f90flags="-mtune=generic" build_ext build
	pip install .
fi
