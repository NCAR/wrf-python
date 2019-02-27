cd ..\fortran\build_help
CALL gfortran -o sizes -fopenmp omp_sizes.f90
CALL python sub_sizes.py

cd ..
CALL gfortran -E ompgen.F90 -cpp -fopenmp -o omp.f90
cd ..

CALL python setup.py clean --all
    
IF %PROCESSOR_ARCHITECTURE% == AMD64 (
    CALL python setup.py config_fc --f90flags="-O2 -mtune=generic -fopenmp" build_ext --libraries="gomp" build --compiler=msvc --fcompiler=gnu95       
) ELSE (
    CALL python setup.py config_fc --f90flags="-O2 -mtune=generic -fopenmp -mincoming-stack-boundary=2" build_ext --libraries="gomp" build --compiler=msvc --fcompiler=gnu95
)

CALL pip install .

