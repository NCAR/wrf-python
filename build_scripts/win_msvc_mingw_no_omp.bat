cd ../fortran
CALL gfortran -E ompgen.F90 -cpp -o omp.f90
cd ..

CALL python setup.py clean --all

IF %PROCESSOR_ARCHITECTURE% == AMD64 (
    CALL python setup.py config_fc --f90flags="-O2 -mtune=generic" build --compiler=mingw32 --fcompiler=gnu95       
) ELSE (
    CALL python setup.py config_fc --f90flags="-O2 -mtune=generic -mincoming-stack-boundary=2" build --compiler=msvc --fcompiler=gnu95
)
    CALL pip install .
