cd ../fortran
gfortran -E ompgen.F90 -cpp -o omp.f90
REM Wildcards not working on Windows for some reason
f2py -m _wrffortran -h wrffortran.pyf --overwrite-signature wrf_constants.f90 wrf_testfunc.f90 wrf_user.f90 rip_cape.f90 wrf_cloud_fracf.f90 wrf_fctt.f90 wrf_user_dbz.f90 wrf_relhl.f90 calc_uh.f90 wrf_user_latlon_routines.f90 wrf_pvo.f90 eqthecalc.f90 wrf_rip_phys_routines.f90 wrf_pw.f90 wrf_vinterp.f90 wrf_wind.f90 omp.f90
cd ..
python setup.py clean --all

%PROCESSOR_ARCHITECTURE% == AMD64 ( 
    python setup.py config_fc --f90flags="-O2 -mtune=generic -mfpmath=sse -msse2" build --compiler=mingw32 --fcompiler=gnu95 
) else (
    python setup.py config_fc --f90flags="-O2 -mtune=generic -mfpmath=sse -msse2 -mincoming-stack-boundary=2" build --compiler=mingw32 --fcompiler=gnu95
)
python setup.py install --single-version-externally-managed --record=record.txt
