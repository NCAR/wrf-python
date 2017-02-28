
%ARCH% == 64 ( 
    %PYTHON% setup.py config_fc --f90flags="-O2 -mtune=generic -mfpmath=sse -msse2" build --compiler=mingw32 --fcompiler=gnu95 
) else (
    %PYTHON% setup.py config_fc --f90flags="-O2 -mtune=generic -mfpmath=sse -msse2 -mincoming-stack-boundary=2" build --compiler=mingw32 --fcompiler=gnu95
)
%PYTHON% setup.py install --single-version-externally-managed --record=record.txt

