#!/bin/bash

if [ `uname` == Darwin ]; then
    LDFLAGS="$LDFLAGS -undefined dynamic_lookup -bundle"
fi

$PYTHON setup.py config_fc --f90flags="-mtune=generic -mfpmath=sse" build 
$PYTHON setup.py install --single-version-externally-managed --record=record.txt


