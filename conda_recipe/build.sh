#!/bin/sh

if [ `uname` == Darwin ]; then
    LDFLAGS="$LDFLAGS -undefined dynamic_lookup -bundle"
fi

if [`uname` == Darwin] || [`uname` == Linux]; then
    pip install .
else
    pip install --global-option build --global-option --compiler=mingw32 --global-option --fcompiler=gnu95
fi


