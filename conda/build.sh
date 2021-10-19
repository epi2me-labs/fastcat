#!/bin/bash

NAME=fastcat

export EXTRA_CFLAGS="-I$PREFIX/include"
export EXTRA_LDFLAGS="-L$PREFIX/lib"
export EXTRA_LIBS="-ldl"

OS=$(uname)
if [[ "$OS" == "Darwin" ]]; then
    echo "Setting Darwin args"
    export ARGP=${PREFIX}/lib/libargp.a
    export EXTRA_CFLAGS="${EXTRA_CFLAGS} -isysroot ${CONDA_BUILD_SYSROOT} -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
fi

make clean $NAME

mkdir -p $PREFIX/bin
cp $NAME $PREFIX/bin && chmod +x $PREFIX/bin/$NAME
