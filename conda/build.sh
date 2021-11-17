#!/bin/bash

NAME=fastcat

export HTS_CONF_ARGS="--prefix=${PREFIX} --enable-libcurl --with-libdeflate --enable-plugins --enable-gcs --enable-s3"
export EXTRA_CFLAGS="-I$PREFIX/include"
export EXTRA_LDFLAGS="-L$PREFIX/lib"
export EXTRA_LIBS="-ldl -ldeflate"

OS=$(uname)
if [[ "$OS" == "Darwin" ]]; then
    echo "Setting Darwin args"
    export ARGP=${PREFIX}/lib/libargp.a
    export EXTRA_CFLAGS="${EXTRA_CFLAGS} -isysroot ${CONDA_BUILD_SYSROOT} -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
fi

make clean

mkdir -p $PREFIX/bin
for binary in fastcat bamstats; do
    make $binary
    cp $binary $PREFIX/bin && chmod +x $PREFIX/bin/$binary
done
