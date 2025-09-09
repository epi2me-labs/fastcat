#!/bin/bash

NAME=fastcat


# don't enable libdeflate -- seems to cause hangs when used with threaded decompression
#export HTS_CONF_ARGS="--prefix=${PREFIX} --enable-libcurl --enable-plugins --enable-gcs --enable-s3"
# ignore that, just link to htslib from bioconda
export EXTRA_CFLAGS="-I$PREFIX/include"
export STATIC_HTSLIB=""
export EXTRA_LDFLAGS="-L$PREFIX/lib"
export EXTRA_LIBS="-ldl -lhts"

OS=$(uname)
if [[ "$OS" == "Darwin" ]]; then
    echo "Setting Darwin args"
    export ARGP=${PREFIX}/lib/libargp.a
    export EXTRA_CFLAGS="${EXTRA_CFLAGS} -isysroot ${CONDA_BUILD_SYSROOT} -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
fi

make clean clean_htslib

mkdir -p $PREFIX/bin
for binary in fastcat fastlint bamstats bamindex; do
    make $binary
    cp $binary $PREFIX/bin && chmod +x $PREFIX/bin/$binary
done
