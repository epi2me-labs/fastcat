#!/bin/bash

NAME=fastcat

# use bundled htslib build for patched APIs
export EXTRA_CFLAGS="-I$PREFIX/include"
export STATIC_HTSLIB="third_party/htslib/libhts.a"
export EXTRA_LDFLAGS="-L$PREFIX/lib"
export EXTRA_LIBS="-ldl -ldeflate"
export HTS_CONF_ARGS="${HTS_CONF_ARGS} --with-libdeflate --enable-libcurl --enable-plugins --enable-gcs --enable-s3"
export HTS_BUILD_FLAGS="CPPFLAGS=\"${EXTRA_CFLAGS}\" LDFLAGS=\"${EXTRA_LDFLAGS}\""

OS=$(uname)
if [[ "$OS" == "Darwin" ]]; then
    echo "Setting Darwin args"
    export ARGP=${PREFIX}/lib/libargp.a
    export EXTRA_CFLAGS="${EXTRA_CFLAGS} -isysroot ${CONDA_BUILD_SYSROOT} -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
fi

make clean_all
make

mkdir -p $PREFIX/bin
for binary in fastcat; do
    cp build/$binary $PREFIX/bin/$binary
    chmod +x $PREFIX/bin/$binary
done
