#!/usr/bin/env bash

set -euxo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 {build|test|check}" >&2
    exit 2
fi

export PATH=/ucrt64/bin:/usr/bin:/bin:/c/Program\ Files/Git/bin:"$PATH"
export CC=/ucrt64/bin/gcc
export MAKEFLAGS="-j$(nproc) --output-sync=target --no-print-directory"

case "$1" in
    build)
        echo "=== Building native Windows fastcat ==="
        make --trace CC="$CC" "$PWD/build/fastcat"
        ;;
    test)
        echo "=== Running native Windows unit tests ==="
        make --trace CC="$CC" test_unit
        ;;
    check)
        "$0" build
        ./build/fastcat fastq --help
        ./build/fastcat bam --help
        ./build/fastcat qc --help
        ./build/fastcat index --help
        "$0" test
        ;;
    *)
        echo "Unknown Windows build phase: $1" >&2
        exit 2
        ;;
esac
