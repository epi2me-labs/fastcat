#!/usr/bin/env bash

set -euxo pipefail

export PATH=/ucrt64/bin:/usr/bin:/bin:/c/Program\ Files/Git/bin:"$PATH"
export CC=/ucrt64/bin/gcc
export MAKEFLAGS="-j$(nproc) --output-sync=target --no-print-directory"

echo "=== Packaging Windows Desktop bundle ==="
make --trace CC="$CC" BUILD_MODE=desktop DESKTOP_TARGET=windows-x86_64 desktop_bundle

shopt -s nullglob
archives=(build/desktop/*-windows-x86_64.tar.gz)
test "${#archives[@]}" -eq 1
archive="${archives[0]}"
test -f "$archive"
echo "=== Desktop archive: $archive ==="
tar -tzf "$archive"
for file in fastcat.exe nextflow_schema.json README.md LICENSE SHA256SUMS; do
    test -f build/desktop-stage/fastcat-*-windows-x86_64/"$file"
done
