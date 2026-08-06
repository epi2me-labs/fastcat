#!/usr/bin/env bash
set -euo pipefail

: "${OSX_SIGN_IDENTITY:?OSX_SIGN_IDENTITY must name the Developer ID certificate}"

bundle_dir=${1:?usage: sign_macos.sh BUNDLE_DIRECTORY}
binary="$bundle_dir/fastcat"

echo "Signing $binary with $OSX_SIGN_IDENTITY"
codesign --force --options runtime --timestamp --sign "$OSX_SIGN_IDENTITY" "$binary"
echo "Verifying macOS code signature"
codesign --verify --deep --strict --verbose=2 "$binary"

# A standalone executable cannot be stapled. The CI job submits a ZIP wrapper
# to Apple so the notarization ticket can be associated with this signature.
