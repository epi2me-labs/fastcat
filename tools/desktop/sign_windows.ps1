param(
    [Parameter(Mandatory = $true)]
    [string]$BundleDirectory
)

$ErrorActionPreference = 'Stop'

if (-not $env:CODESIGNTHUMBPRINT) {
    throw 'CODESIGNTHUMBPRINT must identify the Windows signing certificate.'
}
if (-not $env:SIGNTOOL_PATH) {
    throw 'SIGNTOOL_PATH must point to signtool.exe.'
}

$binary = Join-Path $BundleDirectory 'fastcat.exe'
Write-Host "Signing $binary"
& $env:SIGNTOOL_PATH sign /sha1 $env:CODESIGNTHUMBPRINT /sm /fd SHA256 /tr http://timestamp.digicert.com /td SHA256 $binary
if ($LASTEXITCODE -ne 0) { throw 'Authenticode signing failed.' }
Write-Host 'Verifying Authenticode signature'
& $env:SIGNTOOL_PATH verify /pa /v $binary
if ($LASTEXITCODE -ne 0) { throw 'Authenticode verification failed.' }
