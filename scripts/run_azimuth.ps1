# Run Azimuth PBMC annotation (install check + h5ad prep + R script)
# Usage: powershell -ExecutionPolicy Bypass -File scripts/run_azimuth.ps1

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
Set-Location $ProjectRoot

$Rscript = "C:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe"
if (-not (Test-Path $Rscript)) {
    $Rscript = (Get-Command Rscript -ErrorAction SilentlyContinue).Source
}
if (-not $Rscript) {
    Write-Error "Rscript not found. Install R from https://cran.r-project.org/"
}

# Rtools g++ for presto (if Azimuth not yet installed)
$env:PATH = "C:\rtools45\usr\bin;C:\rtools45\x86_64-w64-mingw32.static.posix\bin;" + $env:PATH

Write-Host "=== Step 1/3: Check / install R packages ===" -ForegroundColor Cyan
& $Rscript scripts/install_r_packages.R
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "=== Step 2/3: Prepare Azimuth-ready h5ad (Python) ===" -ForegroundColor Cyan
python scripts/prepare_azimuth_h5ad.py
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "=== Step 3/3: Run Azimuth annotation (R) ===" -ForegroundColor Cyan
& $Rscript scripts/azimuth_annotation.R
exit $LASTEXITCODE
