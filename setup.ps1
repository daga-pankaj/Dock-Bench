# Install everything the app needs, all inside this directory:
#   .env-tools\  micromamba prefix holding smina (+ its Open Babel runtime)
#   .venv\       Python venv holding FastAPI, uvicorn and RDKit
#
# Safe to re-run; existing installs are reused.
#
# Native Windows is CPU-only: gnina (the GPU-capable engine) has no official
# Windows build. For GPU docking on a Windows box, run this app inside WSL2
# and use setup.sh there instead - see README.md.

$ErrorActionPreference = "Stop"

$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$Tools = Join-Path $Here ".env-tools"
$Bin = Join-Path $Here ".bin"
$Mamba = Join-Path $Bin "micromamba.exe"
$Venv = Join-Path $Here ".venv"

function Say($msg) {
    Write-Host ""
    Write-Host "==> $msg" -ForegroundColor Cyan
}

# --- smina via micromamba ---------------------------------------------------
$SminaExe = Join-Path $Tools "Library\bin\smina.exe"
if (Test-Path $SminaExe) {
    Say "smina already installed at $SminaExe"
} else {
    if (-not (Test-Path $Mamba)) {
        Say "downloading micromamba (win-64)"
        New-Item -ItemType Directory -Force -Path $Bin | Out-Null
        Invoke-WebRequest -UseBasicParsing `
            -Uri "https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-win-64.exe" `
            -OutFile $Mamba
    }
    Say "installing smina from conda-forge into $Tools"
    & $Mamba create -y -q -p $Tools -c conda-forge smina
    if ($LASTEXITCODE -ne 0) { throw "micromamba install failed" }
}
& $SminaExe --version | Select-Object -First 1

# --- gnina: no official Windows build upstream ------------------------------
Say "gnina has no Windows build - skipping. smina above covers CPU docking; use WSL2 + setup.sh for GPU."

# --- python side -------------------------------------------------------------
$PyExe = if ($env:PYTHON) { $env:PYTHON } else { "python" }
$VenvPython = Join-Path $Venv "Scripts\python.exe"
if (-not (Test-Path $VenvPython)) {
    Say "creating venv with $PyExe"
    & $PyExe -m venv $Venv
    if ($LASTEXITCODE -ne 0) { throw "venv creation failed" }
}
Say "installing Python requirements"
& $VenvPython -m pip install -q --upgrade pip
& $VenvPython -m pip install -q -r (Join-Path $Here "requirements.txt")

# --- frontend viewer bundle --------------------------------------------------
if (-not (Get-Command node -ErrorAction SilentlyContinue)) {
    throw "node/npm not found - required to build the Mol* viewer bundle (frontend/app.js)"
}
$Frontend = Join-Path $Here "frontend"
Say "installing frontend build dependencies (Mol*, React, esbuild)"
Push-Location $Frontend
try {
    npm install --no-audit --no-fund -q
    if ($LASTEXITCODE -ne 0) { throw "npm install failed" }
    Say "bundling the viewer (frontend/src/app.js -> frontend/app.js)"
    npm run build
    if ($LASTEXITCODE -ne 0) { throw "npm run build failed" }
} finally {
    Pop-Location
}

$Component = Join-Path $Here "streamlit_app\molstar_component\frontend"
Say "building the Streamlit app's Mol* viewer component"
Push-Location $Component
try {
    npm install --no-audit --no-fund -q
    if ($LASTEXITCODE -ne 0) { throw "npm install failed" }
    npm run build
    if ($LASTEXITCODE -ne 0) { throw "npm run build failed" }
} finally {
    Pop-Location
}

Say "done - start the app with:  .\run_streamlit.ps1  (or .\run.ps1 for the original web page)"
