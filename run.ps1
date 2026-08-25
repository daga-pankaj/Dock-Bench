# Start the docking web app on http://127.0.0.1:8000
$ErrorActionPreference = "Stop"

$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$Port = if ($env:PORT) { $env:PORT } else { "8000" }
$VenvPython = Join-Path $Here ".venv\Scripts\python.exe"

if (-not (Test-Path $VenvPython)) {
    Write-Error "no .venv found - run .\setup.ps1 first"
    exit 1
}

Set-Location $Here
& $VenvPython -u -m uvicorn backend.main:app --host 127.0.0.1 --port $Port @args
