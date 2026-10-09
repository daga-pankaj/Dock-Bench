# Start the Streamlit app on http://127.0.0.1:8501 (set $env:HOST = "0.0.0.0"
# to share it on the network - and set DOCKING_AUTH_PASS first).
$ErrorActionPreference = "Stop"

$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$Port = if ($env:PORT) { $env:PORT } else { "8501" }
$BindHost = if ($env:HOST) { $env:HOST } else { "127.0.0.1" }
$VenvPython = Join-Path $Here ".venv\Scripts\python.exe"

if (-not (Test-Path $VenvPython)) {
    Write-Error "no .venv found - run .\setup.ps1 first"
    exit 1
}

Set-Location $Here
& $VenvPython -m streamlit run streamlit_app/app.py --server.port $Port --server.address $BindHost @args
