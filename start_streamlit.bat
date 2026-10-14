@echo off
rem Set up (first run only) and start the Streamlit docking app on Windows.
rem Double-click this file, or run it from cmd. Extra arguments are passed
rem through to streamlit, e.g.  start_streamlit.bat --server.headless true
setlocal

cd /d "%~dp0"

rem --- pick a Python 3 interpreter --------------------------------------------
rem Plain "python" on this machine can resolve to MGLTools' Python 2.7, so use
rem the py launcher unless PYTHON is already set.
if not defined PYTHON (
    for /f "delims=" %%i in ('py -3.12 -c "import sys; print(sys.executable)" 2^>nul') do set "PYTHON=%%i"
)
if not defined PYTHON (
    for /f "delims=" %%i in ('py -3 -c "import sys; print(sys.executable)" 2^>nul') do set "PYTHON=%%i"
)
if not defined PYTHON (
    echo No Python 3 found. Install Python 3.12 or set PYTHON to python.exe.
    goto :fail
)
echo Using Python: %PYTHON%

rem --- remove a venv built on macOS/Linux (bin\ instead of Scripts\) ----------
if exist ".venv\pyvenv.cfg" if not exist ".venv\Scripts\python.exe" (
    echo.
    echo .venv was created on another OS and cannot run on Windows.
    choice /m "Delete .venv and rebuild it"
    if errorlevel 2 goto :fail
    rmdir /s /q ".venv"
)

rem --- remove a macOS/Linux symlink copied over as a plain file ---------------
rem micromamba aborts with "Target prefix already exists and is not a folder".
if exist ".env-tools" if not exist ".env-tools\" (
    echo Removing stale .env-tools file ^(symlink copied from another OS^).
    del /f /q ".env-tools"
)

rem --- run setup if the venv or smina is missing ------------------------------
set "NEED_SETUP="
if not exist ".venv\Scripts\python.exe" set "NEED_SETUP=1"
if not exist ".env-tools\Library\bin\smina.exe" set "NEED_SETUP=1"
if defined NEED_SETUP (
    echo.
    echo Running setup.ps1 ...
    powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0setup.ps1"
    if errorlevel 1 (
        echo setup.ps1 failed - see the messages above.
        goto :fail
    )
)

rem --- start the app ----------------------------------------------------------
echo.
echo Starting Streamlit at http://127.0.0.1:8501 ...
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_streamlit.ps1" %*
if errorlevel 1 goto :fail
endlocal
exit /b 0

:fail
echo.
pause
endlocal
exit /b 1
