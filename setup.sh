#!/usr/bin/env bash
# Install everything the app needs, all inside this directory:
#   .env-tools/  micromamba prefix holding smina (+ its Open Babel runtime)
#   .venv/       Python venv holding FastAPI, uvicorn and RDKit
#
# Safe to re-run; existing installs are reused.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS="$HERE/.env-tools"
BIN="$HERE/.bin"
MAMBA="$BIN/micromamba"
VENV="$HERE/.venv"

say() { printf '\n\033[1;34m==>\033[0m %s\n' "$1"; }

# --- platform ---------------------------------------------------------------
case "$(uname -s)-$(uname -m)" in
  Darwin-arm64)  PLATFORM=osx-arm64 ;;
  Darwin-x86_64) PLATFORM=osx-64 ;;
  Linux-x86_64)  PLATFORM=linux-64 ;;
  Linux-aarch64) PLATFORM=linux-aarch64 ;;
  *) echo "unsupported platform: $(uname -s)-$(uname -m)" >&2; exit 1 ;;
esac

# --- smina via micromamba ---------------------------------------------------
if [ -x "$TOOLS/bin/smina" ]; then
  say "smina already installed at $TOOLS/bin/smina"
else
  if [ ! -x "$MAMBA" ]; then
    say "downloading micromamba ($PLATFORM)"
    mkdir -p "$BIN"
    curl -Ls "https://micro.mamba.run/api/micromamba/$PLATFORM/latest" \
      | tar -xj -C "$BIN" --strip-components=1 bin/micromamba
    chmod +x "$MAMBA"
  fi
  say "installing smina from conda-forge into $TOOLS"
  "$MAMBA" create -y -q -p "$TOOLS" -c conda-forge smina
fi

"$TOOLS/bin/smina" --version | head -1

# --- python side ------------------------------------------------------------
PY="${PYTHON:-python3}"
if [ ! -x "$VENV/bin/python" ]; then
  say "creating venv with $PY"
  "$PY" -m venv "$VENV"
fi
say "installing Python requirements"
"$VENV/bin/pip" install -q --upgrade pip
"$VENV/bin/pip" install -q -r "$HERE/requirements.txt"

say "done — start the app with:  ./run.sh"
