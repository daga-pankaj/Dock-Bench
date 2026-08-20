#!/usr/bin/env bash
# Start the docking web app on http://127.0.0.1:8000
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORT="${PORT:-8000}"

if [ ! -x "$HERE/.venv/bin/python" ]; then
  echo "no .venv found — run ./setup.sh first" >&2
  exit 1
fi

cd "$HERE"
exec "$HERE/.venv/bin/python" -m uvicorn backend.main:app --host 127.0.0.1 --port "$PORT" "$@"
