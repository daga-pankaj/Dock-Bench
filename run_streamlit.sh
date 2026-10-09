#!/usr/bin/env bash
# Start the Streamlit app. Defaults to http://127.0.0.1:8501 (this machine
# only); set HOST=0.0.0.0 to let colleagues on the network open
# http://<this-machine>:8501 - and set DOCKING_AUTH_PASS first.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORT="${PORT:-8501}"
HOST="${HOST:-127.0.0.1}"

if [ ! -x "$HERE/.venv/bin/python" ]; then
  echo "no .venv found — run ./setup.sh first" >&2
  exit 1
fi

cd "$HERE"
exec "$HERE/.venv/bin/python" -m streamlit run streamlit_app/app.py \
  --server.port "$PORT" --server.address "$HOST" "$@"
