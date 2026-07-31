#!/bin/bash
# Runs the test suite. No install needed: sets PYTHONPATH to the source tree.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

PYTHONPATH=. python3 -m pytest tests/ -v "$@"
