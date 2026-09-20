#!/usr/bin/env bash
set -u
cd "$(dirname "$0")/.."
status=0

python -m coverage erase
python -m coverage run -m unittest discover tests || status=1
PYTHON="python -m coverage run -a" ./scripts/run_doctests.sh || status=1
python -m coverage report || status=1
python -m coverage xml

exit $status
