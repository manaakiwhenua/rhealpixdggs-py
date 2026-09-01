#!/usr/bin/env bash
# Run every unit test. Works from any directory.
set -u
cd "$(dirname "$0")/.."

python -m unittest discover tests
