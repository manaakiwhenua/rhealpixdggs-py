#!/usr/bin/env bash
# Run every doctest in the project. Exits nonzero if ANY file fails, not
# just the last one, so CI can rely on the exit code. Works from any
# directory.
set -u
cd "$(dirname "$0")/.."
status=0

# the .rst files with their own doctest code
for f in docs/source/introduction.rst docs/source/isolatitude.rst
do
  python -m doctest "$f" || status=1
done

# test the other python files
for f in rhealpixdggs/*.py
do
  python -m doctest "$f" || status=1
done

exit $status
