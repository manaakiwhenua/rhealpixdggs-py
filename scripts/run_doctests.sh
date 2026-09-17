#!/usr/bin/env bash
# Run every doctest in the project. Exits nonzero if ANY file fails, not
# just the last one, so CI can rely on the exit code. Works from any
# directory.
set -u
cd "$(dirname "$0")/.."
status=0

# the .rst files with their own doctest code; harmonics needs the optional
# ducc0, as does the module of that name below
have_ducc0=1
python -c "import ducc0" 2>/dev/null || have_ducc0=0
for f in docs/source/introduction.rst docs/source/isolatitude.rst docs/source/harmonics.rst docs/source/export.rst docs/source/zoneset.rst
do
  if [ "$f" = docs/source/harmonics.rst ] && [ $have_ducc0 = 0 ]
  then
    echo "skipping $f: ducc0 is not installed"
    continue
  fi
  python -m doctest "$f" || status=1
done

# test the other python files
for f in rhealpixdggs/*.py
do
  if [ "$f" = rhealpixdggs/harmonics.py ] && [ $have_ducc0 = 0 ]
  then
    echo "skipping $f: ducc0 is not installed"
    continue
  fi
  python -m doctest "$f" || status=1
done

exit $status
