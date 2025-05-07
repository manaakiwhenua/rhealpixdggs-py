#!/usr/bin/env bash
## CTRL-c to quit
trap exit int
## various safety features
set -o nounset -o errexit -o pipefail

echo "Make / empty build directory (see build/log for output)"
rm -rf build || true
mkdir build || true

echo "Install build dependencies (in a conda prefix with micromamba).."
micromamba create --quiet -p ./build/env python==3.11 rattler-build 
eval "$(micromamba shell hook --shell bash)"
micromamba activate ./build/env

echo "Build conda package.."
rattler-build build --recipe meta_remote_source.yaml --output-dir=build/conda_package

echo "Install locally.."
micromamba install --quiet ./build/conda_package::rhealpixdggs 

echo "test.."
echo "import rhealpixdggs ; print(rhealpixdggs.__file__)" | python
python -m unittest discover ../tests/
python -m doctest ../docs/source/introduction.rst ../rhealpixdggs/*.py
