#!/usr/bin/env bash
## CTRL-c to quit
trap exit int
## various safety features
set -o nounset -o errexit -o pipefail

echo "Make / empty build directory (see build/log for output)"
rm -rf build || true
mkdir build || true

echo "Install build dependencies (in a conda prefix with micromamba).."
micromamba create -p ./build/env python==3.11 python-build rattler-build 
eval "$(micromamba shell hook --shell bash)"
micromamba activate ./build/env

echo "Build python dist.."
python -m build .. --outdir=build/python_dist

echo "Determine path to python dist and its checksum"
export PYTHON_DIST_ARCHIVE=$(realpath build/python_dist/rhealpixdggs-*.tar.gz)
export PYTHON_DIST_SHA26=$(sha256sum "$PYTHON_DIST_ARCHIVE" | awk '{print $1}')

echo "Build conda package.."
### rattler-build build --recipe meta.yaml --output-dir=build/conda_package
micromamba install conda-build
conda-build --override-channels --channel=conda-forge --output-folder=build/conda_package meta.yaml

echo "Install locally.."
micromamba install ./build/conda_package::rhealpixdggs 

echo "test.."
echo "import rhealpixdggs ; print(rhealpixdggs.__file__)" | python
python -m unittest discover ../tests/
python -m doctest ../docs/source/introduction.rst ../rhealpixdggs/*.py
