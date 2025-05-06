#!/usr/bin/env bash
## CTRL-c to quit
trap exit int
## various safety features
set -o nounset -o errexit -o pipefail


echo "Make / empty build directory (see build/log for output)"
rm -rf build || true
mkdir build || true

## I use a conda environment and micromamba to obtain the build dependcies
echo "Install build dependencies.."
micromamba create -p ./env
eval "$(micromamba shell hook --shell bash)"
micromamba activate ./env
micromamba install python==3.11 python-build rattler-build 

echo "Build python dist.."
python -m build .. --outdir=build/python_dist

echo "Build conda recipe with correct local checksum and directory."
PYTHON_DIST_ARCHIVE=$(realpath build/python_dist/rhealpixdggs-*.tar.gz)
PYTHON_DIST_SHA26=$(sha256sum "$PYTHON_DIST_ARCHIVE" | awk '{print $1}')
sed "s|PYTHON_DIST_SHA26|$PYTHON_DIST_SHA26|;s|PYTHON_DIST_ARCHIVE|$PYTHON_DIST_ARCHIVE|" conda_recipe.yaml > build/conda_recipe_build.yaml

echo "Build conda package.."
rattler-build build --recipe build/conda_recipe_build.yaml --output-dir=build/conda_package

echo "Install locally.."
micromamba install ./build/conda_package::rhealpixdggs 

echo "test.."
echo "import rhealpixdggs ; print(rhealpixdggs.__file__)" | python
python -m unittest discover ../tests/
python -m doctest ../docs/source/introduction.rst ../rhealpixdggs/*.py
