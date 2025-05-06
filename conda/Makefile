.SHELL := bash
.ONESHELL: # multiple commands are run in the same shell
.SHELLFLAGS := -euc
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules
.RECIPEPREFIX = > # commands start with >, not TAB

.PHONY: conda_package
conda_package: 
>
> echo "Install non-python dependencies"
> micromamba install python==3.11
> micromamba install rattler-build
> pip install --upgrade pip
> pip install build
>
> echo "Products in build/"
> rm -rf build || true
> mkdir build || true
>
> echo "Build python dist"
> python -m build .. --outdir=build/python_dist
>
> echo "Build conda recipe with correct local checksum and directory"
> PYTHON_DIST_SHA26=$$(sha256sum build/python_dist/rhealpixdggs-*.tar.gz | awk '{print $$1}')
> PYTHON_DIST_DIR=$$(realpath  build/python_dist/rhealpixdggs-0.5.6.tar.gz)
> sed "s|PYTHON_DIST_SHA26|$$PYTHON_DIST_SHA26|;s|PYTHON_DIST_DIR|$$PYTHON_DIST_DIR|" conda_recipe.yaml > build/conda_recipe_build.yaml
>
> echo "Build conda package"
> rattler-build build --recipe build/conda_recipe_build.yaml --output-dir=build/conda_package
>
> echo "Install locally"
> micromamba install ./build/conda_package::rhealpixdggs
>
> echo "test"
> echo "import rhealpixdggs ; print(rhealpixdggs.__file__)" | python
> python -m unittest discover ../tests/
> python -m doctest ../docs/source/introduction.rst ../rhealpixdggs/*.py
