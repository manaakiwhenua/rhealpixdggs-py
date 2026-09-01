# rHEALPixDGGS

[![Tests](https://github.com/manaakiwhenua/rhealpixdggs-py/actions/workflows/tests.yml/badge.svg)](https://github.com/manaakiwhenua/rhealpixdggs-py/actions/workflows/tests.yml)
[![manaakiwhenua-standards](https://github.com/manaakiwhenua/rhealpixdggs-py/workflows/manaakiwhenua-standards/badge.svg)](https://github.com/manaakiwhenua/manaakiwhenua-standards)

## Introduction

rHEALPixDGGS is a Python package that implements the rHEALPix Discrete Global Grid System (DGGS).

### Release Notes

This package was originally authored by Alexander Raichev in 2013 and is
actively maintained. See `CHANGES.rst` for the release history, including
explicit callouts of any breaking changes.

### Error Budget

rHEALPix is a theoretically equal-area projection: all cells at a given
resolution have the same area by construction. In practice, floating-point
arithmetic introduces a small rounding error when computing cell areas. The
`RHEALPixDGGS.area_error_budget()` method returns these values analytically
for the default WGS84-based DGGS (`WGS84_003`, `N_side=3`,
`max_areal_resolution=1 m²`):

| Resolution | Cell area (m²) | Cell area (km²) | Abs. tolerance (m²) |
|------------|----------------|-----------------|---------------------|
| 0          | 8.501 × 10¹³   | 85,010,937      | 0.189               |
| 1          | 9.446 × 10¹²   | 9,445,660       | 0.021               |
| 2          | 1.050 × 10¹²   | 1,049,518       | 0.0023              |
| 3          | 1.166 × 10¹¹   | 116,613         | 2.59 × 10⁻⁴         |
| 4          | 1.296 × 10¹⁰   | 12,957          | 2.88 × 10⁻⁵         |
| 5          | 1.440 × 10⁹    | 1,440           | 3.20 × 10⁻⁶         |
| 6          | 1.600 × 10⁸    | 160             | 3.55 × 10⁻⁷         |
| 7          | 1.777 × 10⁷    | 17.8            | 3.95 × 10⁻⁸         |
| 8          | 1.975 × 10⁶    | 1.97            | 4.39 × 10⁻⁹         |
| 9          | 2.194 × 10⁵    | 0.219           | 4.87 × 10⁻¹⁰        |
| 10         | 2.438 × 10⁴    | 0.0244          | 5.41 × 10⁻¹¹        |
| 11         | 2,709          | 0.00271         | 6.02 × 10⁻¹²        |
| 12         | 301            | 3.01 × 10⁻⁴     | 6.68 × 10⁻¹³        |
| 13         | 33.4           | 3.34 × 10⁻⁵     | 7.43 × 10⁻¹⁴        |
| 14         | 3.72           | 3.72 × 10⁻⁶     | 8.25 × 10⁻¹⁵        |
| 15         | 0.413          | 4.13 × 10⁻⁷     | 9.17 × 10⁻¹⁶        |

The **absolute tolerance** answers the question: "how different can two
computed cell areas be before we should conclude they are genuinely different,
rather than just floating-point noise?" For example, at resolution 9 the
tolerance is ~5 × 10⁻¹⁰ m² — if two cells' computed areas differ by less than
half a nanometre squared, they should be treated as equal-area.

The relative tolerance is constant at 2.22 × 10⁻¹⁵ (10 × machine epsilon)
for all resolutions, reflecting the fixed number of floating-point operations
in the area formula. The absolute tolerance shrinks proportionally with cell
area as resolution increases. Resolution 15 is the finest level supported by
the default configuration (cells of approximately 0.41 m²).

### Requirements

Dependencies are declared in `pyproject.toml` (this project is managed
with [Poetry](https://python-poetry.org/); development dependencies are
declared there too, under `tool.poetry.group.dev.dependencies`):

- [NumPy >=2.0](https://www.numpy.org/) Base N-dimensional array package
- [SciPy >=1.11](https://www.scipy.org/) Fundamental library for scientific computing
- [Pyproj >=3.6](https://pyproj4.github.io/pyproj/) Python interface to the PROJ cartographic library
- [Shapely >=2.1](https://shapely.readthedocs.io/) Manipulation and analysis of planar GEOS geometries

### Installation

This package is available on PyPI, the Python Package Index from where it can be installed as follows:

    pip install rhealpixdggs

rHEALPixDGGS is also available for download from the GitHub repository <https://github.com/manaakiwhenua/rhealpixdggs-py> from where the latest version can be cloned.

You can install from source using Poetry in a virtual environment (MacOS and Linux, assuming Poetry is already installed on your system):

    python3 -m venv rhealpixdggs
    source rhealpixdggs/bin/activate
    python -m pip install --upgrade pip
    poetry install

Or on Windows:

    python3 -m venv rhealpixdggs
    rhealpixdggs\Scripts\activate
    python -m pip install --upgrade pip
    poetry install

For development:

    poetry shell
    poetry install

### Tests

The files in the `tests` directory test the rHEALPixDGGS modules. These files are plain `unittest` files (the Python testing framework contained within the standard distribution), and the examples embedded in docstrings are tested with `doctest` (also part of the standard distribution).

Two UNIX shell scripts are included in this repository to run all unit and doc tests:

- `scripts/run_doctests.sh`
- `scripts/run_unittests.sh`

Running the command `python tests/test_<foo>.py` performs a sequence of automated tests of `<foo>.py`.

For example, `tests/test_cell.py` automatically tests `cell.py`.

If you update a module, then update its test file to test the changes you made!

Test early, test often, test automatically!

There are a couple of files in the main package directory that can be used to run all tests, starting `run_...`.

### Documentation

Documentation can be found at:

- [The rHEALPix Discrete Global Grid System](https://datastore.landcareresearch.co.nz/dataset/rhealpix-discrete-global-grid-system) - The rHEALPix Discrete Global Grid System
- `docs/build/latex/rHEALPixDGGS.pdf` - The rHEALPixDGGS manual
- `docs/build/html/index.html` - The rHEALPixDGGS manual in HTML format

The latter two documents are generated automatically from the source code of the `rhealpixdggs` package modules.
To automatically build these yourself, install the Python package [Sphinx](http://sphinx-doc.org/) (but do not run `sphinx-quickstart`, because the make file `Makefile` and the configuration file `docs/source/conf.py` already exist) and then from the `docs` directory run the command `make latexpdf` to make the PDF documentation or `make html` to make the HTML documentation.
For the PDF documentation, you might also need to install [LaTeX](http://www.latex-project.org/).

The `source` and `build` directories contain all the Sphinx source and build files, respectively.

### License

This code is dual-licensed. You may choose to use it under either:

- [GNU Lesser General Public License v3.0](http://www.gnu.org/licenses/lgpl-3.0.html) (see `LICENSE-LGPL`)
- [MIT License](https://opensource.org/licenses/MIT) (see `LICENSE-MIT`)

See the file `LICENSE` for more information about the dual licensing arrangement.

This dual licensing approach allows the software to be used in contexts where LGPL may not be permitted (such as certain corporate environments or standards bodies like the OGC) while still maintaining a copyleft option for those who prefer it.

### Releasing

See `RELEASING.md` for the release process. Note that working with this
repository's packaging (`poetry install`, `poetry build`, ...) requires
[Poetry](https://python-poetry.org/) version 2.0 or later.

A **conda package** `rhealpixdggs` is also maintained at
[conda-forge](https://github.com/conda-forge/rhealpixdggs-feedstock).

### Contact

*Maintainer*:  
**Richard Law**  
[Bioeconomy Science Institute](https://www.bioeconomyscience.co.nz/)  
<lawr@landcareresearch.co.nz>  
  
*Release Manager*:  
**Dr Nicholas J. Car**  
[SURROUND Australia Pty Ltd](https://surround.com)  
<nicholas.car@surroundaustralia.com>  
  
*Original author*:  
**Alexander Raichev**  
<https://raichev.net/>  
<alex@raichev.net>  
  
*Former Maintainer*:  
**Robert Gibb**  
[Manaaki Whenua – Landcare Research](https://www.landcareresearch.co.nz/)  
<gibbr@landcareresearch.co.nz>
