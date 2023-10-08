************
rHEALPixDGGS
************

Introduction
============
rHEALPixDGGS is a Python package that implements the rHEALPix Discrete Global Grid System (DGGS).

Release Notes
-------------
This package was originally authored in 2013 and has had only minor code updates since then.

0.5.4 - current
^^^^^^^^^^^^^^^
Code unchanged from 0.5.3 other than updating to Python 3.11

Refer to file CHANGES.rst for a more detailed history of changes.

Requirements 
-------------
* ``requirements.txt`` - all the module requirements for operation
    - `NumPy >=1.7 <http://www.numpy.org/>`_ Base N-dimensional array package
    - `SciPy >=0.12 <http://www.scipy.org/>`_ Fundamental library for scientific computing
    - `Matplotlib >=1.2.1 <http://matplotlib.org/>`_ Comprehensive 2D Plotting
    - `Pyproj >=1.9.3 <http://code.google.com/p/pyproj/>`_ Python interface to the PROJ.4 cartographic library
* ``requirements.dev.txt`` - packages needed for developing this package

Installation
--------------
This package is available on PyPI, the Python Package Index from where it can be installed as follows:

::

    pip install rhealpixdggs

rHEALPixDGGS is also available for download from the github repository `<https://github.com/manaakiwhenua/rhealpixdggs-py>`_ from where the latest version can be cloned.
  
Tests
------
The files in the ``tests`` directory test the rHEALPixDGGS modules. These files are plain ``unittest`` files (the Python testing framework contained within the standard distribution). Tests for examples in documents need the ``doctest`` module installed (see ``requirements.dev.txt``).

Two UNIX shell scripts are included in this repository to run all unit and doc tests:

* ``run_doctests.sh``
* ``run_unittests.sh``

Running the command ``python tests/test_<foo>.py`` performs a sequence of automated tests of ``<foo>.py``.

For example, ``tests/test_distortion.py`` automatically tests ``distortion.py``.

If you update a module, then update its test file to test the changes you made!

Test early, test often, test automatically!

There are a couple of files in the main package directory that can be used to run all tests, starting ``run_...``.

Documentation
--------------
Documentation can be found at:

- `The rHEALPix Discrete Global Grid System <https://datastore.landcareresearch.co.nz/dataset/rhealpix-discrete-global-grid-system>`_ - The rHEALPix Discrete Global Grid System
- ``docs/build/latex/rHEALPixDGGS.pdf`` - The rHEALPixDGGS manual
- ``docs/build/html/index.html`` - The rHEALPixDGGS manual in HTML format

The latter two documents are generated automatically from the source code of the ``rhealpixdggs`` package modules.
To automatically build these yourself, install the Python package `Sphinx <http://sphinx-doc.org/>`_ (but do not run ``sphinx-quickstart``, because the make file ``Makefile`` and the configuration file ``docs/source/conf.py`` already exist) and then from the ``docs`` directory run the command ``make latexpdf`` to make the PDF documentation or ``make html`` to make the HTML documentation.
For the PDF documentation, you might also need to install `LaTeX <http://www.latex-project.org/>`_.

The ``source`` and ``build`` directories contain all the Sphinx source and build files, respectively.  

License
-------
This code is licensed under the `GNU Lesser General Public License v3.0, <http://www.gnu.org/licenses/lgpl-3.0.html>`_. See the file ``LICENSE`` for a copy of the deed.

Contact
-------
| *Maintainer*:
| **Robert Gibb**
| `Manaaki Whenua – Landcare Research <https://www.landcareresearch.co.nz/>`_
| `Gibbr@landcareresearch.co.nz <mailto:Gibbr@landcareresearch.co.nz>`_
|
| *Release Manager*:
| **Dr Nicholas J. Car**
| `SURROUND Australia Pty Ltd <https://surround.com>`_
| `nicholas.car@surroundaustralia.com <mailto:nicholas.car@surroundaustralia.com>`_
|
| *Original author*:
| **Alexander Raichev**
| `<https://raichev.net/>`_
| `alex@raichev.net <mailto:alex@raichev.net>`_
