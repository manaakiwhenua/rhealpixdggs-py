Introduction
============
rHEALPixDGGS is a Python package that implements the rHEALPix Discrete Global Grid System (DGGS).

Release Note
------------
This package was originally authored in 2013 and has had only minor code updates since then. It is currently - July 2020 - at Release 0.5.1.

Many of the 0.5.1 release tests fail but mostly for trivial reasons. A 0.5.2 release is planned shortly (within July, 2020) that will see many of the trivial test issues fixed.

Requirements 
-------------
* ``requirements.txt`` - all the module requirements for operation
    - `NumPy >=1.7 <http://www.numpy.org/>`_ Base N-dimensional array package
    - `SciPy >=0.12 <http://www.scipy.org/>`_ Fundamental library for scientific computing
    - `Matplotlib >=1.2.1 <http://matplotlib.org/>`_ Comprehensive 2D Plotting
    - `Pyproj >=1.9.3 <http://code.google.com/p/pyproj/>`_ Python interface to the PROJ.4 cartographic library
* ``requirements.opt.txt`` - Optional comprehensive mathematics package needed only for a few optional graphics methods.
    - `Sage >=5.10 <http://www.sagemath.org>`_ - To use the optional graphics methods, start a Sage notebook session and import/attach the Python module that contains the methods. For examples, see the Sage worksheet ``tests/test_rhealpixdggs.sws``.
* ``requirements.dev.txt`` - packages for developing this package

Installation
--------------
This package is available on PyPI, the Python Package Index from where it can be installed as follows:

::

    pip install rhealpixdggs

rHEALPixDGGS is also available for download from the github repository `<https://github.com/manaakiwhenua/rhealpixdggs-py>`_ from where the latest version can be cloned.
  
Tests
------
The files in the ``tests`` directory test the rHEALPixDGGS modules. These files are plain unittest files (the Python testing framework contained within the standard distribution) but, movinf forward, `pytest <https://docs.pytest.org/>`_ will be used.

Running the command ``python tests/test_<foo>.py`` performs a sequence of automated tests of ``<foo>.py``.

For example, ``tests/test_distortion.py`` automatically tests ``distortion.py``.

If you update a module, then update its test file to test the changes you made!

Test early, test often, test automatically!

The ``.sws`` files in the ``tests`` directory are `Sage <http://www.sagemath.org>`_ worksheets.
They are not automated tests, but rather supplementary graphical tests.
To run these, install Sage, install Pyproj in Sage (download the Pyproj source code, change to the Pyproj directory, start a Sage shell via ``sage -sh``, then type ``python setup.py build``, then ``python setup.py install``; if that doesn't work, try again but first start in superuser mode via ``sudo su``), start up a Sage notebook session, and open the worksheets.

There are a couple of files in the main package directory that can be used to run all tests, starting ``run_...``.

Documentation
--------------
Documentation can be found at:

- `The rHEALPix Discrete Global Grid System <https://datastore.landcareresearch.co.nz/dataset/rhealpix-discrete-global-grid-system>`_ - The rHEALPix Discrete Global Grid System
- ``docs/build/latex/rHEALPixDGGS.pdf`` - The rHEALPixDGGS manual (the main paper)
- ``docs/build/html/index.html`` - The rHEALPixDGGS manual in HTML format

The latter two documents are generated automatically from the source code of the ``rhealpixdggs`` package modules.
To automatically build these yourself, install the Python package `Sphinx <http://sphinx-doc.org/>`_ (but do not run ``sphinx-quickstart``, because the make file ``Makefile`` and the configuration file ``docs/source/conf.py`` already exist) and then from the ``docs`` directory run the command ``make latexpdf`` to make the PDF documentation or ``make html`` to make the HTML documentation.
For the PDF documentation, you might also need to install `LaTeX <http://www.latex-project.org/>`_.

The ``source`` and ``build`` directories contain all the Sphinx source and build files, respectively.  

License
-------
This code is licensed under the `GNU General Public License, v3 <https://www.gnu.org/licenses/gpl-3.0.html>`_. See the file ``LICENSE`` for a copy of the deed.

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