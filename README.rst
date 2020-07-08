Introduction
============
rHEALPixDGGS is a Python package that implements the rHEALPix Discrete Global Grid System (DGGS).

Requirements 
-------------
- `Python >=3.3 <http://python.org/>`_ 
- `NumPy >=1.7 <http://www.numpy.org/>`_ Base N-dimensional array package 
- `SciPy >=0.12 <http://www.scipy.org/>`_ Fundamental library for scientific computing 
- `Matplotlib >=1.2.1 <http://matplotlib.org/>`_ Comprehensive 2D Plotting
- `Pyproj >=1.9.3 <http://code.google.com/p/pyproj/>`_
  Python interface to the PROJ.4 cartographic library
- `Sage >=5.10 <http://www.sagemath.org>`_
  (Optional) Comprehensive mathematics package. 
  Needed only for a few optional graphics methods.
  To use the optional graphics methods, start a Sage notebook session and import/attach the Python module that contains the methods.
  For examples, see the Sage worksheet ``tests/test_rhealpix_dggs.sws``.

Installation
--------------
This package is available on PyPI, the Python Package Index from where it can be installed as follows:

::

    pip install rhealpixdggs

rHEALPixDGGS is available for download from the github repository `<https://github.com/rggibb/rhealpixdggs-py>`_ from where the latest version can be cloned.
  
Tests
------
The files in the ``rhealpix_dggs/tests`` directory test the rHEALPixDGGS modules. 
Running the command ``python tests/test_<foo>.py`` from the ``rhealpix_dggs`` directory performs a sequence of automated tests of ``<foo>.py``.
For example, ``tests/test_distortion.py`` automatically tests ``distortion.py``.
If you update a module, then update its test file to test the changes you made!
Test early, test often, test automatically!

The ``.sws`` files in the ``tests`` directory are `Sage <http://www.sagemath.org>`_ worksheets.
They are not automated tests, but rather supplementary graphical tests.
To run these, install Sage, install Pyproj in Sage (download the Pyproj source code, change to the Pyproj directory, start a Sage shell via ``sage -sh``, then type ``python setup.py build``, then ``python setup.py install``; if that doesn't work, try again but first start in superuser mode via ``sudo su``), start up a Sage notebook session, and open the worksheets.

Run ``python run_all_tests.py`` from the ``rhealpix_dggs`` directory to run all automated tests (but not the Sage tests) and doctest all the modules of ``rhealpix_dggs``.

Documentation
--------------
Documentation can be found at::

- ``docs/rhealpix_dggs_preprint.pdf``
  Introduction to the rHEALPix discrete global grid system 
- ``docs/build/latex/rHEALPixDGGS.pdf`` 
  The rHEALPixDGGS manual
- ``docs/build/html/index.html`` 
  The rHEALPixDGGS manual in HTML format

The latter two documents are generated automatically from the source code of the ``rhealpix_dggs`` package modules.
To automatically build these yourself, install the Python package `Sphinx <http://sphinx-doc.org/>`_ (but do not run ``sphinx-quickstart``, because the make file ``Makefile`` and the configuration file ``docs/source/conf.py`` already exist) and then from the ``docs`` directory run the command ``make latexpdf`` to make the PDF documentation or ``make html`` to make the HTML documentation.
For the PDF documentation, you might also need to install `LaTeX <http://www.latex-project.org/>`_.

The ``source`` and ``build`` directories contain all the Sphinx source and build files, respectively.  
