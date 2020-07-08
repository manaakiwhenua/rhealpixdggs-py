Introduction
============
SCENZ-Grid is a geographic analysis system currently under development by Landcare Research, New Zealand.
scenzgrid-py is a collection of Python modules that implements the rHEALPix discrete global grid system (DGGS), an instance of which underlies SCENZ-Grid.
The latest stable version of scenzgrid-py is 0.4 which is located in the directory ``<current> = scenzgrid-py/0.4``.
All the source code and documentation lies there.

Dependencies 
-------------
- `Python 2.7.x <http://python.org/>`_ 
- `numpy and scipy <http://www.scipy.org/>`_
  Third-party Python modules for scientific computation.
- `pyproj <http://code.google.com/p/pyproj/>`_
  Third-party Python module that interfaces with the PROJ.4 cartographic library.
  Needed only for ``distortion.py``.
- `Sage <http://www.sagemath.org>`_
  (Optional) Third-party Python package for symbolic mathematics. 
  Needed only for a few optional graphics methods.
  Instead of installing Sage locally, you can use it `online <http://sagenb.org/>`_.
  To use the optional graphics methods, start a Sage notebook session and import/attach the Python module that contains the methods.
  For examples, see the Sage worksheet ``<current>/docs/figures.sws``.
  
Installation
--------------
The modules of scenzgrid-py are open source and available for download at Landcare's git repository `<http://code.scenzgrid.org/index.php/p/scenzgrid-py/>`_ and can be cloned via the command ``git clone git@code.scenzgrid.org:scenzgrid-py.git``.
  
Tests
------
The files in the ``<current>/src/tests`` directory test the scenzgrid-py modules. 
Running ``test_<foo>.py`` from the command line performs a sequence of automated tests of ``<foo>.py``.
For example, ``test_distortion.py`` automatically tests ``distortion.py``.
If you update a module, then update its test file to test the changes you made!
Test early, test often, test automatically!

The ``.sws`` files in the tests directory are `Sage <http://www.sagemath.org>`_ worksheets.
They are not automated tests, but rather supplementary graphical tests.
To run these, install Sage, start up a Sage notebook session, and open the worksheets.

Use ``<current>/src/run_all_tests.py`` to run all automated tests (but not the Sage tests) and doctest all the modules.
(Doctesting is testing the examples in the docstrings of the modules.)

Documentation
--------------
Documentation can be found at::

- ``<current>/docs/rhealpix_dggs_preprint.pdf``
  Introduction to the rHEALPix discrete global grid system underlying scenzgrid-py
- ``<current>/docs/build/latex/scenzgrid-py.pdf`` 
  scenzgrid-py manual
- ``<current>/docs/build/html/index.html`` 
  scenzgrid-py manual in HTML format

The latter two documents are generated automatically from the source code of the scenzgrid-py modules lying in the ``<current>/src`` directory.
To automatically build these yourself, install the Python module `Sphinx <http://sphinx-doc.org/>`_ (but do not run ``sphinx-quickstart``, because the make file ``Makefile`` and the configuration file ``source/conf.py`` already exist) and then run the command ``make latexpdf`` to make the PDF documentation or ``make html`` to make the HTML documentation.
For the PDF documentation, you might also need to install `LaTeX <http://www.latex-project.org/>`_.

The ``source`` and ``build`` directories contain all the Sphinx source and build files, respectively.  