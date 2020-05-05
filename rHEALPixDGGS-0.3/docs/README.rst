scenzgrid-dggs is documented in these files:

- ``scenzgrid-py/rhealpix_dggs_paper/rhealpix_dggs.pdf``
  Introduction to the rHEALPix discrete global grid system underlying scenzgrid-dggs
- ``scenzgrid-py/0.3/docs/build/latex/scenzgrid-py.pdf`` 
  scenzgrid-dggs manual
- ``scenzgrid-py/0.3/docs/build/html/index.html`` 
  scenzgrid-dggs manual in HTML format

The first document is produced via LaTeX from the source file ``scenzgrid-py/rhealpix_dggs_paper/rhealpix_dggs.tex`` and the bibliography file ``scenzgrid-py/scenzgrid.bib``.

The latter two documents are generated automatically from the source code of the scenzgrid-dggs modules lying in the ``scenzgrid-py/0.3/src`` directory.
To automatically build the latter two documents yourself, install the Python module `Sphinx <http://sphinx-doc.org/>`_ (but do not run ``sphinx-quickstart``, because the make file ``Makefile`` and the configuration file ``source/conf.py`` already exist) and then run the command ``make latexpdf`` to make the PDF documentation or ``make html`` to make the HTML documentation.
For the PDF documentation, you might also need to install `LaTeX <http://www.latex-project.org/>`_, but i'm not sure.

The ``source`` and ``build`` directories contain all the Sphinx source and build files, respectively.  
The rest of the files in the current directory are unrelated to Sphinx.
