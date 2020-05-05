Introduction
============
SCENZ-Grid is a geographic analysis system currently under development by Landcare Research, New Zealand.
scenzgrid-dggs is a collection of Python modules that implements SCENZ-Grid's discrete global grid system (DGGS), which is an instance of the rHEALPix DGGS.

This documentation assumes you are familiar with the rHEALPix DGGS as described in [GRS2013]_ and familiar with basic Python 2.7.x usage as described in, say, `The Python Tutorial <http://docs.python.org/2/tutorial/>`_.

.. [GRS2013] Robert Gibb, Alexander Raichev, Michael Speth, `The rHEALPix discrete global grid system <http://code.scenzgrid.org/index.php/p/scenzgrid-py/source/tree/master/rhealpix_dggs_paper/rhealpix_dggs.pdf>`_, in preparation, 2013.

Dependencies
---------------
- `Python 2.7 <http://python.org/>`_ 
- `numpy and scipy <http://www.scipy.org/>`_
  Third-party Python modules for scientific computation.
- `Sage <http://www.sagemath.org>`_
  (Optional) Third-party Python package for symbolic mathematics. 
  Needed only for a few optional graphics methods.
  
Installation
--------------
The modules of scenzgrid-dggs are open source and available for download at Landcare's git repository `<http://code.scenzgrid.org/index.php/p/scenzgrid-py/>`_ and can be cloned via the command ``git clone git@code.scenzgrid.org:scenzgrid-py.git``.

Usage
------
To use the scenzgrid-dggs modules, start a Python session in the directory where you downloaded the modules and import the modules.
Here are some examples.
For a list of all methods available, see the application programming interface (API) in the following chapters.

Using ``projections.py``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
The ``projections`` module implements the HEALPix and rHEALPix map projections on oblate ellipsoids of revolution whose authalic sphere is the unit sphere.
Because of this restriction, it is not so useful by itself.
Still, here are a few examples.  

Import all the classes, methods, and constants from the module ::

    >>> from projections import *

Project some points of an oblate ellipsoid of revolution using the HEALPix and (1, 2)-rHEALPix projections.

    >>> e = 0.2  # Eccentricity of ellipsoid
    >>> p = (0, 60)
    >>> q = healpix_ellipsoid(*p, e=e, degrees=True); print q
    (0.27646508139326409, 1.0618632447907124)
    >>> print healpix_ellipsoid(*q, e=e, degrees=True, inverse=True), p
    (6.3611093629270335e-15, 59.999997971669899) (0, 60)
    >>> q = rhealpix_ellipsoid(0, 60, e=e, north=1, south=2, degrees=True); print q
    (-0.27646508139326409, 1.0618632447907124)
    >>> print rhealpix_ellipsoid(*q, e=e, north=1, south=2, degrees=True, inverse=True), p
    (6.3611093629270335e-15, 59.999997971669899) (0, 60)

Using ``grids.py``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
The ``grids`` module implements the rHEALPix DGGS and various operations thereupon.
It depends upon the ``projections`` module.

Import all the classes, methods, and constants from the module ::

    >>> from grids import *

Create the (0, 0)-rHEALPix DGGS based upon the WGS84 ellipsoid. 
Use degrees for angular measurements ::

    >>> E = Earth(ellps='WGS84', degrees=True, south=0, north=0)
    >>> print E
    ____Earth model____
    lengths measured in meters and angles measured in degrees
    ellipsoid: WGS84
        major radius: 6378137
        flattening factor: 0.00335281066475
        minor radius: 6356752.31425
        eccentricity: 0.0818191908426
        authalic sphere radius: 6371007.18092
    central meridian: 0
    north pole square position: 0
    south pole square position: 0
    max resolution: 1.0
    max cell level: 15

Pick a (longitude-latitude) point on the ellipsoid and find the level 1 cell that contains it ::

    >>> p = (0, 15)
    >>> c = E.cell_from_point(1, p, surface='ellipsoid'); print c
    Q0
    
Find the planar (edge) neighbors of this cell ::

    >>> for (direction, cell) in c.neighbors('plane').items():
    ...     print direction, cell 
    down Q3
    right Q1
    up N2
    left P2

Compute the ellipsoidal shape and ellipsoidal nuclei of these cells ::

    >>> cells = c.neighbors().values()
    >>> for cell in cells:
    ...     print cell, cell.ellipsoidal_shape(), cell.nucleus_and_vertices(surface='ellipsoid')[0]
    Q3 rectangle (14.999999999999998, 0.0)
    Q1 rectangle (45.0, 26.490118738229611)
    N2 dart (0.0, 58.528017480415983)
    P2 rectangle (-14.999999999999998, 26.490118738229611)