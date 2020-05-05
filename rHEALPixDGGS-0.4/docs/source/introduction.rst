Introduction
============
SCENZ-Grid is a geographic analysis system currently under development by Landcare Research, New Zealand.
scenzgrid-dggs is a collection of Python modules that implements SCENZ-Grid's discrete global grid system (DGGS), which is an instance of the rHEALPix DGGS.

This documentation assumes you are familiar with the rHEALPix DGGS as described in [GRS2013]_ and familiar with basic Python 2.7 usage as described in, say, `The Python Tutorial <http://docs.python.org/2/tutorial/>`_.

.. [GRS2013] Robert Gibb, Alexander Raichev, Michael Speth, `The rHEALPix discrete global grid system <http://code.scenzgrid.org/index.php/p/scenzgrid-py/source/tree/master/rhealpix_dggs_paper/rhealpix_dggs.pdf>`_, in preparation, 2013.

Dependencies
---------------
- `Python 2.7 <http://python.org/>`_ 
- `numpy and scipy <http://www.scipy.org/>`_
  Third-party Python modules for scientific computation.
- `pyproj <http://code.google.com/p/pyproj/>`_
  Third-party Python module that interfaces with the `PROJ.4 cartographic C library <http://trac.osgeo.org/proj/>`_.
  Needed only for when using PROJ.4 projections, e.g. for distortion calculations with ``distortion.py``.
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

Using ``projection_wrapper.py``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
The ``projection_wrapper`` module implements a wrapper for map projections of ellipsoids of revolution (but not of general triaxial ellipsoids) defined in pj_healpix, pj_rhealpix, and PROJ.4.
For brevity hereafter, the word 'ellipsoid' abbreviates 'ellipsoid of revolution'.

Import all the classes, methods, and constants from the module::

    >>> from projection_wrapper import *
    
Create an ellipsoid to work on (which can be a sphere if you wish), say, an ellipsoid with major radius 5 and eccentricity 0.8::

    >>> from ellipsoids import *
    >>> ellps_1 = Ellipsoid(a=5, e=0.8)
    >>> print ellps_1
    ellipsoid:
        R_A = 4.32200117119
        a = 5
        b = 3.0
        e = 0.8
        f = 0.4
        lat_0 = 0
        lon_0 = 0
        radians = False
        sphere = False

The names of the ellipsoid attributes agree with the names of the `PROJ.4 ellipsoid parameters <http://trac.osgeo.org/proj/wiki/GenParms>`_.

By default, angles are measured in degrees.
If you prefer radians, then do::

    >>> ellps_2 = Ellipsoid(a=5, e=0.8, radians=True)
    >>> print ellps_2
    ellipsoid:
        R_A = 4.32200117119
        a = 5
        b = 3.0
        e = 0.8
        f = 0.4
        lat_0 = 0
        lon_0 = 0
        radians = True
        sphere = False

Ellipsoid instances are parametrized by geographic longitude and latitude with the central meridian at ``lon_0`` and the parallel of origin at ``lat_0``. 

Project some points of the ellipsoid using the HEALPix and rHEALPix projections::

    >>> h = Proj(ellps_1, 'healpix')
    >>> rh = Proj(ellps_1, 'rhealpix', north_square=1, south_square=2)
    >>> print h(0, 60)
    (0.0, 3.351278550178026)
    >>> print rh(0, 60)
    (0.0, 3.351278550178026)
    >>> print h(0, 70)
    (0.864006732389895, 4.2584985144432679)
    >>> print rh(0, 70)
    (-0.86400673238989445, 4.2584985144432679)
    
Some common ellipsoids are predefined as constants.

    >>> print UNIT_SPHERE
    ellipsoid:
        R = 1
        R_A = 1
        a = 1
        b = 1
        e = 0
        f = 0
        lat_0 = 0
        lon_0 = 0
        radians = False
        sphere = True
    >>> print WGS84_ELLIPSOID
    ellipsoid:
        R_A = 6374581.4671
        a = 6378137.0
        b = 6356752.31414
        e = 0.0578063088401
        f = 0.00335281068118
        lat_0 = 0
        lon_0 = 0
        radians = False
        sphere = False
    >>> print WGS84_ELLIPSOID_RADIANS
    ellipsoid:
        R_A = 6374581.4671
        a = 6378137.0
        b = 6356752.31414
        e = 0.0578063088401
        f = 0.00335281068118
        lat_0 = 0
        lon_0 = 0
        radians = True
        sphere = False

Using ``dggs_rhealpix.py``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
The ``dggs_rhealpix`` module implements the rHEALPix DGGS and various operations thereupon.
It depends upon the ``projection_wrapper`` module.

Import all the classes, methods, and constants from the module ::

    >>> from dggs_rhealpix import *

Create the (0, 0)-rHEALPix DGGS with N_side=3 that is based upon the WGS84 ellipsoid::

    >>> from ellipsoids import WGS84_ELLIPSOID
    >>> E = WGS84_ELLIPSOID
    >>> rdggs = RHEALPixDGGS(ellipsoid=E, north_square=0, south_square=0, N_side=3)
    >>> print rdggs
    rHEALPix DGGS:
        N_side = 3
        north_square = 0
        south_square = 0
        max_areal_res = 1
        max_level = 15
        ellipsoid:
            R_A = 6374581.4671
            a = 6378137.0
            b = 6356752.31414
            e = 0.0578063088401
            f = 0.00335281068118
            lat_0 = 0
            lon_0 = 0
            radians = False
            sphere = False
    
Pick a (longitude-latitude) point on the ellipsoid and find the level 1 cell that contains it ::

    >>> p = (0, 15)
    >>> c = rdggs.cell_from_point(1, p, plane=False); print c
    Q0
    
Find the ellipsoidal (edge) neighbors of this cell ::

    >>> for (direction, cell) in c.neighbors(plane=False).items():
    ...     print direction, cell 
    west P2
    east Q1
    north N2
    south Q3

Find the planar (edge) neighbors of this cell ::

    >>> for (direction, cell) in c.neighbors('plane').items():
    ...     print direction, cell 
    down Q3
    right Q1
    up N2
    left P2

Find all the level 1 cells intersecting the longitude-latitude aligned ellipsoidal quadrangle with given northwest and southeast corners ::

    >>> nw = (0, 45)
    >>> se = (90, 0)
    >>> cells = rdggs.cells_from_region(1, nw, se, plane=False)
    >>> for row in cells:
    ...     print [str(cell) for cell in row]
    ['N2', 'N1', 'N0']
    ['Q0', 'Q1', 'Q2', 'R0']
    ['Q3', 'Q4', 'Q5', 'R3']

Compute the ellipsoidal shape and ellipsoidal nuclei of these cells ::

    >>> for row in cells:
    ...     for cell in row:
    ...         print cell, cell.ellipsoidal_shape(), cell.nucleus_and_vertices(plane=False)[0]
    N2 dart (5.0888874903416268e-14, 58.470677829627355)
    N1 skew_quad (45.000000000000036, 58.470677829627355)
    N0 dart (89.999999999999957, 58.470677829627363)
    Q0 quad (14.999999999999998, 26.438744923100096)
    Q1 quad (45.0, 26.438744923100096)
    Q2 quad (74.999999999999986, 26.438744923100096)
    R0 quad (105.0, 26.438744923100096)
    Q3 quad (14.999999999999998, 3.560649871414923e-15)
    Q4 quad (45.0, 3.560649871414923e-15)
    Q5 quad (74.999999999999986, 3.560649871414923e-15)
    R3 quad (105.0, 3.560649871414923e-15)
    
Create the (0, 0)-rHEALPix DGGS with N_side = 3 that is based on the WGS84 ellipsoid.
Orient the DGGS so that the planar origin (0, 0) is on Auckland, New Zealand::

    >>> p = (174, -37)  # Approximate Auckland lon-lat coordinates
    >>> from projection_wrapper import *
    >>> E = Ellipsoid(a=WGS84_A, f=WGS84_F, radians=False, lon_0=p[0], lat_0=p[1])
    >>> rdggs = RHEALPixDGGS(E, N_side=3, north_square=0, south_square=0)
    >>> print rdggs
    rHEALPix DGGS:
        N_side = 3
        north_square = 0
        south_square = 0
        max_areal_res = 1
        max_level = 15
        ellipsoid:
            R_A = 6374581.4671
            a = 6378137.0
            b = 6356752.31414
            e = 0.0578063088401
            f = 0.00335281068118
            lat_0 = -37
            lon_0 = 174
            radians = False
            sphere = False
    >>> print rdggs.cell_from_point(1, p, plane=False)
    Q3

Using ``distortion.py``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
The ``distortion`` module computes distortions of map projections created via ``projection_wrapper`` and was used to produce the tables in [GRS2013]_.
The module is not necessary for manipulating the rHEALPix DGGS.

Import all the classes, methods, and constants from the module::

    >>> from distortion import *

Import the WGS84 ellipsoid, define the rHEALPix projection on the ellipsoid, and compute linear distortion and areal distortion for the projection at a given point::

    >>> from ellipsoids import WGS84_ELLIPSOID
    >>> E = WGS84_ELLIPSOID
    >>> from projection_wrapper import Proj
    >>> f = Proj(ellipsoid=E, proj='rhealpix', north_square=0, south_square=0)
    >>> p = (0, 30)
    >>> print distortion(f, *p)  # First entry of output is f(*p)
    ((0.0, 3748655.115049511), 6.976075406932126, 1.1295629192757011, 1.1780969079183845)   
     
Sample 100 random points from the WGS84 ellipsoid and compute the sample minimum, sample maximum, sample median, sample mean, and sample standard deviation of the linear and area distortion functions of the rHEALPix projection of those points::

    >>> sample = [E.random_point() for i in range(100)]
    >>> print distortion_stats(f, sample)[1]  # doctest: +SKIP
    [[14.918, 13.82, 0.124, 48.086, 8.623], [1.345, 0.371, 1.002, 2.375, 1.163], [1.178, 0.0, 1.178, 1.178, 1.178]]
    
Do the same for 100 points chosen at random from an rHEALPix DGGS dart cell::

    >>> from dggs_rhealpix import *
    >>> rdggs = WGS84_003
    >>> cell = rdggs.cell(['N', 6])
    >>> sample = [cell.random_point(plane=False) for i in range(100)]
    >>> print distortion_stats(f, sample)[1]
    [[38.907, 6.044, 28.189, 48.963, 38.701], [2.013, 0.226, 1.644, 2.415, 1.991], [1.178, 0.0, 1.178, 1.178, 1.178]]