Introduction
============
rHEALPixDGGS is a Python 3 package that implements the rHEALPix discrete global grid system (DGGS).
This documentation assumes you are familiar with the rHEALPix DGGS as described in [GRS2013]_ and familiar with basic Python 3.3 usage as described in `The Python Tutorial <http://docs.python.org/3/tutorial/>`_.

.. [GRS2013] Robert Gibb, Alexander Raichev, Michael Speth, `The rHEALPix discrete global grid system <https://github.com/manaakiwhenua/rhealpixdggs-2013-gibb-raichev-speth/blob/master/rhealpix_dggs_preprint.pdf>`_, in preparation, 2013.

Requirements
---------------
- `Python >=3.11 <http://python.org/>`_
- `NumPy >=1.25.2 <http://www.numpy.org/>`_ Base N-dimensional array package
- `SciPy >=1.11.2 <http://www.scipy.org/>`_ Fundamental library for scientific computing
- `Matplotlib >=3.7.2 <http://matplotlib.org/>`_ Comprehensive 2D Plotting
- `Pyproj >=3.6 <http://code.google.com/p/pyproj/>`_
  Python interface to the PROJ.4 cartographic library

Installation
--------------
The package is available on PyPI, the Package Index from where it can be installed as follows:

pip install rhealpixdggs

rHEALPixDGGS is also available for download from Landcare Research's github repository `<https://github.com/manaakiwhenua/rhealpixdggs-py>`_ from wher the latest version can be cloned.

Usage
------
To use rHEALPixDGGS after installing it, start a Python session in the directory where you downloaded the modules and import the modules.
Here are some examples.
For a list of all methods available, see the application programming interface (API) in the following chapters.

Using the ``ellipsoids`` and ``projection_wrapper`` Modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The ``ellipsoids`` module implements functions and constants dealing with ellipsoids of revolution (which include spheres but not general triaxial ellipsoids).
For brevity hereafter, the word 'ellipsoid' abbreviates 'ellipsoid of revolution'.

The module ``projection_wrapper`` implements a wrapper for the map projections of ellipsoids defined in ``pj_healpix``, ``pj_rhealpix``, and Pyproj.

Currently ``projection_wrapper`` uses the HEALPix and rHEALPix projections
defined in ``pj_healpix`` and ``pj_rhealpix`` and *not* the buggy versions  defined in Pyproj 1.9.3 as ``PJ_healpix.c``.
Alternatively, you can download a corrected version of ``PJ_healpix.c`` from
`trac.osgeo.org/proj/changeset/2378 <http://trac.osgeo.org/proj/changeset/2378>`_, rebuild Pyproj with it, and use it in ``dggs`` by editing the ``HOMEMADE_PROJECTIONS`` line in ``projection_wrapper``.

Import all the classes, methods, and constants from the module::

    >>> from rhealpixdggs.projection_wrapper import *
    >>> from rhealpixdggs.ellipsoids import *

Create an ellipsoid, say, one with major radius 5 and eccentricity 0.8::

    >>> ellps_1 = Ellipsoid(a=5, e=0.8)
    >>> print(ellps_1)
    ellipsoid:
        R_A = 4.322001171188888
        a = 5
        b = 2.999999999999999
        e = 0.8
        f = 0.4
        lat_0 = 0
        lon_0 = 0
        radians = False
        sphere = False

The names of the ellipsoid attributes agree with the names of the `PROJ.4 ellipsoid parameters <http://trac.osgeo.org/proj/wiki/GenParms>`_.
For example, ``R_A`` is the authalic radius of the ellipsoid, the radius of the sphere that has the same area as the ellipsoid.

By default, angles are measured in degrees.
If you prefer radians, then do::

    >>> ellps_2 = Ellipsoid(a=5, e=0.8, radians=True)
    >>> print(ellps_2)
    ellipsoid:
        R_A = 4.322001171188888
        a = 5
        b = 2.999999999999999
        e = 0.8
        f = 0.4
        lat_0 = 0
        lon_0 = 0
        radians = True
        sphere = False

Some common ellipsoids are predefined as constants.

    >>> print(UNIT_SPHERE)
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
    >>> print(WGS84_ELLIPSOID)
    ellipsoid:
        R_A = 6371007.180918476
        a = 6378137.0
        b = 6356752.314245179
        e = 0.081819190842621
        f = 0.003352810664747
        lat_0 = 0
        lon_0 = 0
        radians = False
        sphere = False
    >>> print(WGS84_ELLIPSOID_RADIANS)
    ellipsoid:
        R_A = 6371007.180918476
        a = 6378137.0
        b = 6356752.314245179
        e = 0.081819190842621
        f = 0.003352810664747
        lat_0 = 0
        lon_0 = 0
        radians = True
        sphere = False

Ellipsoid instances are parametrized by geographic longitude and latitude with the central meridian at ``lon_0`` and the parallel of origin at ``lat_0``.

Project some points of the ellipsoid using the HEALPix and rHEALPix projections::

    >>> from numpy.testing import assert_allclose
    >>> h = Projection(ellps_1, 'healpix')
    >>> rh = Projection(ellps_1, 'rhealpix', north_square=1, south_square=2)
    >>> assert_allclose(h(0, 60), (0.0, 3.35127855017803), rtol=1e-14, atol=0) == None
    True
    >>> assert_allclose(rh(0, 60), (0.0, 3.35127855017803), rtol=1e-14, atol=0) == None
    True
    >>> assert_allclose(h(0, 70), (0.864006732389895, 4.258498514443268), rtol=1e-14, atol=0) == None
    True
    >>> assert_allclose(rh(0, 70), (-0.864006732389895, 4.258498514443268), rtol=1e-14, atol=0) == None
    True

Using the ``dggs`` Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The module ``dggs`` implements the rHEALPix DGGS and various operations thereupon.

Import all the classes, methods, and constants from the module ::

    >>> from rhealpixdggs.dggs import *

Create the (0, 0)-rHEALPix DGGS with N_side=3 that is based upon the WGS84 ellipsoid::

    >>> from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
    >>> E = WGS84_ELLIPSOID
    >>> rdggs = RHEALPixDGGS(ellipsoid=E, north_square=0, south_square=0, N_side=3)
    >>> print(rdggs)
    rHEALPix DGGS:
        N_side = 3
        north_square = 0
        south_square = 0
        max_areal_resolution = 1
        max_resolution = 15
        ellipsoid:
            R_A = 6371007.180918476
            a = 6378137.0
            b = 6356752.314245179
            e = 0.08181919084262149
            f = 0.0033528106647474805
            lat_0 = 0
            lon_0 = 0
            radians = False
            sphere = False

Some common rHEALPix DGGSs are predefined as constants::

    >>> print(UNIT_003)
    rHEALPix DGGS:
        N_side = 3
        north_square = 0
        south_square = 0
        max_areal_resolution = 1
        max_resolution = 1
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
    >>> print(WGS84_003)
    rHEALPix DGGS:
        N_side = 3
        north_square = 0
        south_square = 0
        max_areal_resolution = 1
        max_resolution = 15
        ellipsoid:
            R_A = 6371007.180918476
            a = 6378137.0
            b = 6356752.314245179
            e = 0.08181919084262149
            f = 0.0033528106647474805
            lat_0 = 0
            lon_0 = 0
            radians = False
            sphere = False
    >>> print(UNIT_003_RADIANS)
    rHEALPix DGGS:
        N_side = 3
        north_square = 0
        south_square = 0
        max_areal_resolution = 1
        max_resolution = 1
        ellipsoid:
            R = 1
            R_A = 1
            a = 1
            b = 1
            e = 0
            f = 0
            lat_0 = 0
            lon_0 = 0
            radians = True
            sphere = True

Pick a (longitude-latitude) point on the ellipsoid and find the level 1 cell that contains it ::

    >>> p = (0, 15)
    >>> c = rdggs.cell_from_point(1, p, plane=False); print(c)
    Q0

Find the ellipsoidal (edge) neighbors of this cell ::

    >>> for (direction, cell) in sorted(c.neighbors(plane=False).items()):
    ...     print(direction, cell)
    east Q1
    north N2
    south Q3
    west P2

Find the planar (edge) neighbors of this cell ::

    >>> for (direction, cell) in sorted(c.neighbors('plane').items()):
    ...     print(direction, cell)
    down Q3
    left P2
    right Q1
    up N2

Find all the level 1 cells intersecting the longitude-latitude aligned ellipsoidal quadrangle with given northwest and southeast corners ::

    >>> nw = (0, 45)
    >>> se = (90, 0)
    >>> cells = rdggs.cells_from_region(1, nw, se, plane=False)
    >>> for row in cells:
    ...     print([str(cell) for cell in row])
    ['N2', 'N1', 'N0']
    ['Q0', 'Q1', 'Q2', 'R0']
    ['Q3', 'Q4', 'Q5', 'R3']

Compute the ellipsoidal shape and ellipsoidal nuclei of these cells ::

    >>> expected_results = [
    ...    [
    ...        (0.0, 58.52801748206219),
    ...        (44.999999999999964, 58.52801748206219),
    ...        (90.0, 58.52801748206219)
    ...    ], [
    ...        (14.999999999999998, 26.490118751439734),
    ...        (45.0, 26.490118751439734),
    ...        (74.99999999999999, 26.490118751439734),
    ...        (105.00000000000001, 26.490118751439734)
    ...    ], [
    ...        (14.999999999999998, 0),
    ...        (45.0, 0),
    ...        (74.99999999999999, 0),
    ...        (105.00000000000001, 0)
    ...    ]]
    >>> for i, row in enumerate(cells):
    ...     for j, cell in enumerate(row):
    ...         print(cell, cell.ellipsoidal_shape(), assert_allclose(cell.nucleus(plane=False), expected_results[i][j], rtol=1e-15, atol=0) == None)
    N2 dart True
    N1 skew_quad True
    N0 dart True
    Q0 quad True
    Q1 quad True
    Q2 quad True
    R0 quad True
    Q3 quad True
    Q4 quad True
    Q5 quad True
    R3 quad True

Create the (0, 0)-rHEALPix DGGS with N_side = 3 that is based on the WGS84 ellipsoid.
Orient the DGGS so that the planar origin (0, 0) is on Auckland, New Zealand::

    >>> p = (174, -37)  # Approximate Auckland lon-lat coordinates
    >>> from rhealpixdggs.projection_wrapper import *
    >>> E = Ellipsoid(a=WGS84_A, f=WGS84_F, radians=False, lon_0=p[0], lat_0=p[1])
    >>> rdggs = RHEALPixDGGS(E, N_side=3, north_square=0, south_square=0)
    >>> print(rdggs)
    rHEALPix DGGS:
        N_side = 3
        north_square = 0
        south_square = 0
        max_areal_resolution = 1
        max_resolution = 15
        ellipsoid:
            R_A = 6371007.180918476
            a = 6378137.0
            b = 6356752.314245179
            e = 0.08181919084262149
            f = 0.0033528106647474805
            lat_0 = -37
            lon_0 = 174
            radians = False
            sphere = False
    >>> print(rdggs.cell_from_point(1, p, plane=False))
    Q3
