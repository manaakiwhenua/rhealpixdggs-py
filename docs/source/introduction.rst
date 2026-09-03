Introduction
============
rHEALPixDGGS is a Python 3 package that implements the rHEALPix discrete global grid system (DGGS).

.. figure:: images/planar_grid.*
   :alt: The (0, 0)-rHEALPix planar grid: six resolution 0 cells labeled
         N, O, P, Q, R, S unfolded from a cube, each subdivided into a
         3-by-3 resolution 1 sub-grid, with cell P's children labeled
         with their SUID digits 0 through 8.
   :width: 85%
   :align: center

   The (0, 0)-rHEALPix planar grid: the six resolution 0 cells unfolded
   from a cube, each subdivided into its resolution 1 sub-grid. Cell P's
   children are labeled with their SUID digits, so e.g. the cell in its
   lower-right corner is P8.

.. figure:: images/ellipsoidal_cells.*
   :alt: The resolution 1 ellipsoidal cells of the (0, 0)-rHEALPix DGGS on
         WGS84, drawn in longitude-latitude coordinates and colored by
         their resolution 0 cell.
   :width: 100%

   The same resolution 1 cells mapped back onto the WGS84 ellipsoid, in
   longitude-latitude coordinates: quad cells in the equatorial band;
   dart, skew-quad, and cap cells in the polar regions.

.. figure:: images/hierarchy.*
   :alt: Four nested squares illustrating the cell hierarchy: cell N
         subdivided into nine children, its center child N4 subdivided
         in turn, then N44, then N444.
   :width: 75%
   :align: center

   The cell hierarchy: each cell subdivides into an N_side-by-N_side
   grid of children (here N_side = 3), and a cell's address (its SUID)
   is its parent's address plus one digit -- N contains N4 contains N44
   contains N444. An address prefix is always an ancestor, which is
   what makes the predicates ``contains_cell``, ``within`` and the
   ``compact_cells`` function pure string operations.

This documentation assumes you are familiar with the rHEALPix DGGS as described in [Gibb2016]_ (with full mathematical detail in the preprint [GRS2013]_) and with basic Python usage as described in `The Python Tutorial <https://docs.python.org/3/tutorial/>`_.

Requirements
---------------
- `Python >=3.11 <https://www.python.org/>`_
- `NumPy >=2.0 <https://www.numpy.org/>`_ Base N-dimensional array package
- `SciPy >=1.11 <https://www.scipy.org/>`_ Fundamental library for scientific computing
- `Pyproj >=3.6 <https://pyproj4.github.io/pyproj/>`_
  Python interface to the PROJ cartographic library
- `Shapely >=2.1 <https://shapely.readthedocs.io/>`_ Manipulation and analysis of planar GEOS geometries

Installation
--------------
The package is available on PyPI, the Package Index from where it can be installed as follows:

pip install rhealpixdggs

rHEALPixDGGS is also available for download from the Bioeconomy Science Institute's GitHub repository `<https://github.com/manaakiwhenua/rhealpixdggs-py>`_ from where the latest version can be cloned.

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

``projection_wrapper`` uses the HEALPix and rHEALPix projections
defined in ``pj_healpix`` and ``pj_rhealpix``, this package's own pure-Python
reference implementations of the projection formulas in [Gibb2016]_, rather
than the C implementations in the PROJ library (reachable through Pyproj).
This keeps all grid computations identical regardless of the installed PROJ
version. To use PROJ's implementations instead, edit the
``HOMEMADE_PROJECTIONS`` line in ``projection_wrapper``.

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

Ellipsoid instances are parametrized by geographic longitude and latitude with the central meridian at ``lon_0``. (``lat_0`` is kept for PROJ name compatibility and must be 0: recentring in latitude is not a rotation of the ellipsoid and is not supported.)

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
    ...         print(cell, cell.ellipsoidal_shape, assert_allclose(cell.nucleus(plane=False), expected_results[i][j], rtol=1e-15, atol=0) == None)
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
Rotate the DGGS about the polar axis so that New Zealand sits in the middle of
an equatorial face: face edges and polar dart cells lie on the meridians
``lon_0 + k * 90``, so a ``lon_0`` 45 degrees west of Auckland keeps them clear
of the country. (Recentring in latitude is not supported; ``Ellipsoid`` rejects
a nonzero ``lat_0`` because a latitude shift is not a rotation of the ellipsoid
and would break the equal-area property.)::

    >>> p = (174, -37)  # Approximate Auckland lon-lat coordinates
    >>> from rhealpixdggs.projection_wrapper import *
    >>> E = Ellipsoid(a=WGS84_A, f=WGS84_F, radians=False, lon_0=p[0] - 45)
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
            lat_0 = 0
            lon_0 = 129
            radians = False
            sphere = False
    >>> print(rdggs.cell_from_point(1, p, plane=False))
    Q7

References
----------
.. [Gibb2016] Robert Gibb, `The rHEALPix discrete global grid system <https://doi.org/10.1088/1755-1315/34/1/012012>`__, IOP Conference Series: Earth and Environmental Science 34, 012012, 2016.
.. [GRS2013] Robert Gibb, Alexander Raichev, Michael Speth, `The rHEALPix discrete global grid system <https://github.com/manaakiwhenua/rhealpixdggs-2013-gibb-raichev-speth/blob/master/rhealpix_dggs_preprint.pdf>`__, preprint, 2013. Gives the full mathematical detail of the system this package implements.
.. [BS2018] David Bowater, Emmanuel Stefanakis, `The rHEALPix Discrete Global Grid System: considerations for Canada <https://doi.org/10.1139/geomat-2018-0008>`__, Geomatica 72(1), 27-37, 2018. Discusses cell shape and orientation in the polar region and how to choose ``lon_0`` and the polar squares for a region of interest.

Related reading:

- Mark R. Calabretta and Boudewijn F. Roukema, `Mapping on the HEALPix grid <https://doi.org/10.1111/j.1365-2966.2007.12297.x>`__, Monthly Notices of the Royal Astronomical Society 381(2), 865--872, 2007. Defines the HEALPix projection that rHEALPix rearranges.
- Charles F. F. Karney, `On auxiliary latitudes <https://doi.org/10.48550/arXiv.2212.05818>`__, 2023. Source of the authalic latitude power series used in ``utils.auth_lat``.
