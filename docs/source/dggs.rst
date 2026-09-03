The dggs Module
===========================

.. figure:: images/polar_squares.*
   :alt: Two unfoldings of the rHEALPix cube into the plane, drawn with
         coastlines. In the first, with north_square=0 and
         south_square=0, the polar squares N and S sit above and below
         the leftmost equatorial square O. In the second, with
         north_square=1 and south_square=2, N sits above P and S sits
         below Q.
   :width: 90%
   :align: center

   The ``north_square`` and ``south_square`` parameters choose which of
   the four equatorial squares the polar squares unfold above and
   below. Both grids cover the ellipsoid identically -- the choice only
   changes where the polar faces land in the planar projection.

.. figure:: images/recentred.*
   :alt: A world map in longitude-latitude coordinates showing the
         resolution 1 cell boundaries of a DGGS rotated so that New
         Zealand sits in the middle of face Q, with the cell containing
         Auckland filled and a star at Auckland itself.
   :width: 100%

   A DGGS rotated about the polar axis by constructing its ellipsoid
   with ``lon_0=129``. The meridian ``lon_0`` becomes the planar
   ``x = 0`` line, which is the shared edge of faces P and Q, and the
   other face edges and the polar dart meridians follow at
   ``lon_0 + k * 90``. Choosing ``lon_0`` 45° west of Auckland therefore
   puts New Zealand in the middle of face Q, as far from those edges as
   possible. Cells of the rotated grid can straddle the antimeridian, as
   the filled cell does here.

   Rotating about the polar axis (``lon_0``) and choosing where the
   polar squares unfold (``north_square``, ``south_square``) are the
   only ways to reorient an rHEALPix grid; they are exactly the
   transformations that keep every cell equal-area and every nucleus on
   a parallel of latitude. Recentring in latitude is not one of them:
   shifting latitude is not a rotation of an ellipsoid of revolution, so
   the resulting grid is neither equal-area nor coherent near the poles,
   and ``Ellipsoid`` rejects a nonzero ``lat_0`` (issue #93). A region
   that straddles the equatorial/polar boundary at about 41.9° latitude
   keeps that seam, but a good ``lon_0`` still helps: face edges and
   the polar dart cells lie on the meridians ``lon_0 + k * 90``, and the
   skew quad cells between them are north-aligned midway between those
   meridians. Choosing ``lon_0`` 45° away from a region's central
   meridian therefore puts the region mid-face with north-aligned
   cells on both sides of the seam [BS2018]_; the two figures below show
   this for New Zealand and for Canada.

.. figure:: images/recentred_nz.*
   :alt: Two maps of New Zealand side by side with resolution 4 cells
         coloured by ellipsoidal shape. On the left, the standard grid,
         where the cells south of the seam lean towards the dart
         meridian at 180 degrees just east of the country. On the
         right, the grid rotated to lon_0 of 129 degrees, where those
         cells are north-aligned and no dart meridian is in view.
   :width: 100%

   Resolution 4 cells over New Zealand coloured by ellipsoidal shape,
   on the standard grid (left) and on the grid rotated to ``lon_0=129``
   (right). The dotted line is the equatorial/polar seam at 41.9°S,
   which no rotation can move: cells north of it are quads, cells south
   of it are skew quads. On the standard grid the dart meridian (and
   face edge) at 180° runs just east of the country and the South
   Island's skew quads lean towards it. On the rotated grid the nearest
   dart meridians are at 129°E and 141°W, New Zealand sits midway
   between them, and its skew quads are north-aligned.

.. figure:: images/recentred_canada.*
   :alt: Two maps of Canada side by side with resolution 3 cells
         coloured by ellipsoidal shape. On the left, the standard grid,
         with a dart meridian at 90 degrees west running through the
         middle of the country. On the right, the grid rotated to lon_0
         of minus 50 degrees, with the dart meridians at 50 and 140
         degrees west at the country's edges.
   :width: 100%

   Resolution 3 cells over Canada, on the standard grid (left) and on
   the grid rotated to ``lon_0=-50`` (right), following [BS2018]_.
   Almost all of Canada lies north of the seam, so its cells are skew
   quads and darts whatever ``lon_0`` is. On the standard grid the dart
   meridian at 90°W runs through the middle of the country; the
   rotation moves the dart meridians to 50°W and 140°W, at its eastern
   and western edges, and the skew quads in between lean progressively
   from north-west to north-east across the country.

.. currentmodule:: rhealpixdggs.dggs

.. autosummary::

    RHEALPixDGGS

.. automodule:: rhealpixdggs.dggs
    :members:
    :undoc-members:
    :show-inheritance:
