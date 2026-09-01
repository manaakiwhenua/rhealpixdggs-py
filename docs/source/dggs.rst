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
         resolution 1 cell boundaries of a DGGS recentred on the
         Auckland meridian, with the cell containing Auckland filled
         and a star at Auckland itself.
   :width: 100%

   A DGGS recentred on the meridian of Auckland by constructing its
   ellipsoid with ``lon_0=174``. Changing ``lon_0`` rotates the whole
   grid about the polar axis so that the chosen meridian becomes the
   planar ``x = 0`` line -- which by construction is the shared edge of
   faces P and Q, not the middle of a face, so Auckland lands in the
   cell that owns that edge. Cells of the rotated grid can straddle the
   antimeridian. (The ellipsoid also accepts ``lat_0``, but that shifts
   the grid in the *planar* y direction, which is not a rotation of the
   ellipsoid: the polar cells of such a grid are geographically
   incoherent, so recentring is best done in longitude only.)

.. currentmodule:: rhealpixdggs.dggs

.. autosummary::

    RHEALPixDGGS

.. automodule:: rhealpixdggs.dggs
    :members:
    :undoc-members:
    :show-inheritance:
