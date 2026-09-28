The rhp_wrappers Module
===========================

.. figure:: images/wrappers_nz.*
   :alt: Two maps of New Zealand. Left: a polygon roughly enclosing the
         country, with the resolution 4 cells whose centroids fall inside
         it drawn filled. Right: a three-point line from the southwest of
         the South Island to the east of the North Island, with the
         sequence of resolution 4 cells it passes through drawn filled.
   :width: 100%

   The H3-style wrappers over New Zealand: ``polyfill`` (left) returns
   the cells whose centroids fall inside a polygon; ``linetrace`` (right)
   returns the cells a line passes through, in order. The color change
   partway down the South Island is the boundary between the equatorial
   R cells (quads) and the south polar S cells (skew quads).

.. figure:: images/wrappers_cap_trace.*
   :alt: View from above the north pole showing the nine resolution 3
         cells around it, with a linestring curving around the pole and
         the cells it passes through -- including the circular cap cell
         N444 -- drawn filled.
   :width: 80%
   :align: center

   ``linetrace`` across the north polar cap, viewed from above the pole
   (Arctic coastlines in grey: northern Greenland at lower left,
   Svalbard below, Franz Josef Land at lower right; the resolution 4
   sub-grid is drawn as a faint outline). A linestring's
   segments are straight in longitude-latitude coordinates, so at these
   latitudes each leg curves *around* the pole rather than hopping
   across it, passing through the circular cap cell N444 and the dart
   and skew-quad cells surrounding it -- the cell shapes for which no
   longitude-latitude rectangle is a valid approximation.

.. figure:: images/rings.*
   :alt: Left: a planar nine-by-nine grid with a center cell and its
         first three cell rings shaded in distinct colors, containing 8,
         16 and 24 cells. Right: an orthographic view centered on a cube
         corner, where the first ring around cell N00 has only 7 cells.
   :width: 100%

   ``cell_ring``: away from the cube's corners (left, around Q44) the
   ring at distance *k* has the familiar 8k cells. Around a cell at a
   cube corner (right, around N00, corner marked) the rings are smaller
   -- the k=1 ring has 7 cells -- because only three cells meet at the
   corner, so the ring spans three cube faces without a fourth,
   diagonally-opposite cell.

.. figure:: images/wrap_antimeridian.*
   :alt: Two world maps tracing the segment from longitude 179 to
         longitude -179 at latitude 10. Top, with wrap_antimeridian
         False: the traced cells span the whole map through longitude
         zero. Bottom, with wrap_antimeridian True: only two cells, one
         on each side of the antimeridian, are traced.
   :width: 100%

   The ``wrap_antimeridian`` flag on ``linetrace``. Coordinates are
   treated as plate carree, so the segment (179, 10) to (-179, 10)
   literally crosses the whole plane through longitude 0 (top, the
   default). With ``wrap_antimeridian=True`` (bottom) the segment is
   interpreted as taking the short way across the antimeridian. Users
   who split their geometries at the antimeridian beforehand never need
   the flag.

Line semantics
--------------

Two points do not by themselves say which line joins them, and the three
answers give different cells. ``line`` names the one wanted:
``"geodesic"`` for the shortest path on the ellipsoid, ``"plate_carree"``
for a line straight in longitude and latitude, and ``"plane"`` for one
straight in the rHEALPix plane. The older ``plane`` flag still selects
between the last two.

A linestring digitised in a geographic CRS usually means geodesic
segments between its vertices, which is what ``line="geodesic"`` traces,
exactly and without densifying. ``wrap_antimeridian`` does not apply to a
geodesic, which already takes the short way, and an antipodal pair raises
``ValueError`` because no single geodesic between them is shortest.

.. figure:: images/line_semantics.*
   :alt: Nine panels. The top row shows one pair of points near 47
         degrees north joined three ways, drawn on the ellipsoid seen
         from above the pole, in plate carree, and in the north polar
         square of the plane. The middle row shows the cells each
         reading traces, in the plane. The bottom row shows the same
         cells back on the ellipsoid.
   :width: 100%

   Each line is straight in exactly one of the three spaces, and the
   panels are those spaces. The plate carree line stays flat near 47
   degrees north while the geodesic and the planar line arc over the
   pole. All three trace sixteen cells here, but not the same sixteen.

.. figure:: images/line_semantics_faces.*
   :alt: Nine panels laid out as the previous figure, for a line from
         Sao Paulo to Berlin. The planar reading leaves the grid part
         way and its cells fall in two disconnected groups.
   :width: 100%

   Once a line crosses faces the readings separate further. The
   geodesic and the plate carree line stay close; the straight planar
   line does not, and part of it leaves the grid's cross-shaped planar
   image altogether, so the cells it returns are not a connected path.
   ``cells_from_line`` warns when that happens.

Attributing length to cells
---------------------------

:meth:`~rhealpixdggs.dggs.RHEALPixDGGS.line_crossings` takes the same
arguments and returns, per cell in the order ``cells_from_line`` gives
them, the fraction of the line's length on the ellipsoid lying in that
cell, so a feature's length can be shared out among the cells it crosses.
A cell visited twice gets an entry per visit, and the fractions sum to 1.

.. currentmodule:: rhealpixdggs.rhp_wrappers

.. autosummary::

    cell_area
    cell_ring
    geo_to_rhp
    k_ring
    linetrace
    polyfill
    rhp_get_base_cell
    rhp_get_resolution
    rhp_is_valid
    rhp_to_center_child
    rhp_to_geo
    rhp_to_geo_boundary
    rhp_to_parent

.. automodule:: rhealpixdggs.rhp_wrappers
    :members:
    :undoc-members:
    :show-inheritance:
