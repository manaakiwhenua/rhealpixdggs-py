0.6.1
^^^^^
**Breaking change (harmless):** removed the ``RhealPolygon`` stub class
(``__init__`` only, no other methods, no references or tests anywhere).
Follow-up from a whole-repository code review; see ``CODE_REVIEW.md``.

0.6.0
^^^^^
**Breaking change:** ``Cell.ellipsoidal_shape`` is now a ``cached_property``
instead of a method. Replace any calls to ``cell.ellipsoidal_shape()`` with
``cell.ellipsoidal_shape``.

Performance improvements (issue #7): ``RHEALPixDGGS.healpix()`` and
``RHEALPixDGGS.rhealpix()`` now cache their ``Projection`` objects so the
construction cost is paid only once per DGGS instance per region.
``Cell.ellipsoidal_shape`` is computed once per ``Cell`` instance and cached,
eliminating repeated work in ``vertices()``, ``boundary()``, ``neighbors()``,
and other hot paths. ``Cell.boundary(plane=False)`` now short-circuits to
``vertices()`` for quad and cap cells, whose edges are already well-represented
by their four vertices on the ellipsoid.

0.5.16
^^^^^^
Added ``RHEALPixDGGS.area_error_budget()``: returns the theoretical equal cell
area and floating-point tolerance for equal-area testing at each resolution,
derived analytically from machine epsilon (issue #19).

0.5.15
^^^^^^

Fixed ``get_finest_containing_cell`` returning ``None`` for polygons inside cap
cells (issue #13). Cap cell boundaries cannot be represented as simple polygons
in geographic coordinates, so containment is now tested in the rHEALPix plane
for those cells.

0.5.14
^^^^^^
Added ``compact_cells`` to ``conversion``: takes a single-resolution iterable of
cell IDs and returns a fully compacted ``set``, recursing until no further merges
are possible. Added ``recursive`` keyword argument to ``compress_order_cells``
(default ``False``) for the same behaviour with sorted-list output. Updated
``polyfill`` to use ``compact_cells`` when ``compress=True``.

0.5.13
^^^^^^
Normalise floats and correct return types (rhp_wrappers).

0.5.12
^^^^^^
Fixes bug (off-by-one error) in `RHEALPixDGGS.rhp_to_center_child` function in the wrapper API.

0.5.11
^^^^^^
Additional doctests added wrapper API; and changes made for packaging this repository for Conda.

0.5.10
^^^^^^
Changes to suport conda packaging and non-Linux build environments

0.5.9
^^^^^
Changes to suport conda packaging and non-Linux build environments

0.5.8
^^^^^
Releaxes dependency constraints on NumPy and SciPy and other packages to allow for the latest compatible versions of these packages to be used.

0.5.7
^^^^^
Re-release of 0.5.6 due to broken PyPI upload. This version is identical to 0.5.6.

0.5.6
^^^^^
Added support for a "linetrace" algorithm. This is a new algorithm that allows for the representation of lines on the sphere using a series of cells. Functions: `RHEALPixDGGS.cells_from_line` and `linetrace` in the wrapper API.

Implementation of cells_from_line is incomplete. It finds the correct path along cells and across cube faces, but only if there's no cap cell involved. Lines traversing or starting/ending in a cap cell may find the correct neighbour cell to continue to, but there's no guarantee. It's best to avoid having a line intersect with a cap cell (e.g. by increasing the resolution) until those bugs have been fixed.

As with polyfill, the new features can't handle lines that cross the antimeridian. Cell vertices that touch it are resolved correctly in most cases.

0.5.5
^^^^^
Added a series of wrapper functions that expose some of the rHEALPix functionality in an API roughly
in line with what Uber's H3 provides.

0.5.4
^^^^^
Code unchanged from 0.5.3 other than updating to Python 3.11

0.5.3
^^^^^
- v0.5.3: 0.5.1 with:
-- tests that pass
-- Sage removed
-- Test fixes where just relating to rounding errors.

- dggs.py:
-- RG, 2020-09-08: Issue #6 Added optional region="none" arg to rhealpix projection calls, and
                            forced region to cell.region() in cell.vertex() and cell.boundary()
-- RG, 2020-07-31: Issue #5 Moved plot_cells to GRS2013 to remove sage dependence
-- Robert Gibb (RG), 2020-07-13: Issue #1 Multiple tests fail due to rounding errors

- ellipsoids.py:
-- Robert Gibb (RG), 2020-07-13: Issue #1 Multiple tests fail due to rounding errors

- pj_healpix.py:
-- RG, 2020-09-08: Issue #6 In in_healpix_image added +-eps to the extreme corner vertices
                            added calling function abbrev to error statements                            
-- RG, 2020-07-31: Issue #5 Moved healpix_diagram to GRS2013 to remove sage dependence
-- Robert Gibb (RG), 2020-07-13: Issue #1 Multiple tests fail due to rounding errors

- pj_rhealpix.py:
-- RG, 2020-09-08: Issue #6 Added optional region="none" arg to all projection calls, and
                            used it to enforce region behaviour in calls to combine_triangles
                            added calling function abbrev to error statements      
-- RG, 2020-07-31: Issue #5 Moved rhealpix_diagram to GRS2013 to remove sage dependence
-- Robert Gibb (RG), 2020-07-13: Issue #1 Multiple tests fail due to rounding errors

- projection_wrapper.py:
-- RG, 2020-09-08: Issue #6 Added optional region="none" arg to all projection calls
-- Robert Gibb (RG), 2020-07-13: Issue #1 Multiple tests fail due to rounding errors

- utils.py:
-- Robert Gibb (RG), 2020-07-13: Issue #1 Multiple tests fail due to rounding errors                                                  

0.5.1
^^^^^
- v0.5.1, 2013-07-29: Port to Python 3.7
- Rounded floats in doctests.

0.5
^^^
- v0.5, 2013-07-26: Port to Python 3.3.

- dggs.py:
-- AR, 2013-03-14: Fixed bug in nw_vertex().
-- AR, 2013-01-16: Changed the string keyword 'surface' to a boolean keyword 'plane'.
-- AR, 2013-01-14: Added intersects_meridian(), cell_latitudes(), cells_from_meridian(), cells_from_parallel(), cells_from_region().
-- AR, 2013-01-01: Added ellipsoidal functionality to neighbor() and neighbors().
-- AR, 2012-12-19: Tested all the methods and added examples.
-- AR, 2012-12-10: Corrected centroid() and moved some methods from graphics.py to here.

- pj_healpix.py:
-- AR, 2013-03-05: In in_healpix_image() increased eps to 1e-10 to decrease out-of-bounds errors i was getting when drawing figures.

Older releases
^^^^^^^^^^^^^^
For releases 0.3, 0.4 & 0.5, see the ``_old`` folder within the `0.5.1 release <https://github.com/manaakiwhenua/rhealpixdggs-py/tree/0.5.1/_old>`_ on GitHub
