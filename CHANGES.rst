0.6.1
^^^^^
**Breaking change:** ``rhp_wrappers.cell_ring()``/``k_ring()`` are rebuilt on
a breadth-first search over ``Cell.neighbor()``/``diagonal_neighbor()``
steps, growing the ring outward one shell at a time, rather than jumping
directly to a computed "start corner" and walking its 4 sides. The previous
approach was badly wrong for any ring touching more than 2 resolution 0
cube faces (documented as a known limitation, issue #60): e.g. a distance-1
ring around a cell sitting at a genuine cube corner used to return several
cells that aren't even neighbors of the centre cell, plus a duplicate. The
new approach handles any number of faces uniformly and correctly, including
genuine 3-valent cube corners, where a ring has fewer than the usual 8*k
cells for an interior ring rather than an incorrect or duplicate one.
Consequences of the rewrite:

- Both functions no longer take a ``verbose`` parameter. It existed to print
  a warning about the now-fixed >2-face limitation; there's nothing left to
  warn about.
- Requesting a ring beyond the grid's actual graph diameter at that
  resolution (the largest number of edge/corner hops between any two cells)
  now correctly returns an empty list, rather than an arbitrary guess at an
  "opposite" cell -- there is no cell at that exact distance once the whole
  grid has already been covered by closer rings.
- The specific cell order within a ring has changed (same cells, differently
  sequenced): rings now always list cells in a fixed clockwise order
  starting from "up", whatever that direction happens to be for the centre
  cell, rather than the previous "walk the perimeter starting from the
  upper-left corner" convention (which didn't generalize to corner cases
  anyway).

Added a ``Cell.diagonal_neighbor(direction)`` method (``direction`` one of
``'up_left'``, ``'up_right'``, ``'down_left'``, ``'down_right'``), returning
this cell's corner-touching-only planar neighbor, or ``None`` if this cell
sits at one of the 8 genuine cube corners in that direction (exactly 3 cells
meet there, not 4, so there's no distinct diagonal cell to return). Also
added a set of DE-9IM-style topological predicates between two cells,
possibly of different resolutions: ``equals``, ``contains_cell``, ``within``,
``covers``, ``covered_by``, ``touches``, and ``disjoint`` (issue #55).
``contains_cell`` is named to avoid confusion with the pre-existing
``contains()``, which tests point-in-cell containment, not cell-in-cell.
Both the new ``diagonal_neighbor`` primitive and ``touches`` (built on it,
for the case where neither cell is an ancestor of the other) were validated
against real vertex/boundary geometry across many DGGS configurations
(``N_side`` 2 through 5, all 16 ``(north_square, south_square)``
combinations) before relying on them here.

Fixed ``Cell.subcell()`` and ``Cell.__le__`` (and so ``<``/``<=``/``>``/``>=``,
and anything built on them, such as ``RHEALPixDGGS.interval()``) comparing a
comma-joined **string** rendering of a cell's suid rather than the suid
tuple itself. This happened to work for the default ``N_side=3``, where
every digit is a single character, but for ``N_side >= 4`` some digits are
multi-character (e.g. 10-15 for ``N_side=4``), and the string comparison
could misfire: e.g. ``('N', 15).subcell(('N', 1))`` incorrectly returned
``True`` (``"N,15"`` starts with ``"N,1"`` as a string, even though digit 15
isn't a descendant of digit 1), and same-resolution siblings like
``('N', 9)``/``('N', 10)`` could compare in the wrong order (``"9" > "1"``
as strings, despite 9 < 10 numerically). Both now compare the suid tuple
directly (issue #71).

**Breaking change:** ``Cell``'s constructor and ``Cell.overlaps()`` validated
their arguments with bare ``assert`` statements, which ``python -O``/
``PYTHONOPTIMIZE`` compiles out entirely (silently disabling the checks),
and which raised ``AssertionError`` for what is really invalid input.
Replaced with explicit ``TypeError``/``ValueError`` matching the kind of
problem: ``TypeError`` if ``suid`` isn't a list or tuple, ``ValueError``
for a bad length, an invalid first character, an out-of-range digit, or
calling ``overlaps()`` on an empty cell. If you were catching
``AssertionError`` from these, catch ``TypeError``/``ValueError`` instead
(issue #54).

Fixed ``Cell.neighbors(plane=False)`` for dart and skew_quad cells, which
computed east/west (and, for darts, southeast/southwest/etc.) by
temporarily reassigning ``self.rdggs.ellipsoid.lon_0`` -- a singleton
shared by every ``Cell`` built from the same ``rdggs`` -- to recentre the
prime meridian on this cell's nucleus, then restoring it afterwards. This
was not thread-safe (concurrent use of the same ``rdggs`` could observe or
clobber the temporarily-shifted value) and not exception-safe (an
exception between the mutation and the restore left the shared ellipsoid
permanently corrupted). It also turned out not to be correct in every
case: confirmed by testing across many DGGS configurations, when the
ellipsoid's own ``lon_0`` happens to be near +-180 degrees, the old
approach could assign "east" and "west" (and their compound directions)
backwards for the same physical cell, purely as a function of where the
prime meridian was drawn -- which should never affect compass-direction
labelling. Replaced with a comparison of each neighbor's longitude
relative to this cell's own nucleus, wrapped to avoid the antimeridian
issue directly (``rhealpixdggs.utils.wrap_longitude``), without touching
any shared state at all (issue #53).

**Breaking change:** ``healpix_sphere_inverse``, ``healpix_ellipsoid_inverse``,
``rhealpix_sphere_inverse``, and ``rhealpix_ellipsoid_inverse``
(``pj_healpix``/``pj_rhealpix``) now raise ``ValueError`` for out-of-bounds
input coordinates, instead of printing an error message and returning a
sentinel: ``float("inf")`` for ``healpix_sphere_inverse``, or bare ``None``
for the other three. The ``None`` case was already effectively a crash for
any caller going through the public ``healpix()``/``rhealpix()`` projection
factories -- their closures do ``array(...)`` then unpack the result, which
raises a confusing, unrelated-looking ``TypeError: iteration over a 0-d
array`` a few frames away from the actual problem. If you were catching
these failures by checking for ``None``/``inf``, catch ``ValueError``
instead (issue #52).

Fixed ``compact_cells``/``compress_order_cells`` (``conversion``) hardcoding
9 as the size of "a complete group of siblings" -- correct only for the
default ``N_side=3``. Under an ``N_side=2`` DGGS (e.g. ``WGS84_002``, 4
children per cell), a genuinely complete sibling group silently never
compressed at all. Both functions now take an ``N_side`` keyword argument
(default 3, so existing callers are unaffected) and use ``N_side ** 2``.
``polyfill()``'s ``compress=True`` path now passes its ``dggs``
parameter's actual ``N_side`` through (issue #51).

Fixed ``Cell.boundary(plane=False)`` silently returning only 4 points for
quad/cap cells regardless of ``n``, instead of the documented ``4*n - 4``
(a regression from 0.6.0's short-circuit to ``vertices()``). The
short-circuit now only applies at ``n=2``, where it's provably equivalent;
``n>2`` on quad/cap cells falls through to the same per-point-projected
algorithm used for other shapes, so it is correct again but not faster
than before 0.6.0 for that case (issue #49).

Performance (issue #7, #62): ``Projection.__call__`` was rebuilding its
underlying projection callable (re-importing its module and, for homemade
projections, recomputing the authalic radius) on every single point
projected, rather than once -- now cached. No behavior change -- same
inputs still produce identical outputs, verified against the full test
suite -- only less redundant setup work per point.

Dropped the ``matplotlib`` dependency (issue #64). It was only ever used
for a point-in-polygon check inside ``in_healpix_image``/
``in_rhealpix_image`` -- a vestige of plotting code that was moved out of
this package back in 0.5.3 -- and is now replaced with an equivalent,
cached ``shapely`` check (``shapely`` was already a required dependency,
used elsewhere in ``conversion.py``). No behavior change, verified
against the full test suite including every doctest that pins exact
boundary cases.

Combined, these two caching/dependency fixes measured on a real-world
repro (``Cell.boundary(10, plane=False)`` on a dart-shaped cell, 2000
iterations, issue #7's own reported case): better than 2x faster
(3.62s -> ~1.5-2.1s depending on inputs), and it fixes the correctness
issue at the same time -- not just unchanged despite dropping a
dependency, faster.

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
