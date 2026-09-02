0.7.1
^^^^^
Adopted mypy (issue #109, first three stages): every function in the
package is now annotated and checked in CI (``check_untyped_defs`` plus
``disallow_untyped_defs``). Small behavioural edges tightened along the
way: comparing a ``Cell`` against a non-``Cell`` returns ``False``
instead of raising ``AttributeError``; ``Cell()`` without a DGGS raises
``TypeError`` immediately (it could never construct a usable cell); and
several operations that would crash on the empty cell now fail with an
explicit assertion. ``Cell.suid`` and
``Cell.resolution`` gained honest attribute types
(``tuple[str | int, ...]`` and ``int | None`` — None only for the empty
cell). Annotations that claimed more than the code
delivers were corrected — e.g. ``cell_from_point``, ``cell_from_region``,
``get_finest_containing_cell``, ``rhp_get_resolution`` and
``rhp_get_base_cell`` are now declared ``... | None``, matching their
long-documented behaviour, and ``pj_rhealpix.triangle`` is declared to
return ``None`` for equatorial points as its docstring always said. Three
edge-case behaviour fixes fell out of the same pass:
``RHEALPixDGGS.minimal_cover`` now skips points that lie outside the
planar grid instead of collecting a spurious ``None`` entry,
``cells_from_region`` returns ``[]`` when any corner of a planar region
falls outside the grid instead of raising ``AttributeError``, and
``cell_from_region`` likewise returns ``None`` instead of raising when a
corner is off-grid.
**Breaking change:** removed the plotting-era point generators, none of
which has had a caller in any tree in this repository's history (issue
#110): ``Ellipsoid.get_points`` (which has also returned ``[]`` for any
input for the whole of the project's Python 3 life, due to a bytes-vs-str
comparison), ``Ellipsoid.lattice``, ``Ellipsoid.graticule``,
``Ellipsoid.meridian``, ``Ellipsoid.parallel``,
``pj_healpix.healpix_vertices`` and ``pj_rhealpix.rhealpix_vertices``.
They predate the DGGS layer; cell geometry (``Cell.ul_vertex``,
``Cell.width``, ``Cell.boundary``, ``RHEALPixDGGS.cell_boundaries``) is
the supported way to draw the grid and the projection image, as
``docs/make_figures.py`` demonstrates.

0.7.0
^^^^^
Added ``RHEALPixDGGS.cell_boundaries(cells, n=2, plane=True)``: boundary
points for a whole set of cells at once, projecting every shared boundary
point once instead of once per adjacent cell (issue #87, proposed in
issue #7). For contiguous blocks of cells this roughly halves the
projection work of rendering cell outlines, and it makes adjacent cells'
copies of their shared points identical floats rather than two
independently computed values. Each cell's entry agrees with its own
``boundary(n, plane)`` in count, order, and coordinates. ``Cell`` is now
hashable (consistent with its equality: cells with equal addresses from
the same DGGS hash alike), so cells can be dictionary keys and set
members.


**Breaking change:** rewrote ``RHEALPixDGGS.cells_from_line`` (the engine
behind ``rhp_wrappers.linetrace``) to be exact for every cell shape,
replacing an edge-walking algorithm that modelled each cell as a
straight-edged quadrilateral -- wrong for polar cells (cap cells are
circles in longitude-latitude coordinates; dart and skew-quad edges are
curves), forbidden from ever revisiting a cell (correct polar paths do),
and equipped with a fail-safe that silently jumped to the end cell,
skipping arbitrarily many cells. Confirmed failures included skipping a
cell a segment crossed straight through the middle of, and polar traces
missing 20 cells at once. The new implementation sweeps the segment's
planar image for cell-edge crossings: the crossing parameters are exact
(piece boundaries -- region and polar-triangle changes and the longitude
wrap -- are known in closed form, and every crossing is bracketed and
solved to machine precision), and each inter-crossing interval midpoint
is located with the shape-exact ``cell_from_point``. Validated against
an independent dense point-location oracle on hundreds of random
segments (equatorial, polar, region-crossing, planar): zero
disagreements, and the sweep catches sliver crossings as short as 1/500
of a cell that sampling misses. Consequences:

- The cap-cell limitation is gone: ``tests``' known-failure case now has
  its correct expected sequence (confirmed leg-by-leg by dense point
  location) and passes as a regular test. ``LINETRACE_WARNING`` and the
  warning emission are removed from ``rhp_wrappers.linetrace``
  (``verbose`` still controls the malformed-geometry warnings).
- Some previously returned sequences change: cells the old algorithm
  skipped are now included (e.g. a cell a segment crossed through the
  middle of), and cells it invented are gone. Cells met in a single
  point only (a segment passing exactly through a cell corner) are not
  included.
- Longitude-latitude segments are documented as straight in coordinate
  space (plate carree, matching the shapely input type and GeoJSON RFC
  7946's planar reading): segments spanning more than half a turn of
  longitude run the long way around by default rather than wrapping.
  Splitting inputs at the antimeridian, as RFC 7946 prescribes, avoids
  the ambiguity; the new ``wrap_antimeridian`` keyword (on
  ``cells_from_line`` and ``linetrace``, default False) traces such
  segments the short way instead, equivalent to splitting. This
  replaces the previous blanket "cannot handle the antimeridian"
  caveat, and ``polyfill``'s antimeridian TODO is likewise resolved by
  documenting the split-first requirement.
- **Breaking change (harmless):** removed
  ``RHEALPixDGGS.antimeridian_check_and_flip``, which existed solely to
  patch cell edges for the old algorithm's shapely intersection tests.


Adopted the organisation's repository-standards check
(``manaakiwhenua-standards``) as a GitHub Actions workflow, with its
status badge in the README (issue #84). The check requires a non-empty
``README.md``, so the README is converted from reStructuredText to
Markdown (``README.rst`` -> ``README.md``, with ``pyproject.toml``'s
``readme`` field updated; PyPI renders Markdown long descriptions
natively -- verified in the built wheel's metadata). The README's
"Release Notes" section, which duplicated part of ``CHANGES.rst`` inline
and claimed the package "has had only minor code updates" since 2013, is
replaced by a pointer to ``CHANGES.rst``.

Packaging and release-process fixes (issue #58): the ``pyproject.toml``
license field now declares the actual dual license,
``LGPL-3.0-or-later OR MIT`` -- it previously declared
``GPL-3.0-or-later``, which misrepresented the project as *more*
restrictively licensed than it is (the kind of error automated license
scanners act on); built wheels now carry the correct PEP 639
``License-Expression`` plus all three license files. Robert Gibb joins
the declared authors and the current maintainer is declared in the
``maintainers`` field. Added GitHub Actions CI running the unit tests
and doctests on Python 3.11 through 3.14 plus a docs build, with
matching Python version classifiers and a README badge; fixed
``run_doctests.sh`` to exit nonzero when any file's doctests fail rather
than only the last one's (a latent flaw that would have let CI pass on
failing doctests). Replaced the outdated ``pypi.txt`` (it described a
``setup.py``/``twine`` release flow the project no longer has) with
``RELEASING.md``, and noted that working with this repository's
packaging requires Poetry >= 2.0. Fixed the sdist sweeping transient
files under ``docs/`` past gitignore (Poetry's include patterns take
priority); the figure generator's download cache now lives at the
repository root. Enforcement of ``black`` formatting is deliberately
deferred (issue #82).

``RELEASING.md``'s checklist is now carried out by ``scripts/release.py``,
which stages the release (preflight checks, version bump, tests, build,
tag, upload, GitHub release) so that each step that is visible outside
the machine has to be asked for by name. Beyond the checklist it verifies
the built wheel's metadata and that the sdist has not swept up build
leftovers. The two test runners move alongside it: ``run_unittests.sh``
and ``run_doctests.sh`` are now ``scripts/run_unittests.sh`` and
``scripts/run_doctests.sh``, leaving the repository root to project
metadata, and both now work from any directory rather than only from the
repository root.


Documentation fixes (issue #57): the Sphinx manual now includes pages for
the ``cell``, ``conversion``, and ``rhp_wrappers`` modules, which were
missing from its toctree (so the entire H3-style wrapper API, among other
things, never appeared in the built manual); ``docs/source/conf.py`` now
reads the release number from ``pyproject.toml`` instead of pinning a
stale copy (it said 0.5.4); the tracked copies of stock Sphinx theme
assets under ``docs/source/_static`` -- which overrode the theme's own,
newer versions at build time -- are removed; and the README's
requirements section now matches ``pyproject.toml`` (it referenced
``requirements.txt`` files that no longer exist and a ``numpy<2`` bound
that contradicted the actual ``numpy>=2.0`` requirement), alongside
fixes to its broken ``python install`` command and an example naming a
test file that doesn't exist.

Fixed ``Cell.centroid(plane=False)`` for quad cells, whose centroid
latitude was computed as the midpoint of the cell's two edge latitudes
rather than the area-weighted mean latitude the centroid is defined as
(the average value of latitude over the cell). Latitude is a nonlinear
function of planar y, so the two differ -- by up to ~0.63 degrees for
resolution 1 quad cells, shrinking quadratically with cell size, and
vanishing for cells symmetric about the equator. The defect traces back
to the founding paper's section 7 summary table, whose quad-cell entry
contradicts the paper's own integral definition of the centroid directly
above it; dart and skew_quad cells always used the integral and were
unaffected, as were all centroid longitudes and everything planar. Quad
centroid latitudes are now computed by (effectively exact) fixed-order
Gauss-Legendre quadrature of the mean-latitude integral, validated
against adaptive quadrature at a different meridian and an independent
equal-area Monte Carlo estimate. Affects anything consuming ellipsoidal
quad centroids, notably ``rhp_wrappers.rhp_to_geo(..., plane=False)``
and ``rhp_wrappers.polyfill(..., plane=False)`` membership decisions for
quad cells whose true and previously-reported centroids straddle the
polygon boundary. One observable side effect: an exactly-symmetric
cell's centroid latitude may now come out as floating-point dust (e.g.
6e-16) rather than exactly 0 (issue #75).

Test-coverage additions (issue #56), no library behavior changes:
``utils.auth_lat()``'s large-flattening branch (``f > 1/150``, unreachable
by any predefined ellipsoid) is now validated against the numerically
integrated defining integral of authalic latitude, alongside the power
series branch; ``projection_wrapper`` gained functional tests
(forward/inverse round trips for homemade and PROJ-backed projections,
``lon_0``/``lat_0`` recentring, projection keyword pass-through -- note
that ``lat_0`` round trips cleanly without the resulting grid being
sound, see the known limitation at the end of this entry);
``Cell.area()``, ``Cell.color()``, ``Cell.region_overlaps()``, and
``RHEALPixDGGS.area_error_budget()`` gained their first tests (including
checking that every resolution's cells sum to the authalic sphere's
surface area); and ``Cell.centroid()`` gained its first active test (the
previous one was a disabled stub), validating the dart/skew_quad
integration against an independent fixed-seed Monte Carlo estimate --
which flagged that quad cells' centroid latitudes are computed
incorrectly, now tracked separately as issue #75. Also renamed
``tests/test_PJ_healpix.c.py`` (whose embedded extra ``.c`` suffix made
``unittest`` discovery silently skip its 4 pyproj/PROJ cross-checks) to
``tests/test_pj_healpix_c_reference.py``, and repaired it: it had rotted
unrunnable (removed ``scipy`` re-exports of numpy functions, a missing
``itertools.product`` import, scalar inputs to
``scipy.spatial.distance.euclidean``), and its ellipsoidal inverse
round-trip tolerances are now geographic (pole/antimeridian aware) and
two-tier -- a precision check at WGS84 eccentricity, a documented
smoke-test bound at the file's deliberately extreme ``e = 0.8``, where
PROJ's own low-order inverse authalic series is only good to ~2 degrees.

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

**Known limitation, not fixed in this release:** an ellipsoid with a
nonzero ``lat_0`` does not give a coherent DGGS. ``Projection.__call__``
recentres by translating the input coordinates, which is a true rotation
of the ellipsoid for ``lon_0`` but not for ``lat_0``: a latitude pushed
past a pole is reflected to the antipodal latitude at the *same*
longitude, instead of continuing over the pole and down the far side.
Any ``lat_0 != 0`` therefore folds a band of width ``|lat_0|`` around one
pole, which collides the two poles onto a single planar point, leaves
cells covering two disjoint geographic patches (one near each pole), and
destroys the equal-area property -- measured at 5.1x the uniform share
for the worst resolution 1 cell under ``lon_0=174, lat_0=-37``. Forward
and inverse round trips still succeed, because the reflection is applied
symmetrically in both directions, which is why the defect stayed
invisible. None of the predefined ellipsoids uses a nonzero ``lat_0``,
so only user-constructed ellipsoids are affected, and recentring by
``lon_0`` alone is sound. Note that the ``dggs`` module docstring's
Auckland example does construct its ellipsoid with ``lat_0=-37``; use
``lon_0`` only. This is not a regression: the same recentring appears in
rHEALPixDGGS 0.4, the earliest version of the code in this repository's
history. Tracked as issue #93, milestoned for 0.8.0.

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
