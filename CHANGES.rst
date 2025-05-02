0.5.6
^^^^^
Added support for a "linetrace" algorithm. This is a new algorithm that allows for the representation of lines on the sphere using a series of cells. Functions: `RHEALPixDGGS.cells_from_line` and `linetrace` in the wrapper API.

Implementation of cells_from_line is incomplete. It finds the correct path along cells and across cube faces, but only if there's no cap cell involved. Lines traversing or starting/ending in a cap cell may find the correct neighbour cell to continue to, but there's no guarantee. It's best to avoid having a line intersect with a cap cell (e.g. by increasing the resolution) until those bugs have been fixed.

Just like with polyfill, the new features can't handle lines that cross the antimeridian. Cell vertices that touch it are resolved correctly in most cases when checking whether an edge intersects the line, however.

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
