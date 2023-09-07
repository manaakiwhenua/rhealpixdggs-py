=== 0.5.1 ===
- v0.5.1, 2013-07-29: Rounded floats in doctests.

=== 0.5 ===
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
