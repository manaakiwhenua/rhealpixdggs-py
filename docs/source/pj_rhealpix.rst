The pj_rhealpix Module
===============================

.. currentmodule:: rhealpixdggs.pj_rhealpix

.. autosummary::

    combine_triangles
    in_rhealpix_image
    rhealpix
    rhealpix_ellipsoid
    rhealpix_ellipsoid_inverse
    rhealpix_sphere
    rhealpix_sphere_inverse
    triangle

Array inputs
------------

The callable returned by ``rhealpix`` accepts numpy arrays as well as floats:
pass two arrays of a common shape and every point is projected (or
inverted) in one array pass, returning a pair of float64 arrays. Scalar
calls are unchanged. If any point of an inverse batch lies outside the
projection's image the whole call raises ``ValueError``, naming the count
and the first offending point.

.. automodule:: rhealpixdggs.pj_rhealpix
    :members:
    :undoc-members:
    :show-inheritance:
