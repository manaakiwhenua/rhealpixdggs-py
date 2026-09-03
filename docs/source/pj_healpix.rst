The pj_healpix Module
===============================

.. currentmodule:: rhealpixdggs.pj_healpix

.. autosummary::

    healpix
    healpix_ellipsoid
    healpix_ellipsoid_inverse
    healpix_sphere
    healpix_sphere_inverse
    in_healpix_image

Array inputs
------------

The callable returned by ``healpix`` accepts numpy arrays as well as floats:
pass two arrays of a common shape and every point is projected (or
inverted) in one array pass, returning a pair of float64 arrays. Scalar
calls are unchanged. If any point of an inverse batch lies outside the
projection's image the whole call raises ``ValueError``, naming the count
and the first offending point.

.. automodule:: rhealpixdggs.pj_healpix
    :members:
    :undoc-members:
    :show-inheritance:
