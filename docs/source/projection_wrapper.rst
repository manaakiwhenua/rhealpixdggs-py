The projection_wrapper Module
===============================

.. currentmodule:: rhealpixdggs.projection_wrapper

.. autosummary::

    Projection

Array inputs
------------

``Projection.__call__`` accepts numpy arrays as well as floats: the
``lon_0``/``lat_0`` recentring and the longitude and latitude wrapping
apply elementwise and the underlying projection runs once on the whole
batch, returning a pair of float64 arrays. ``RHEALPixDGGS.rhealpix`` and
``RHEALPixDGGS.healpix`` pass arrays straight through.

.. automodule:: rhealpixdggs.projection_wrapper
    :members:
    :undoc-members:
    :show-inheritance:
