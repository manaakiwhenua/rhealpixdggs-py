The conversion Module
===========================

.. figure:: images/compaction.*
   :alt: Two maps of the same polygon over New Zealand filled with
         cells. Left: 518 cells, all at resolution 5. Right: the same
         coverage after compact_cells, 182 cells of mixed resolutions,
         with large low-resolution cells in the interior and small
         high-resolution cells along the edges.
   :width: 100%

   ``compact_cells``: a uniform resolution 5 ``polyfill`` of a polygon
   (left) and the same set of cells after compaction (right), where any
   nine siblings all present in the set are replaced by their parent,
   repeatedly. The interior collapses to coarse cells while the edge
   detail is preserved exactly -- here 518 cells become 182. The change
   of texture partway down is the boundary between the equatorial Q
   cells (longitude-latitude rectangles) and the south polar S cells
   (skew quads and darts, which are not).

.. currentmodule:: rhealpixdggs.conversion

.. autosummary::

    CellZoneFromPoly
    compact_cells
    compress_order_cells
    get_finest_containing_cell

.. automodule:: rhealpixdggs.conversion
    :members:
    :undoc-members:
    :show-inheritance:
