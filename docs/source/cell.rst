The cell Module
===========================

.. figure:: images/cell_shapes.*
   :alt: Four orthographic globe views, one per ellipsoidal cell shape.
         A quad is a longitude-latitude rectangle on an equatorial face;
         a cap is the round cell containing a pole; a dart is a
         three-cornered cell whose planar square has a corner at the
         projected pole; a skew quad is a four-cornered polar cell with
         curved edges.
   :width: 100%

   The four ellipsoidal cell shapes (each drawn filled, with its
   neighborhood for context). Cells on the equatorial faces O--R are
   *quads*: genuine longitude-latitude rectangles. Cells on the polar
   faces N and S are *caps* (the cell containing the pole), *darts*
   (three-cornered cells) or *skew quads* (four-cornered cells whose
   edges are not meridians and parallels). A cell's shape is reported by
   :attr:`Cell.ellipsoidal_shape`.

.. figure:: images/predicates.*
   :alt: Four three-by-three planar grids of the children of cell P4,
         illustrating touches with a shared edge, touches at a corner
         only, within/contains, and disjoint.
   :width: 100%

   The spatial predicates on the planar grid, illustrated with the
   children of cell P4: ``touches`` is true for neighbors that share an
   edge *and* for cells that meet at a corner only; ``within`` /
   ``contains_cell`` relate a cell to its ancestors and descendants;
   ``disjoint`` is everything else.

.. figure:: images/predicates_polar.*
   :alt: Four orthographic views. First three, centered on the north
         pole: the cap N4 touching skew quad N1 along an edge, the cap
         touching dart N0 at a single point, and the resolution 2 cap
         N44 within the resolution 1 cap N4. Fourth, centered on a cube
         corner: cells N0, Q2 and R0 all touching one another along
         edges.
   :width: 100%

   The same predicates among polar cell shapes. The cap N4 *touches*
   its edge-neighbor skew quads (here N1) and also touches the darts
   (here N0), which meet it at a single point (marked). Caps nest like
   any other cells: N44 is *within* N4. At a cube corner (right) only
   three cells meet, and every pair of them shares an edge -- unlike a
   planar grid corner, no pair touches at a corner only.

.. figure:: images/cube_corner.*
   :alt: Orthographic view of the globe centered on a corner of the
         rHEALPix cube, where exactly three cells meet: N0, its neighbor
         Q2 above, and its neighbor R0 to the left. A dot marks the
         corner and an annotation notes that N0.diagonal_neighbor
         with argument up_left is None.
   :width: 80%
   :align: center

   A corner of the underlying cube, on the ellipsoid. Away from the
   cube's corners every cell corner is shared by four cells, but at the
   eight cube corners exactly *three* cells meet, so one diagonal
   neighbor does not exist: here ``N0.neighbor('up')`` is Q2 and
   ``N0.neighbor('left')`` is R0, but ``N0.diagonal_neighbor('up_left')``
   is ``None``.

.. figure:: images/cube_corner_cube.*
   :alt: An isometric drawing of the rHEALPix cube showing faces N, Q
         and R, each subdivided into its nine resolution 1 cells. The
         children N0, Q2 and R0 are shaded and meet at the front corner
         vertex, which is marked with a dot.
   :width: 70%
   :align: center

   The same corner drawn on the cube itself. The three shaded cells N0,
   Q2 and R0 each contribute one corner to the shared vertex, and each
   pair of them shares an edge along one of the three creases -- there
   is simply no fourth, diagonally-opposite cell.

.. currentmodule:: rhealpixdggs.cell

.. autosummary::

    Cell

.. automodule:: rhealpixdggs.cell
    :members:
    :undoc-members:
    :show-inheritance:
