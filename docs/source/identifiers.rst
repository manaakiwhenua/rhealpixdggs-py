Cell Identifiers
===========================

A cell can be named in several ways. Each names a cell only together with
its grid: the same identifier names different cells in grids with a
different ``N_side`` or different polar squares, so data exchanged by
identifier should say which grid it uses, for instance by its
:attr:`~rhealpixdggs.dggs.RHEALPixDGGS.uri`.

Representations
---------------

.. list-table::
   :widths: 22 38 40
   :header-rows: 1

   * - Representation
     - Example (``WGS84_003_OGC``)
     - Notes
   * - Index string
     - ``'R553110'``
     - A face letter and one digit per resolution, for ``N_side`` 2 and 3
       only. The registered rHEALPix DGGRS's text identifier.
   * - suid
     - ``('R', 5, 5, 3, 1, 1, 0)``
     - The same as a tuple, for every grid.
   * - Position
     - ``(6, 1080, 2847)``
     - Level, row and column among all cells of that level in the planar
       grid, counted from the top left.
   * - uint64
     - ``3458765673461713695``
     - The position packed as ``level << 59 | row << 30 | column``. The
       registered rHEALPix DGGRS's integer identifier, for every grid.
   * - Level-row-col text
     - ``'6-438-B1F'``
     - The position in hexadecimal, the text form of OGC 2D tile matrix
       sets. Not an identifier the registered DGGRS declares.

Converting between them
-----------------------

For one cell, each representation is an attribute of the
:class:`~rhealpixdggs.cell.Cell`, and the grid builds the cell back:

.. list-table::
   :widths: 22 39 39
   :header-rows: 1

   * - Representation
     - From a cell
     - To a cell
   * - Index string
     - ``str(cell)``
     - :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.parse_index`, then
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell`
   * - suid
     - ``cell.suid``
     - :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell`
   * - Position
     - :attr:`~rhealpixdggs.cell.Cell.level_row_col`
     - :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell_from_level_row_col`
   * - uint64
     - :attr:`~rhealpixdggs.cell.Cell.uint64`
     - :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell_from_uint64`
   * - Level-row-col text
     - :attr:`~rhealpixdggs.cell.Cell.level_row_col_text`
     - :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell_from_level_row_col_text`

For many cells, the grid converts arrays without building cells, in under
a microsecond per cell for the integer forms. Positions are the hub, three
int64 arrays of level, row and column:

- index strings and positions:
  :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.level_row_col_from_indices`,
  :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.indices_from_level_row_col`;
- positions and uint64, for every grid:
  :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.uint64_from_level_row_col`,
  :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.level_row_col_from_uint64`;
- index strings and uint64 directly:
  :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.uint64_from_indices`,
  :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.indices_from_uint64`;
- index strings and level-row-col text:
  :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.level_row_col_text_from_indices`,
  :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.indices_from_level_row_col_text`.

:func:`~rhealpixdggs.rhp_wrappers.rhp_to_uint64` and
:func:`~rhealpixdggs.rhp_wrappers.rhp_from_uint64` convert single index
strings, returning None for invalid input like the other wrappers.

Sub-zone orders
---------------

Data laid out by position rather than by identifier, as in OGC API - DGGS
and DGGS-JSON, depends on the order of a cell's subcells. The orders take
their names from OGC's DGGRS definitions: ``"mortonCurve"`` (the order of
the digits), ``"scanline"`` (rows across the whole cell), ``"hilbertCurve"``
(``N_side`` a power of 2) and ``"spiralFromCenter"`` (odd ``N_side``).
:attr:`~rhealpixdggs.dggs.RHEALPixDGGS.sub_zone_order` is the order a grid
declares, scanline for the registered grid and Morton for every other, and
the default of :meth:`~rhealpixdggs.cell.Cell.subcells`,
:meth:`~rhealpixdggs.cell.Cell.sub_zone_index` (a cell's position among an
ancestor's subcells) and :meth:`~rhealpixdggs.cell.Cell.sub_zone_at` (the
subcell at a position).

The identifiers module
----------------------

.. automodule:: rhealpixdggs.identifiers
    :members: CELLS0, INDEX_STRING_N_SIDES, SubZoneOrder, SUB_ZONE_ORDERS
