Standards conformance
=====================

rHEALPixDGGS implements code to both define an rHEALPix DGGS Reference
System and to perform topological queries on its identifiers. The
project's roadmap is for version 1.0 to be fully compliant with
`OGC Abstract Specification Topic 21 version 2.0 (OGC 20-040r3)
<https://docs.ogc.org/as/20-040r3/20-040r3.html>`_ /
ISO 19170-1:2020, *Discrete Global Grid Systems — Part 1: Core
Reference System and Operations, and Equal Area Earth Reference
System*.

The conformance target is the specification's strongest conformance
class, **Equal-Area Earth DGGS**, which comprises the Common
Spatio-temporal Classes tests A.2, A.3 and A.5, all DGGS Core tests
(A.6–A.19) and all Equal-Area Earth Reference System tests
(A.20–A.29) of the specification's Annex A. Every test in Annex A is
of type *Basic*, assessed by inspection of the DGGS specification's
documentation — so this page is both the project's self-assessment
and, as the remaining rows are completed, the conformance
documentation itself.

Scope notes:

- The purely temporal tests A.1 and A.4 belong to the separate
  *Spatio-temporal DGGS Core* conformance class and are not required
  for a spatial DGGS; rHEALPixDGGS does not claim them.
- The interoperation modules deliberately stop at the library
  boundary: rHEALPixDGGS provides query execution and standard-format
  serialization, while transport protocols (HTTP and friends, e.g. an
  `OGC API — DGGS <https://ogcapi.ogc.org/dggs/>`_ service) belong to
  software built on top of it.

Each row below records a status:

- **Met** — the substance is implemented and stated here or in the
  linked API documentation.
- **Partial** — implemented in substance, with named work remaining.
- **Planned** — not yet implemented; the release target is named.

Common Spatio-temporal Classes
------------------------------

.. list-table::
   :widths: 16 34 10 40
   :header-rows: 1

   * - Test
     - What it requires
     - Status
     - Evidence and plan
   * - A.2 zone geometry
     - Zonal geometry and topology conforming to the common data
       model (zones with boundaries, neighborhoods and hierarchy).
     - Partial
     - :meth:`~rhealpixdggs.cell.Cell.boundary`,
       :meth:`~rhealpixdggs.cell.Cell.vertices`,
       :meth:`~rhealpixdggs.cell.Cell.neighbors`; the explicit mapping
       onto the model's classes is to be written (v1.0).
   * - A.3 spatial location
     - Locations identified by reference to zones (ISO 19112 style).
     - Partial
     - :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell_from_point`,
       :func:`~rhealpixdggs.rhp_wrappers.geo_to_rhp`; mapping
       documentation to be written (v1.0).
   * - A.5 RS using zonal identifiers
     - A reference system whose zones carry unique identifiers, a
       parent/child hierarchy, an extent and a representative
       position.
     - Partial
     - The identifiers are specified in `Zonal identifiers`_ below.
       :func:`~rhealpixdggs.rhp_wrappers.rhp_to_parent` and
       :meth:`~rhealpixdggs.cell.Cell.subcells` give the hierarchy,
       :meth:`~rhealpixdggs.cell.Cell.boundary` the extent and
       :meth:`~rhealpixdggs.cell.Cell.centroid` the representative
       position (A.27; the :meth:`~rhealpixdggs.cell.Cell.nucleus` is
       the indexing point, interior by construction). The mapping onto
       the model's classes is to be written (v1.0).

DGGS Core — reference system
----------------------------

.. list-table::
   :widths: 16 34 10 40
   :header-rows: 1

   * - Test
     - What it requires
     - Status
     - Evidence and plan
   * - A.6 harmonized model
     - The reference system's definition conforms to the Core data
       model (grids, refinement strategy, globe geometry).
     - Partial
     - All defining parameters live on
       :class:`~rhealpixdggs.dggs.RHEALPixDGGS` (ellipsoid, N_side,
       north_square, south_square, max_resolution); the formal
       declaration against the model's classes is to be written
       (v1.0).
   * - A.7 CRS
     - A defined coordinate reference system, with coordinate-epoch
       provisions.
     - Met
     - :attr:`~rhealpixdggs.dggs.RHEALPixDGGS.crs` is the projected CRS
       in WKT2 and PROJJSON; `Coordinate reference system`_ below
       gives the definition and the coordinate-epoch statement.
   * - A.8 global domain
     - A specified global domain and its dimensionality.
     - Met
     - The domain is the entire surface of the ellipsoid of
       revolution; spatial dimension 2, no temporal dimension.
   * - A.9 domain completeness
     - The level 0 grid covers the entire global domain.
     - Met
     - The six resolution 0 cells (N, O, P, Q, R, S) partition the
       plane image of the projection, whose inverse covers the whole
       ellipsoid.
   * - A.10 location uniqueness
     - Every location lies in exactly one level 0 cell.
     - Met
     - Cells own their upper and left boundaries;
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell_from_point` is
       deterministic and regression-tested at cell edges.
   * - A.11 simple cells
     - Cells at every level are geometrically simple.
     - Met
     - All cells are simple polygons; on the ellipsoid they take one
       of four shapes (quad, cap, dart, skew quad; see
       :attr:`~rhealpixdggs.cell.Cell.ellipsoidal_shape`).
   * - A.12 direct position
     - Every zone is assigned a direct position inside its boundary.
     - Met
     - The direct position is the centroid (A.27),
       :meth:`~rhealpixdggs.cell.Cell.centroid`, which lies inside its
       cell for every cell shape (tested in ``tests/test_cell.py``); the
       indexing point :meth:`~rhealpixdggs.cell.Cell.nucleus` is interior
       by construction.
   * - A.13 address
     - Every zone has a globally unique identifier, structured on a
       recognized indexing method.
     - Met
     - `Zonal identifiers`_ below gives the grammar, the validity
       rules and the codespace that makes an identifier globally
       unique. Both zone identifier reference systems the
       OGC-registered rHEALPix DGGRS declares are supported: the text
       one (index strings) and the uint64 one
       (:attr:`~rhealpixdggs.cell.Cell.uint64`,
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell_from_uint64`).
   * - A.14 discrete global grid
     - Grids defined as the aggregation of all cells at one level.
     - Met
     - :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.grid` enumerates the
       cells of any resolution.
   * - A.15 grid sequence
     - The grid hierarchy is ordered by increasing refinement.
     - Met
     - Resolutions are the integers 0 through max_resolution;
       every API is indexed by them.

DGGS Core — functions
---------------------

.. list-table::
   :widths: 16 34 10 40
   :header-rows: 1

   * - Test
     - What it requires
     - Status
     - Evidence and plan
   * - A.16 quantization
     - Defined operations assigning external data to cells, conforming
       to the quantization model.
     - Partial
     - Vector quantization exists:
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell_from_point`,
       :func:`~rhealpixdggs.rhp_wrappers.polyfill`,
       :func:`~rhealpixdggs.rhp_wrappers.linetrace`,
       :class:`~rhealpixdggs.conversion.CellZoneFromPoly`,
       :func:`~rhealpixdggs.conversion.compact_cells`. The mapping of
       these onto the model's quantization roles is to be written
       (v1.0).
   * - A.17 zone queries
     - The full DE-9IM query interface between zones: contains,
       within, touches, disjoint, equals, intersects, overlaps,
       crosses, plus distance and relative position.
     - Met
     - Every operation and attribute of Table 53 is implemented. The
       DE-9IM and one-dimensional operations are methods of
       :class:`~rhealpixdggs.cell.Cell`:
       :meth:`~rhealpixdggs.cell.Cell.equals`,
       :meth:`~rhealpixdggs.cell.Cell.contains_cell`,
       :meth:`~rhealpixdggs.cell.Cell.within`,
       :meth:`~rhealpixdggs.cell.Cell.touches`,
       :meth:`~rhealpixdggs.cell.Cell.disjoint`,
       :meth:`~rhealpixdggs.cell.Cell.intersects`,
       :meth:`~rhealpixdggs.cell.Cell.overlaps` and
       :meth:`~rhealpixdggs.cell.Cell.crosses` (the last two are always
       False between two cells of one hierarchy, which nest, touch or
       are disjoint; :meth:`~rhealpixdggs.cell.Cell.region_overlaps`
       is the cell-against-region form), the DE-9IM pattern test
       :meth:`~rhealpixdggs.cell.Cell.relate`,
       :meth:`~rhealpixdggs.cell.Cell.distance` (the infimum distance
       between the cells: 0 when they intersect, otherwise the least
       distance between their boundaries, Euclidean in the plane or
       geodesic on the ellipsoid) with
       :meth:`~rhealpixdggs.cell.Cell.within_distance`, and
       :meth:`~rhealpixdggs.cell.Cell.relative_position` returning the
       Table 54 enumeration :class:`~rhealpixdggs.cell.RelativePosition`
       with :meth:`~rhealpixdggs.cell.Cell.relate_position`. The
       ``projectTo`` vector of the spec selects a dimension; for this
       two-dimensional grid ``relative_position`` takes its spatial part
       as a planar ``direction``. The hierarchy predicates are
       :meth:`~rhealpixdggs.cell.Cell.parent_of`,
       :meth:`~rhealpixdggs.cell.Cell.child_of` and
       :meth:`~rhealpixdggs.cell.Cell.sibling_of`, with the spec's
       ``inheritID`` flag. The operations that return zone sets are
       methods of :class:`~rhealpixdggs.zoneset.ZoneSet`, with
       zone-to-zone delegates on ``Cell``: ``union``, ``intersection``,
       ``difference`` and ``sym_difference`` by exact hierarchical
       algebra under the ``rangeRefine`` resolution filter, ``buffer``
       from ``distance``, and ``parent``, ``child`` and ``sibling`` with
       ``levels`` and ``inheritID`` per clause 8.3.3, whose worked
       examples are the tests. The query attributes ``boundary``,
       ``convexHull`` and ``boundaryType`` are ``ZoneSet.boundary``,
       ``ZoneSet.convex_hull`` and ``ZoneSet.boundary_type``, the last a
       value of the EA_BoundaryType code list (Table 76,
       :class:`~rhealpixdggs.zoneset.BoundaryType`); every rHEALPix edge
       is a ``projectedLine``. Beyond the mandate, the DE-9IM predicates
       are also defined between two zone sets read as regions. See
       :doc:`zoneset`.
   * - A.18 interoperation query
     - Operations to read, interpret and execute external data
       queries.
     - Partial
     - External geometries (Shapely/WKT) are accepted as queries by
       :func:`~rhealpixdggs.rhp_wrappers.polyfill` and
       :func:`~rhealpixdggs.rhp_wrappers.linetrace`; the documented
       mapping onto the model's query interface is to be written
       (v1.0). Transport protocols are out of library scope (see
       Scope notes).
   * - A.19 interoperation broadcast
     - Operations translating query results into standard data
       formats for delivery.
     - Met
     - :mod:`rhealpixdggs.export` translates cells and query results
       into the ``TranslationType`` codes ``toGeoJSON``
       (:func:`~rhealpixdggs.export.to_geojson`, an RFC 7946
       FeatureCollection: longitude-latitude order, counter-clockwise
       rings, antimeridian-straddling cells split into MultiPolygons,
       cap and pole-touching cells closed as section 3.1.9 recommends;
       every feature valid and the features of a resolution tiling the
       longitude-latitude rectangle, verified over whole resolutions in
       ``tests/test_export.py``) and ``toASCII``
       (:func:`~rhealpixdggs.export.to_csv`, delimited text with the
       same attributes), publishing the centroid of A.27 as each cell's
       position; :func:`~rhealpixdggs.export.geometries` gives the
       shapely geometries for OGR-backed writers. ``toGML`` follows in
       v1.0 with the application schema (#103); ``toHDF``,
       ``toJSON-LD``, ``toNetCDF`` and ``toXML`` are not planned. The
       ``broadcastResult`` operation of Table 60 (transport) is outside
       the library's scope, per the note above.

Equal-Area Earth Reference System
---------------------------------

.. list-table::
   :widths: 16 34 10 40
   :header-rows: 1

   * - Test
     - What it requires
     - Status
     - Evidence and plan
   * - A.20 harmonized model
     - The EAERS definition conforms to the Equal-Area data model.
     - Partial
     - Substance as per A.6; the formal declaration (including the
       grid-constraint list with cellEquiSized first) is to be
       written (v1.0).
   * - A.21 global domain
     - The domain is the whole Earth-model surface and the equal-size
       cell constraint is declared.
     - Partial
     - The domain is the whole ellipsoid surface (A.8); the formal
       cellEquiSized declaration rides with A.20 (v1.0).
   * - A.22 initial tessellation
     - An initial grid of equal-area cells produced by mapping the
       faces of a base polyhedron to the Earth model.
     - Met
     - The six cube faces map to the six resolution 0 cells through
       the rHEALPix projection; their areas are exactly equal by the
       projection's equal-area property.
   * - A.23 tessellation sequence
     - Tessellation operations generating progressively finer grids.
     - Met
     - Each cell subdivides into N_side² congruent planar children;
       :meth:`~rhealpixdggs.cell.Cell.subcells`.
   * - A.24 sequence maximum
     - A declared limit on refinement iterations tied to the cell-area
       error budget.
     - Met
     - ``RHEALPixDGGS.max_resolution`` is the limit (15 for
       ``WGS84_003``); identifiers beyond it are rejected. `Area error
       budget and sequence maximum`_ below derives the limit and shows
       every grid up to it is within the budget.
   * - A.25 global area preservation
     - Domain completeness and position uniqueness preserved at every
       level (dynamic-datum form).
     - Met
     - WGS 84 is a dynamic datum; completeness and uniqueness hold at
       every level by the subdivision construction (A.9, A.10 apply
       per level), and at every coordinate epoch, because the grid is
       fixed to the ellipsoid (`Coordinate reference system`_).
   * - A.26 cells are simple 2D polygons
     - Every EAERS cell is a simple polygon.
     - Met
     - As A.11; the boundary-curve types per cell shape are to be
       described alongside the GML work (v1.0).
   * - A.27 direct position is the centroid
     - Each cell's direct position is its centroid, the geodesic
       center of surface area.
     - Met
     - Requirement 27: "the DirectPosition of an EA_Cell [is] the
       centroid of the EA_cell, computed as the geodesic centre of
       surface area", and its test "requires verification that the
       attribute EA_Zone(centroid) returns a direct position on the
       surface of the cell". :meth:`~rhealpixdggs.cell.Cell.centroid`
       is the published position: the area-weighted mean of longitude
       and latitude over the cell, integrated in the equal-area plane
       (where equal planar area is equal surface area) with a fixed
       Gauss-Legendre rule per cell shape, exact to about 1e-14 degrees;
       :func:`~rhealpixdggs.rhp_wrappers.rhp_to_geo` and
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.centroids` report it, the
       exports of A.19 publish it, and ``tests/test_cell.py`` verifies
       it lies inside its cell for every shape. The nucleus (the centre
       of the planar square) remains the indexing point and the point on
       the isolatitude rings.
   * - A.28 area error budget
     - A declared cell-area error budget of at most 1% per grid.
     - Met
     - Declared per grid in `Area error budget and sequence maximum`_
       below, in the dynamic-datum form: at most 6 × 10⁻⁸ for
       ``WGS84_003``, computed by
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.area_error_budget`.
   * - A.29 equal-area cells
     - Cells in each grid are equal-area within the declared budget.
     - Met
     - Planar cell areas are exactly equal at each resolution, and
       the projection is equal-area onto the ellipsoid; verified by
       the area tests in the test suite (budgets published under
       A.28).

Zonal identifiers
-----------------

Requirement 13 asks that every zone in every grid have a globally unique
identifier, and Requirement 5 that identifiers follow the common data
model, in which a ``ZonalIdentifier`` is a label or code that names one
zone. This section specifies the identifiers; :doc:`identifiers`
describes working with them.

The text identifier of a cell, its index string, is a face letter
followed by one digit per refinement level. In ABNF (:rfc:`5234`):

.. code-block:: abnf

   index-string = face *digit     ; at most max_resolution digits
   face         = "N" / "O" / "P" / "Q" / "R" / "S"
   digit        = %x30-38         ; "0" to "8" for N_side 3
                                  ; %x30-33, "0" to "3", for N_side 2

The letter names a resolution 0 cell: ``N`` and ``S`` the polar faces,
``O`` to ``R`` the equatorial faces eastward from ``lon_0``. Each digit
picks one of the ``N_side²`` children of the cell named so far, numbered
from the top left along rows, so ``N12`` is child 2 of child 1 of ``N``.
The number of digits is the cell's resolution, from 0 to
``max_resolution``, 15 for ``WGS84_003`` (Requirement 24, below). Index
strings are defined for ``N_side`` 2 and 3, whose digits are single
characters. The sequence form below is the same grammar over integers,
a face letter followed by at most ``max_resolution`` integers from 0 to
``N_side² - 1``, and every grid has it; a grid with a larger ``N_side``
names its cells by the sequence and uint64 forms.

A string is an identifier of a grid exactly when it matches the grammar
for that grid's ``N_side`` with at most ``max_resolution`` digits: upper
case, with no whitespace, separator or prefix.
:meth:`~rhealpixdggs.dggs.RHEALPixDGGS.parse_index` returns None for any
other string, :func:`~rhealpixdggs.rhp_wrappers.rhp_is_valid` False, and
the array methods that take index strings reject it. The regular
expression the grammar gives for ``WGS84_003`` agrees with the parser::

    >>> import re
    >>> from rhealpixdggs.dggs import WGS84_003
    >>> index_string = re.compile(r"[NOPQRS][0-8]{0,15}")
    >>> for s in ["N", "N12", "S001450634", "n12", "N9", "N 1", "N" + "0" * 16]:
    ...     print(repr(s), bool(index_string.fullmatch(s)), WGS84_003.parse_index(s) is not None)
    'N' True True
    'N12' True True
    'S001450634' True True
    'n12' False False
    'N9' False False
    'N 1' False False
    'N0000000000000000' False False

The index string is the canonical text form, and each cell has exactly
one: the string records the path from the face to the cell, so distinct
cells have distinct strings, and a string's prefixes name the cell's
ancestors. The sequence form, the suid ``('N', 1, 2)``, is the same
identifier as a tuple of its symbols, and
:meth:`~rhealpixdggs.dggs.RHEALPixDGGS.format_index` and
:meth:`~rhealpixdggs.dggs.RHEALPixDGGS.parse_index` convert between the
two without loss. The uint64 form is the integer identifier the OGC
registration declares beside the text one, the cell's level, row and
column packed into one integer (:attr:`~rhealpixdggs.cell.Cell.uint64`);
it is unique because a position is unique within its level. The
level-row-col text of OGC tile matrix sets is a position, not an
identifier the reference system declares.

Within a grid, uniqueness is therefore structural. Across grids the same
string names different cells: ``N12`` lies at a different place in a grid
with another ``lon_0``, another ``N_side``, other polar squares or
another ellipsoid. An identifier is globally unique together with its
codespace, the URI of the DGGRS definition that fixes the grid, which
:attr:`~rhealpixdggs.dggs.RHEALPixDGGS.codespace` gives. For the
OGC-registered grid the URI is OGC's. For any other grid the authority is
whoever defines the grid, so the URI is given when the grid is
constructed, and should resolve to the grid's definition in the DGGRS
definition schema of OGC API - DGGS, which
:meth:`~rhealpixdggs.dggs.RHEALPixDGGS.to_dggrs_definition` writes for
any grid; for the registered grid the result agrees with the definition
OGC publishes. A grid given no URI has no codespace, and an export that
needs one, such as GML, requires it::

    >>> from rhealpixdggs.dggs import WGS84_003_OGC, WGS84_003, RHEALPixDGGS
    >>> WGS84_003_OGC.codespace
    'https://www.opengis.net/def/dggrs/OGC/1.0/rHEALPix'
    >>> WGS84_003.codespace is None
    True
    >>> RHEALPixDGGS(N_side=2, uri="https://example.org/dggrs/n2").codespace
    'https://example.org/dggrs/n2'

Data exchanged by identifier should carry the codespace, as the
``codeSpace`` attribute of a GML identifier does, and data meant for
other software should use the registered grid, whose URI other software
can resolve.

Coordinate reference system
---------------------------

Requirement 7 asks a DGGS reference system to define a CRS and to meet
the coordinate-epoch provisions of ``MD_ReferenceSystem``. The CRS of a
grid is the projected CRS of the rHEALPix projection with the grid's
central meridian and polar squares, in metres, on the grid's ellipsoid:
for the WGS84 grids the base CRS is WGS 84, EPSG:4326, and for any
other ellipsoid a geographic CRS of that ellipsoid.
:attr:`~rhealpixdggs.dggs.RHEALPixDGGS.crs` gives the CRS as a
``pyproj.CRS``, from which ``to_wkt()`` and ``to_json()`` write the WKT2
and PROJJSON definitions; the DGGRS definition of the previous section
carries the PROJJSON. For the registered grid::

    >>> from rhealpixdggs.dggs import WGS84_003_OGC
    >>> print(WGS84_003_OGC.crs.to_wkt(pretty=True))
    PROJCRS["rHEALPix (WGS 84, lon_0=50, north_square=0, south_square=0)",
        BASEGEOGCRS["WGS 84",
            ENSEMBLE["World Geodetic System 1984 ensemble",
                MEMBER["World Geodetic System 1984 (Transit)"],
                MEMBER["World Geodetic System 1984 (G730)"],
                MEMBER["World Geodetic System 1984 (G873)"],
                MEMBER["World Geodetic System 1984 (G1150)"],
                MEMBER["World Geodetic System 1984 (G1674)"],
                MEMBER["World Geodetic System 1984 (G1762)"],
                MEMBER["World Geodetic System 1984 (G2139)"],
                MEMBER["World Geodetic System 1984 (G2296)"],
                ELLIPSOID["WGS 84",6378137,298.257223563,
                    LENGTHUNIT["metre",1]],
                ENSEMBLEACCURACY[2.0]],
            PRIMEM["Greenwich",0,
                ANGLEUNIT["degree",0.0174532925199433]],
            ID["EPSG",4326]],
        CONVERSION["rHEALPix lon_0=50, north_square=0, south_square=0",
            METHOD["PROJ rhealpix"],
            PARAMETER["lon_0",50,
                ANGLEUNIT["degree",0.0174532925199433,
                    ID["EPSG",9122]]],
            PARAMETER["north_square",0,
                SCALEUNIT["unity",1,
                    ID["EPSG",9201]]],
            PARAMETER["south_square",0,
                SCALEUNIT["unity",1,
                    ID["EPSG",9201]]]],
        CS[Cartesian,2],
            AXIS["(E)",east,
                ORDER[1],
                LENGTHUNIT["metre",1,
                    ID["EPSG",9001]]],
            AXIS["(N)",north,
                ORDER[2],
                LENGTHUNIT["metre",1,
                    ID["EPSG",9001]]]]

The conversion is PROJ's ``rhealpix`` method, so the definition is one
other software can apply: every planar coordinate this package computes,
cell corners and nuclei included, is a coordinate in this CRS, and
PROJ's own ``rhealpix`` reproduces the package's projection to within a
few micrometres over the whole globe, which the test suite checks for
each kind of ellipsoid. The projection is defined on the authalic sphere
of the ellipsoid, so the plane's extent is ``2 π R_A`` by ``3 π R_A / 2``
in metres. The angular unit a grid is built with (degrees or radians)
is a convention of the package's interface, not of the CRS, whose
geographic coordinates are in degrees.

WGS 84 is a datum ensemble, EPSG:4326 with an ensemble accuracy of 2 m,
and its realisations are dynamic reference frames, so ISO 19111:2019
asks that coordinates referenced to it carry a coordinate epoch. The
grid is defined on the ellipsoid and has no epoch: a zone is a fixed
region of the ellipsoid, the same at every epoch, so the sequence of
grids preserves domain completeness and position uniqueness at any epoch
(A.25) and the error budget above holds at any epoch (A.28). What an
epoch changes is which zone a point on the ground falls in, because the
ground moves relative to the frame with its tectonic plate, by up to
about 7 cm a year. A decade of that motion, 0.7 m, is the width of a
resolution 15 cell; a century, 7 m, the width of a resolution 13 cell.
The package therefore takes coordinates as given, in whatever
realisation and epoch they have, and makes no epoch transformation.
Data assigned to zones at a resolution where such displacements are a
fraction of a cell should be transformed to a common frame and epoch
first, which PROJ does between ITRF realisations with a plate-motion
model, and the frame and epoch recorded with the identifiers and their
codespace.

Area error budget and sequence maximum
--------------------------------------

Requirement 28 asks an equal-area reference system to declare, for each
of its grids, an error budget of at most 1 %: the largest ratio by which
a cell's area may differ from the theoretical cell area, the surface
area of the Earth model divided by the number of cells. Requirement 24
asks for a limit on the sequence of grids that keeps every grid within
the budget. This section is that declaration, with ``WGS84_003`` as the
worked case. WGS84 is a dynamic datum, so the dynamic-datum form of
Requirement 28 applies: one budget value per grid.

The rHEALPix projection is equal-area from the ellipsoid to the plane,
where the six resolution 0 cells are squares of side ``π R_A / 2`` and
each cell has ``N_side²`` congruent children. The resolution ``r`` grid
has ``6 N_side^(2r)`` cells, so the theoretical cell area is

.. math::

   A_r = \frac{4 \pi R_A^2}{6\, N_{side}^{2r}},

and because the projection is equal-area, every cell's area is also
``A_r``: the ratio of Requirement 28 is 1 at every resolution in exact
arithmetic. :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell_area` returns
``A_r`` and :meth:`~rhealpixdggs.cell.Cell.area` reports it for every
cell.

Of the error sources Topic 21 lists, two apply to a closed-form system in
double precision: data-type precision and error propagation. The datum,
its realisation and the Earth model are the CRS's (A.7), and the area
another program measures from a cell's boundary points carries that
program's error. Two floating-point effects remain. The computed value
of ``A_r`` passes through at most 10 rounding steps, so carries at most
``10 ε = 2.2 × 10⁻¹⁵`` relative error at every resolution;
:meth:`~rhealpixdggs.dggs.RHEALPixDGGS.area_error_budget` derives the
count and gives the resulting absolute tolerance per grid, the
tolerance for testing two computed areas for equality. And the corners
of a cell, as the plane holds them, lie within a few units in the last
place of the plane's extent ``π R_A``, under ``δ = 10⁻⁸ m``, of their
exact positions, so the area enclosed by a cell's corners differs from
``A_r`` by at most ``4 δ / w_r``, where ``w_r`` is the cell width. The
declared budget for the resolution ``r`` grid is the sum,
``10 ε + 4 δ / w_r``. For ``WGS84_003``::

    >>> import sys
    >>> from rhealpixdggs.dggs import WGS84_003
    >>> budget = WGS84_003.area_error_budget()
    >>> for r, b in budget.items():
    ...     declared = 10 * sys.float_info.epsilon + 4e-8 / WGS84_003.cell_width(r)
    ...     print(f"{r:2d}  {b['cell_area_m2']:10.4e}  {b['abs_tolerance']:8.2e}  {declared:8.2e}")
     0  8.5011e+13  1.89e-01  6.22e-15
     1  9.4457e+12  2.10e-02  1.42e-14
     2  1.0495e+12  2.33e-03  3.82e-14
     3  1.1661e+11  2.59e-04  1.10e-13
     4  1.2957e+10  2.88e-05  3.26e-13
     5  1.4397e+09  3.20e-06  9.73e-13
     6  1.5996e+08  3.55e-07  2.92e-12
     7  1.7774e+07  3.95e-08  8.74e-12
     8  1.9749e+06  4.39e-09  2.62e-11
     9  2.1943e+05  4.87e-10  7.87e-11
    10  2.4381e+04  5.41e-11  2.36e-10
    11  2.7090e+03  6.02e-12  7.08e-10
    12  3.0100e+02  6.68e-13  2.12e-09
    13  3.3444e+01  7.43e-14  6.37e-09
    14  3.7160e+00  8.25e-15  1.91e-08
    15  4.1289e-01  9.17e-16  5.74e-08

The columns are the resolution, the cell area ``A_r`` in m², the
absolute tolerance on its computed value in m², and the declared
budget as a ratio. The budget grows by ``N_side`` per resolution,
because the corner displacement is fixed in metres while the cell
shrinks, and reaches 1 % at a cell width of about 4 µm: resolution 26
for ``N_side = 3``. This is the limit Requirement 24 refers to, and the
sequence stops well short of it. ``max_resolution`` is the first
resolution whose cells are no larger than ``max_areal_resolution``,
one square metre by default:

.. math::

   \text{max\_resolution} = \left\lceil
   \frac{\log(A_0 / \text{max\_areal\_resolution})}{2 \log N_{side}}
   \right\rceil,

15 for ``N_side = 3`` (cells of 0.41 m²) and 24 for ``N_side = 2``
(0.30 m²). Identifiers longer than ``max_resolution + 1`` digits are
rejected. At resolution 15 the budget is 6 × 10⁻⁸, five orders of
magnitude inside the 1 % that Requirement 28 allows, and the binding
limit on the sequence is the area floor. A ``max_areal_resolution``
below about 10⁻¹¹ m² would carry the sequence past the budget.

Planned work
------------

The rows above that are not yet Met are v1.0 work: GML export with its
application schema
(`#103 <https://github.com/manaakiwhenua/rhealpixdggs-py/issues/103>`_),
and the formal reference-system declaration with the quantization and
interoperation mapping documentation
(`#104 <https://github.com/manaakiwhenua/rhealpixdggs-py/issues/104>`_).

Progress is tracked in the project's `GitHub milestones
<https://github.com/manaakiwhenua/rhealpixdggs-py/milestones>`_.
