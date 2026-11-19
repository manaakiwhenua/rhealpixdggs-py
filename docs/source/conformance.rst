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

- The purely temporal tests A.1 and A.4 are not part of the
  Equal-Area Earth DGGS class; they are listed below as not claimed.
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
- **Not claimed** — outside the conformance class the project targets.

Common Spatio-temporal Classes
------------------------------

.. list-table::
   :widths: 16 34 10 40
   :header-rows: 1

   * - Test
     - What it requires
     - Status
     - Evidence and plan
   * - A.1 temporal geometry
     - Temporal geometry and topology (instants, intervals and their
       ordering) conforming to the common data model.
     - Not claimed
     - rHEALPixDGGS currently has no temporal dimension. The test
       belongs to the spatio-temporal conformance classes only.
   * - A.2 zone geometry
     - A zone has a geometry of declared spatial, temporal and
       topological dimension, and a topology of the same dimension.
     - Met
     - A zone is a :class:`~rhealpixdggs.cell.Cell`, a simple spatial
       geometry with the topology of a two-dimensional face. `Zone
       geometry and topology`_ below maps it onto the model's classes
       and declares the dimensions.
   * - A.3 spatial location
     - A location is a place named by an identifier (ISO 19112 style),
       with an extent, a point inside it, and links to the locations
       above and below it.
     - Met
     - Every :class:`~rhealpixdggs.cell.Cell` is such a location.
       `Spatial location`_ below maps it onto the model's class.
   * - A.4 temporal RS using period identifiers
     - A temporal reference system whose periods carry identifiers,
       conforming to the common data model.
     - Not claimed
     - As A.1: currently no temporal dimension, and outside the
       targeted class.
   * - A.5 RS using zonal identifiers
     - A reference system made of zones and zone classes: each zone
       has an identifier, an extent, a representative position and
       parent/child links, and each zone class a refinement level.
     - Met
     - The reference system is an
       :class:`~rhealpixdggs.dggs.RHEALPixDGGS`, a zone a
       :class:`~rhealpixdggs.cell.Cell` and a zone class a resolution.
       `Reference system using zonal identifiers`_ below maps them
       onto the model's classes.

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
     - Met
     - All defining parameters live on
       :class:`~rhealpixdggs.dggs.RHEALPixDGGS` (ellipsoid, N_side,
       north_square, south_square, max_resolution). `Reference system
       declaration`_ below states them as a DGG_ReferenceSystem and
       maps the model's classes.
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
       revolution: spatial dimension 3, topological dimension 2,
       temporal dimension 0 (`Zone geometry and topology`_).
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
     - Met
     - `Quantization`_ below declares the three roles the library
       supports, ``asCoordinates``, ``asTags`` and ``asDataCells``,
       names the operations for each and maps the classes of the
       model. The library assigns identifiers and holds no data: the
       observed values stay with the caller.
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
       shapely geometries for OGR-backed writers; and ``toGML``
       (:func:`~rhealpixdggs.export.to_gml`, GML 3.2 under the
       application schema published at ``schemas/zone.xsd``, with the
       identifier in the grid's codespace, the boundary curve types of
       A.26 and the geometry in EPSG:4326 or, exactly, in the grid's
       own CRS carried in the document). ``toHDF``, ``toJSON-LD``,
       ``toNetCDF`` and ``toXML`` are not planned. The
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
     - Met
     - `Equal-area reference system declaration`_ below extends the
       declaration of A.6 and maps the model's classes, with
       ``cellEquiSized`` the first grid constraint.
   * - A.21 global domain
     - The domain is the whole Earth-model surface and the equal-size
       cell constraint is declared.
     - Met
     - The domain is the whole ellipsoid surface (A.8), and
       ``cellEquiSized`` is declared as the first grid constraint
       (`Reference system declaration`_).
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
     - As A.11. The curve type of every edge, from the
       ``EA_BoundaryType`` code list, is
       :func:`~rhealpixdggs.export.boundary_types` and is written with
       each zone in the GML export: quad cells are bounded by parallels
       (``smallCircle``) and meridians (``greatEllipse``); a cap by one
       parallel; dart and skew quad cells by parallels, by meridians
       through the pole for even ``N_side``, and by images of straight
       planar lines (``projectedLine``), as :doc:`export` describes.
       The test suite checks the classification against the sampled
       edges over whole resolutions. Parallels and meridians are exact
       as segments in EPSG:4326; a projected line is exact only as a
       segment in the grid's CRS, which the planar GML form writes.
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

Zone geometry and topology
--------------------------

Requirement 2 asks that zones conform to the classes of clause 7.2.2
(Tables 17–22). Those classes say what kind of geometric object a
zone is and what dimensions it has.

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - Model class
     - In rHEALPixDGGS
   * - ZoneGeometry
     - :class:`~rhealpixdggs.cell.Cell`. Every cell has
       ``spatialDimension`` 3, ``topologicalDimension`` 2 and
       ``temporalDimension`` 0: it is a region of the surface of an
       ellipsoid embedded in three-dimensional space, as Table 19
       notes and Table 64 requires of an Equal-Area Earth DGGS.
   * - ZoneSimpleGeometry
     - Every cell: a spatial geometry invariant over time. Its shape
       is given by :meth:`~rhealpixdggs.cell.Cell.vertices` and
       :meth:`~rhealpixdggs.cell.Cell.boundary`, on the ellipsoid or
       in the plane of the projection.
   * - ZoneSimpleTopology, ZoneTopology
     - The topology of a cell is that of a two-dimensional face, with
       an interior, a boundary and an exterior. The zone queries of
       A.17 are defined on these, and
       :meth:`~rhealpixdggs.cell.Cell.relate` exposes them directly as
       a DE-9IM matrix.
   * - ZoneCompoundGeometry, ZoneCompoundTopology
     - Not instantiated. These pair a spatial zone with a temporal
       one, and rHEALPixDGGS currently has no temporal dimension.

Spatial location
----------------

Requirement 3 asks that spatial locations conform to the LocationS
class of clause 7.3.1 (Table 23), a specialisation of the ISO 19112
location: a place referred to by an identifier instead of by
coordinates. Every :class:`~rhealpixdggs.cell.Cell` is a LocationS.
As the table notes, these locations are defined by the reference
system itself and not by an independent authority.

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - LocationS element
     - In rHEALPixDGGS
   * - identifier
     - The cell's index string, ``str(cell)``, specified in `Zonal
       identifiers`_.
   * - extent
     - :meth:`~rhealpixdggs.cell.Cell.boundary`.
   * - coordinateTuple
     - :meth:`~rhealpixdggs.cell.Cell.centroid`, the direct position
       of A.12 and A.27, which lies inside the extent.
   * - parentOf, childOf
     - :meth:`~rhealpixdggs.cell.Cell.subcells` gives the locations a
       cell is parent of, and :meth:`~rhealpixdggs.cell.Cell.parent`
       the location it is child of. Resolution 0 cells have no parent
       and cells at ``max_resolution`` have no children (asking for
       them raises ``ValueError``), which is why the model makes both
       links conditional.
   * - Constraint: dimension > 0
     - A cell has topological dimension 2 (`Zone geometry and
       topology`_).

Reference system using zonal identifiers
----------------------------------------

Requirement 5 asks that the reference system conform to the classes
of clause 7.3.3 (Tables 28–39). They describe four things: the
reference system, the identifier, the zone and the zone class. Each
comes as a primitive, which is spatial or temporal, and as a compound
or complex of one spatial and one temporal primitive. rHEALPixDGGS
instantiates the spatial primitives.

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - Model class
     - In rHEALPixDGGS
   * - RSUsingZonalIdentifiers, RSUsingZonalIdentifiersPrimitive
     - An :class:`~rhealpixdggs.dggs.RHEALPixDGGS` instance, such as
       ``WGS84_003``. It comprises its zone classes, the resolutions 0
       to ``max_resolution``, and its zones, every cell of every
       resolution.
   * - ZonalIdentifier, ZonalIdentifierPrimitive
     - The index string of a cell, a geographic identifier specified
       in `Zonal identifiers`_. The 64-bit integer
       :attr:`~rhealpixdggs.cell.Cell.uint64` is a second encoding of
       the same identifier.
   * - Zone, ZonePrimitive
     - :class:`~rhealpixdggs.cell.Cell`, which is the LocationS of
       `Spatial location`_. Its ``identifier`` is ``str(cell)``, its
       ``extent`` :meth:`~rhealpixdggs.cell.Cell.boundary`, its
       ``representativePosition``
       :meth:`~rhealpixdggs.cell.Cell.centroid`, and its ``parentOf``
       and ``childOf`` links :meth:`~rhealpixdggs.cell.Cell.subcells`
       and :meth:`~rhealpixdggs.cell.Cell.parent`.
   * - ZoneClass, ZoneClassPrimitive
     - A resolution. Its ``refinementLevel`` is the integer itself,
       which every cell carries as ``Cell.resolution``; its zones are
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.grid`; the classes it is
       child of and parent of are the resolutions one less and one
       more. Resolution 0 has no class above it and
       ``max_resolution`` none below.
   * - RSUsingZonalIdentifiersCompound, ZonalIdentifierComplex,
       ZoneCompound, ZoneClassComplex
     - Not instantiated. Each pairs a spatial primitive with a
       temporal one, and rHEALPixDGGS currently has no temporal
       dimension.

Reference system declaration
----------------------------

Requirement 6 asks that the reference system be defined in the terms
of clause 8.2 (Tables 40–46). An
:class:`~rhealpixdggs.dggs.RHEALPixDGGS` is a DGG_ReferenceSystem
(Table 41). The library defines a family of them, one for each choice
of ellipsoid, ``N_side``, ``north_square`` and ``south_square``; the
declaration below holds for every member, with the values for
``WGS84_003`` given as the example.

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - Attribute
     - Value
   * - DGGS_Grids
     - ``cellEquiSized``, then ``cellAxisAligned``, in that order of
       priority (Table 45). Cells of one resolution have equal area
       (A.29), and their edges are parallel to the axes of the grid's
       projected CRS.
   * - DGGS_Refinement
     - ``nestedChildCell`` (Table 46): the children of a cell tile it
       exactly, so their combined boundary is the parent's. When
       ``N_side`` is odd the centre child also shares its parent's
       nucleus. That is not declared as ``centredChildCell``, which is
       defined on the representative position, the centroid (A.27),
       and the two centroids differ.
   * - DGGS_RefinementRatio
     - ``N_side`` squared, the same at every level: 9 for
       ``WGS84_003``.
   * - DGGS_RefSys: CRS
     - :attr:`~rhealpixdggs.dggs.RHEALPixDGGS.crs`, with the
       coordinate-epoch statement of `Coordinate reference system`_.
   * - DGGS_RefSys: GLOBE
     - The surface of the grid's ellipsoid of revolution,
       ``RHEALPixDGGS.ellipsoid``: spatial dimension 3, topological
       dimension 2, temporal dimension 0 (`Zone geometry and
       topology`_).
   * - DGGS_RefSys: MDRS
     - The description of the whole system that
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.to_dggrs_definition`
       returns.
   * - DGGS_RefSys: ZIRS
     - :attr:`~rhealpixdggs.dggs.RHEALPixDGGS.uri`, the identifier of
       the grid's zonal identifier reference system and the codespace
       of its identifiers (`Zonal identifiers`_). A grid constructed
       without a URI, ``WGS84_003`` among them, has none.

The remaining classes of the model are these.

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - Model class
     - In rHEALPixDGGS
   * - Cell
     - :class:`~rhealpixdggs.cell.Cell`, the unit of geometry that
       goes with a zone. It has the dimensions of the globe.
   * - DiscreteGlobalGrid
     - The cells of one resolution,
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.grid`, paired with the
       zone class of that resolution (`Reference system using zonal
       identifiers`_). Its ``refinementRatio`` is ``N_side`` squared.
   * - GlobeGeometry
     - The ellipsoid surface. It comprises every cell, and each
       resolution's grid is a tessellation of it (A.9, A.10, A.25).
   * - DGG_ReferenceSystemType, DGG_GridConstraint,
       DGG_RefinementStrategy
     - The union and code lists whose values are declared above.

Quantization
------------

Requirement 16 asks for operations that assign data from external
sources to cells, in the terms of clause 8.3.2 (Tables 47–51).
rHEALPixDGGS is not a data service. It stores no data and returns
none; it is a component for someone building a system that does. Its
part of quantization is to turn the geometry of external data into
zonal identifiers, and the values assigned to those identifiers stay
with the caller.

Clause 8.3.1 lets a specification describe one or more of the six
roles a cell can play. The library supports three.

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Role
     - Operations
   * - ``asCoordinates``
     - A coordinate tuple becomes the identifier of the cell that
       contains it, at a resolution the caller chooses to reflect the
       precision of the position:
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell_from_point`, and
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cells_from_points` for
       arrays. To keep the structure of a line,
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cells_from_line` returns
       the ordered cells a segment passes through, and
       :func:`~rhealpixdggs.conversion.cells_from_linestring` those
       of a Shapely line. The wrappers
       :func:`~rhealpixdggs.rhp_wrappers.geo_to_rhp` and
       :func:`~rhealpixdggs.rhp_wrappers.linetrace` are the same
       operations for a latitude-longitude pair and a line.
   * - ``asTags``
     - An object is tagged with the smallest cell that wholly
       encloses it, whose resolution indicates the object's extent:
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cell_from_region` for a
       rectangle or a longitude-latitude quadrangle, and
       :func:`~rhealpixdggs.conversion.get_finest_containing_cell`
       for a polygon.
   * - ``asDataCells``
     - The cells that take a value from a feature are chosen by their
       geometry.
       :func:`~rhealpixdggs.conversion.cells_from_polygon`
       returns the cells of a Shapely polygon under one of three
       rules: the cell's centroid is inside it, the whole cell is
       inside it, or the cell meets it. It can return them compacted,
       and :func:`~rhealpixdggs.rhp_wrappers.polyfill` is its wrapper.
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.cells_from_region`
       covers a rectangle or quadrangle and
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.minimal_cover` a set of
       points. :func:`~rhealpixdggs.conversion.compact_cells` merges
       complete groups of siblings into their parents.
   * - ``asDataTiles``, ``asGraphicCells``, ``asGraphicTiles``
     - Not provided. Clipping features to cells and rendering them
       are left to software built on the library.

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - Model class
     - In rHEALPixDGGS
   * - QuantizationFunctions
     - The operations above. Its ``cellContext`` is the list
       ``asCoordinates``, ``asTags``, ``asDataCells``.
   * - DataAssignmentProcess
     - One call to one of those operations on an external geometry.
       Its ``cellContext`` is the role of the operation called.
   * - QuantizationContext
     - The ``role`` is fixed by the operation, as the table above
       lists.
   * - DGG_Observation
     - The ``identifier`` values are the index strings an operation
       returns. The ``field`` values are the caller's: the record is
       a mapping from index string to value that the caller builds
       and keeps. :func:`~rhealpixdggs.harmonics.to_ring_order` and
       :func:`~rhealpixdggs.harmonics.analysis` take such a mapping.
   * - Quantization
     - The code list of roles, of which three are supported.

Equal-area reference system declaration
---------------------------------------

Requirement 20 asks that an equal-area Earth reference system be
defined in the terms of clause 9.1 (Tables 63–76). These extend the
classes of `Reference system declaration`_, and an
:class:`~rhealpixdggs.dggs.RHEALPixDGGS` is an EA_DGG_ReferenceSystem
(Table 63) with the declaration given there: ``cellEquiSized`` is the
first grid constraint, as the table requires. Its CRS is a projected
CRS on a geodetic base, which serves both as the base CRS of Table 65
and as the projected CRS that Table 67 associates with the polyhedron.

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - Model class
     - In rHEALPixDGGS
   * - EA_GlobeGeometry
     - The surface of the grid's ellipsoid, with ``spatialDimension``
       3 and ``topologicalDimension`` 2 as Table 64 requires.
   * - EA_SurfaceInterpolation
     - ``elliptical``, or ``spherical`` for a grid on a sphere.
   * - EA_BaseUnitPolyhedron
     - A cube, with ``numFaces`` 6 and ``numEdges`` 12. The plane
       image of the rHEALPix projection is a net of the cube: four
       squares in a row for the equatorial faces, with the polar faces
       attached above and below the squares numbered ``north_square``
       and ``south_square``. The central meridian of the ellipsoid
       fixes its orientation.
   * - EA_PolyhedralTessellation
     - The rHEALPix projection,
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.rhealpix`. It is
       equal-area, so the six faces map to six cells of equal area
       that cover the ellipsoid once (A.22, A.9, A.10).
   * - EA_InitialDiscreteGlobalGrid
     - The six cells N, O, P, Q, R and S,
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.grid` at resolution 0,
       with ``refinementLevel`` 0.
   * - EA_DiscreteGlobalGridTessellation
     - The division of every cell into ``N_side`` squared equal planar
       squares, :meth:`~rhealpixdggs.cell.Cell.subcells` (A.23).
   * - EA_RefinedDiscreteGlobalGrid
     - ``grid`` at each resolution from 1 to ``max_resolution``, its
       ``refinementLevel`` one more than its parent's. The
       ``domainOfValidity`` is global: the grid is fixed to the
       ellipsoid (`Coordinate reference system`_).
   * - EA_DiscreteGlobalGrid
     - Either of the two above. Its ``dggsAxis`` is the curve that
       visits the cells of a resolution in the order of their
       identifiers, which :meth:`~rhealpixdggs.cell.Cell.successor`
       and :meth:`~rhealpixdggs.cell.Cell.predecessor` step along. Its
       ``dggsSurfaceInterpolation`` is that of the globe.
   * - EA_Zone
     - :class:`~rhealpixdggs.cell.Cell` as a location (`Spatial
       location`_), with its ``identifier`` and its ``centroid``,
       :meth:`~rhealpixdggs.cell.Cell.centroid` (A.27).
   * - EA_BoundaryType, EA_BoundaryData
     - The curve type of each edge,
       :func:`~rhealpixdggs.export.boundary_types` (A.26), and the
       curves themselves as the GML export writes them
       (:doc:`export`).

An EA_Cell is a :class:`~rhealpixdggs.cell.Cell` on the ellipsoid. Its
attributes (Table 74) are these.

.. list-table::
   :widths: 28 72
   :header-rows: 1

   * - EA_Cell attribute
     - In rHEALPixDGGS
   * - area
     - :meth:`~rhealpixdggs.cell.Cell.area` with ``plane=False``.
   * - perimeter
     - :meth:`~rhealpixdggs.cell.Cell.perimeter` with ``plane=False``.
   * - numEdges
     - :meth:`~rhealpixdggs.cell.Cell.num_edges` with ``plane=False``:
       4 for a quad or skew quad, 3 for a dart and 1 for a cap.
   * - nodes
     - :meth:`~rhealpixdggs.cell.Cell.vertices` with ``plane=False``,
       clockwise from the north-west. A dart has three, given by
       ``trim_dart=True``; the four points returned for a cap divide
       its one edge.
   * - boundaryTypes, boundaryData
     - As EA_BoundaryType and EA_BoundaryData above, one entry per
       edge.
   * - errorBudget
     - One value for each resolution,
       :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.area_error_budget`
       (A.28).
   * - maxLevel
     - ``max_resolution``, one value for the grid (A.24).

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
(A.25) and the error budget below is the budget of the grid on the
ellipsoid, which holds at any epoch (A.28). Topic 21 asks that a budget
on WGS 84 cater for the approximation in the WGS 84 definition. The
realisations of the ensemble differ by near-rigid transformations, so
the approximation is a positional uncertainty, 2 m, and not an area
one: what it affects is which zone a point on the ground falls in, and
it exceeds a cell width from resolution 14 (2.1 m). The coordinate
epoch affects the same thing, because the ground moves relative to the
frame with its tectonic plate, by up to about 7 cm a year: a decade of
that motion, 0.7 m, is the width of a resolution 15 cell, and a century,
7 m, the width of a resolution 13 cell. The package therefore takes
coordinates as given, in whatever realisation and epoch they have, and
makes no transformation between them. Data assigned to zones at a
resolution where such displacements are a fraction of a cell should be
referenced to a named realisation, not the ensemble, and transformed to
a common frame and epoch first, which PROJ does between ITRF
realisations with a plate-motion model; the frame and epoch belong with
the identifiers and their codespace.
:meth:`~rhealpixdggs.dggs.RHEALPixDGGS.crs_on` gives the grid's CRS on
a named realisation, ``crs_on("EPSG:9755")`` for WGS 84 (G2139), for
data that is referenced that way.

As an example of the epoch's effect, take a point in GDA2020
(EPSG:7844), the static Australian datum, fixed to its plate at epoch
2020.0, and observe it in ITRF2014 (EPSG:9000) in 2010 and again in
2030. The plate carries the point 1.16 m between the observations, and
the two fall in different resolution 15 cells; transformed back to the
static frame, which removes the plate's motion, both fall in the one
cell the ground point occupies::

    >>> import pyproj
    >>> from rhealpixdggs.dggs import WGS84_003_OGC
    >>> to_itrf = pyproj.Transformer.from_crs("EPSG:7844", "EPSG:9000", always_xy=True)
    >>> to_gda = pyproj.Transformer.from_crs("EPSG:9000", "EPSG:7844", always_xy=True)
    >>> lon, lat = 149.0999973, -35.3000031
    >>> observed = {epoch: to_itrf.transform(lon, lat, tt=epoch)[:2] for epoch in (2010.0, 2030.0)}
    >>> [str(WGS84_003_OGC.cell_from_point(15, p, plane=False)) for p in observed.values()]
    ['R665234534347420', 'R665234534347183']
    >>> round(pyproj.Geod(ellps="WGS84").inv(*observed[2010.0], *observed[2030.0])[2], 2)
    1.16
    >>> fixed = [to_gda.transform(*p, tt=epoch)[:2] for epoch, p in observed.items()]
    >>> [str(WGS84_003_OGC.cell_from_point(15, p, plane=False)) for p in fixed]
    ['R665234534347186', 'R665234534347186']

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

The rows above marked Partial or Planned are v1.0 work: the formal
reference-system declaration with the quantization and
interoperation mapping documentation
(`#104 <https://github.com/manaakiwhenua/rhealpixdggs-py/issues/104>`_).

Progress is tracked in the project's `GitHub milestones
<https://github.com/manaakiwhenua/rhealpixdggs-py/milestones>`_.
