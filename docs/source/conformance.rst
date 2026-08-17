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
     - Cell suids, :func:`~rhealpixdggs.rhp_wrappers.rhp_to_parent`,
       :meth:`~rhealpixdggs.cell.Cell.subcells`,
       :meth:`~rhealpixdggs.cell.Cell.nucleus`,
       :meth:`~rhealpixdggs.cell.Cell.centroid`; a formal identifier
       grammar (valid strings, canonical form) is to be published
       (v1.0).

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
     - Partial
     - :class:`~rhealpixdggs.ellipsoids.Ellipsoid` plus the rHEALPix
       projection define the CRS operationally; WKT2/PROJJSON
       definitions and a dynamic-datum epoch statement are to be
       published (v1.0).
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
     - :meth:`~rhealpixdggs.cell.Cell.nucleus` is interior by
       construction (see also A.27 for the equal-area refinement).
   * - A.13 address
     - Every zone has a globally unique identifier, structured on a
       recognized indexing method.
     - Met
     - Suids are hierarchy-based (resolution 0 letter plus one base-
       N_side² digit per refinement); uniqueness is structural. A
       formal grammar is to be published with A.5 (v1.0).
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
     - Partial
     - :meth:`~rhealpixdggs.cell.Cell.touches`,
       :meth:`~rhealpixdggs.cell.Cell.within`,
       :meth:`~rhealpixdggs.cell.Cell.contains_cell` and
       :meth:`~rhealpixdggs.cell.Cell.disjoint` are implemented.
       Remaining (v0.9.0): equals, intersects, crosses, distance and
       relativePosition, and realigning
       :meth:`~rhealpixdggs.cell.Cell.overlaps` with DE-9IM semantics
       (a breaking change: its current meaning is
       containment-in-either-direction).
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
     - Planned
     - GeoJSON FeatureCollection export of zones and query results,
       plus delimited-text export (v0.9.0); GML export with a
       published application schema (v1.0).

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
     - Partial
     - ``RHEALPixDGGS.max_resolution`` (default 15) and
       ``max_areal_resolution`` exist; the derivation tying the limit
       to the area error budget is to be documented (v1.0).
   * - A.25 global area preservation
     - Domain completeness and position uniqueness preserved at every
       level (dynamic-datum form).
     - Met
     - WGS84 is a dynamic datum; completeness and uniqueness hold at
       every level by the subdivision construction (A.9, A.10 apply
       per level).
   * - A.26 cells are simple 2D polygons
     - Every EAERS cell is a simple polygon.
     - Met
     - As A.11; the boundary-curve types per cell shape are to be
       described alongside the GML work (v1.0).
   * - A.27 direct position is the centroid
     - Each cell's direct position is its centroid, the geodesic
       center of surface area.
     - Partial
     - :meth:`~rhealpixdggs.cell.Cell.centroid` computes area-true
       centroids for all four cell shapes (integration is performed
       in the equal-area plane); adopting it as the published
       representative position, alongside the nucleus used for
       indexing, is scheduled (v0.9.0).
   * - A.28 area error budget
     - A declared cell-area error budget of at most 1% per grid.
     - Partial
     - :meth:`~rhealpixdggs.dggs.RHEALPixDGGS.area_error_budget`
       computes the budget; the per-resolution table is to be
       published on this page (v1.0).
   * - A.29 equal-area cells
     - Cells in each grid are equal-area within the declared budget.
     - Met
     - Planar cell areas are exactly equal at each resolution, and
       the projection is equal-area onto the ellipsoid; verified by
       the area tests in the test suite (budgets published under
       A.28).

Planned work by release
-----------------------

- **v0.9.0** — complete the zone-query interface (equals, intersects,
  crosses, distance, relativePosition; DE-9IM semantics for
  ``overlaps``), adopt the centroid as the published representative
  position, and add GeoJSON and delimited-text export of zones and
  query results.
- **v1.0** — finish this page's remaining rows: the formal reference
  system declaration, CRS definitions (WKT2/PROJJSON) with epoch
  statement, the zonal identifier grammar, the published error-budget
  and maximum-resolution derivation, the quantization and
  interoperation mapping documentation, and GML export with its
  application schema.

Progress is tracked in the project's GitHub milestones and issues.
