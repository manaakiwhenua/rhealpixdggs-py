Export
======

Cells, and the sets of cells that queries such as ``polyfill`` return, can
be handed to other software as GeoJSON or as delimited text. This is the
"broadcast" half of interoperation in OGC Topic 21 v2.0 (requirement 19),
whose ``TranslationType`` code list names the formats: ``toGeoJSON`` is
:func:`~rhealpixdggs.export.to_geojson` and ``toASCII`` is
:func:`~rhealpixdggs.export.to_csv`. The library stops at the encoding, a
dictionary and a string; writing files and serving them over a protocol
are the caller's business, as the :doc:`conformance` page's scope note
says.

From a query to a file
----------------------

``polyfill_array`` returns the cells of a region as an array of index
strings, which is what the export functions take (they also accept plain
strings and ``Cell`` objects)::

    >>> import json
    >>> from rhealpixdggs import export
    >>> from rhealpixdggs.dggs import WGS84_003
    >>> cells = [str(c) for c in WGS84_003.cell(['P', 4]).subcells()]
    >>> fc = export.to_geojson(cells)
    >>> fc['type'], len(fc['features'])
    ('FeatureCollection', 9)
    >>> feature = fc['features'][0]
    >>> feature['id'], feature['geometry']['type']
    ('P40', 'Polygon')
    >>> feature['properties']['shape'], feature['properties']['resolution']
    ('quad', 2)
    >>> text = json.dumps(fc)  # then open(path, 'w').write(text)

The same cells as delimited text, one line per cell under a header::

    >>> print(export.to_csv(cells[:2]), end='')
    index,resolution,shape,area,centroid_lon,centroid_lat,nucleus_lon,nucleus_lat
    P40,2,quad,1049517740173.0219,-55.0,8.565250163965553,-55.0,8.557311559978533
    P41,2,quad,1049517740173.0219,-45.00000000000001,8.565250163965553,-45.00000000000001,8.557311559978533

Both carry the same attributes, as flat scalars so that GeoJSON readers
such as QGIS, OGR and geopandas show them as columns:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Property
     - Meaning
   * - ``index``
     - The cell's index string, also the feature's ``id``.
   * - ``resolution``
     - The cell's resolution.
   * - ``shape``
     - :attr:`~rhealpixdggs.cell.Cell.ellipsoidal_shape`: ``quad``,
       ``cap``, ``dart`` or ``skew_quad``.
   * - ``area``
     - The cell's surface area on the ellipsoid, in the ellipsoid's
       length units squared (square metres on the WGS84 grids); the same
       for every cell of a resolution.
   * - ``centroid_lon``, ``centroid_lat``
     - The cell's direct position, :meth:`~rhealpixdggs.cell.Cell.centroid`
       (see :doc:`conformance`, A.27), in degrees.
   * - ``nucleus_lon``, ``nucleus_lat``
     - The cell's indexing point, :meth:`~rhealpixdggs.cell.Cell.nucleus`,
       in degrees.

:func:`~rhealpixdggs.export.cell_table` returns these attributes as arrays,
and :func:`~rhealpixdggs.export.geometries` the geometries as shapely
objects, which is the route to geopandas and to every format OGR writes::

    >>> geoms = export.geometries(cells)
    >>> geoms[0].geom_type, round(geoms[0].area, 6)  # square degrees
    ('Polygon', 86.282871)

How the GeoJSON meets RFC 7946
------------------------------

- **Coordinates** are longitude, latitude in degrees. A grid whose
  ellipsoid works in radians is converted. RFC 7946 declares WGS84; a grid
  on another ellipsoid exports that ellipsoid's geodetic coordinates, which
  is fine for a plot but is not a datum shift.
- **Rings are counter-clockwise** (the right-hand rule), starting at the
  cell's north-west vertex. Quad cells have their four vertices exactly,
  since their edges are parallels and meridians. Dart and skew quad cells
  have curved edges, sampled at ``n`` points per edge (default 6): straight
  segments between samples miss the true edge by up to about half a degree
  at resolution 1 and by less at every finer resolution, shrinking with the
  cell. Raise ``n`` for coarse polar cells drawn at large scale.
- **The antimeridian.** A cell straddling longitude 180 becomes a
  ``MultiPolygon`` of two parts, the eastern hemisphere's first, each cut
  exactly at 180 or -180 (section 3.1.9). Which cells straddle it depends
  on the grid: with the default ``lon_0`` only polar cells do, on a rotated
  grid equatorial cells can too. A cell whose edge lies on the
  antimeridian is a single polygon.
- **The poles.** A cap cell, whose boundary is a single parallel, is closed
  over the pole with points at longitudes -180 and 180 on that parallel and
  at the pole, as the RFC recommends. On a grid with even ``N_side`` the
  cells around a pole have the pole as a vertex; that vertex is written as
  a segment along the pole latitude between the cell's two meridian edges.
- Every exported polygon is valid and the polygons of one resolution tile
  the longitude-latitude rectangle exactly; the tests check both over whole
  resolutions on three grids.

Invalid index strings raise ``ValueError`` naming them; an export must not
drop cells silently. Duplicates are kept and the input order is preserved.

GML export, the ``toGML`` code, is planned for v1.0 alongside the published
application schema (issue #103). The other codes of the list (HDF, JSON-LD,
NetCDF, XML) are not planned.

Module reference
----------------

.. automodule:: rhealpixdggs.export
    :members:
    :undoc-members:
    :show-inheritance:
