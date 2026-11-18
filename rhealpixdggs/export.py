"""
Export of cells and query results to standard formats: the interoperation
"broadcast" of OGC Topic 21 v2.0 (requirement 19), whose TranslationType
code list names the formats. This module provides ``toGeoJSON``
(``to_geojson``) and ``toASCII`` (``to_csv``); the library returns the
encodings as a dictionary and a string, and writing files or serving them
is the caller's business.

GeoJSON follows RFC 7946: longitude-latitude order in degrees, exterior
rings counter-clockwise, a cell straddling the antimeridian split into a
MultiPolygon with its eastern part first, and a polar cap cell closed over
the pole with points at longitudes -180 and 180 (section 3.1.9). Quad
cells have their four vertices; the curved edges of dart and skew quad
cells are sampled at `n` points per edge, and a cap's boundary parallel at
``4 * n - 4`` longitudes. RFC 7946 declares WGS84; a grid on another
ellipsoid exports that ellipsoid's geodetic coordinates, and a grid in
radians is converted to degrees.

Each feature's ``id`` is the cell's index string and its properties are
flat scalars, the same columns the CSV has: ``index``, ``resolution``,
``shape`` (``Cell.ellipsoidal_shape``), ``area`` (the cell's surface area
in the ellipsoid's length units squared, metres on WGS84),
``centroid_lon``, ``centroid_lat`` (the cell's direct position,
``Cell.centroid``) and ``nucleus_lon``, ``nucleus_lat`` (its indexing
point). ``geometries`` returns the same geometries as shapely objects, for
geopandas and other OGR-backed writers.

GML follows GML 3.2 and the application schema ``schemas/zone.xsd`` that
ships with the package and is published at ``GML_SCHEMA_LOCATION``:
``to_gml`` writes a ``ZoneCollection`` of ``Zone`` features with the same
attributes, each edge's curve type from OGC Topic 21's ``EA_BoundaryType``
code list (``boundary_types``), and the geometry in EPSG:4326 with
latitude first, as GML's axis order for that CRS requires.
"""

import csv
import io
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from math import isclose
from typing import Any, NamedTuple

import numpy as np
import shapely
from shapely.affinity import translate
from shapely.geometry import MultiPolygon, Polygon, box

from rhealpixdggs.cell import Cell
from rhealpixdggs.dggs import WGS84_003, RHEALPixDGGS
from rhealpixdggs.utils import FloatArray

PROPERTIES = (
    "index",
    "resolution",
    "shape",
    "area",
    "centroid_lon",
    "centroid_lat",
    "nucleus_lon",
    "nucleus_lat",
)


class CellTable(NamedTuple):
    """
    The attributes of a sequence of cells, as ``cell_table`` returns them,
    one entry per cell in input order: the `indices` (index strings),
    `resolution`, `shape` (``Cell.ellipsoidal_shape``), `area` (surface
    area in the ellipsoid's length units squared), and the `centroid` and
    `nucleus` as ``(cells, 2)`` arrays of longitude and latitude in
    degrees.
    """

    indices: np.ndarray
    resolution: np.ndarray
    shape: np.ndarray
    area: FloatArray
    centroid: FloatArray
    nucleus: FloatArray


def _index_strings(indices: Iterable[str | Cell]) -> list[str]:
    return [str(index) for index in indices]


def _check_valid(dggs: RHEALPixDGGS, indices: list[str]) -> np.ndarray:
    """The shapes of `indices`; raise ValueError naming invalid indices."""
    shapes = dggs.shapes(indices)
    invalid = [
        index for index, shape in zip(indices, shapes, strict=True) if shape == ""
    ]
    if invalid:
        shown = ", ".join(repr(index) for index in invalid[:10])
        more = f" and {len(invalid) - 10} more" if len(invalid) > 10 else ""
        raise ValueError(f"invalid cell indices: {shown}{more}")
    return shapes


def _degrees(dggs: RHEALPixDGGS, values: FloatArray) -> FloatArray:
    return np.degrees(values) if dggs.ellipsoid.radians else values


def cell_table(
    indices: Iterable[str | Cell], dggs: RHEALPixDGGS = WGS84_003
) -> CellTable:
    """
    Return the attributes of the cells `indices` (index strings or
    ``Cell`` objects, read as cells of `dggs`) as a ``CellTable``. Raise
    ValueError if any index is invalid.

    EXAMPLES::

        >>> table = cell_table(['P0', 'N'])
        >>> table.shape.tolist(), table.resolution.tolist()
        (['quad', 'cap'], [1, 0])
        >>> table.centroid.round(3).tolist()
        [[-75.0, 26.79], [-180.0, 90.0]]

    """
    strings = _index_strings(indices)
    shapes = _check_valid(dggs, strings)
    resolution = np.array(
        [dggs.index_resolution(index) for index in strings], dtype=np.int64
    )
    area = np.array(
        [dggs.cell_area(int(r), plane=False) for r in resolution], dtype=np.float64
    )
    return CellTable(
        indices=np.array(strings, dtype=str),
        resolution=resolution,
        shape=shapes,
        area=area,
        centroid=_degrees(dggs, dggs.centroids(strings)),
        nucleus=_degrees(dggs, dggs.nuclei(strings)),
    )


def _normalise(dggs: RHEALPixDGGS, rings: FloatArray) -> FloatArray:
    """
    Rings in degrees with longitudes within about 1e-9 degrees of the
    antimeridian snapped to exactly -180 or 180, so that the wrap logic
    sees a cell edge on the antimeridian as lying on it, and latitudes
    clipped to the poles.
    """
    rings = _degrees(dggs, rings).copy()
    lon, lat = rings[..., 0], rings[..., 1]
    at_edge = np.abs(np.abs(lon) - 180) <= 1e-9
    rings[..., 0] = np.where(at_edge, np.copysign(180.0, lon), lon)
    at_pole = np.abs(np.abs(lat) - 90) <= 1e-9
    rings[..., 1] = np.where(at_pole, np.copysign(90.0, lat), np.clip(lat, -90.0, 90.0))
    return rings


def _expand_poles(ring: FloatArray) -> FloatArray:
    """
    A ring with a vertex at a pole, whose longitude there is arbitrary,
    with that vertex replaced by two points on the pole latitude at the
    longitudes of its neighbours: the two edges meeting at a pole are
    meridians, so this is the RFC 7946 rendering of a polygon touching
    the pole. Rings without a pole vertex are returned unchanged.
    """
    at_pole = np.abs(ring[:, 1]) == 90.0
    if not at_pole.any():
        return ring
    points: list[tuple[float, float]] = []
    m = len(ring)
    for k in range(m):
        if not at_pole[k]:
            points.append((float(ring[k, 0]), float(ring[k, 1])))
            continue
        pole = float(ring[k, 1])
        points.append((float(ring[(k - 1) % m, 0]), pole))
        points.append((float(ring[(k + 1) % m, 0]), pole))
    return np.array(points)


def _unwrap(ring: FloatArray) -> tuple[FloatArray, bool]:
    """
    A ring whose longitudes jump by more than half a turn, unwrapped
    eastwards past 180, and whether it genuinely crosses the antimeridian.
    A cell edge on the antimeridian never causes the jump: the boundary
    methods give it the sign that keeps the ring's span under half a turn.
    """
    lon = ring[:, 0]
    if lon.max() - lon.min() <= 180:
        return ring, False
    ring = ring.copy()
    ring[:, 0] = np.where(lon < 0, lon + 360, lon)
    return ring, bool(ring[:, 0].max() > 180)


def _counter_clockwise(rings: FloatArray) -> FloatArray:
    """Reverse the clockwise rings, keeping the first point first."""
    return np.roll(rings[:, ::-1, :], 1, axis=1)


def _split_antimeridian(polygon: Polygon) -> MultiPolygon:
    """
    Split a polygon whose longitudes have been unwrapped eastwards past
    180 into its eastern part (longitudes up to 180) and its western part
    (shifted back by a full turn), both counter-clockwise.
    """
    east = polygon.intersection(box(-180.0, -90.0, 180.0, 90.0))
    west = translate(polygon.intersection(box(180.0, -90.0, 540.0, 90.0)), xoff=-360.0)
    parts = [
        part
        for geometry in (east, west)
        for part in shapely.get_parts(geometry)
        if isinstance(part, Polygon) and not part.is_empty
    ]
    parts = [
        shapely.transform(
            part, lambda c: np.column_stack([np.clip(c[:, 0], -180, 180), c[:, 1]])
        )
        for part in parts
    ]
    return MultiPolygon([shapely.orient_polygons(part) for part in parts])


def _cap_polygon(ring: FloatArray, north: bool) -> Polygon:
    """
    The RFC 7946 polygon of a cap cell from the samples of its boundary
    parallel: the parallel from longitude -180 to 180, closed over the
    pole along the antimeridian's two sides, counter-clockwise.
    """
    latitude = float(ring[0, 1])
    inner = sorted(float(lon) for lon in ring[:, 0] if abs(lon) < 180)
    pole = 90.0 if north else -90.0
    parallel = [
        (-180.0, latitude),
        *((lon, latitude) for lon in inner),
        (180.0, latitude),
    ]
    if north:
        coords = [*parallel, (180.0, pole), (-180.0, pole), parallel[0]]
    else:
        coords = [(-180.0, pole), (180.0, pole), *parallel[::-1], (-180.0, pole)]
    return Polygon(coords)


def geometries(
    indices: Iterable[str | Cell], dggs: RHEALPixDGGS = WGS84_003, n: int = 6
) -> np.ndarray:
    """
    Return the cells `indices` (index strings or ``Cell`` objects, read as
    cells of `dggs`) as an object array of shapely geometries in
    longitude-latitude degrees, in input order: a ``Polygon`` for a cell,
    or a ``MultiPolygon`` of an eastern and a western part for a cell
    straddling the antimeridian. Exterior rings are counter-clockwise.
    Quad cells have their four vertices; dart and skew quad cells are
    sampled at `n` (at least 2) points per curved edge, and a cap's
    boundary parallel at ``4 * n - 4`` longitudes before it is closed over
    the pole. Raise ValueError if any index is invalid.

    EXAMPLES::

        >>> quad, cap, crossing = geometries(['P4', 'N', 'N46'])
        >>> quad.geom_type, len(quad.exterior.coords)
        ('Polygon', 5)
        >>> cap.geom_type, cap.bounds
        ('Polygon', (-180.0, 41.93785391016014, 180.0, 90.0))
        >>> crossing.geom_type, len(crossing.geoms)
        ('MultiPolygon', 2)

    """
    strings = _index_strings(indices)
    shapes = _check_valid(dggs, strings)
    n = max(n, 2)
    result = np.empty(len(strings), dtype=object)
    quad = shapes == "quad"
    cap = shapes == "cap"
    curved = ~quad & ~cap
    for mask, points in ((quad, 2), (curved, n)):
        if not mask.any():
            continue
        rings = _normalise(dggs, dggs.boundary_array(np.array(strings)[mask], n=points))
        # Rings with a vertex at a pole (the cells around a pole of a grid
        # with even N_side) change length when the pole is expanded, so
        # they are built one by one; the rest go through in one array.
        touches_pole = (np.abs(rings[:, :, 1]) == 90.0).any(axis=1)
        polygons = np.empty(len(rings), dtype=object)
        if (~touches_pole).any():
            plain = rings[~touches_pole]
            lon = plain[:, :, 0]
            wraps = lon.max(axis=1) - lon.min(axis=1) > 180
            plain[wraps, :, 0] = np.where(lon[wraps] < 0, lon[wraps] + 360, lon[wraps])
            crosses = wraps & (plain[:, :, 0].max(axis=1) > 180)
            built = np.asarray(
                shapely.polygons(_counter_clockwise(plain)), dtype=object
            )
            for k in np.flatnonzero(crosses):
                built[k] = _split_antimeridian(built[k])
            polygons[~touches_pole] = built
        for k in np.flatnonzero(touches_pole):
            ring, crosses_meridian = _unwrap(_expand_poles(rings[k]))
            polygon = Polygon(np.roll(ring[::-1], 1, axis=0))
            polygons[k] = _split_antimeridian(polygon) if crosses_meridian else polygon
        result[mask] = polygons
    if cap.any():
        rings = _normalise(dggs, dggs.boundary_array(np.array(strings)[cap], n=n))
        for k, ring in zip(np.flatnonzero(cap), rings, strict=True):
            result[k] = _cap_polygon(ring, north=strings[k][0] == "N")
    return result


def _coordinates(geometry: Any) -> list[Any]:
    """GeoJSON coordinate lists, of Python floats, for a Polygon or MultiPolygon."""
    if geometry.geom_type == "Polygon":
        return [[[float(x), float(y)] for x, y in geometry.exterior.coords]]
    return [_coordinates(part) for part in geometry.geoms]


def to_geojson(
    indices: Iterable[str | Cell], dggs: RHEALPixDGGS = WGS84_003, n: int = 6
) -> dict[str, Any]:
    """
    Return the cells `indices` (index strings or ``Cell`` objects, read as
    cells of `dggs`) as a GeoJSON FeatureCollection dictionary of builtin
    Python types, ready for ``json.dumps``: one Feature per cell in input
    order, its ``id`` the index string, its geometry as ``geometries``
    builds it with `n` points per curved edge, and the flat properties
    ``index``, ``resolution``, ``shape``, ``area``, ``centroid_lon``,
    ``centroid_lat``, ``nucleus_lon`` and ``nucleus_lat``. Raise
    ValueError if any index is invalid.

    EXAMPLES::

        >>> feature = to_geojson(['P0'])['features'][0]
        >>> feature['id'], feature['geometry']['type']
        ('P0', 'Polygon')
        >>> sorted(feature['properties'])
        ['area', 'centroid_lat', 'centroid_lon', 'index', 'nucleus_lat', 'nucleus_lon', 'resolution', 'shape']
        >>> feature['properties']['shape'], feature['properties']['resolution']
        ('quad', 1)

    """
    strings = _index_strings(indices)
    table = cell_table(strings, dggs)
    geoms = geometries(strings, dggs, n)
    features = []
    for k, index in enumerate(strings):
        features.append(
            {
                "type": "Feature",
                "id": index,
                "geometry": {
                    "type": geoms[k].geom_type,
                    "coordinates": _coordinates(geoms[k]),
                },
                "properties": {
                    "index": index,
                    "resolution": int(table.resolution[k]),
                    "shape": str(table.shape[k]),
                    "area": float(table.area[k]),
                    "centroid_lon": float(table.centroid[k, 0]),
                    "centroid_lat": float(table.centroid[k, 1]),
                    "nucleus_lon": float(table.nucleus[k, 0]),
                    "nucleus_lat": float(table.nucleus[k, 1]),
                },
            }
        )
    return {"type": "FeatureCollection", "features": features}


def to_csv(
    indices: Iterable[str | Cell], dggs: RHEALPixDGGS = WGS84_003, delimiter: str = ","
) -> str:
    """
    Return the attributes of the cells `indices` (index strings or
    ``Cell`` objects, read as cells of `dggs`) as delimited text: a header
    line then one line per cell in input order, with the columns
    ``index``, ``resolution``, ``shape``, ``area``, ``centroid_lon``,
    ``centroid_lat``, ``nucleus_lon``, ``nucleus_lat`` separated by
    `delimiter`, lines ending in a newline, floats written so that
    ``float()`` reads them back exactly. Raise ValueError if any index is
    invalid.

    EXAMPLES::

        >>> print(to_csv(['P0', 'N']), end='')
        index,resolution,shape,area,centroid_lon,centroid_lat,nucleus_lon,nucleus_lat
        P0,1,quad,9445659661557.2,-74.99999999999999,26.790327229701848,-74.99999999999999,26.490118751439734
        N,0,cap,85010936954014.78,-180.0,90.0,-180.0,90.0

    """
    table = cell_table(indices, dggs)
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter=delimiter, lineterminator="\n")
    writer.writerow(PROPERTIES)
    for k in range(len(table.indices)):
        writer.writerow(
            [
                str(table.indices[k]),
                int(table.resolution[k]),
                str(table.shape[k]),
                float(table.area[k]),
                float(table.centroid[k, 0]),
                float(table.centroid[k, 1]),
                float(table.nucleus[k, 0]),
                float(table.nucleus[k, 1]),
            ]
        )
    return buffer.getvalue()


# GML -----------------------------------------------------------------------

#: The namespace of the GML application schema, ``schemas/zone.xsd``.
GML_NAMESPACE = "https://manaakiwhenua.github.io/rhealpixdggs-py/schemas/zone"
#: Where the application schema is published.
GML_SCHEMA_LOCATION = "https://manaakiwhenua.github.io/rhealpixdggs-py/schemas/zone.xsd"
_GML = "http://www.opengis.net/gml/3.2"
_XSI = "http://www.w3.org/2001/XMLSchema-instance"
_XLINK = "http://www.w3.org/1999/xlink"
_PLANAR_SRS = "#rhealpix"
_QUAD_EDGES = ("smallCircle", "greatEllipse", "smallCircle", "greatEllipse")


def boundary_types(
    index: str | Cell, dggs: RHEALPixDGGS = WGS84_003
) -> tuple[str, ...]:
    """
    The curve type of each edge of the cell `index` (an index string or
    ``Cell``, read as a cell of `dggs`) on the ellipsoid, from the
    ``EA_BoundaryType`` code list of OGC Topic 21: one entry per planar
    edge clockwise from the upper-left vertex (top, right, bottom, left),
    or one entry for a cap cell, whose four edges lie on one parallel.

    A quad cell's edges are parallels (``smallCircle``) and meridians
    (``greatEllipse``). The diagonals of a polar square cut it into four
    sectors, and within a sector the lines parallel to its rim are
    parallels of latitude, so a polar cell's edge is a ``smallCircle``
    when it runs parallel to the rim of its sector (cell edges never cross
    a diagonal), a ``greatEllipse`` when its line passes through the pole
    (the square's centre lines, cell edges only for even ``N_side``), and
    otherwise the image of a straight planar line, ``projectedLine``.
    Raise ValueError for an invalid index.

    EXAMPLES::

        >>> boundary_types('P0')
        ('smallCircle', 'greatEllipse', 'smallCircle', 'greatEllipse')
        >>> boundary_types('N')
        ('smallCircle',)
        >>> boundary_types('N1')
        ('smallCircle', 'projectedLine', 'smallCircle', 'projectedLine')
        >>> from rhealpixdggs.dggs import WGS84_002
        >>> boundary_types('N3', WGS84_002)
        ('greatEllipse', 'smallCircle', 'smallCircle', 'greatEllipse')

    """
    cell = index if isinstance(index, Cell) else dggs.cell(dggs.parse_index(index))
    shape = cell.ellipsoidal_shape
    if shape == "quad":
        return _QUAD_EDGES
    if shape == "cap":
        return ("smallCircle",)
    (x0, x1), (y0, y1) = cell.xy_range()
    (fx0, fx1), (fy0, fy1) = cell.rdggs.cell((cell.suid[0],)).xy_range()
    xc, yc = (fx0 + fx1) / 2, (fy0 + fy1) / 2
    tol = (fx1 - fx0) * 1e-9

    def edge(along: float, across: float) -> str:
        # `along` is the edge's offset from the pole across its direction,
        # `across` its midpoint's offset along it.
        if isclose(along, 0, abs_tol=tol):
            return "greatEllipse"
        return "smallCircle" if abs(along) > abs(across) else "projectedLine"

    mx, my = (x0 + x1) / 2 - xc, (y0 + y1) / 2 - yc
    return (edge(y1 - yc, mx), edge(x1 - xc, my), edge(y0 - yc, mx), edge(x0 - xc, my))


def _pos(lon: float, lat: float, planar: bool = False) -> str:
    """A GML position: easting northing in the plane, else latitude first as
    EPSG:4326's axis order has it."""
    if planar:
        return f"{float(lon)!r} {float(lat)!r}"
    return f"{float(lat)!r} {float(lon)!r}"


def _pos_list(ring: Any, planar: bool = False) -> str:
    return " ".join(_pos(x, y, planar) for x, y in ring.coords)


def _gml_crs(parent: ET.Element, dggs: RHEALPixDGGS) -> None:
    """The grid's CRS as an inline gml:ProjectedCRS with gml:id "rhealpix":
    the conversion's parameters, with the base CRS and the Cartesian CS by
    EPSG URN. Raise ValueError unless the base CRS has an EPSG code."""
    base = dggs.crs.geodetic_crs
    authority = base.to_authority() if base is not None else None
    if authority is None:
        raise ValueError(
            "planar GML needs a grid on WGS84, whose base CRS has an EPSG code"
        )
    gml, xlink = f"{{{_GML}}}", f"{{{_XLINK}}}"
    operation = dggs.crs.coordinate_operation
    assert operation is not None
    crs = ET.SubElement(parent, f"{gml}ProjectedCRS", {f"{gml}id": "rhealpix"})
    ET.SubElement(crs, f"{gml}identifier", {"codeSpace": GML_NAMESPACE}).text = (
        dggs.crs.name
    )
    ET.SubElement(crs, f"{gml}scope").text = (
        "The planar coordinates of the rHEALPix grid."
    )
    conversion = ET.SubElement(
        ET.SubElement(crs, f"{gml}conversion"),
        f"{gml}Conversion",
        {f"{gml}id": "rhealpix.conversion"},
    )
    ET.SubElement(conversion, f"{gml}identifier", {"codeSpace": GML_NAMESPACE}).text = (
        operation.name
    )
    ET.SubElement(conversion, f"{gml}scope").text = (
        "The rHEALPix projection with this grid's parameters."
    )
    method = ET.SubElement(
        ET.SubElement(conversion, f"{gml}method"),
        f"{gml}OperationMethod",
        {f"{gml}id": "rhealpix.method"},
    )
    ET.SubElement(
        method, f"{gml}identifier", {"codeSpace": "https://proj.org"}
    ).text = operation.method_name
    ET.SubElement(method, f"{gml}formula").text = (
        "The rHEALPix projection, https://proj.org/operations/projections/rhealpix.html"
    )
    ET.SubElement(method, f"{gml}sourceDimensions").text = "2"
    ET.SubElement(method, f"{gml}targetDimensions").text = "2"
    units = {
        "degree": "urn:ogc:def:uom:EPSG::9102",
        "unity": "urn:ogc:def:uom:EPSG::9201",
    }
    for param in operation.params:
        definition = ET.SubElement(
            ET.SubElement(method, f"{gml}parameter"),
            f"{gml}OperationParameter",
            {f"{gml}id": f"rhealpix.{param.name}"},
        )
        ET.SubElement(
            definition, f"{gml}identifier", {"codeSpace": "https://proj.org"}
        ).text = param.name
    for param in operation.params:
        value = ET.SubElement(
            ET.SubElement(conversion, f"{gml}parameterValue"), f"{gml}ParameterValue"
        )
        ET.SubElement(value, f"{gml}value", {"uom": units[param.unit_name]}).text = (
            format(param.value, ".15g")
        )
        ET.SubElement(
            value,
            f"{gml}operationParameter",
            {f"{xlink}href": f"#rhealpix.{param.name}"},
        )
    ET.SubElement(
        crs,
        f"{gml}baseGeodeticCRS",
        {f"{xlink}href": f"urn:ogc:def:crs:{authority[0]}::{authority[1]}"},
    )
    ET.SubElement(
        crs, f"{gml}cartesianCS", {f"{xlink}href": "urn:ogc:def:cs:EPSG::4400"}
    )


def _gml_point(
    parent: ET.Element,
    gml_id: str,
    lon: float,
    lat: float,
    srs: dict[str, str],
    planar: bool = False,
) -> None:
    point = ET.SubElement(parent, f"{{{_GML}}}Point", {f"{{{_GML}}}id": gml_id, **srs})
    ET.SubElement(point, f"{{{_GML}}}pos").text = _pos(lon, lat, planar)


def _gml_polygon(
    parent: ET.Element,
    gml_id: str,
    polygon: Any,
    srs: dict[str, str],
    planar: bool = False,
) -> None:
    element = ET.SubElement(
        parent, f"{{{_GML}}}Polygon", {f"{{{_GML}}}id": gml_id, **srs}
    )
    exterior = ET.SubElement(element, f"{{{_GML}}}exterior")
    ring = ET.SubElement(exterior, f"{{{_GML}}}LinearRing")
    ET.SubElement(ring, f"{{{_GML}}}posList", {"srsDimension": "2"}).text = _pos_list(
        polygon.exterior, planar
    )


def to_gml(
    indices: Iterable[str | Cell],
    dggs: RHEALPixDGGS = WGS84_003,
    n: int = 6,
    codespace: str | None = None,
    planar: bool = False,
) -> str:
    """
    Return the cells `indices` (index strings or ``Cell`` objects, read as
    cells of `dggs`) as a GML 3.2 document: a ``ZoneCollection`` of the
    application schema at ``GML_SCHEMA_LOCATION`` with one ``Zone`` per
    cell in input order. Each zone's ``gml:identifier`` is the index
    string in `codespace`, by default the grid's
    ``codespace``, the URI of its DGGRS definition; a grid with none
    needs `codespace` given, since the identifier is unique only with it.
    The zone carries the resolution, shape, area, the curve type of each
    edge (``boundary_types``), the centroid as its position, the nucleus,
    and the geometry ``geometries`` builds with `n` points per curved
    edge: a ``gml:Polygon``, or a ``gml:MultiSurface`` of two for a cell
    straddling the antimeridian. Geometry is written latitude first, as
    GML requires for ``urn:ogc:def:crs:EPSG::4326``, the ``srsName`` of a
    grid on WGS84; a grid on another ellipsoid writes its own geodetic
    coordinates without an ``srsName``. With `planar`, the geometry is
    instead the cell's exact square in the grid's own CRS, which the
    document carries as a ``gml:ProjectedCRS`` that the geometries refer
    to by ``srsName="#rhealpix"``; this needs a grid on WGS84. Raise
    ValueError if any index is invalid, no codespace is known, or a
    planar document is asked of another ellipsoid.

    EXAMPLES::

        >>> from rhealpixdggs.dggs import WGS84_003_OGC
        >>> text = to_gml(['P0', 'N'], WGS84_003_OGC)
        >>> print(text[:text.index('>') + 1])
        <?xml version='1.0' encoding='utf-8'?>
        >>> import xml.etree.ElementTree as ET
        >>> root = ET.fromstring(text)
        >>> zones = root.findall('{%s}member/{%s}Zone' % (GML_NAMESPACE, GML_NAMESPACE))
        >>> [z.get('{http://www.opengis.net/gml/3.2}id') for z in zones]
        ['P0', 'N']
        >>> zones[0].find('{http://www.opengis.net/gml/3.2}identifier').get('codeSpace')
        'https://www.opengis.net/def/dggrs/OGC/1.0/rHEALPix'
        >>> [b.text for b in zones[1].iter('{%s}boundaryType' % GML_NAMESPACE)]
        ['smallCircle']

    """
    if codespace is None:
        codespace = dggs.codespace
    if codespace is None:
        raise ValueError(
            "the grid has no codespace: give it a uri, or pass codespace, so the "
            "zone identifiers are unique"
        )
    strings = _index_strings(indices)
    table = cell_table(strings, dggs)
    if planar:
        srs = {"srsName": _PLANAR_SRS}
        geoms = np.array(
            [
                Polygon(dggs.cell(dggs.parse_index(index)).vertices(plane=True)[::-1])
                for index in strings
            ],
            dtype=object,
        )
        centroid = dggs.centroids(strings, plane=True)
        nucleus = dggs.nuclei(strings, plane=True)
    else:
        geoms = geometries(strings, dggs, n)
        centroid, nucleus = table.centroid, table.nucleus
        base = dggs.crs.geodetic_crs
        authority = base.to_authority() if base is not None else None
        srs = (
            {"srsName": f"urn:ogc:def:crs:{authority[0]}::{authority[1]}"}
            if authority
            else {}
        )
    ET.register_namespace("gml", _GML)
    ET.register_namespace("rhp", GML_NAMESPACE)
    ET.register_namespace("xsi", _XSI)
    ET.register_namespace("xlink", _XLINK)
    rhp = f"{{{GML_NAMESPACE}}}"
    gml = f"{{{_GML}}}"
    root = ET.Element(
        f"{rhp}ZoneCollection",
        {
            f"{gml}id": "zones",
            f"{{{_XSI}}}schemaLocation": f"{GML_NAMESPACE} {GML_SCHEMA_LOCATION}",
        },
    )
    ET.SubElement(root, f"{rhp}dggrs").text = codespace
    if planar:
        _gml_crs(ET.SubElement(root, f"{rhp}crs"), dggs)
    for k, index in enumerate(strings):
        zone = ET.SubElement(
            ET.SubElement(root, f"{rhp}member"), f"{rhp}Zone", {f"{gml}id": index}
        )
        ET.SubElement(zone, f"{gml}identifier", {"codeSpace": codespace}).text = index
        ET.SubElement(zone, f"{rhp}resolution").text = str(int(table.resolution[k]))
        ET.SubElement(zone, f"{rhp}shape").text = str(table.shape[k])
        ET.SubElement(zone, f"{rhp}area", {"uom": "m2"}).text = repr(
            float(table.area[k])
        )
        for curve in boundary_types(index, dggs):
            ET.SubElement(zone, f"{rhp}boundaryType").text = curve
        _gml_point(
            ET.SubElement(zone, f"{rhp}position"),
            f"{index}.position",
            centroid[k, 0],
            centroid[k, 1],
            srs,
            planar,
        )
        _gml_point(
            ET.SubElement(zone, f"{rhp}nucleus"),
            f"{index}.nucleus",
            nucleus[k, 0],
            nucleus[k, 1],
            srs,
            planar,
        )
        geometry = ET.SubElement(zone, f"{rhp}geometry")
        if geoms[k].geom_type == "Polygon":
            _gml_polygon(geometry, f"{index}.geometry", geoms[k], srs, planar)
        else:
            multi = ET.SubElement(
                geometry, f"{gml}MultiSurface", {f"{gml}id": f"{index}.geometry", **srs}
            )
            for j, part in enumerate(geoms[k].geoms):
                _gml_polygon(
                    ET.SubElement(multi, f"{gml}surfaceMember"),
                    f"{index}.geometry.{j}",
                    part,
                    {},
                )
    ET.indent(root)
    return ET.tostring(root, encoding="unicode", xml_declaration=True) + "\n"
