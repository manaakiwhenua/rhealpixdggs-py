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
"""

import csv
import io
from collections.abc import Iterable
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
    invalid = [index for index, shape in zip(indices, shapes) if shape == ""]
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
    resolution = np.array([len(index) - 1 for index in strings], dtype=np.int64)
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
    eastwards past 180, and whether it genuinely crosses the antimeridian
    (a ring whose east edge lies on it does not).
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
        for k, ring in zip(np.flatnonzero(cap), rings):
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
