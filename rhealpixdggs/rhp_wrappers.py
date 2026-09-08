from typing import Any, Literal
from warnings import warn

import numpy as np
import shapely
from shapely import is_valid_reason
from shapely.affinity import translate
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon, box

# ======== Messages and constants ======== #
# Pre-defined DGGS with WGS84 ellipsoid, coordinates in degrees, n == 3 to subdivide
# cell sides, and both N and S polar cube face attached to O equatorial cube face:
# N
# O P Q R
# S
import rhealpixdggs.dggs as rhp_dggs

# List of resolution 0 cell addresses (i.e. cube faces)
from rhealpixdggs.cell import CELLS0, Cell
from rhealpixdggs.dggs import WGS84_003, RHEALPixDGGS

# Warnings
PARENT_RESOLUTION_WARNING = "WARNING: You requested a parent resolution that is higher than the cell resolution. Returning the cell address itself."
CHILD_RESOLUTION_WARNING = "WARNING: You requested a child resolution that is lower than the cell resolution. Returning the cell address itself."
CELL_CENTRE_WARNING = "WARNING: You requested a centre cell for a DGGS that has an even number of cells on a side. Returning None."
POLYFILL_GEOMETRY_WARNING = "WARNING: Empty or missing geometry, unsupported geometry type (not Polygon or MultiPolygon), or geometry with no area. Returning None."
LINETRACE_GEOMETRY_WARNING = "WARNING: Empty or missing line geometry, unsupported line type (not LineString or MultiLineString), or line with no length. Returning None."


# ======== Main API ======== #


def geo_to_rhp(
    lat: float,
    lng: float,
    resolution: int,
    plane: bool = True,
    dggs: RHEALPixDGGS = WGS84_003,
) -> str | None:
    """
    Turn a latitute and longitude (in degrees) into an rHEALPix cell address at
    the requested resolution.

    Mostly passes through the parameters to the function turning coordinate points
    into cells, but converts the address tuple from the resulting cell into a
    string.

    Returns None if no cell matching the coordinates is found.

    EXAMPLES::
        >>> geo_to_rhp(-43.738058, -176.260506, 9, plane=False)
        'S001450634'
        >>> geo_to_rhp(-43.738379,-176.258807, 9, plane=False)
        'S001450635'
    """
    # Get the grid cell corresponding to the coordinates
    cell = dggs.cell_from_point(resolution, (lng, lat), plane)

    # Bail out if there's no matching cell
    if cell is None:
        return None

    # Return the cell ID after converting int digits to str
    return str(cell)


def rhp_to_geo(
    rhpindex: str,
    geo_json: bool = True,
    plane: bool = True,
    dggs: RHEALPixDGGS = WGS84_003,
) -> tuple[float, float] | None:
    """
    Look up the centroid (in degrees) of the cell identified by rhpindex.

    Returns None if the cell index is invalid.

    If geojson is requested as the output format:
        - Will return a (longitude, latitude) coordinate pair.

    if geojson is NOT requested as the output format:
        - Will return a (latitude, longitude) coordinate pair in order to be consistent with
          h3 coordinate ordering.

    EXAMPLES::

        >>> rhp_to_geo('S001450634', True, False)
        (-176.2606635452476, -43.7365450535837)
        >>> rhp_to_geo('S001450635', True, False)
        (-176.25592420875034, -43.7365450535837)
        >>> rhp_to_geo('NotACellId', True, False)
    """
    # Stop early if the cell index is invalid
    if not rhp_is_valid(rhpindex, dggs):
        return None

    # Grab cell centroid matching rhpindex string
    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    cell = dggs.cell(suid)
    centroid = cell.centroid(plane=plane)

    # rhealpix coordinates come out natively as lng/lat, h3 ones as lat/lng
    if not geo_json:
        # lng/lat -> lat/lng to make it consistent with h3
        centroid = centroid[::-1]

    return (float(centroid[0]), float(centroid[1]))


def rhp_to_parent(
    rhpindex: str,
    res: int | None = None,
    verbose: bool = True,
    dggs: RHEALPixDGGS = WGS84_003,
) -> str | None:
    """
    Returns parent of rhpindex at resolution res (immediate parent if res == None).

    Returns None if the cell index is invalid.

    EXAMPLES::

        >>> rhp_to_parent('S001450634', res=0)
        'S'
        >>> rhp_to_parent('S001450634')
        'S00145063'
        >>> rhp_to_parent('INVALID')
    """
    # Stop early if the cell index is invalid
    if not rhp_is_valid(rhpindex, dggs):
        return None

    # Top-level cells are their own parent, regardless of the requested resolution (by convention)
    child_res = len(rhpindex) - 1
    if child_res < 1:
        return rhpindex

    # res == None returns the first address up (by convention)
    elif res is None:
        return rhpindex[:-1]

    # Handle mismatch between cell resolution and requested parent resolution
    elif res > child_res:
        if verbose:
            warn(PARENT_RESOLUTION_WARNING)
        return rhpindex

    # Standard case (including child_res == res)
    else:
        return rhpindex[: res + 1]


def rhp_to_center_child(
    rhpindex: str,
    res: int | None = None,
    verbose: bool = True,
    dggs: RHEALPixDGGS = WGS84_003,
) -> str | None:
    """
    Returns central child of rhpindex at resolution res (immediate central
    child if res == None).

    Returns None if the cell index is invalid.

    Returns None if the DGGS has an even number of cells on a side.

    EXAMPLES::

        >>> rhp_to_center_child('S001450634')
        'S0014506344'
        >>> rhp_to_center_child('S001450634', res=13)
        'S0014506344444'
        >>> rhp_to_center_child('S001450634', res=9)
        'S001450634'
        >>> rhp_to_center_child('S001450634', res=10)
        'S0014506344'
        >>> rhp_to_center_child('INVALID')
    """
    # Stop early if the cell index is invalid
    if not rhp_is_valid(rhpindex, dggs):
        return None

    # DGGSs with even numbers of cells on a side never have a cell at the centre
    if (dggs.N_side % 2) == 0:
        if verbose:
            warn(CELL_CENTRE_WARNING)
        return None

    # Handle mismatch between cell resolution and requested child resolution
    parent_res = len(rhpindex) - 1
    if res is not None and res < parent_res:
        if verbose:
            warn(CHILD_RESOLUTION_WARNING)
        return rhpindex

    # Standard case (including parent_res == res)
    else:
        # res == None returns the central child from one level down (by convention)
        added_levels = 1 if res is None else res - parent_res

        # Derive index of centre child and append that to rhpindex
        # NOTE: only works for odd values of N_side
        c_index = int((dggs.N_side**2 - 1) / 2)

        # Append the required number of child digits to cell index
        child_index = rhpindex + "".join(str(c_index) for _ in range(added_levels))

        return child_index


def rhp_to_geo_boundary(
    rhpindex: str,
    geo_json: bool = True,
    plane: bool = True,
    dggs: RHEALPixDGGS = WGS84_003,
) -> tuple[tuple[float, float], ...] | None:
    """
    Extract the corner coordinates of a cell at a given cell ID and returns them as
    a tuple of coordinate pairs (in degrees).

    Returns None if the cell index is invalid.

    If geojson is requested as the output format:
        - Will return (longitude, latitude) coordinate pairs.
        - Will repeat the first vertex and insert it at the end if geojson is requested as
          the output format.

    If geojson is NOT requested as the output format:
        - Will return (latitude, longitude) coordinate pairs in order to be consistent with
          rHEALPix coordinate ordering.

    EXAMPLES::

        >>> rhp_to_geo_boundary('S001450634', False, False)
        ((-43.73395872598705, -176.26086040756147), (-43.73395872598705, -176.25612132062557), (-43.73913136381169, -176.26046658591815), (-43.73913136381169, -176.2652061719943))
        >>> rhp_to_geo_boundary('S001450634', True, False)
        ((-176.26086040756147, -43.73395872598705), (-176.25612132062557, -43.73395872598705), (-176.26046658591815, -43.73913136381169), (-176.2652061719943, -43.73913136381169), (-176.26086040756147, -43.73395872598705))
        >>> rhp_to_geo_boundary('INVALID')
    """
    # Stop early if the cell index is invalid
    if not rhp_is_valid(rhpindex, dggs):
        return None

    # Grab the cell vertices (includes non-corner point in darts if plane == False)
    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    cell = dggs.cell(suid)
    verts = tuple(cell.vertices(plane=plane))

    # rhealpix coordinates come out natively as lng/lat, h3 ones as lat/lng
    # Neither has the repeated vertex that geo_json wants so it's inserted here when needed
    if not geo_json:
        # lng/lat -> lat/lng to make it consistent with h3
        verts = tuple(v[::-1] for v in verts)
    else:
        # last point same as first
        verts += (verts[0],)

    return tuple((float(v[0]), float(v[1])) for v in verts)


def rhp_get_resolution(rhpindex: str, dggs: RHEALPixDGGS = WGS84_003) -> int | None:
    """
    Returns the resolution of a given cell index (or None if invalid).

    EXAMPLES::

        >>> rhp_get_resolution('S001450634')
        9
        >>> rhp_get_resolution('INVALID')
    """
    if not rhp_is_valid(rhpindex, dggs):
        return None

    return len(rhpindex) - 1


def rhp_get_base_cell(rhpindex: str, dggs: RHEALPixDGGS = WGS84_003) -> str | None:
    """
    Returns the resolution 0 cell id of a given cell index (or None if invalid).

    EXAMPLES::
        >>> rhp_get_base_cell('S001450634')
        'S'
        >>> rhp_get_base_cell('INVALID')
    """
    if not rhp_is_valid(rhpindex, dggs):
        return None

    return rhpindex[0]


def rhp_is_valid(rhpindex: str, dggs: RHEALPixDGGS = WGS84_003) -> bool:
    """
    Checks if the given cell address is valid within the DGGS

    EXAMPLES::
        >>> rhp_is_valid('S001450634')
        True
        >>> rhp_is_valid('N')
        True
        >>> rhp_is_valid('INVALID')
        False
        >>> rhp_is_valid(None)
        False
        >>> rhp_is_valid(123456)
        Traceback (most recent call last):
        ...
        TypeError: object of type 'int' has no len()
    """
    # Empty strings are invalid
    if rhpindex is None or len(rhpindex) == 0:
        return False

    # Addresses that don't start with the resolution 0 face are invalid
    if rhpindex[0] not in CELLS0:
        return False

    # Addresses that have digits out of range are invalid
    num_subcells = dggs.N_side**2
    for d in rhpindex[1:]:
        if not d.isdigit() or (int(d) >= num_subcells):
            return False

    # Passed all checks - must be the real thing
    return True


def cell_area(
    rhpindex: str,
    unit: Literal["km^2", "m^2"] = "km^2",
    plane: bool = True,
    dggs: RHEALPixDGGS = WGS84_003,
) -> float | None:
    """
    Returns the area of a cell in the requested unit (or None if rhpindex is invalid).

    TODO: investigate use case where unit is 'rads^2'

    EXAMPLES::

        >>> cell_area('S001450634', unit='m^2', plane=False)
        219428.08748562282
        >>> cell_area('S001450635', unit='m^2', plane=False)
        219428.08748562282
        >>> cell_area('S00145063', unit='km^2', plane=False)
        1.9748527873706059
        >>> cell_area('INVALID', unit='km^2', plane=False)
    """
    if not rhp_is_valid(rhpindex, dggs):
        return None

    # Grab cell area in native unit (m^2)
    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    cell = dggs.cell(suid)
    area = float(cell.area(plane=plane))

    # Scale area if needed
    if unit == "km^2":
        area = area / 10**6

    return float(area)


def cell_ring(
    rhpindex: str, k: int = 1, dggs: RHEALPixDGGS = WGS84_003
) -> list[str] | None:
    """
    Returns the ring of cell indices around rhpindex at distance k (in graph
    distance -- shortest number of edge/corner hops -- not straight-line
    distance), or None if rhpindex is invalid.

    Also returns None if k < 0.

    Returns [ rhpindex ] if k == 0 (by convention).

    Built by growing the ring outward one shell at a time from `neighbor()`
    (edge-adjacent) and `diagonal_neighbor()` (corner-adjacent) steps, so
    this is correct near a cube corner too: exactly 3 cells meet at each of
    the cube's 8 corners rather than 4, so a ring passing near one has fewer
    cells there than the usual 8*k for an interior ring, rather than an
    incorrect or duplicate cell.

    EXAMPLES::

        >>> cell_ring('S001450634')
        ['S001450631', 'S001450632', 'S001450635', 'S001450638', 'S001450637', 'S001450636', 'S001450633', 'S001450630']
        >>> cell_ring('S001450634', k=2)
        ['S001450607', 'S001450608', 'S001450606', 'S001450616', 'S001450640', 'S001450643', 'S001450646', 'S001450670', 'S001450662', 'S001450661', 'S001450660', 'S001442882', 'S001442858', 'S001442855', 'S001442852', 'S001442828']
        >>> cell_ring('S001450634', k=0)
        ['S001450634']
        >>> cell_ring('S', k=1)
        ['O', 'P', 'Q', 'R']
        >>> cell_ring('S', k=-1)
    """
    if not rhp_is_valid(rhpindex, dggs) or (k < 0):
        return None

    if k == 0:
        return [rhpindex]

    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    center = dggs.cell(suid)
    return [str(cell) for cell in _rings_up_to(center, k)[k]]


def k_ring(
    rhpindex: str, k: int = 1, dggs: RHEALPixDGGS = WGS84_003
) -> list[str] | None:
    """
    Returns the k-ring (the centre cell together with every ring out to
    distance k, in graph distance -- see cell_ring()) of cell indices around
    rhpindex, or None if rhpindex is invalid.

    Also returns None if k < 0.

    EXAMPLES::

        >>> k_ring('S001450634')
        ['S001450634', 'S001450631', 'S001450632', 'S001450635', 'S001450638', 'S001450637', 'S001450636', 'S001450633', 'S001450630']
        >>> k_ring('S001450634', k=2)
        ['S001450634', 'S001450631', 'S001450632', 'S001450635', 'S001450638', 'S001450637', 'S001450636', 'S001450633', 'S001450630', 'S001450607', 'S001450608', 'S001450606', 'S001450616', 'S001450640', 'S001450643', 'S001450646', 'S001450670', 'S001450662', 'S001450661', 'S001450660', 'S001442882', 'S001442858', 'S001442855', 'S001442852', 'S001442828']
        >>> k_ring('S001450634', k=0)
        ['S001450634']
        >>> k_ring('S001450634', k=-1)
        >>> k_ring('INVALID')

    """
    if not rhp_is_valid(rhpindex, dggs) or (k < 0):
        return None

    if k == 0:
        return [rhpindex]

    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    center = dggs.cell(suid)
    return [str(cell) for ring in _rings_up_to(center, k) for cell in ring]


def polyfill(
    geometry: Polygon | MultiPolygon,
    res: int,
    plane: bool = True,
    compress: bool = False,
    verbose: bool = False,
    dggs: RHEALPixDGGS = WGS84_003,
    containment: Literal["center", "full", "overlapping"] = "center",
    min_res: int = 0,
) -> set[str] | None:
    """
    Turns the area contained in a shapely polygon or multipolygon into a set of cell
    indices at the requested resolution. Which cells count is chosen by
    `containment`, following H3's modes:

    - ``"center"`` (default): a cell is included if its centroid is inside the
      geometry defined by the boundaries and holes.
    - ``"full"``: a cell is included if the whole cell lies within the geometry
      (its boundary may touch the geometry's).
    - ``"overlapping"``: a cell is included if any part of it meets the
      geometry, so the result covers the geometry.

    ``full`` and ``overlapping`` test the cell's boundary polygon. In the plane
    that polygon is exact. On the ellipsoid, equatorial cells are
    longitude-latitude aligned and exact too; polar cells have curved edges,
    which are represented by ``6`` points per edge (``Cell.boundary(n=6)``), so
    a decision there is made against that piecewise-linear approximation. A
    cap cell, which spans every longitude, is fully inside only a geometry
    that covers the whole latitude band from the cap's boundary to the pole,
    and overlaps any geometry reaching that latitude.

    The candidate cells are those of the geometry's bounding box
    (``RHEALPixDGGS.cells_in_box``), and every test runs on all candidates
    at once: centroids from ``RHEALPixDGGS.centroids``, boundaries from
    ``RHEALPixDGGS.boundary_array``, predicates from shapely's vectorised
    functions.

    Returns an empty set if no cell centroids fall within the input geometry.

    Returns None if the geom_type field in the input geometry is anything other than
    'Polygon' or 'MultiPolygon'.

    Returns None if the geometry is empty, or if it has no area.

    Returns None if no cells match the geometry for some reason.

    Returns None if the geometry is invalid in other ways, e.g. if a point on a hole
    boundary is outside the exterior boundary of its polygon, or if two polygons in a
    multipolygon overlap.

    Polygon edges are straight in the input coordinate space (shapely
    geometries are planar) and do not wrap around the antimeridian; a
    region straddling it must be split into a MultiPolygon along the
    antimeridian by the caller first, as is standard for planar GIS
    geometry.

    The fill descends the hierarchy from resolution `min_res`: a cell wholly
    inside the geometry is kept whole, a cell missing it is dropped, and only
    the cells the boundary passes through are split, down to resolution
    `res`, where the remaining cells are decided by `containment`. So the
    work grows with the length of the boundary, not the area. With
    `compress` = False every kept coarse cell is expanded into its
    resolution-`res` descendants, giving exactly the cells a cell-by-cell
    test at `res` would give. With `compress` = True the coarse cells are
    returned as they are, a mixed-resolution cover with no cell coarser than
    `min_res` and no complete group of siblings left unmerged; so every cell
    of the result has exactly one ancestor at resolution `min_res`, which
    suits partitioning by that resolution. "Wholly inside" is judged by the
    cell's boundary polygon for ``full`` and ``overlapping`` and, for
    ``center``, by the cell's longitude-latitude bounding box, which every
    descendant's centroid lies within, so the centroid rule is kept exactly.

    ``polyfill_array`` returns the same cells as a sorted numpy string array,
    which costs less than half the memory per cell for large fills.

    EXAMPLES::
        >>> from shapely import Polygon
        >>> coords = ((0., 0.), (0., 1.), (1., 1.), (1., 0.), (0., 0.))
        >>> polygon = Polygon(coords)
        >>> polyfill(polygon, res=1, plane=False)
        set()
        >>> polyfill(polygon, res=2, plane=False)
        set()
        >>> polyfill(polygon, res=3, plane=False)
        set()
        >>> polyfill(polygon, res=4, plane=False)
        {'Q3330'}
        >>> sorted(polyfill(polygon, res=5, plane=False))
        ['Q33303', 'Q33304', 'Q33305', 'Q33306', 'Q33307', 'Q33308', 'Q33330', 'Q33331', 'Q33332']
        >>> sorted(polyfill(polygon, res=6, plane=False))
        ['Q333033', 'Q333034', 'Q333035', 'Q333036', 'Q333037', 'Q333038', 'Q333043', 'Q333044', 'Q333045', 'Q333046', 'Q333047', 'Q333048', 'Q333053', 'Q333054', 'Q333056', 'Q333057', 'Q333060', 'Q333061', 'Q333062', 'Q333063', 'Q333064', 'Q333065', 'Q333066', 'Q333067', 'Q333068', 'Q333070', 'Q333071', 'Q333072', 'Q333073', 'Q333074', 'Q333075', 'Q333076', 'Q333077', 'Q333078', 'Q333080', 'Q333081', 'Q333083', 'Q333084', 'Q333086', 'Q333087', 'Q333300', 'Q333301', 'Q333302', 'Q333303', 'Q333304', 'Q333305', 'Q333306', 'Q333307', 'Q333308', 'Q333310', 'Q333311', 'Q333312', 'Q333313', 'Q333314', 'Q333315', 'Q333316', 'Q333317', 'Q333318', 'Q333320', 'Q333321', 'Q333323', 'Q333324', 'Q333326', 'Q333327', 'Q333330', 'Q333331', 'Q333332', 'Q333333', 'Q333334', 'Q333335', 'Q333340', 'Q333341', 'Q333342', 'Q333343', 'Q333344', 'Q333345', 'Q333350', 'Q333351', 'Q333353', 'Q333354']

        Four of those res-6 groups have all 9 siblings present (Q33306x, Q33307x,
        Q33330x, Q33331x), so compress=True collapses them one level to res-5 parents.
        The remaining partial groups stay at res-6, giving a mixed-resolution result:

        >>> sorted(c for c in polyfill(polygon, res=6, plane=False, compress=True) if len(c) < 7)
        ['Q33306', 'Q33307', 'Q33330', 'Q33331']

        At res=7 the same four cells are again fully covered (all 81 of each cell's
        res-7 grandchildren fall inside the polygon), so compaction runs two levels:
        res-7 siblings collapse to res-6 parents, those in turn collapse to res-5:

        >>> result7 = polyfill(polygon, res=7, plane=False, compress=True)
        >>> 'Q33306' in result7  # res-5 cell present despite res=7 request
        True
        >>> 'Q333060' not in result7  # intermediate res-6 cell absorbed
        True
        >>> 'Q3330600' not in result7  # original res-7 cell absorbed
        True
    """
    found = _polyfill_arrays(
        geometry, res, plane, verbose, dggs, containment, compress, min_res
    )
    if found is None:
        return None
    return set(_index_strings(dggs, *found).tolist())


def polyfill_array(
    geometry: Polygon | MultiPolygon,
    res: int,
    plane: bool = True,
    compress: bool = False,
    verbose: bool = False,
    dggs: RHEALPixDGGS = WGS84_003,
    containment: Literal["center", "full", "overlapping"] = "center",
    min_res: int = 0,
) -> np.ndarray | None:
    """
    ``polyfill`` returning the cell indices as a numpy string array instead of
    a set: the same cells for the same arguments, sorted (by base cell, then
    digit by digit, the order ``Cell`` objects sort in; with `compress` a
    coarse cell comes where its descendants would) and without duplicates,
    or ``None`` where ``polyfill`` returns ``None``. An empty array means no
    cell qualified.

    The array costs about ``4 * (res + 2)`` bytes per cell against roughly a
    hundred for a set entry, and feeds array-based callers such as
    ``RHEALPixDGGS.boundary_array``, ``RHEALPixDGGS.centroids`` or a data
    frame column directly. Like the rest of the index-string API, this is
    for DGGSs with single-character digits, ``N_side`` 2 or 3.

    EXAMPLES::
        >>> from shapely import Polygon
        >>> polygon = Polygon(((0., 0.), (0., 1.), (1., 1.), (1., 0.), (0., 0.)))
        >>> polyfill_array(polygon, res=5, plane=False).tolist()
        ['Q33303', 'Q33304', 'Q33305', 'Q33306', 'Q33307', 'Q33308', 'Q33330', 'Q33331', 'Q33332']
        >>> polyfill_array(polygon, res=2, plane=False).size
        0
        >>> polyfill_array(Polygon(), res=2, plane=False) is None
        True
    """
    found = _polyfill_arrays(
        geometry, res, plane, verbose, dggs, containment, compress, min_res
    )
    if found is None:
        return None
    return _index_strings(dggs, *found)


def _index_strings(
    dggs: RHEALPixDGGS, face: np.ndarray, digits: np.ndarray, resolution: np.ndarray
) -> np.ndarray:
    """
    The index strings of cells given as base cell codes, digit rows and
    per-cell resolutions, in input order, as one array wide enough for the
    deepest.
    """
    width = int(resolution.max()) + 1 if len(resolution) else 1
    out = np.empty(len(face), dtype=f"<U{width}")
    for r in np.unique(resolution):
        rows = resolution == r
        out[rows] = dggs._format_indices(face[rows], digits[rows], int(r))
    return out


def _polyfill_arrays(
    geometry: Polygon | MultiPolygon,
    res: int,
    plane: bool,
    verbose: bool,
    dggs: RHEALPixDGGS,
    containment: str,
    compress: bool,
    min_res: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """
    The cells ``polyfill`` selects, as ``(face, digits, resolution)`` arrays
    in the form of ``RHEALPixDGGS._parse_indices`` (digits zero-padded to
    `res` columns), sorted by base cell and then digit by digit and without
    duplicates; ``None`` for a malformed geometry.

    The fill descends from resolution `min_res`. At each level a cell wholly
    inside the geometry is accepted, a cell disjoint from it dropped, and the
    rest split into their children, until `res`, where the leaves are decided
    by `containment`. Accepted cells are kept as they are when `compress`,
    else expanded to their descendants at `res`. Every cell is held as one
    int64 key: its base cell followed by its digits as a base N**2 number,
    zero-padded to `res` digits, so that the keys' numeric order is the order
    of base cell then digit by digit and a cover, having no cell that is an
    ancestor of another, has no ties.
    """
    # Stop early if the geometry is malformed
    if _malformed_geometry(geometry):
        if verbose:
            message = is_valid_reason(geometry)
            if not message or message == "Valid Geometry":
                warn(POLYFILL_GEOMETRY_WARNING)
            else:
                warn(str.format("WARNING: {0}. Returning None.", message))

        return None

    if containment not in ("center", "full", "overlapping"):
        raise ValueError(
            f"containment must be 'center', 'full' or 'overlapping', not {containment!r}"
        )
    if not 0 <= min_res <= res:
        raise ValueError(f"min_res must be between 0 and res={res}, not {min_res!r}")
    N2 = dggs.N_side**2
    if 6 * N2**res >= 2**62:
        raise ValueError(f"resolution {res} is too deep for N_side={dggs.N_side}")

    # Extract list of polygons from geometry: Polygon needs to be wrapped in
    # one, MultiPolygon has it stashed in a property
    if isinstance(geometry, Polygon):
        geoms = [geometry]
    else:
        geoms = list(geometry.geoms)

    chunk = rhp_dggs._LATTICE_CHUNK
    keys_out: list[np.ndarray] = []
    res_out: list[np.ndarray] = []

    def emit(keys: np.ndarray, level: int) -> None:
        keys_out.append(keys)
        res_out.append(np.full(len(keys), level, dtype=np.int8))

    # Above the leaves a cell is accepted or rejected only when it is clear
    # of the geometry's boundary by this margin (about a centimetre), so
    # that a cell edge coinciding with a polygon edge, whose coordinates
    # differ by an ulp between a parent and its children, is left to the
    # exact test at the leaves.
    margin = 1e-9 * (dggs.cell_width(0) if plane else dggs.ellipsoid.pi() / 2)
    for geom in geoms:
        shapely.prepare(geom)
        minx, miny, maxx, maxy = geom.bounds
        boxes = dggs._candidate_boxes((minx, maxy), (maxx, miny), plane)
        # The frontier at min_res: the lattice cells of the bounding box.
        frontier = [
            (face, digits[:, :min_res])
            for face, digits in dggs._lattice_cells(boxes, min_res)
        ]
        for level in range(min_res, res + 1):
            if not frontier:
                break
            face = np.concatenate([f for f, _ in frontier])
            digits = np.concatenate([d for _, d in frontier])
            frontier = []
            scale = N2 ** (res - level)  # descendants at res per cell here
            # A cell accepted for the flat fill is expanded into `scale`
            # keys; below a chunk per cell that is done here, above it the
            # cell is split instead, so no expansion exceeds a chunk.
            expand_here = compress or scale <= chunk
            for start in range(0, len(face), chunk):
                f = face[start : start + chunk]
                d = digits[start : start + chunk]
                resolution = np.full(len(f), level)
                if level == res:
                    if containment == "center":
                        keep = _cells_with_centroid_inside(
                            geom, dggs, plane, f, d, resolution
                        )
                    else:
                        keep = _cells_within_or_overlapping(
                            geom, dggs, plane, f, d, resolution, containment == "full"
                        )
                    if keep.any():
                        emit(_keys(f[keep], d[keep], level, res, N2), res)
                    continue
                accept, reject = _classify(
                    geom, dggs, plane, f, d, resolution, containment, margin
                )
                if not expand_here:
                    accept = np.zeros_like(accept)
                if accept.any():
                    base = _keys(f[accept], d[accept], level, res, N2)
                    if compress:
                        emit(base, level)
                    else:
                        # Expand in groups whose total stays within a chunk.
                        per_group = max(1, chunk // scale)
                        for g0 in range(0, len(base), per_group):
                            group = base[g0 : g0 + per_group]
                            emit(
                                (
                                    group[:, None] + np.arange(scale, dtype=np.int64)
                                ).ravel(),
                                res,
                            )
                split = ~accept & ~reject
                if split.any():
                    frontier.append(_children(f[split], d[split], N2))

    width = max(res, 1)
    if not keys_out:
        empty = np.empty(0, dtype=np.int64)
        return empty, np.empty((0, width), dtype=np.int16), empty.astype(np.int8)
    keys = np.concatenate(keys_out)
    levels = np.concatenate(res_out)
    if compress:
        keys, levels = _merge_sibling_groups(keys, levels, res, min_res, N2)
    else:
        keys = np.unique(keys)  # a cell met by two parts appears twice
        levels = np.full(len(keys), res, dtype=np.int8)
    return _decode(keys, levels, res, N2)


def _keys(
    face: np.ndarray, digits: np.ndarray, level: int, res: int, N2: int
) -> np.ndarray:
    """The res-padded int64 keys of resolution-`level` cells."""
    weights = N2 ** np.arange(res - 1, res - 1 - level, -1, dtype=np.int64)
    return face.astype(np.int64) * N2**res + (
        digits[:, :level].astype(np.int64) * weights
    ).sum(axis=1)


def _children(
    face: np.ndarray, digits: np.ndarray, N2: int
) -> tuple[np.ndarray, np.ndarray]:
    """The N2 children of each cell: a digit column appended, cells repeated."""
    n = len(face)
    child = np.tile(np.arange(N2, dtype=digits.dtype), n)
    face = np.repeat(face, N2)
    digits = np.column_stack([np.repeat(digits, N2, axis=0), child])
    return face, digits


def _classify(
    geom: Polygon,
    dggs: RHEALPixDGGS,
    plane: bool,
    face: np.ndarray,
    digits: np.ndarray,
    resolution: np.ndarray,
    containment: str,
    margin: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    For cells above the leaf resolution: `accept` marks cells wholly inside
    `geom` and clear of its boundary by `margin`, so that every descendant
    qualifies under `containment`, and `reject` cells farther than `margin`
    from it, so that none does. For ``center`` the test is on the cell's
    longitude-latitude bounding box (planar square in the plane), which
    contains every descendant's centroid; for ``full`` and ``overlapping`` on
    the cell's boundary polygon, as the leaf test is. Cells whose box cannot
    be trusted (spanning more than half a turn) are neither.
    """
    if containment == "center":
        if plane:
            x, y, width, _ = dggs._index_geometry(face, digits, resolution)
            boxes = shapely.box(x, y - width, x + width, y)
            trusted = np.ones(len(face), dtype=bool)
        else:
            corners = dggs._boundary_array(face, digits, resolution, 2, False)
            lo = corners.min(axis=1)
            hi = corners.max(axis=1)
            trusted = hi[:, 0] - lo[:, 0] <= dggs.ellipsoid.pi()
            boxes = shapely.box(lo[:, 0], lo[:, 1], hi[:, 0], hi[:, 1])
        boundary = geom.boundary
        shapely.prepare(boundary)
        inside = np.asarray(shapely.contains(geom, boxes))
        near = np.asarray(shapely.dwithin(boundary, boxes, margin))
        close = np.asarray(shapely.dwithin(geom, boxes, margin))
        return inside & ~near & trusted, ~close & trusted
    accept = _cells_within_or_overlapping(
        geom, dggs, plane, face, digits, resolution, True, margin
    )
    meets = _cells_within_or_overlapping(
        geom, dggs, plane, face, digits, resolution, False, margin
    )
    return accept, ~meets


def _merge_sibling_groups(
    keys: np.ndarray, levels: np.ndarray, res: int, min_res: int, N2: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Deduplicate a cover and merge every complete group of N2 siblings into
    its parent, repeatedly, never above `min_res`; return the keys sorted.
    """
    pairs = np.unique(np.column_stack([keys, levels.astype(np.int64)]), axis=0)
    keys, levels = pairs[:, 0], pairs[:, 1]
    for level in range(res, min_res, -1):
        at = levels == level
        if not at.any():
            continue
        stride = N2 ** (res - level + 1)
        parents = (keys[at] // stride) * stride
        unique_parents, counts = np.unique(parents, return_counts=True)
        complete = counts == N2
        if not complete.any():
            continue
        merged = np.isin(parents, unique_parents[complete])
        keep = np.ones(len(keys), dtype=bool)
        keep[np.flatnonzero(at)[merged]] = False
        keys = np.concatenate([keys[keep], unique_parents[complete]])
        levels = np.concatenate(
            [levels[keep], np.full(int(complete.sum()), level - 1, dtype=np.int64)]
        )
    order = np.argsort(keys, kind="stable")
    return keys[order], levels[order].astype(np.int8)


def _decode(
    keys: np.ndarray, levels: np.ndarray, res: int, N2: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Base cell codes and zero-padded int16 digit rows from res-padded keys."""
    width = max(res, 1)
    face = keys // N2**res
    digits = np.zeros((len(keys), width), dtype=np.int16)
    # One digit column at a time, so that no intermediate is wider than the
    # keys themselves.
    for k in range(res):
        digits[:, k] = (keys // N2 ** (res - 1 - k)) % N2
    return face, digits, levels.astype(np.int8)


def _cells_with_centroid_inside(
    geom: Polygon,
    dggs: RHEALPixDGGS,
    plane: bool,
    face: np.ndarray,
    digits: np.ndarray,
    resolution: np.ndarray,
) -> np.ndarray:
    """
    For each cell, given as the arrays ``RHEALPixDGGS._parse_indices``
    produces, whether its centroid lies inside `geom`.

    In the plane the centroid is the nucleus. On the ellipsoid the centroid
    is the mean of longitude and latitude over the cell, so it lies within
    the cell's longitude-latitude bounding box: a cell whose box is wholly
    inside `geom` is in, one whose box is disjoint from it is out, and only
    the cells whose box meets the boundary need their centroid computed.
    Boxes spanning more than half a turn of longitude (cells straddling the
    antimeridian) are not trusted and always get the exact test.
    """
    if plane:
        x, y, width, _ = dggs._index_geometry(face, digits, resolution)
        return shapely.contains_xy(geom, x + width / 2, y - width / 2)
    corners = dggs._boundary_array(face, digits, resolution, 2, False)
    lo = corners.min(axis=1)
    hi = corners.max(axis=1)
    trusted = hi[:, 0] - lo[:, 0] <= dggs.ellipsoid.pi()
    keep = np.zeros(len(face), dtype=bool)
    undecided = ~trusted
    if trusted.any():
        boxes = shapely.box(
            lo[trusted, 0], lo[trusted, 1], hi[trusted, 0], hi[trusted, 1]
        )
        inside = shapely.contains(geom, boxes)
        outside = shapely.disjoint(geom, boxes)
        keep[np.flatnonzero(trusted)[inside]] = True
        undecided[np.flatnonzero(trusted)[~inside & ~outside]] = True
    if undecided.any():
        centroids = dggs._centroids(
            face[undecided], digits[undecided], resolution[undecided], False
        )
        keep[undecided] = shapely.contains_xy(geom, centroids[:, 0], centroids[:, 1])
    return keep


def _cells_within_or_overlapping(
    geom: Polygon,
    dggs: RHEALPixDGGS,
    plane: bool,
    face: np.ndarray,
    digits: np.ndarray,
    resolution: np.ndarray,
    full: bool,
    margin: float = 0.0,
) -> np.ndarray:
    """
    For each cell, given as the arrays ``RHEALPixDGGS._parse_indices``
    produces, whether it lies fully within `geom` (`full` = True) or meets
    it at all (`full` = False), tested on the cells' boundary polygons:
    exact in the plane and for equatorial cells, a 6-point-per-edge
    approximation of the curved edges of polar cells. Cap cells and cells
    crossing the antimeridian are handled apart, as their rings are not
    polygons in longitude-latitude coordinates.

    With a positive `margin` the tests are those the hierarchical descent
    needs: "fully within" also requires the cell to be clear of `geom`'s
    boundary by `margin`, and "meets" becomes "within `margin` of `geom`";
    cap cells then count as neither fully within nor clear of it.
    """
    polar = (face == 0) | (face == 5)
    n = 6 if (not plane and polar.any()) else 2
    rings = dggs._boundary_array(face, digits, resolution, n, plane)
    boundary = geom.boundary if margin > 0 else None
    if boundary is not None:
        shapely.prepare(boundary)

    def within(polys: Any) -> np.ndarray:
        inside = np.asarray(shapely.contains(geom, polys))
        if boundary is None:
            return inside
        return inside & ~np.asarray(shapely.dwithin(boundary, polys, margin))

    def meets(g: Polygon, polys: Any) -> np.ndarray:
        if boundary is None:
            return np.asarray(shapely.intersects(g, polys))
        return np.asarray(shapely.dwithin(g, polys, margin))

    if plane:
        polys = shapely.polygons(rings)
        return within(polys) if full else meets(geom, polys)

    keep = np.zeros(len(face), dtype=bool)
    cap = dggs._shape_code(face, digits, resolution) == 1
    half = dggs.ellipsoid.pi()
    # A ring whose longitudes jump between -180 and 180 is unwrapped
    # eastwards past 180. A cell whose east edge lies exactly on the
    # antimeridian then reads as an ordinary polygon within [-180, 180];
    # one that genuinely crosses it reaches beyond 180.
    lon = rings[:, :, 0]
    wraps = ~cap & (lon.max(axis=1) - lon.min(axis=1) > half)
    rings[wraps, :, 0] = np.where(lon[wraps] < 0, lon[wraps] + 2 * half, lon[wraps])
    crosses = wraps & (rings[:, :, 0].max(axis=1) > half * (1 + 1e-12))
    plain = ~cap & ~crosses

    if plain.any():
        polys = shapely.polygons(rings[plain])
        keep[plain] = within(polys) if full else meets(geom, polys)
    if crosses.any() and not full:
        # A crossing cell can only be fully inside a geometry that itself
        # crosses, which callers must have split, so it is never `full`; it
        # overlaps the geometry if its unwrapped ring meets the geometry or
        # the geometry shifted one turn east.
        polys = shapely.polygons(rings[crosses])
        shifted = translate(geom, xoff=2 * half)
        keep[crosses] = meets(geom, polys) | meets(shifted, polys)
    if cap.any():
        _, miny, _, maxy = geom.bounds
        for k in np.flatnonzero(cap):
            if boundary is not None:
                # Neither clear of the boundary nor clear of the geometry:
                # a cap is only ever decided at the leaves.
                keep[k] = not full
                continue
            lat_c = rings[k, 0, 1]
            north = face[k] == 0
            band = (
                box(-half, lat_c, half, half / 2)
                if north
                else box(-half, -half / 2, half, lat_c)
            )
            if full:
                keep[k] = geom.covers(band)
            else:
                keep[k] = maxy >= lat_c if north else miny <= lat_c
    return keep


def linetrace(
    geometry: LineString | MultiLineString,
    res: int,
    plane: bool = True,
    verbose: bool = False,
    dggs: RHEALPixDGGS = WGS84_003,
    wrap_antimeridian: bool = False,
) -> list[str] | None:
    """
    Returns the list of cell indices touched by a shapely linestring or multilinestring
    at the requested resolution. Removes internal sequences of duplicate cells before
    returning result.

    Returns None if the geom_type field in the input geometry is anything other than
    'LineString' or 'MultiLineString'.

    Returns None if the geometry is empty, or if it has no length.

    Returns None if no cells match the geometry for some reason.

    Returns None if the geometry is invalid in other ways, e.g. if a linestring contains
    self intersecting segments.

    A linestring's segments are straight in the input coordinate space
    (shapely geometries are planar): planar rHEALPix coordinates if
    `plane` = True, longitude-latitude coordinates otherwise -- so
    longitude-latitude segments are plate carree lines, not geodesics,
    and by default do not wrap around the antimeridian (a segment from
    longitude 179 to -179 runs the long way around, through longitude 0
    -- the literal planar reading). Splitting inputs at the antimeridian,
    as GeoJSON's RFC 7946 prescribes for data producers, avoids the
    ambiguity entirely; alternatively pass `wrap_antimeridian` = True to
    trace segments spanning more than half a turn of longitude the short
    way, across the antimeridian, which traces the same cells as
    splitting them there. To trace a geodesic, densify it into short
    segments first (e.g. with pyproj.Geod). See
    RHEALPixDGGS.cells_from_line, which traces each segment exactly, for
    the details.

    EXAMPLES::

        >>> from shapely import LineString
        >>> line = LineString([[-176.260506, -43.738058], [-176.258807, -43.738379]])
        >>> linetrace(line, res=9, plane=False)
        ['S001450634', 'S001450635']
    """
    # Stop early if the line geometry is malformed
    if _malformed_lines(geometry):
        if verbose:
            message = is_valid_reason(geometry)
            if not message or message == "Valid Geometry":
                warn(LINETRACE_GEOMETRY_WARNING)
            else:
                warn(str.format("WARNING: {0}. Returning None.", message))

        return None

    # Extract list of linestrings from geometry: LineString needs to be wrapped in
    # one, MultiLineString has it stashed in a property
    if isinstance(geometry, LineString):
        lines = [geometry]
    else:
        lines = list(geometry.geoms)

    cells: list[str] = []
    for linestring in lines:
        # Extract coordinate pairs along the line segments
        coords = zip(linestring.coords, linestring.coords[1:])

        # Walk along line segments
        while (vertex_pair := next(coords, None)) is not None:
            # Extract vertex pair defining line segment in (lng, lat) order
            i, j = vertex_pair

            # Convert line segment to cell ids
            line_cells = dggs.cells_from_line(
                res,
                (i[0], i[1]),
                (j[0], j[1]),
                plane,
                wrap_antimeridian=wrap_antimeridian,
            )

            # Convert cells to string ids and add to collection
            if line_cells:
                cells = cells + [str(cell) for cell in line_cells]

        # Remove duplicates along sequence
        cells = _remove_sequential_duplicates(cells)

    return cells


# ======== Helper functions ======== #


# Fixed clockwise order for a cell's up-to-8 edge/corner neighbors, paired
# with whether that direction is an edge step (neighbor()) or a corner-only
# step (diagonal_neighbor()).
_RING_STEP_DIRECTIONS = [
    ("up", False),
    ("up_right", True),
    ("right", False),
    ("down_right", True),
    ("down", False),
    ("down_left", True),
    ("left", False),
    ("up_left", True),
]


def _ring_step_neighbors(cell: Cell) -> list[Cell]:
    """
    Return this cell's up to 8 distinct edge- and corner-adjacent
    neighbors, in a fixed clockwise order starting from "up". Skips any
    direction with no neighbor: a genuine cube corner, where exactly 3
    cells meet rather than 4 (see Cell.diagonal_neighbor()).
    """
    neighbors = []
    for direction, diagonal in _RING_STEP_DIRECTIONS:
        neighbor = (
            cell.diagonal_neighbor(direction)
            if diagonal
            else cell.neighbor(direction, plane=True)
        )
        if neighbor is not None:
            neighbors.append(neighbor)
    return neighbors


def _rings_up_to(center: Cell, k: int) -> list[list[Cell]]:
    """
    Return `[ring_0, ring_1, ..., ring_k]`, where `ring_0 == [center]` and
    `ring_i` (`i` >= 1) is the list of cells at graph distance exactly `i`
    from `center` -- the shortest number of edge/corner hops, via
    `_ring_step_neighbors()` -- computed by growing the ring outward one
    shell at a time.

    Each ring's cells are visited in a stable order: each cell in the
    previous ring is expanded via `_ring_step_neighbors()`'s fixed
    clockwise order, in previous-ring order, keeping the first occurrence
    of any cell reachable from more than one previous-ring cell. Cells
    already seen in an earlier ring (tracked cumulatively, not just the
    immediately preceding one, so this is a proper breadth-first search)
    are skipped, so a ring near a genuine cube corner naturally comes out
    with fewer cells than the usual 8*i for an interior ring, rather than
    a duplicate or incorrect cell.
    """
    rings = [[center]]
    seen_suids = {center.suid}
    for _ in range(k):
        next_ring = []
        seen_this_ring = set()
        for cell in rings[-1]:
            for neighbor in _ring_step_neighbors(cell):
                if neighbor.suid in seen_suids or neighbor.suid in seen_this_ring:
                    continue
                seen_this_ring.add(neighbor.suid)
                next_ring.append(neighbor)
        seen_suids.update(seen_this_ring)
        rings.append(next_ring)
    return rings


def _malformed_geometry(geometry: Polygon | MultiPolygon) -> bool:
    # Geometry has to have things in it
    if geometry is None or geometry.is_empty:
        return True

    # Geometry needs to be of the correct type
    if geometry.geom_type != "Polygon" and geometry.geom_type != "MultiPolygon":
        return True

    # This catches e.g. self intersecting hulls and holes, or overlapping polygons
    if not geometry.is_valid:
        return True

    # Geometry has to have an area, i.e. not be collapsed to a line
    return geometry.area == 0


def _malformed_lines(lines: LineString | MultiLineString) -> bool:
    # There have to be lines
    if lines is None or lines.is_empty:
        return True

    # Lines need to be of the correct type
    if lines.geom_type != "LineString" and lines.geom_type != "MultiLineString":
        return True

    if not lines.is_valid:
        return True

    # Lines need to have a length, i.e. not be collapsed into points
    return lines.length == 0


def _remove_sequential_duplicates(cells: list[str]) -> list[str]:
    """
    EXAMPLES::

        >>> _remove_sequential_duplicates(['S','S','P','S'])
        ['S', 'P', 'S']
    """
    if not cells:
        return []

    trimmed_cells = []
    prev = None
    for cell in cells:
        if cell != prev:
            trimmed_cells.append(cell)
            prev = cell

    return trimmed_cells
