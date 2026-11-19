from collections.abc import Iterable
from itertools import compress, pairwise
from typing import Any, Literal, cast

import numpy as np
import shapely
from shapely import is_valid_reason
from shapely.affinity import translate
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon, box
from shapely.geometry.base import BaseGeometry

import rhealpixdggs.dggs as rhp_dggs
from rhealpixdggs.cell import Cell
from rhealpixdggs.dggs import WGS84_003, LineSemantics, RHEALPixDGGS
from rhealpixdggs.identifiers import _format_index, _index_width, _parse_index


def cells_from_polygon(
    geometry: Polygon | MultiPolygon,
    resolution: int,
    plane: bool = True,
    containment: Literal["center", "full", "overlapping"] = "center",
    compress: bool = False,
    min_res: int = 0,
    rdggs: RHEALPixDGGS = WGS84_003,
) -> np.ndarray:
    """
    Return the index strings of the cells of the Shapely polygon
    `geometry`, read as cells of `rdggs`, as a sorted numpy string array
    without duplicates.
    The polygon's coordinates are planar `x` and `y` if `plane` =
    True, else longitude and latitude. Its edges are straight in
    those coordinates and do not wrap around the antimeridian: split
    a region that straddles it into a MultiPolygon first.

    `containment` chooses the cells at resolution `resolution`:

    - ``"center"``: the cell's centroid is inside the geometry.
    - ``"full"``: the whole cell lies within the geometry (their
      boundaries may touch).
    - ``"overlapping"``: any part of the cell meets the geometry, so
      the result covers it.

    ``full`` and ``overlapping`` test the cell's boundary polygon,
    which is exact in the plane and for equatorial cells; the curved
    edges of polar cells are represented by 6 points each. A cap
    cell is fully inside only a geometry that covers the whole band
    from the cap's boundary to the pole.

    The fill descends from resolution `min_res`, keeping whole each
    cell wholly inside the geometry, dropping each cell clear of it
    and splitting only the cells its boundary passes through, so the
    work grows with the length of the boundary and not the area.
    With `compress` = False the kept cells are expanded to
    `resolution`. With `compress` = True the result is a
    mixed-resolution cover with no cell coarser than `min_res` and
    no complete group of siblings left unmerged; a coarse cell sorts
    where its descendants would.

    Return an empty array for an empty geometry or one with no area.
    Raise TypeError unless `geometry` is a Polygon or MultiPolygon,
    and ValueError if it is invalid, if `containment` is not one of
    the three rules or if `min_res` is not between 0 and
    `resolution`.

    EXAMPLES::

        >>> from shapely import Polygon
        >>> polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
        >>> cells_from_polygon(polygon, 4, plane=False).tolist()
        ['Q3330']
        >>> len(cells_from_polygon(polygon, 6, plane=False))
        80
        >>> len(cells_from_polygon(polygon, 6, plane=False, compress=True))
        48

    """
    rdggs._require_index_strings()
    if containment not in ("center", "full", "overlapping"):
        raise ValueError(
            f"containment must be 'center', 'full' or 'overlapping', not {containment!r}"
        )
    if not 0 <= min_res <= resolution:
        raise ValueError(
            f"min_res must be between 0 and res={resolution}, not {min_res!r}"
        )
    if not _usable(geometry, (Polygon, MultiPolygon)) or geometry.area == 0:
        return np.empty(0, dtype=f"<U{resolution + 1}")
    return _index_strings(
        rdggs,
        *_fill_arrays(
            geometry, resolution, plane, rdggs, containment, compress, min_res
        ),
    )


def cells_from_linestring(
    geometry: LineString | MultiLineString,
    resolution: int,
    plane: bool = True,
    wrap_antimeridian: bool = False,
    line: LineSemantics | None = None,
    rdggs: RHEALPixDGGS = WGS84_003,
) -> list[str]:
    """
    Return the index strings of the resolution `resolution` cells that
    the Shapely line `geometry` passes through, in order along it, a
    cell appearing again each time the line re-enters it, read as cells
    of `rdggs`. Each segment
    is traced by ``RHEALPixDGGS.cells_from_line``, which describes `plane`,
    `wrap_antimeridian` and `line`: by default a segment is straight
    in the input coordinates, planar `x` and `y` if `plane` = True,
    else longitude and latitude, where it does not wrap around the
    antimeridian.

    Return an empty list for an empty geometry. Raise TypeError unless
    `geometry` is a LineString or MultiLineString, and ValueError if
    it is invalid.

    EXAMPLES::

        >>> from shapely import LineString
        >>> line = LineString([(-176.260506, -43.738058), (-176.258807, -43.738379)])
        >>> cells_from_linestring(line, 9, plane=False)
        ['S001450634', 'S001450635']

    """
    rdggs._require_index_strings()
    if not _usable(geometry, (LineString, MultiLineString)):
        return []
    lines = [geometry] if isinstance(geometry, LineString) else list(geometry.geoms)
    cells: list[str] = []
    for linestring in lines:
        for start, end in pairwise(linestring.coords):
            cells += [
                str(cell)
                for cell in rdggs.cells_from_line(
                    resolution,
                    (start[0], start[1]),
                    (end[0], end[1]),
                    plane,
                    wrap_antimeridian=wrap_antimeridian,
                    line=line,
                )
            ]
        cells = _remove_sequential_duplicates(cells)
    return cells


def get_finest_containing_cell(
    polygon: Polygon | MultiPolygon, rdggs: RHEALPixDGGS = WGS84_003
) -> Cell | None:
    """
    Finds the finest DGGS Cell containing a given cartesian polygon
    """
    # Cap cells span all longitudes at the pole, so their boundary cannot be
    # represented as a simple polygon in geographic coordinates; containment
    # is tested in the rHEALPix plane instead. Project the polygon there
    # once, each part's ring in one array call. (MultiPolygon has no
    # .exterior; check each part individually.)
    parts = list(polygon.geoms) if hasattr(polygon, "geoms") else [polygon]
    planar_parts = []
    for part in parts:
        coords = np.asarray(part.exterior.coords)
        x, y = rdggs.rhealpix(coords[:, 0], coords[:, 1])
        planar_parts.append(Polygon(np.column_stack([x, y])))

    def _cell_contains(cell: Cell, poly: Polygon | MultiPolygon) -> bool:
        if cell.ellipsoidal_shape == "cap":
            plane_cell = Polygon(cell.vertices(plane=True))
            return all(plane_cell.contains(part) for part in planar_parts)
        return Polygon(cell.vertices(plane=False)).contains(poly)

    def _get_finest_cell(
        polygon: Polygon | MultiPolygon,
        suid: tuple[str | int, ...],
        rdggs: RHEALPixDGGS = rdggs,
    ) -> Cell | None:
        parent_cell = Cell(rdggs=rdggs, suid=suid)
        # get the children cells and polygons for these cells
        children_cells = list(parent_cell.subcells())
        # function and truth list for multipolygon / polygon (polygon) contained within multipolygon / polygon (cell)
        truth = [_cell_contains(cell, polygon) for cell in children_cells]
        # if we get something back, check the next level lower
        returned_cells = list(compress(children_cells, truth))
        if returned_cells:
            finest = _get_finest_cell(polygon, returned_cells[0].suid, rdggs)
        else:
            if _cell_contains(parent_cell, polygon):
                finest = parent_cell
            else:
                finest = None
        return finest

    for suid in [tuple(x) for x in ["N", "O", "P", "Q", "R", "S"]]:
        finest = _get_finest_cell(polygon, suid, rdggs)
        if finest is not None:
            return finest
    return None


def _parent_index(index: str, N_side: int) -> str:
    """
    The index string of the parent of the cell with index string `index` in
    a grid with `N_side`, or the empty string for a resolution 0 cell.
    """
    suid = _parse_index(index, N_side)
    if suid is None:
        raise ValueError(f"invalid cell index {index!r} for N_side={N_side}")
    return _format_index(suid[:-1], N_side) if len(suid) > 1 else ""


def compress_order_cells(
    cells: list[str], recursive: bool = False, N_side: int = 3
) -> list[str]:
    """
    Compress and sort a collection of cells.

    Each complete group of ``N_side ** 2`` siblings is replaced by their
    parent cell. ``N_side`` defaults to 3 (the default for `RHEALPixDGGS`
    and most presets, e.g. ``WGS84_003``); pass the ``N_side`` of the DGGS
    the cells actually belong to if it's not 3 (e.g. 2 for ``WGS84_002``),
    otherwise a complete sibling group will never be recognised as complete
    and nothing will compress.
    By default only one level of compression is performed. Pass
    ``recursive=True`` to repeat until no further compression is possible.
    The result is always deduplicated and sorted alphanumerically.
    """
    import re

    _index_width(N_side)

    def alphanum_sort(lst: list[str]) -> list[str]:
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
        return sorted(lst, key=alphanum_key)

    def _compress_once(cell_set: set[str]) -> set[str]:
        upper_cells: dict[str, list[str]] = {}
        for cell in cell_set:
            upper_cells.setdefault(_parent_index(cell, N_side), []).append(cell)
        result: set[str] = set()
        for k, v in upper_cells.items():
            if len(v) == N_side**2:
                result.add(k)
            else:
                result.update(v)
        return result

    cell_set = set(cells)
    if recursive:
        while True:
            next_set = _compress_once(cell_set)
            if next_set == cell_set:
                break
            cell_set = next_set
    else:
        cell_set = _compress_once(cell_set)
    return alphanum_sort(list(cell_set))


def compact_cells(cells: Iterable[str], N_side: int = 3) -> set[str]:
    """
    Compact a set of rHEALPix DGGS cells by repeatedly merging complete
    groups of ``N_side ** 2`` siblings into their parent until no further
    merging is possible. ``N_side`` defaults to 3 (the default for
    `RHEALPixDGGS` and most presets, e.g. ``WGS84_003``); pass the
    ``N_side`` of the DGGS the cells actually belong to if it's not 3
    (e.g. 2 for ``WGS84_002``), otherwise a complete sibling group will
    never be recognised as complete and nothing will compact.

    All input cells must be at the same resolution (equal string length).
    Raises ValueError if the input contains cells at mixed resolutions.
    Returns a set; if the caller needs ordering, sort the result separately.
    """
    _index_width(N_side)
    cell_set = set(cells)
    if not cell_set:
        return cell_set

    lengths = {len(c) for c in cell_set}
    if len(lengths) > 1:
        raise ValueError(
            f"All input cells must be at the same resolution; "
            f"got cells with string lengths {sorted(lengths)}"
        )

    while True:
        upper: dict[str, list[str]] = {}
        for cell in cell_set:
            upper.setdefault(_parent_index(cell, N_side), []).append(cell)

        next_set: set[str] = set()
        for parent, children in upper.items():
            if len(children) == N_side**2:
                next_set.add(parent)
            else:
                next_set.update(children)

        if next_set == cell_set:
            return cell_set
        cell_set = next_set


def _usable(geometry: Any, kinds: tuple[type, ...]) -> bool:
    """
    Raise TypeError unless `geometry` is one of `kinds` and ValueError if
    it is invalid; return whether it has anything to fill or trace.
    """
    if not isinstance(geometry, kinds):
        names = " or ".join(kind.__name__ for kind in kinds)
        raise TypeError(f"geometry must be a {names}, not {type(geometry).__name__}")
    shape = cast(BaseGeometry, geometry)
    if not shape.is_valid:
        raise ValueError(f"invalid geometry: {is_valid_reason(shape)}")
    return not shape.is_empty


def _remove_sequential_duplicates(cells: list[str]) -> list[str]:
    """
    EXAMPLES::

        >>> _remove_sequential_duplicates(['S','S','P','S'])
        ['S', 'P', 'S']
    """
    trimmed: list[str] = []
    for cell in cells:
        if not trimmed or cell != trimmed[-1]:
            trimmed.append(cell)
    return trimmed


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


def _fill_arrays(
    geometry: Polygon | MultiPolygon,
    res: int,
    plane: bool,
    dggs: RHEALPixDGGS,
    containment: str,
    compress: bool,
    min_res: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    The cells ``cells_from_polygon`` selects, as ``(face, digits, resolution)`` arrays
    in the form of ``RHEALPixDGGS._parse_indices`` (digits zero-padded to
    `res` columns), sorted by base cell and then digit by digit and without
    duplicates.

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
    dggs._require_index_strings()
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
    the cell's boundary polygon where that is exact (the plane, equatorial
    cells) and otherwise on its longitude-latitude bounding box, which
    contains the cell, so that no decision here contradicts the leaf test.
    Cells whose box cannot be trusted (spanning more than half a turn) are
    neither.
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
        reject = ~close & trusted
        if not trusted.all():
            # A cell crossing the antimeridian has a box once its western
            # longitudes are moved a turn east, to be tested against the
            # geometry there and a turn east of there.
            rows = np.flatnonzero(~trusted)
            turn = 2 * dggs.ellipsoid.pi()
            lon = corners[rows, :, 0]
            lon = np.where(lon < 0, lon + turn, lon)
            west, east = lon.min(axis=1), lon.max(axis=1)
            unwrapped = shapely.box(west, lo[rows, 1], east, hi[rows, 1])
            shifted = translate(geom, xoff=turn)
            clear = ~np.asarray(shapely.dwithin(geom, unwrapped, margin)) & ~np.asarray(
                shapely.dwithin(shifted, unwrapped, margin)
            )
            reject[rows] = clear & (east - west <= turn / 2)
        return inside & ~near & trusted, reject
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
    cap cells then count as neither fully within nor clear of it, and polar
    cells are judged by their longitude-latitude bounding box rather than
    their 6-point polygon. The box contains the whole cell, as a polar
    cell's curved edges are monotone in longitude and latitude, whereas the
    polygon of a coarse cell can miss the true edge by a degree, so a
    decision on it could drop or accept descendants the leaf test would
    decide the other way.
    """
    polar = (face == 0) | (face == 5)
    # Always 6 points per edge on the ellipsoid, not only when the batch
    # has polar cells: the lattice the points are projected from has pitch
    # width / (n - 1), so a cell's coordinates depend on n in the last bit,
    # and a decision about a cell must not depend on which cells it is
    # tested with.
    n = 2 if plane else 6
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

    def cells(mask: np.ndarray) -> Any:
        polys = shapely.polygons(rings[mask])
        if boundary is not None and (polar & mask).any():
            # Above the leaves, a polar cell is its bounding box.
            box_rows = polar[mask]
            polys[box_rows] = shapely.envelope(polys[box_rows])
        return polys

    if plain.any():
        polys = cells(plain)
        keep[plain] = within(polys) if full else meets(geom, polys)
    if crosses.any() and not full:
        # A crossing cell can only be fully inside a geometry that itself
        # crosses, which callers must have split, so it is never `full`; it
        # overlaps the geometry if its unwrapped ring meets the geometry or
        # the geometry shifted one turn east.
        polys = cells(crosses)
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
