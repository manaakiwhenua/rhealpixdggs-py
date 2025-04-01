from typing import Literal, Union
from warnings import warn
from shapely import Point, Polygon, MultiPolygon, is_valid_reason

from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs.cell import Cell
from rhealpixdggs.conversion import compress_order_cells

# ======== Messages and constants ======== #


# Pre-defined DGGS with WGS84 ellipsoid, coordinates in degrees and n == 3 to subdivide cell sides
from rhealpixdggs.dggs import WGS84_003

# List of resolution 0 cell addresses (i.e. cube faces)
from rhealpixdggs.cell import CELLS0

PARENT_RESOLUTION_WARNING = "WARNING: You requested a parent resolution that is higher than the cell resolution. Returning the cell address itself."
CHILD_RESOLUTION_WARNING = "WARNING: You requested a child resolution that is lower than the cell resolution. Returning the cell address itself."
CELL_CENTRE_WARNING = "WARNING: You requested a centre cell for a DGGS that has an even number of cells on a side. Returning None."
CELL_RING_WARNING = "WARNING: Implementation of cell rings is incomplete. Requesting a {0} ring that involves more than two resolution 0 cube faces will return unexpected results."
POLYFILL_GEOMETRY_WARNING = "WARNING: Empty or missing geometry, unsupported geometry type (not Polygon or MultiPolygon), or geometry with no area. Returning None."

# ======== Main API ======== #


def geo_to_rhp(
    lat: float,
    lng: float,
    resolution: int,
    plane: bool = True,
    dggs: RHEALPixDGGS = WGS84_003,
) -> str:
    """
    Turn a latitute and longitude (in degrees) into an rHEALPix cell address at
    the requested resolution.

    Mostly passes through the parameters to the function turning coordinate points
    into cells, but converts the address tuple from the resulting cell into a
    string.

    Returns None if no cell matching the coordinates is found.
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
) -> tuple[float, float]:
    """
    Look up the centroid (in degrees) of the cell identified by rhpindex.

    Returns None if the cell index is invalid.

    If geojson is requested as the output format:
        - Will return a (longitude, latitude) coordinate pair.

    if geojson is NOT requested as the output format:
        - Will return a (latitude, longitude) coordinate pair in order to be consistent with
          h3 coordinate ordering.
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

    return centroid


def rhp_to_parent(
    rhpindex: str, res: int = None, verbose: bool = True, dggs: RHEALPixDGGS = WGS84_003
) -> str:
    """
    Returns parent of rhpindex at resolution res (immediate parent if res == None).

    Returns None if the cell index is invalid.
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
    rhpindex: str, res: int = None, verbose: bool = True, dggs: RHEALPixDGGS = WGS84_003
) -> str:
    """
    Returns central child of rhpindex at resolution res (immediate central
    child if res == None).

    Returns None if the cell index is invalid.

    Returns None if the DGGS has an even number of cells on a side.
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
    parent_res = len(rhpindex)
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
        child_index = rhpindex + "".join(str(c_index) for _ in range(0, added_levels))

        return child_index


def rhp_to_geo_boundary(
    rhpindex: str,
    geo_json: bool = True,
    plane: bool = True,
    dggs: RHEALPixDGGS = WGS84_003,
) -> tuple[tuple[float, float]]:
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

    return verts


def rhp_get_resolution(rhpindex: str, dggs: RHEALPixDGGS = WGS84_003) -> int:
    """
    Returns the resolution of a given cell index (or None if invalid).
    """
    if not rhp_is_valid(rhpindex, dggs):
        return None

    return len(rhpindex) - 1


def rhp_get_base_cell(rhpindex: str, dggs: RHEALPixDGGS = WGS84_003) -> str:
    """
    Returns the resolution 0 cell id of a given cell index (or None if invalid).
    """
    if not rhp_is_valid(rhpindex, dggs):
        return None

    return rhpindex[0]


def rhp_is_valid(rhpindex: str, dggs: RHEALPixDGGS = WGS84_003) -> bool:
    """
    Checks if the given cell address is valid within the DGGS
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
    plane=True,
    dggs: RHEALPixDGGS = WGS84_003,
) -> float:
    """
    Returns the area of a cell in the requested unit (or None if rhpindex is invalid).

    TODO: investigate use case where unit is 'rads^2'
    """
    if not rhp_is_valid(rhpindex, dggs):
        return None

    # Grab cell area in native unit (m^2)
    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    cell = dggs.cell(suid)
    area = cell.area(plane=plane)

    # Scale area if needed
    if unit == "km^2":
        area = area / 10**6

    return area


def cell_ring(
    rhpindex: str, k: int = 1, verbose: bool = True, dggs: RHEALPixDGGS = WGS84_003
) -> list[str]:
    """
    Returns the ring of cell indices around rhpindex at distance k, or None if rhpindex
    is invalid.

    Also returns None if k < 0.

    Returns [ rhpindex ] if k == 0 (by convention).

    Returns the four neighbouring faces at resolution 0 if k > 0 and cell resolution is 0
    (by convention).
    """
    if verbose:
        warn(str.format(CELL_RING_WARNING, "cell"))

    if not rhp_is_valid(rhpindex, dggs) or (k < 0):
        return None

    # A cell ring at distance 0 just consists of the cell itself
    if k == 0:
        return [rhpindex]

    # Grab the centre cell
    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    cell = dggs.cell(suid)

    # Maximum ring distance from centre cell
    half_circle = 2 * cell.N_side ** rhp_get_resolution(rhpindex)

    # Just return the opposite cell if k is beyond what the resolution can do
    if k > half_circle:
        cell = _mirror_cell_on_cube(cell)
        return [str(cell)]

    # Init the ring and directions
    ring = []
    directions = ["right", "down", "left", "up"]

    # Top-level cells (cube faces) are special
    if len(rhpindex) == 1:
        for direction in directions:
            ring.append(cell.neighbor(direction).suid[0])

    # Start in the upper left corner of the ring: it's k times left and k times up
    else:
        # Mapping to detect direction changes
        direction_inverse = {
            "right": "left",
            "down": "up",
            "left": "right",
            "up": "down",
        }

        # Initialise iteration parameters
        k_eff, max_steps, cell = _cell_ring_setup(cell, half_circle / 2, k)

        # We're done if k_eff takes us all the way to the opposite cell (shouldn't happen at this point...)
        if k_eff < 1:
            ring.append(str(cell))

        # We have to do the full walk around the ring
        else:
            # Set starting point
            cell, direction, n_steps = _find_cell_ring_start(
                cell, k_eff, max_steps, directions, direction_inverse
            )

            # Walk around the ring one side at a time and collect cell addresses
            for _ in range(0, len(directions)):
                step = 0
                while step < n_steps:
                    # Add index to ring, take a step
                    ring.append(str(cell))
                    next = cell.neighbor(direction)

                    # Looking back not being the same as looking ahead means we need to realign
                    if next.neighbor(direction_inverse[direction]) != cell:
                        direction = direction_inverse[_neighbor_direction(next, cell)]

                    # Take the step
                    cell = next
                    step = step + 1

                # Prepare walking direction for next ring side
                if n_steps == 2 * k_eff:
                    direction = directions[
                        (directions.index(direction) + 1) % len(directions)
                    ]

                # Reset number of steps along a side
                n_steps = max_steps

    return ring


def k_ring(
    rhpindex: str, k: int = 1, verbose: bool = True, dggs: RHEALPixDGGS = WGS84_003
) -> list[str]:
    """
    Returns the k-ring of cell indices around rhpindex at distance k (or None if rhpindex is invalid).

    Also returns None if k < 0.
    """
    if verbose:
        warn(str.format(CELL_RING_WARNING, "k"))

    if not rhp_is_valid(rhpindex, dggs) or (k < 0):
        return None

    # A k-ring at distance 0 just consists of the centre cell itself
    if k == 0:
        return [rhpindex]

    distance = min(2 * dggs.N_side ** rhp_get_resolution(rhpindex), k)
    kring = [rhpindex]

    for d in range(1, distance + 1):
        kring = kring + cell_ring(rhpindex, d, verbose=False, dggs=dggs)

    return kring


def polyfill(
    geometry: Union[Polygon, MultiPolygon],
    res: int,
    plane: bool = True,
    compress: bool = False,
    verbose: bool = False,
    dggs: RHEALPixDGGS = WGS84_003,
) -> set[str]:
    """
    Turn the area contained in a shapely polygon or multipolygon into a set of cell
    indices at the requested resolution. A cell index is included if its centroid is
    inside the geometry defined by the boundaries and holes.

    Returns an empty set if no cell centroids fall within the input geometry.

    Returns None if the geom_type field in the input geometry is anything other than
    'Polygon' or 'MultiPolygon'.

    Returns None if the geometry is empty, or if it has no area.

    Returns None if no cells match the geometry for some reason.

    Returns None if the geometry is invalid in other ways, e.g. if a point on a hole
    boundary is outside the exterior boundary of its polygon, or if two polygons in a
    multipolygon overlap.

    TODO: decide what to do with the antimeridian (if anything)
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

    # Extract list of polygons from geometry: Polygon needs to be wrapped in
    # one, MultiPolygon has it stashed in a property
    if geometry.geom_type == "Polygon":
        geoms = [geometry]
    else:
        geoms = geometry.geoms

    # Collect cells in regions of interest
    cells = set()
    for geom in geoms:
        # Region of interest is the bounding box around the geometry
        bbox = geom.bounds

        # rhealpixdggs wants nw and se corners of region of interest
        nw = (bbox[0], bbox[3])
        se = (bbox[2], bbox[1])

        # Cells in bounding box at requested resolution
        roi_cells = dggs.cells_from_region(res, nw, se, plane)

        if roi_cells:
            # Flatten list of lists of cells in bbox
            roi_cells = [cell for nested_list in roi_cells for cell in nested_list]

            # Check each cell against geometry, add to results if inside polygon
            for cell in roi_cells:
                if geom.contains(Point(cell.centroid(plane))):
                    cells.add(str(cell))

    # Merge cells inside polygon into larger ones where possible
    if compress:
        cells = set(compress_order_cells(cells))

    return cells


def linetrace(
    geometry, res: int, plane: bool = True, dggs: RHEALPixDGGS = WGS84_003
) -> list[str]:
    raise NotImplementedError()


# ======== Helper functions ======== #


def _neighbor_direction(cell: Cell, neighbor: Cell) -> str:
    n_dict = cell.neighbors()
    for dir in n_dict:
        if n_dict[dir] == neighbor:
            return dir

    return None


def _mirror_cell_on_cube(cell: Cell) -> Cell:
    # Cube faces map to their opposites
    face_mapping = {"N": "S", "S": "N", "O": "Q", "P": "R", "Q": "O", "R": "P"}
    transformed_suid = [face_mapping[cell.suid[0]]]

    # Skip the numerical digits part for top-level cells
    if len(cell.suid) > 1:
        # Transform row or column indices depending on region and rearrange into pairs
        region = cell.region()
        rowcolidxs = cell.suid_rowcol()
        rowidxs = rowcolidxs[0][1:]
        colidxs = rowcolidxs[1][1:]
        transformed_idxs = [
            cell.N_side - idx - 1
            for idx in (rowidxs if region == "equatorial" else colidxs)
        ]
        rowcols = (
            zip(transformed_idxs, colidxs)
            if region == "equatorial"
            else zip(rowidxs, transformed_idxs)
        )

        # Reapply transformed indices to suid digits
        for row, col in rowcols:
            transformed_suid.append(cell.N_side * row + col)

    return Cell(cell.rdggs, transformed_suid)


def _cell_ring_setup(
    cell: Cell, quarter_circle: int, k: int
) -> tuple[int, int, Cell, bool]:
    # Cell ring distance farther than the hemisphere equator requires mirroring
    if k > quarter_circle:
        k_eff = max(2 * quarter_circle - k, 0)
        starting_cell = _mirror_cell_on_cube(cell)
    else:
        k_eff = k
        starting_cell = cell

    # Cell ring distance taking k beyond resolution requires clamping
    if 2 * k_eff > quarter_circle:
        max_steps = quarter_circle
    else:
        max_steps = 2 * k_eff

    return (k_eff, max_steps, starting_cell)


def _find_cell_ring_start(
    cell: Cell,
    k: int,
    max_steps: int,
    directions: list[str],
    direction_inverse: dict[str, str],
) -> tuple[Cell, str, int]:
    # Always start by going left
    dir_idx = directions.index("left")

    # Work your way to the starting point one (left, up) pair of steps at a time
    steps_from_start = -1
    num_edges = 0
    d = 0
    while d < k:
        # Prep for later
        d = d + 1

        # One step to local left
        direction = directions[dir_idx]
        next = cell.neighbor(direction)

        # Detect edge crossing on cube
        if cell.suid[0] != next.suid[0]:
            num_edges = num_edges + 1
            # Looking back not being the same as looking ahead means we need to realign as well
            if next.neighbor(direction_inverse[direction]) != cell:
                dir_idx = directions.index(
                    direction_inverse[_neighbor_direction(next, cell)]
                )

        # Take the step
        cell = next

        # One step to local up
        direction = directions[(dir_idx + 1) % len(directions)]
        next = cell.neighbor(direction)

        # Detect edge crossing on cube
        if cell.suid[0] != next.suid[0]:
            num_edges = num_edges + 1

            # Looking back not being the same as looking ahead means we need to realign as well
            if next.neighbor(direction_inverse[direction]) != cell:
                dir_idx = (
                    directions.index(direction_inverse[_neighbor_direction(next, cell)])
                    - 1
                ) % len(directions)

            # Handle unexpected corner and end the loop
            if num_edges > 1:
                dir_idx = (dir_idx - 1) % len(directions)
                steps_from_start = d
                d = k

        # Take the step
        cell = next

    # Initialise walking direction and side length
    direction = direction_inverse[directions[dir_idx]]
    if steps_from_start >= 0:
        n_steps = min(k + steps_from_start - 1, max_steps)  # TODO: this is wrong
        local_up = directions[(directions.index(direction) - 1) % len(directions)]
        for _ in range(0, k - steps_from_start):
            next = cell.neighbor(local_up)
            cell = next
    else:
        n_steps = max_steps

    return (cell, direction, n_steps)


def _malformed_geometry(geometry: Union[Polygon, MultiPolygon]) -> bool:
    # Geometry has to have things in it
    if geometry is None or geometry.is_empty:
        return True

    # Geometry needs to be of the correct type
    if geometry.geom_type != "Polygon" and geometry.geom_type != "MultiPolygon":
        return True

    if not geometry.is_valid:
        return True

    # Geometry has to have an area, i.e. not be collapsed to a line
    if geometry.area == 0:
        return True

    return False
