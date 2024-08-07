from typing import Literal

# Pre-defined DGGS using WGS84 ellipsoid and n == 3 for cell side subpartitioning
from rhealpixdggs.dggs import WGS84_003

from rhealpixdggs.cell import Cell, CELLS0


def geo_to_rhp(lat: float, lng: float, resolution: int, plane: bool = True) -> str:
    """
    Turn a latitute and longitude (in degrees) into an rHEALPix cell address at
    the requested resolution.

    Uses the predefined WGS84_003 DGGS with the WGS84 ellipsoid and n = 3 to
    subdivide the cell sides.

    Mostly passes through the parameters to the function turning coordinate points
    into cells, but converts the address tuple from the resulting cell into a
    string.

    Returns None if no cell matching the coordinates is found.

    TODO: give the option to select another predefined DGGS, or pass in a custom one
    TODO: give the option to set n to something other than 3
    TODO: give the option to enter the coordinates in radians
    """
    # Get the grid cell corresponding to the coordinates
    cell = WGS84_003.cell_from_point(resolution, (lng, lat), plane)

    # Bail out if there's no matching cell
    if cell is None:
        return None

    # Return the cell ID after converting int digits to str
    return "".join([str(d) for d in cell.suid])


def rhp_to_geo(
    rhpindex: str, geo_json: bool = True, plane: bool = True
) -> tuple[float, float]:
    """
    Look up the centroid (in degrees) of the cell identified by rhpindex.

    Returns None if the cell index is invalid.

    If geojson is requested as the output format:
        - Will return a (longitude, latitude) coordinate pair.

    if geojson is NOT requested as the output format:
        - Will return a (latitude, longitude) coordinate pair in order to be consistent with
          h3 coordinate ordering.

    TODO: give the option to select another predefined DGGS, or pass in a custom one
    TODO: give the option to set n to something other than 3
    TODO: give the option of requesting centroid coordinates in radians
    """
    # Stop early if the cell index is invalid
    if not rhp_is_valid(rhpindex):
        return None

    # Grab cell centroid matching rhpindex string
    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    cell = WGS84_003.cell(suid)
    centroid = cell.centroid(plane=plane)

    # rhealpix coordinates come out natively as lng/lat, h3 ones as lat/lng
    if not geo_json:
        # lng/lat -> lat/lng to make it consistent with h3
        centroid = centroid[::-1]

    return centroid


def rhp_to_parent(rhpindex: str, res: int = None, verbose: bool = True) -> str:
    """
    Returns parent of rhpindex at resolution res (immediate parent if res == None).

    Returns None if the cell index is invalid.
    """
    # Stop early if the cell index is invalid
    if not rhp_is_valid(rhpindex):
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
            print(
                f"Warning: You requested a parent resolution that is higher than the cell resolution. Returning the cell address itself."
            )
        return rhpindex

    # Standard case (including child_res == res)
    else:
        return rhpindex[: res + 1]


def rhp_to_center_child(rhpindex: str, res: int = None, verbose: bool = True) -> str:
    """
    Returns central child of rhpindex at resolution res (immediate central
    child if res == None).

    Returns None if the cell index is invalid.

    TODO: come up with a scheme for even numbers on a side
    """
    # Stop early if the cell index is invalid
    if not rhp_is_valid(rhpindex):
        return None

    # Handle mismatch between cell resolution and requested child resolution
    parent_res = len(rhpindex)
    if res is not None and res < parent_res:
        if verbose:
            print(
                f"Warning: You requested a child resolution that is lower than the cell resolution. Returning the cell address itself."
            )
        return rhpindex

    # Standard case (including parent_res == res)
    else:
        # res == None returns the central child from one level down (by convention)
        added_levels = 1 if res is None else res - parent_res

        # Derive index of centre child and append that to rhpindex
        # NOTE: only works for odd values of N_side
        c_index = int((WGS84_003.N_side**2 - 1) / 2)

        # Append the required number of child digits to cell index
        child_index = rhpindex + "".join(str(c_index) for _ in range(0, added_levels))

        return child_index


def rhp_to_geo_boundary(
    rhpindex: str, geo_json: bool = True, plane: bool = True
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
          h3 coordinate ordering.

    TODO: give the option to select another predefined DGGS, or pass in a custom one
    TODO: give the option to set n to something other than 3
    TODO: give the option of requesting corner coordinates in radians
    """
    # Stop early if the cell index is invalid
    if not rhp_is_valid(rhpindex):
        return None

    # Grab the cell vertices (includes non-corner point in darts if plane == False)
    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    cell = WGS84_003.cell(suid)
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


def rhp_get_resolution(rhpindex: str) -> int:
    """
    Returns the resolution of a given cell index (or None if invalid).
    """
    if not rhp_is_valid(rhpindex):
        return None

    return len(rhpindex) - 1


def rhp_get_base_cell(rhpindex: str) -> str:
    """
    Returns the resolution 0 cell id of a given cell index (or None if invalid).
    """
    if not rhp_is_valid(rhpindex):
        return None

    return rhpindex[0]


def rhp_is_valid(rhpindex: str) -> bool:
    # TODO: call this function in those that include address-based functionality
    #       i.e. anything that accepts a string 'rhpindex' as an argument
    """
    Checks if the given cell address is valid within the DGGS

    TODO: give the option to select another predefined DGGS, or pass in a custom one
    TODO: give the option to set n to something other than 3
    """
    # Empty strings are invalid
    if rhpindex is None or len(rhpindex) == 0:
        return False

    # Addresses that don't start with the resolution 0 face are invalid
    if rhpindex[0] not in CELLS0:
        return False

    # Addresses that have digits out of range are invalid
    num_subcells = WGS84_003.N_side**2
    for d in rhpindex[1:]:
        if not d.isdigit() or (int(d) >= num_subcells):
            return False

    # Passed all checks - must be the real thing
    return True


def cell_area(
    rhpindex: str, unit: Literal["km^2", "m^2"] = "km^2", plane=True
) -> float:
    """
    Returns the area of a cell in the requested unit (or None if rhpindex is invalid).

    TODO: investigate use case where unit is 'rads^2'
    """
    if not rhp_is_valid(rhpindex):
        return None

    # Grab cell area in native unit (m^2)
    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    cell = WGS84_003.cell(suid)
    area = cell.area(plane=plane)

    # Scale area if needed
    if unit == "km^2":
        area = area / 10**6

    return area


def cell_ring(rhpindex: str, k: int = 1) -> list[str]:
    """
    Returns the ring of cell indices around rhpindex at distance k, in clockwise order and
    without duplicates (or None if rhpindex is invalid).

    Also returns None if k < 0.

    Returns [ rhpindex ] if k == 0 (by convention).

    Returns the four neighbouring faces at resolution 0 if k > 0 and cell resolution is 0
    (by convention).
    """
    if not rhp_is_valid(rhpindex) or (k < 0):
        return None

    # A cell ring at distance 0 just consists of the cell itself
    if k == 0:
        return [rhpindex]

    # Grab the centre cell
    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    cell = WGS84_003.cell(suid)

    # Init the ring and directions
    ring = []
    directions = ["right", "down", "left", "up"]

    # Top-level cells (cube faces) are special
    if len(rhpindex) == 1:
        for direction in directions:
            ring.append(cell.neighbor(direction).suid[0])

    # Start in the upper left corner of the ring: it's k times up and k times left
    else:
        # Mapping to detect direction changes
        direction_inverse = {
            "right": "left",
            "down": "up",
            "left": "right",
            "up": "down",
        }
        # Always start by going left
        dir_idx = directions.index("left")

        # Work your way to the starting point
        for _ in range(0, k):
            # One step to local left
            direction = directions[dir_idx]
            next = cell.neighbor(direction)
            # Looking back not being the same as looking ahead means we need to realign
            if next.neighbor(direction_inverse[direction]) != cell:
                dir_idx = directions.index(
                    direction_inverse[_neighbor_direction(next, cell)]
                )
            cell = next

            # One step to local up
            direction = directions[(dir_idx + 1) % len(directions)]
            next = cell.neighbor(direction)
            # Looking back not being the same as looking ahead means we need to realign
            if next.neighbor(direction_inverse[direction]) != cell:
                dir_idx = (
                    directions.index(direction_inverse[_neighbor_direction(next, cell)])
                    - 1
                ) % len(directions)
            cell = next

        # Initialise walking direction
        direction = direction_inverse[directions[dir_idx]]

        # Walk around the ring one side at a time and collect cell addresses
        for _ in range(0, len(directions)):
            for _ in range(0, 2 * k):
                # Add index to ring, take a step
                ring.append("".join([str(d) for d in cell.suid]))
                next = cell.neighbor(direction)

                # Looking back not being the same as looking ahead means we need to realign
                if next.neighbor(direction_inverse[direction]) != cell:
                    direction = direction_inverse[_neighbor_direction(next, cell)]

                cell = next

            # Update walking direction before going around a corner
            direction = directions[(directions.index(direction) + 1) % len(directions)]

    return _eliminate_duplicates(ring)


def _eliminate_duplicates(l: list) -> list:
    """
    Helper function to eliminate duplicates in a list while preserving list order
    """
    if l is not None:
        u = []
        for i in l:
            if i not in u:
                u.append(i)

        return u

    return None


def _neighbor_direction(cell: Cell, neighbor: Cell) -> str:
    n_dict = cell.neighbors()
    for dir in n_dict:
        if n_dict[dir] == neighbor:
            return dir

    return None
