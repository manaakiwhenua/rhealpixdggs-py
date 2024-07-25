# Pre-defined DGGS using WGS84 ellipsoid and n == 3 for cell side subpartitioning
from rhealpixdggs.dggs import WGS84_003


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


def rhp_to_center_child() -> str:  # TODO: function arguments
    pass


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
    if rhpindex[0] not in WGS84_003.cells0:
        return False

    # Addresses that have digits out of range are invalid
    num_subcells = WGS84_003.N_side * WGS84_003.N_side
    for d in rhpindex[1:]:
        if not d.isdigit() or (int(d) >= num_subcells):
            return False

    # Passed all checks - must be the real thing
    return True


def cell_area() -> float:  # TODO: function arguments
    pass
