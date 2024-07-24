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

    TODO: give the option to select another predefined DGGS, or pass in a custom one
    TODO: give the option to set n to something other than 3
    TODO: give the option to enter the coordinates in radians
    """
    # Get the grid cell corresponding to the coordinates
    cell = WGS84_003.cell_from_point(resolution, (lng, lat), plane)

    # Return the cell ID after converting int digits to str
    return "".join([str(d) for d in cell.suid])


def rhp_to_geo(
    rhpindex: str, geo_json: bool = True, plane: bool = True
) -> tuple[float, float]:
    """
    Look up the centroid (in degrees) of the cell identified by rhpindex.

    If geojson is requested as the output format:
        - Will return a (longitude, latitude) coordinate pair.

    if geojson is NOT requested as the output format:
        - Will return a (latitude, longitude) coordinate pair in order to be consistent with
          h3 coordinate ordering.

    TODO: give the option of requesting centroid coordinates in radians
    """
    # Grab cell centroid matching rhpindex string
    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    cell = WGS84_003.cell(suid)
    centroid = cell.centroid(plane=plane)

    # rhealpix coordinates come out natively as lng/lat, h3 ones as lat/lng
    if not geo_json:
        # Swap coordinates
        centroid = centroid[::-1]

    return centroid


def rhp_to_parent(rhpindex: str, res: int = None, verbose: bool = True) -> str:
    """
    Return parent of rhpindex at resolution res (immediate parent if res == None)
    """
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


def rhp_to_geo_boundary(
    rhpindex: str, geo_json: bool = True, plane: bool = True
) -> tuple[tuple[float, float]]:
    """
    Extract the corner coordinates of a cell at a given cell ID and returns them as
    a tuple of coordinate pairs (in degrees).

    If geojson is requested as the output format:
        - Will return (longitude, latitude) coordinate pairs.
        - Will repeat the first vertex and insert it at the end if geojson is requested as
          the output format.

    If geojson is NOT requested as the output format:
        - Will return (latitude, longitude) coordinate pairs in order to be consistent with
          h3 coordinate ordering.

    TODO: give the option of requesting corner coordinates in radians
    """
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
