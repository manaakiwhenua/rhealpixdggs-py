# Pre-defined DGGS using WGS84 ellipsoid and n == 3 for cell side subpartitioning
from rhealpixdggs.dggs import WGS84_003


def geo_to_rhp(lat: float, lng: float, resolution: int, plane: bool = True) -> str:
    # Get the grid cell corresponding to the coordinates
    cell = WGS84_003.cell_from_point(resolution, (lng, lat), plane)

    # Return the cell ID after converting int digits to str
    return "".join([str(d) for d in cell.suid])


def rhp_to_parent(rhpindex: str, res: int) -> str:
    """
    Return parent of rhpindex at resolution res (immediate parent if res == None)
    """
    # Top-level cells are their own parent, regardless of the requested resolution
    # TODO: check if this convention actually applies?
    child_res = len(rhpindex) - 1
    if child_res < 1:
        return rhpindex

    # res == None returns the first address up (by convention)
    elif res is None:
        return rhpindex[:-1]

    # Handle mismatch between cell resolution and requested parent resolution
    elif res > child_res:
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
    # Grab the cell vertices as tuple((lat, lng), (lat, lng), (lat, lng), (lat, lng))
    # (Includes non-corner point in darts if plane == False)
    suid = [int(d) if d.isdigit() else d for d in rhpindex]
    cell = WGS84_003.cell(suid)
    verts = tuple(cell.vertices(plane=plane))

    if geo_json:
        # lat/lng -> lng/lat and last point same as first (pinched from h3 package)
        verts += (verts[0],)
        verts = tuple(v[::-1] for v in verts)

    return verts
