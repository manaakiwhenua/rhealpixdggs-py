"""
Measure how far cheap distance and bearing measures stray from the geodesic
on the WGS84 ellipsoid, for the error bounds the proximity measures state
(issue #170):

- great-circle distance and initial bearing on the authalic sphere, with
  the longitude and the authalic latitude of each point;
- Euclidean distance in the rHEALPix plane between points on the same face,
  which hop counts approximate.

Run from the repository root to print the bounds:

    python docs/proximity_error.py

``sample`` returns the raw errors, for the figure in the docs.
"""

import itertools
import pathlib
import sys

import numpy as np
from pyproj import Geod

sys.path.insert(0, str(pathlib.Path(__file__).parents[1]))

from rhealpixdggs.dggs import WGS84_003_OGC
from rhealpixdggs.ellipsoids import WGS84_A, WGS84_E, WGS84_R_A

GEOD = Geod(ellps="WGS84")
# Distances in metres, log-spaced from 10 m to beyond the antipode.
DISTANCE_RANGE = (1e1, 2.0e7)
BINS = (1e1, 1e3, 1e5, 1e6, 5e6, 1.5e7, 1.9e7, 2.1e7)


def authalic_latitude(lat):
    """The authalic latitude, in degrees, by the exact formula."""
    e = WGS84_E
    s = np.sin(np.radians(lat))
    q = (1 - e**2) * (
        s / (1 - (e * s) ** 2) - np.log((1 - e * s) / (1 + e * s)) / (2 * e)
    )
    qp = 1 - (1 - e**2) / (2 * e) * np.log((1 - e) / (1 + e))
    return np.degrees(np.arcsin(np.clip(q / qp, -1, 1)))


def sphere_inverse(lon1, lat1, lon2, lat2):
    """Great-circle distance and initial bearing on the authalic sphere."""
    b1 = np.radians(authalic_latitude(lat1))
    b2 = np.radians(authalic_latitude(lat2))
    dl = np.radians(np.asarray(lon2) - lon1)
    h = np.sin((b2 - b1) / 2) ** 2 + np.cos(b1) * np.cos(b2) * np.sin(dl / 2) ** 2
    distance = 2 * WGS84_R_A * np.arcsin(np.sqrt(h))
    bearing = np.degrees(
        np.arctan2(
            np.sin(dl) * np.cos(b2),
            np.cos(b1) * np.sin(b2) - np.sin(b1) * np.cos(b2) * np.cos(dl),
        )
    )
    return distance, bearing


def frame_block(x, y, rdggs=WGS84_003_OGC):
    """The (row, column) of the face block of the plane containing (x, y)."""
    side = rdggs.ellipsoid.R_A * np.pi / 2
    left = rdggs.ul_vertex["O"][0]
    top = rdggs.ul_vertex["O"][1] + side
    return np.floor((top - y) / side).astype(int), np.floor((x - left) / side).astype(
        int
    )


def sample(n=400_000, seed=0):
    """
    Random pairs, uniform over the ellipsoid in the first point, uniform in
    direction and log-uniform in distance. Returns a dict of arrays: the
    geodesic distance ``distance`` and ``latitude`` of the first point,
    ``sphere_ratio`` (sphere distance over geodesic), ``bearing_error``
    (degrees), and ``plane_ratio`` (planar distance over geodesic, NaN for
    pairs on different faces) with ``polar`` marking pairs on a polar face.
    """
    rng = np.random.default_rng(seed)
    lon1 = rng.uniform(-180, 180, n)
    lat1 = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    azimuth = rng.uniform(-180, 180, n)
    lo, hi = np.log10(DISTANCE_RANGE)
    lon2, lat2, _ = GEOD.fwd(lon1, lat1, azimuth, 10 ** rng.uniform(lo, hi, n))
    # Past the antipode fwd's path is no longer the shortest, so measure anew.
    azimuth, _, distance = GEOD.inv(lon1, lat1, lon2, lat2)
    sphere_distance, sphere_bearing = sphere_inverse(lon1, lat1, lon2, lat2)

    x1, y1 = map(np.asarray, WGS84_003_OGC.rhealpix(lon1, lat1))
    x2, y2 = map(np.asarray, WGS84_003_OGC.rhealpix(lon2, lat2))
    (r1, c1), (r2, c2) = frame_block(x1, y1), frame_block(x2, y2)
    same_face = (r1 == r2) & (c1 == c2)
    plane_ratio = np.where(same_face, np.hypot(x2 - x1, y2 - y1) / distance, np.nan)

    return {
        "distance": distance,
        "latitude": lat1,
        "sphere_ratio": sphere_distance / distance,
        "bearing_error": (sphere_bearing - azimuth + 180) % 360 - 180,
        "plane_ratio": plane_ratio,
        "polar": r1 != 1,
    }


def main():
    s = sample()
    d = s["distance"]
    sphere = np.abs(s["sphere_ratio"] - 1) * 100
    bearing = np.abs(s["bearing_error"])
    # Local scale of the ellipsoid on the authalic sphere along the equator,
    # whose reciprocal is the scale along the meridian there.
    k = WGS84_A / WGS84_R_A
    print("Authalic sphere against the WGS84 geodesic")
    print(
        f"  analytic: distance {(k - 1) * 100:.3f} %, bearing at short range "
        f"{np.degrees(np.arcsin((k - 1 / k) / (k + 1 / k))):.3f} deg"
    )
    print("  geodesic distance (km)      n   distance %   bearing deg")
    for lo, hi in itertools.pairwise(BINS):
        m = (d >= lo) & (d < hi)
        print(
            f"  {lo / 1e3:8g} - {hi / 1e3:7g} {m.sum():7d} {sphere[m].max():10.3f} "
            f"{bearing[m].max():13.3f}"
        )

    print("Euclidean distance in the plane over the geodesic, same face only")
    print(
        f"  analytic at short range, equatorial: {3 * np.pi / 8 * np.sqrt(5) / 3:.3f} to "
        f"{3 / np.sqrt(5):.3f}, before the ellipsoid's +-{(k - 1) * 100:.3f} %"
    )
    for name, m in (("equatorial", ~s["polar"]), ("polar", s["polar"])):
        ratio = s["plane_ratio"][m]
        ratio = ratio[~np.isnan(ratio)]
        print(f"  {name:10s} {ratio.size:7d}   {ratio.min():.3f} to {ratio.max():.3f}")


if __name__ == "__main__":
    main()
