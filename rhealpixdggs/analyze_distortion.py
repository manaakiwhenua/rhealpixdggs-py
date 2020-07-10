"""
This Python 3.3 module analyzes angular, linear, and areal distortion of various map projections.

CHANGELOG:

- Alexander Raichev (AR), 2013-02-04: Initial version.
- AR, 2013-07-23: Ported to Python 3.3.

"""
# *****************************************************************************
#       Copyright (C) 2013 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import my modules.
from rhealpixdggs.ellipsoids import EMR_SPHERE, WGS84_ELLIPSOID
from rhealpixdggs.projection_wrapper import Proj
from rhealpixdggs.distortion import distortion_stats


if __name__ == "__main__":
    # List of map projections to analyze.
    # They must be implemented in the pyproj third-party Python module or
    # homemade, e.g. csea as implemented in pj_csea.py.
    projections = [
        ("HEALPix", "healpix"),
        ("rHEALPix", "rhealpix"),
        ("Lambert cylindrical equal area", "cea"),
        ("Lambert azimuthal equal area", "laea"),
        ("Quadrilateralized spherical cube", "qsc"),
        ("Cubic Snyder equal area", "csea"),
        ("Icosahedral Snyder equal area", "isea"),
        ("Mercator", "merc"),
        ("Stereographic", "stere"),
        ("Lambert conformal conic", "lcc"),
        ("Equirectangular", "eqc"),
        ("Azimuthal equidistant", "aeqd"),
        ("Winkel triple", "wintri"),
    ]

    N = 30000
    for E in [EMR_SPHERE, WGS84_ELLIPSOID]:
        print(E)
        sample = [E.random_point(phi_min=-89.5, phi_max=89.5) for i in range(N)]
        for name, proj in projections:
            f = Proj(E, proj)
            stats = distortion_stats(f, sample)[1]
            print((name + " (" + proj + "):"))
            print(stats)
