"""
This Python 3.11 code implements ellipsoids of revolution.

- Alexander Raichev (AR), 2012-01-26: Refactored code from release 0.3.

NOTE:

All lengths are measured in meters and all angles are measured in radians
unless indicated otherwise.

By 'ellipsoid' throughout, I mean an ellipsoid of revolution and *not* a general (triaxial) ellipsoid.
Points lying on an ellipsoid are given in geodetic (longitude, latitude) coordinates.
"""
# *****************************************************************************
#       Copyright (C) 2013 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import third-party modules.
from numpy import pi, sqrt, sin, cos, arcsin, arctanh, deg2rad, rad2deg

# Import standard modules.
from random import uniform

# Import my modules.
from rhealpixdggs.utils import my_round, auth_lat, auth_rad

# Parameters of some common ellipsoids.
WGS84_A = 6378137.0
WGS84_F = (
    1 / 298.257222101
)  # ORIGINAL: GRS80 from EPSG:42310 298.257222101, based on WGS84+GRS80
# WGS84_F = 1 / 298.257222100882711  # GRS80 from https://en.wikipedia.org/wiki/World_Geodetic_System
# WGS84_F = 1 / 298.257223563  # new value from EPSG:7030
WGS84_B = WGS84_A * (1 - WGS84_F)
WGS84_E = sqrt(WGS84_F * (1 - WGS84_F))
WGS84_R_A = sqrt(WGS84_A**2 / 2 + WGS84_B**2 / 2 * (arctanh(WGS84_E) / WGS84_E))
R_EM = 6371000  # Earth's mean radius


class Ellipsoid(object):
    """
    Represents an ellipsoid of revolution (possibly a sphere) with a
    geodetic longitude-latitude coordinate frame.

    INSTANCE ATTRIBUTES:

    - `sphere` - True if the ellipsoid is a sphere, and False otherwise.
    - `R` - The radius of the ellipsoid in meters, implying that it is a
      sphere.
    - `a` - Major radius of the ellipsoid in meters.
    - `b` - Minor radius of the ellipsoid in meters.
    - `e` - Eccentricity of the ellipsoid.
    - `f` - Flattening of the ellipsoid.
    - `R_A` - Authalic radius of the ellipsoid in meters.
    - `lon_0` - Central meridian.
    - `lat_0` - Latitude of origin.
    - `radians` - If True, use angles measured in radians for all calculations.
      Use degrees otherwise.
    - `phi_0` - The latitude separating the equatorial region and
      the north polar region in the context of the (r)HEALPix projection.

    Except for phi_0, these attribute names match the names of the
    `PROJ.4 ellipsoid parameters <http://trac.osgeo.org/proj/wiki/GenParms>`_.
    """

    def __init__(
        self,
        R=None,
        a=WGS84_A,
        b=None,
        e=None,
        f=WGS84_F,
        lon_0=0,
        lat_0=0,
        radians=False,
    ):
        self.lon_0 = lon_0
        self.lat_0 = lat_0
        self.radians = radians
        if R is not None:
            # The ellipsoid is a sphere.
            # Override the other geometric parameters.
            self.sphere = True
            self.R = R
            self.a = R
            self.b = R
            self.e = 0
            self.f = 0
            self.R_A = R
        else:
            self.sphere = False
            self.a = a
            if b is not None:
                # Derive the other geometric parameters from a and b.
                self.b = b
                self.e = sqrt(1 - (b / a) ** 2)
                self.f = (a - b) / a
            elif e is not None:
                # Derive the other geometric parameters from a and e.
                self.e = e
                self.b = a * sqrt(1 - e**2)
                self.f = 1 - sqrt(1 - e**2)
            else:
                self.f = f
                self.b = self.a * (1 - f)
                self.e = sqrt(f * (1 - f))
            self.R_A = auth_rad(self.a, self.e)
        self.phi_0 = auth_lat(arcsin(2.0 / 3), e=self.e, radians=True, inverse=True)
        if not self.radians:
            # Convert to degrees.
            self.phi_0 = rad2deg(self.phi_0)

    def __str__(self):
        result = ["ellipsoid:"]
        # result.append('lengths measured in meters')
        for k, v in sorted(self.__dict__.items()):
            if k == "phi_0":
                continue
            if k in {"sphere", "radians"}:
                result.append("    " + k + " = " + str(v))
            else:
                result.append("    " + k + " = " + str(my_round(v, 15)))
        return "\n".join(result)

    def __eq__(self, other):
        if self.a == other.a and self.b == other.b:
            return True
        else:
            return False

    def __ne__(self, other):
        """
        The inequality relation on cells.
        Since Python 3.3 doesn't automatically create reverse relations
        from given ones, i must define this seemingly redundant relation.
        """
        return not self.__eq__(other)

    def pi(self):
        """
        Return pi if `self.radians` = True and 180 otherwise.
        """
        if self.radians:
            return pi
        else:
            return 180.0

    def random_point(self, lam_min=None, lam_max=None, phi_min=None, phi_max=None):
        """
        Return a point (given in geodetic coordinates) sampled uniformly at
        random from the section of this ellipsoid with longitude in the range
        `lam_min <= lam < lam_max` and latitude in the range
        `phi_min <= phi < phi_max`.
        But avoid the poles.

        EXAMPLES::

           >>> E = UNIT_SPHERE
           >>> print(E.random_point()) # doctest: +SKIP
           (-1.0999574573422948, 0.21029104897701129)

        """
        PI = self.pi()
        if lam_min is None:
            lam_min = -PI
        if lam_max is None:
            lam_max = PI
        if phi_min is None:
            phi_min = -PI / 2
        if phi_max is None:
            phi_max = PI / 2
        if not self.radians:
            # Convert to radians.
            lam_min, lam_max, phi_min, phi_max = deg2rad(
                [lam_min, lam_max, phi_min, phi_max]
            )
        # Pick a longitude.
        while True:
            u = uniform(0, 1)
            lam = (lam_max - lam_min) * u + lam_min
            # Don't include lam_max.
            if lam < lam_max:
                # Success.
                break
        # Pick a latitude.
        delta = pi / 360
        while True:
            v = uniform(0, 1)
            if self.sphere:
                phi = arcsin((sin(phi_max) - sin(phi_min)) * v + sin(phi_min))
            else:
                # Sample from the authalic sphere.
                # The map from the ellipsoid to the authalic sphere is
                # an equiareal diffeomorphism.
                # So a uniform distribution on the authalic sphere gives
                # rise to a uniform distribution on the ellipsoid.
                beta0 = auth_lat(phi_min, e=self.e, radians=True)
                beta1 = auth_lat(phi_max, e=self.e, radians=True)
                beta = arcsin((sin(beta1) - sin(beta0)) * v + sin(beta0))
                phi = auth_lat(beta, e=self.e, radians=True, inverse=True)
            # Avoid the poles.
            if abs(phi) <= pi / 2 - delta:
                # Success.
                break
        if not self.radians:
            # Convert back to degrees.
            lam, phi = rad2deg([lam, phi])
        return lam, phi

    def lattice(self, n=90):
        """
        Return a 2n x n square lattice of longitude-latitude points.

        EXAMPLES::

            >>> E = UNIT_SPHERE
            >>> for p in E.lattice(n=3):
            ...     print(p)
            (-150.0, -60.0)
            (-150.0, 0.0)
            (-150.0, 60.0)
            (-90.0, -60.0)
            (-90.0, 0.0)
            (-90.0, 60.0)
            (-30.0, -60.0)
            (-30.0, 0.0)
            (-30.0, 60.0)
            (30.0, -60.0)
            (30.0, 0.0)
            (30.0, 60.0)
            (90.0, -60.0)
            (90.0, 0.0)
            (90.0, 60.0)
            (150.0, -60.0)
            (150.0, 0.0)
            (150.0, 60.0)

        """
        PI = self.pi()
        # Longitudinal and latitudinal spacing between points.
        delta = PI / n
        return [
            (-PI + delta * (0.5 + i), -PI / 2 + delta * (0.5 + j))
            for i in range(2 * n)
            for j in range(n)
        ]

    def meridian(self, lam, n=200):
        """
        Return a list of `n` equispaced longitude-latitude
        points lying along the meridian of longitude `lam`.
        Avoid the poles.
        """
        PI = self.pi()
        delta = PI / n
        return [(lam, -PI / 2 + delta * (0.5 + i)) for i in range(n)]

    def parallel(self, phi, n=200):
        """
        Return a list of `2*n` equispaced longitude-latitude
        points lying along the parallel of latitude `phi`.
        """
        PI = self.pi()
        delta = PI / n
        return [(-PI + delta * (0.5 + i), phi) for i in range(2 * n)]

    def graticule(self, n=400, spacing=None):
        """
        Return a list of longitude-latitude points sampled from a
        longitude-latitude graticule on this ellipsoid with the given
        spacing between meridians and between parallels.
        The number of points on longitude and latitude per pi radians is `n`.
        The spacing should be specified in the angle units used for this
        ellipsoid.
        If `spacing=None`, then a default spacing of pi/16 radians will be set.

        EXAMPLES::

            >>> E = UNIT_SPHERE
            >>> print(len(E.graticule(n=400)))
            25600

        """
        PI = self.pi()
        result = []
        # delta = PI/n
        # Set default spacing.
        if spacing is None:
            spacing = PI / 16
        # Longitude lines.
        lam = -PI
        while lam < PI:
            # result.extend([(lam, -PI/2 + delta*(0.5 + i)) for i in range(n)])
            result.extend(self.meridian(lam, n))
            lam += spacing
        # Latitude lines. Avoid the poles.
        eps = PI / 360
        phi = -PI / 2 + eps
        while phi < PI / 2:
            # result.extend([(-PI + delta*(0.5 + i), phi) for i in range(2*n)])
            result.extend(self.parallel(phi, n))
            phi += spacing
        return result

    def get_points(self, filename):
        """
        Return a list of longitude-latitude points contained in
        the file with filename `filename`.
        Assume the file is a text file containing at most one
        longitude-latitude point per line with the coordinates separated by
        whitespace and angles given in degrees.
        """
        result = []
        for line in open(filename, "rb"):
            if line[0] not in ["-", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
                # Ignore line.
                continue
            else:
                # Split coordinate pair on whitespace.
                p = [float(x) for x in line.split()]
                result.append(p)
        if self.radians:
            # Convert to radians.
            result = [deg2rad(p) for p in result]
        return result

    def xyz(self, lam, phi):
        """
        Given a point on this ellipsoid with longitude-latitude coordinates
        `(lam, phi)`, return the point's 3D rectangular coordinates.

        EXAMPLES::

            >>> E = UNIT_SPHERE
            >>> print(my_round(E.xyz(0, 45), 15))
            (0.707106781186548, 0.0, 0.707106781186548)

        NOTES:: .. Issue #1 was ..
            (0.70710678118654802, 0.0, 0.70710678118654802)

        """
        a = self.a
        e = self.e
        if not self.radians:
            lam, phi = deg2rad([lam, phi])
        # Equals a iff e = 0 (sphere):
        N = a / sqrt(1 - e**2 * sin(phi) ** 2)
        return (
            N * cos(lam) * cos(phi),
            N * sin(lam) * cos(phi),
            N * (1 - e**2) * sin(phi),
        )


# Define some common ellipsoids.
WGS84_ELLIPSOID = Ellipsoid(a=WGS84_A, f=WGS84_F)
WGS84_ELLIPSOID_RADIANS = Ellipsoid(a=WGS84_A, f=WGS84_F, radians=True)
WGS84_ASPHERE = Ellipsoid(R=WGS84_R_A)
WGS84_ASPHERE_RADIANS = Ellipsoid(R=WGS84_R_A, radians=True)
EMR_SPHERE = Ellipsoid(R=R_EM)
EMR_SPHERE_RADIANS = Ellipsoid(R=R_EM, radians=True)
UNIT_SPHERE = Ellipsoid(R=1)
UNIT_SPHERE_RADIANS = Ellipsoid(R=1, radians=True)
