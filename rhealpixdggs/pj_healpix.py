"""
This Python 3.11 module implements the HEALPix map projection as described in [CaRo2007]_.

.. [CaRo2007] Mark R. Calabretta and Boudewijn F. Roukema, Mapping on the healpix grid, Monthly Notices of the Royal Astronomical Society 381 (2007), no. 2, 865--872.

- Alexander Raichev (AR), 2013-01-26: Refactored code from release 0.3.

NOTE:

All lengths are measured in meters and all angles are measured in radians
unless indicated otherwise.
By 'ellipsoid' below, I mean an oblate ellipsoid of revolution.
"""
# *****************************************************************************
#       Copyright (C) 2013 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import third-party modules.
from numpy import pi, floor, sqrt, sin, arcsin, sign, array, deg2rad, rad2deg
from typing import Callable

# Import my modules.
from rhealpixdggs.utils import my_round, auth_lat, auth_rad


def healpix_sphere(lam: float, phi: float) -> tuple[float, float]:
    """
    Compute the signature function of the HEALPix
    projection of the unit sphere.

    INPUT:

    - `lam, phi` - Geodetic longitude-latitude coordinates in radians.
      Assume -pi <= `lam` < pi and -pi/2 <= `phi` <= pi/2.

    EXAMPLES::

        >>> print(healpix_sphere(0, arcsin(2.0/3)) == (0, pi/4))
        True

    """
    phi0 = arcsin(2.0 / 3)
    # Equatorial region.
    if abs(phi) <= phi0:
        x = lam
        y = 3 * pi / 8 * sin(phi)
    # Polar region.
    else:
        sigma = sqrt(3 * (1 - abs(sin(phi))))
        cap_number = floor(2 * lam / pi + 2)
        if cap_number >= 4:
            # Rounding error
            cap_number = 3
        lamc = -3 * pi / 4 + (pi / 2) * cap_number
        x = lamc + (lam - lamc) * sigma
        y = sign(phi) * pi / 4 * (2 - sigma)
    return x, y


def healpix_sphere_inverse(x: float, y: float) -> tuple[float, float]:
    """
    Compute the inverse of the healpix_sphere().

    INPUT:

    - `x, y` - Planar coordinates in meters in the image of the
      HEALPix projection of the unit sphere.

    EXAMPLES::

        >>> print(healpix_sphere_inverse(0, pi/4) == (0, arcsin(2.0/3)))
        True

    """
    # Throw error if input coordinates are out of bounds.
    if not in_healpix_image(x, y):
        print("Error (hsi): input coordinates (%.20f,%.20f) are out of bounds" % (x, y))
        return float("inf"), float("inf")
    y0 = pi / 4
    # Equatorial region.
    if abs(y) <= y0:
        lam = x
        phi = arcsin(8 * y / (3 * pi))
    # Polar region but not the poles.
    elif abs(y) < pi / 2:
        cap_number = floor(2 * x / pi + 2)
        if cap_number >= 4:
            # Rounding error.
            cap_number = 3
        xc = -3 * pi / 4 + (pi / 2) * cap_number
        tau = 2 - 4 * abs(y) / pi
        lam = xc + (x - xc) / tau
        phi = sign(y) * arcsin(1 - tau**2 / 3)
        # Handle rounding errors in longitude.
        if lam < -pi:
            lam = -pi
        elif lam > pi:
            lam = pi
    # Poles.
    else:
        lam = -pi
        phi = sign(y) * pi / 2
    return lam, phi


def healpix_ellipsoid(lam: float, phi: float, e: float = 0) -> tuple[float, float]:
    """
    Compute the signature functions of the HEALPix projection of an oblate
    ellipsoid with eccentricity `e` whose authalic sphere is the unit sphere.
    Works when `e` = 0 (spherical case) too.

    INPUT:

    - `lam, phi` - Geodetic longitude-latitude coordinates in radians.
      Assume -pi <= `lam` < pi and -pi/2 <= `phi` <= pi/2.
    - `e` - Eccentricity of the oblate ellipsoid.

    EXAMPLES::

        >>> print(my_round(healpix_ellipsoid(0, pi/7), 15))
        (0, 0.511157237746422)
        >>> print(my_round(healpix_ellipsoid(0, pi/7, e=0.8), 15))
        (0, 0.268484450857837)

    """
    beta = auth_lat(phi, e, radians=True)
    return healpix_sphere(lam, beta)


def healpix_ellipsoid_inverse(x: float, y: float, e: float = 0) -> tuple[float, float]:
    """
    Compute the inverse of healpix_ellipsoid().

    EXAMPLES::

        >>> p = (0, pi/7)
        >>> q = healpix_ellipsoid(*p)
        >>> print(my_round(healpix_ellipsoid_inverse(*q), 15))
        (0, 0.448798950512828)
        >>> print(my_round(p, 15))
        (0, 0.448798950512828)

    """
    # Throw error if input coordinates are out of bounds.
    if not in_healpix_image(x, y):
        print("Error (hei): input coordinates (%.20f,%.20f) are out of bounds" % (x, y))
        return

    lam, beta = healpix_sphere_inverse(x, y)
    phi = auth_lat(beta, e, radians=True, inverse=True)
    return lam, phi


def in_healpix_image(x: float, y: float) -> bool:
    """
    Return True if and only if `(x, y)` lies in the image of the HEALPix
    projection of the unit sphere.

    EXAMPLES::

        >>> eps = 0     # Test boundary points.
        >>> hp = [
        ... (-pi - eps, pi/4),
        ... (-3*pi/4, pi/2 + eps),
        ... (-pi/2, pi/4 + eps),
        ... (-pi/4, pi/2 + eps),
        ... (0, pi/4 + eps),
        ... (pi/4, pi/2 + eps),
        ... (pi/2, pi/4 + eps),
        ... (3*pi/4, pi/2 + eps),
        ... (pi + eps, pi/4),
        ... (pi + eps,-pi/4),
        ... (3*pi/4,-pi/2 - eps),
        ... (pi/2,-pi/4 - eps),
        ... (pi/4,-pi/2 - eps),
        ... (0,-pi/4 - eps),
        ... (-pi/4,-pi/2 - eps),
        ... (-pi/2,-pi/4 - eps),
        ... (-3*pi/4,-pi/2 - eps),
        ... (-pi - eps,-pi/4)
        ... ]
        >>> for p in hp:
        ...     if not in_healpix_image(*p):
        ...             print('Fail')
        ...
        >>> in_healpix_image(0, 0)
        True
        >>> in_healpix_image(0, pi/4 + 0.1)
        False

    """
    # matplotlib is a third-party module.
    from matplotlib.path import Path

    # Fuzz to slightly expand HEALPix image boundary so that
    # points on the boundary count as lying in the image.
    eps = 1e-10
    vertices = [
        (-pi - eps, pi / 4 + eps),
        (-3 * pi / 4, pi / 2 + eps),
        (-pi / 2, pi / 4 + eps),
        (-pi / 4, pi / 2 + eps),
        (0, pi / 4 + eps),
        (pi / 4, pi / 2 + eps),
        (pi / 2, pi / 4 + eps),
        (3 * pi / 4, pi / 2 + eps),
        (pi + eps, pi / 4 + eps),
        (pi + eps, -pi / 4 - eps),
        (3 * pi / 4, -pi / 2 - eps),
        (pi / 2, -pi / 4 - eps),
        (pi / 4, -pi / 2 - eps),
        (0, -pi / 4 - eps),
        (-pi / 4, -pi / 2 - eps),
        (-pi / 2, -pi / 4 - eps),
        (-3 * pi / 4, -pi / 2 - eps),
        (-pi - eps, -pi / 4 - eps),
    ]
    poly = Path(vertices)
    return bool(poly.contains_point([x, y]))


def healpix_vertices() -> list[tuple[float, float, float]]:
    """
    Return a list of the planar vertices of the HEALPix projection of
    the unit sphere.
    """
    return [
        (pi, pi / 4),
        (3 * pi / 4, pi / 2),
        (pi / 2, pi / 4),
        (pi / 4, pi / 2),
        (0, pi / 4),
        (-pi / 4, pi / 2),
        (-pi / 2, pi / 4),
        (-3 * pi / 4, pi / 2),
        (-pi, pi / 4),
        (-pi, -pi / 4),
        (-3 * pi / 4, -pi / 2),
        (-pi / 2, -pi / 4),
        (-pi / 4, -pi / 2),
        (0, -pi / 4),
        (pi / 4, -pi / 2),
        (pi / 2, -pi / 4),
        (3 * pi / 4, -pi / 2),
        (pi, -pi / 4),
    ]


def healpix(
    a: float = 1, e: float = 0
) -> Callable[[float, float, bool, bool], tuple[float, float]]:
    """
    Return a function object that wraps the HEALPix projection and its inverse
    of an ellipsoid with major radius `a` and eccentricity `e`.

    EXAMPLES::

        >>> f = healpix(a=2, e=0)
        >>> print(my_round(f(0, pi/3, radians=True), 15))
        (0.574951359778215, 2.145747686573111)
        >>> p = (0, 60)
        >>> q = f(*p, radians=False); print(my_round(q, 15))
        (0.574951359778215, 2.145747686573111)
        >>> print(my_round(f(*q, radians=False, inverse=True), 15))
        (6e-15, 59.999999999999986)
        >>> print(my_round(p, 15))
        (0, 60)

    OUTPUT:

    - A function object of the form f(u, v, radians=False, inverse=False).
    """
    R_A = auth_rad(a, e)

    def f(
        u: float, v: float, radians: bool = False, inverse: bool = False
    ) -> tuple[float, float]:
        if not inverse:
            lam, phi = u, v
            if not radians:
                # Convert to radians.
                lam, phi = deg2rad([lam, phi])
            return tuple(R_A * array(healpix_ellipsoid(lam, phi, e=e)))
        else:
            # Scale down to R_A = 1.
            x, y = array((u, v)) / R_A
            lam, phi = array(healpix_ellipsoid_inverse(x, y, e=e))
            if not radians:
                # Convert to degrees.
                lam, phi = rad2deg([lam, phi])
            return lam, phi

    return f
