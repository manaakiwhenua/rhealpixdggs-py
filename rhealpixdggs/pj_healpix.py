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

from typing import cast

import numpy as np
from numpy import arcsin, array, deg2rad, floor, pi, rad2deg, sign, sin, sqrt

# my_round is doctest-only: the doctests use it from the module globals.
from rhealpixdggs.utils import (  # noqa: F401
    FloatArray,
    ProjectionFunction,
    _auth_lat_array,
    auth_lat,
    auth_rad,
    my_round,
)

_PHI0 = arcsin(2.0 / 3)
_IMAGE_EPS = 1e-10


def _as_float_arrays(u: object, v: object) -> tuple[FloatArray, FloatArray]:
    """
    Coerce `u` and `v` to float64 arrays broadcast to a common shape.
    """
    u_arr, v_arr = np.broadcast_arrays(
        np.asarray(u, dtype=np.float64), np.asarray(v, dtype=np.float64)
    )
    return u_arr, v_arr


def _raise_out_of_image(
    x: FloatArray, y: FloatArray, ok: np.ndarray, projection: str
) -> None:
    bad = np.flatnonzero(~ok)
    i = int(bad[0])
    raise ValueError(
        f"{bad.size} of {ok.size} points are not in the image of the "
        f"{projection} projection of the unit sphere; first offending point "
        f"at index {i}: ({x.ravel()[i]:.20f},{y.ravel()[i]:.20f})"
    )


def _cap_number(lam: FloatArray) -> FloatArray:
    """
    The polar cap (0-3) each longitude or planar x falls in, as
    ``healpix_sphere`` and ``healpix_sphere_inverse`` compute it: floor,
    clamped from above only.
    """
    cap = np.floor(2 * lam / pi + 2)
    return np.asarray(np.where(cap >= 4, 3.0, cap), dtype=np.float64)


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


def _healpix_sphere_array(
    lam: FloatArray, phi: FloatArray
) -> tuple[FloatArray, FloatArray]:
    """
    ``healpix_sphere`` applied elementwise to arrays.
    """
    lam, phi = _as_float_arrays(lam, phi)
    x = np.empty_like(lam)
    y = np.empty_like(lam)
    s = np.sin(phi)
    equatorial = np.abs(phi) <= _PHI0
    x[equatorial] = lam[equatorial]
    y[equatorial] = 3 * pi / 8 * s[equatorial]
    polar = ~equatorial
    sigma = np.sqrt(3 * (1 - np.abs(s[polar])))
    lamc = -3 * pi / 4 + (pi / 2) * _cap_number(lam[polar])
    x[polar] = lamc + (lam[polar] - lamc) * sigma
    y[polar] = np.sign(phi[polar]) * pi / 4 * (2 - sigma)
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
    if not in_healpix_image(x, y):
        raise ValueError(
            f"({x:.20f},{y:.20f}) is not a point in the image of the HEALPix "
            "projection of the unit sphere"
        )
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


def _healpix_sphere_inverse_array(
    x: FloatArray, y: FloatArray
) -> tuple[FloatArray, FloatArray]:
    """
    ``healpix_sphere_inverse`` applied elementwise to arrays. Raises
    ``ValueError`` if any point lies outside the image.
    """
    x, y = _as_float_arrays(x, y)
    ok = _in_healpix_image_array(x, y)
    if not ok.all():
        _raise_out_of_image(x, y, ok, "HEALPix")
    lam = np.empty_like(x)
    phi = np.empty_like(x)
    abs_y = np.abs(y)
    equatorial = abs_y <= pi / 4
    pole = abs_y >= pi / 2
    polar = ~equatorial & ~pole
    lam[equatorial] = x[equatorial]
    phi[equatorial] = np.arcsin(8 * y[equatorial] / (3 * pi))
    xc = -3 * pi / 4 + (pi / 2) * _cap_number(x[polar])
    tau = 2 - 4 * abs_y[polar] / pi
    lam[polar] = np.clip(xc + (x[polar] - xc) / tau, -pi, pi)
    phi[polar] = np.sign(y[polar]) * np.arcsin(1 - np.float_power(tau, 2) / 3)
    lam[pole] = -pi
    phi[pole] = np.sign(y[pole]) * pi / 2
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

        >>> print(tuple(x if type(x) is int else x.tolist() for x in my_round(healpix_ellipsoid(0, pi/7), 15)))
        (0, 0.511157237746422)
        >>> print(tuple(x if type(x) is int else x.tolist() for x in my_round(healpix_ellipsoid(0, pi/7, e=0.8), 15)))
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
        >>> print(tuple(x if type(x) is int else x.tolist() for x in my_round(healpix_ellipsoid_inverse(*q), 15)))
        (0, 0.448798950512828)
        >>> print(my_round(p, 15))
        (0, 0.448798950512828)

    """
    # healpix_sphere_inverse() raises ValueError for out-of-bounds (x, y);
    # no need to duplicate that check here.
    lam, beta = healpix_sphere_inverse(x, y)
    phi = auth_lat(beta, e, radians=True, inverse=True)
    return lam, phi


def _healpix_ellipsoid_array(
    lam: FloatArray, phi: FloatArray, e: float = 0
) -> tuple[FloatArray, FloatArray]:
    """
    ``healpix_ellipsoid`` applied elementwise to arrays.
    """
    beta = _auth_lat_array(np.asarray(phi, dtype=np.float64), e, radians=True)
    return _healpix_sphere_array(lam, beta)


def _healpix_ellipsoid_inverse_array(
    x: FloatArray, y: FloatArray, e: float = 0
) -> tuple[FloatArray, FloatArray]:
    """
    ``healpix_ellipsoid_inverse`` applied elementwise to arrays.
    """
    lam, beta = _healpix_sphere_inverse_array(x, y)
    phi = _auth_lat_array(beta, e, radians=True, inverse=True)
    return lam, phi


def in_healpix_image(x: float, y: float) -> bool:
    """
    Return True if and only if `(x, y)` lies in the image of the HEALPix
    projection of the unit sphere, allowing a margin of about 1e-10 outside
    it so that points computed to lie on the boundary are not rejected for
    rounding error.

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
    # Fuzz to slightly expand the image so that points on its boundary
    # count as lying in it.
    eps = 1e-10
    abs_y = abs(y)
    if not (abs(x) < pi + eps and abs_y < pi / 2 + eps):
        return False
    if abs_y < pi / 4 + eps:
        return True
    # Polar triangle with apex (x_c, +-pi/2) and base corners (x_c +- pi/4,
    # +-pi/4): |x - x_c| + |y| <= pi/2.
    cap_number = min(max(int(2 * x / pi + 2), 0), 3)
    x_c = -3 * pi / 4 + (pi / 2) * cap_number
    return abs(x - x_c) + abs_y < pi / 2 + eps


def _in_healpix_image_array(x: FloatArray, y: FloatArray) -> np.ndarray:
    """
    ``in_healpix_image`` applied elementwise to arrays.
    """
    x, y = _as_float_arrays(x, y)
    eps = _IMAGE_EPS
    abs_y = np.abs(y)
    band = (np.abs(x) < pi + eps) & (abs_y < pi / 2 + eps)
    equatorial = abs_y < pi / 4 + eps
    cap_number = np.clip(np.trunc(2 * x / pi + 2), 0, 3)
    x_c = -3 * pi / 4 + (pi / 2) * cap_number
    triangle = np.abs(x - x_c) + abs_y < pi / 2 + eps
    return np.asarray(band & (equatorial | triangle))


def healpix(a: float = 1, e: float = 0) -> ProjectionFunction:
    """
    Return a function object that wraps the HEALPix projection and its inverse
    of an ellipsoid with major radius `a` and eccentricity `e`.

    EXAMPLES::

        >>> f = healpix(a=2, e=0)
        >>> print(tuple(x.tolist() for x in my_round(f(0, pi/3, radians=True), 15)))
        (0.574951359778215, 2.145747686573111)
        >>> p = (0, 60)
        >>> q = f(*p, radians=False); print(tuple(x.tolist() for x in my_round(q, 15)))
        (0.574951359778215, 2.145747686573111)
        >>> print(tuple(x.tolist() for x in my_round(f(*q, radians=False, inverse=True), 15)))
        (6e-15, 59.999999999999986)
        >>> print(my_round(p, 15))
        (0, 60)

    OUTPUT:

    - A function object of the form f(u, v, radians=False, inverse=False).
      `u` and `v` may be floats or numpy arrays of a common shape; arrays
      come back as a pair of float64 arrays.
    """
    R_A = auth_rad(a, e)

    def f_array(
        u: float | FloatArray, v: float | FloatArray, radians: bool, inverse: bool
    ) -> tuple[FloatArray, FloatArray]:
        u, v = _as_float_arrays(u, v)
        if not inverse:
            if not radians:
                u, v = deg2rad(u), deg2rad(v)
            x, y = _healpix_ellipsoid_array(u, v, e=e)
            return R_A * x, R_A * y
        lam, phi = _healpix_ellipsoid_inverse_array(u / R_A, v / R_A, e=e)
        if not radians:
            lam, phi = rad2deg(lam), rad2deg(phi)
        return lam, phi

    def f(
        u: float | FloatArray,
        v: float | FloatArray,
        radians: bool = False,
        inverse: bool = False,
    ) -> tuple[float, float] | tuple[FloatArray, FloatArray]:
        if isinstance(u, np.ndarray) or isinstance(v, np.ndarray):
            return f_array(u, v, radians, inverse)
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

    return cast(ProjectionFunction, f)
