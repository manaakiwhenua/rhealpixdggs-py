"""
This Python 3.11 module implements several helper functions for coding map projections.

- Alexander Raichev (AR), 2012-01-26: Refactored code from release 0.3.

NOTE:

All lengths are measured in meters and all angles are measured in radians
unless indicated otherwise.
"""

# *****************************************************************************
#       Copyright (C) 2012 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from math import asin, copysign, log, pi, sin, sqrt
from typing import Any, Protocol, overload

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


class ProjectionFunction(Protocol):
    """
    The callable a homemade projection factory returns: ``f(u, v,
    radians=False, inverse=False)`` mapping a pair of floats to a pair of
    floats, or a pair of equal-shape float arrays to a pair of arrays.
    """

    @overload
    def __call__(
        self, u: float, v: float, radians: bool = ..., inverse: bool = ...
    ) -> tuple[float, float]: ...

    @overload
    def __call__(
        self, u: FloatArray, v: FloatArray, radians: bool = ..., inverse: bool = ...
    ) -> tuple[FloatArray, FloatArray]: ...


def my_round(x: Any, digits: int = 0) -> Any:
    """
    Round the floating point number or list/tuple of floating point
    numbers to ``digits`` number of digits.
    Calls Python's ``round()`` function.

    EXAMPLES::

        >>> print(my_round(1./7, 6))
        0.142857
        >>> print(my_round((1./3, 1./7), 6))
        (0.333333, 0.142857)

    """
    try:
        result = round(x, digits)
    except TypeError:
        result = [my_round(xx, digits) for xx in x]
        if isinstance(x, tuple):
            result = tuple(result)
    return result


def wrap_longitude(lam: float, radians: bool = False) -> float:
    """
    Given a point p on the unit circle at angle `lam` from the positive
    x-axis, return its angle theta in the range -pi <= theta < pi.
    If `radians` = True, then `lam` and the output are given in radians.
    Otherwise, they are given in degrees.

    EXAMPLES::

        >>> wrap_longitude(2*pi + pi, radians=True)
        -3.141592653589793

    """
    if radians:
        half_range = pi
    else:
        half_range = 180

    if lam < -half_range or lam >= half_range:
        result = lam % (2 * half_range)
        if result >= half_range:
            result = result - 2 * half_range
    else:
        result = lam
    return result


def wrap_latitude(phi: float, radians: bool = False) -> float:
    """
    Given a point p on the unit circle at angle `phi` from the positive x-axis,
    if p lies in the right half of the circle, then return its angle that lies
    in the interval [-pi/2, pi/2].
    If p lies in the left half of the circle, then reflect it through the
    origin, and return the angle of the reflected point that lies in the
    interval [-pi/2, pi/2].
    If `radians` = True, then `phi` and the output are given in radians.
    Otherwise, they are given in degrees.

    EXAMPLES::

        >>> wrap_latitude(45.0, radians=False)
        45.0
        >>> wrap_latitude(-45.0, radians=False)
        -45.0
        >>> wrap_latitude(90.0, radians=False)
        90.0
        >>> wrap_latitude(-90.0, radians=False)
        -90.0
        >>> wrap_latitude(135.0, radians=False)
        -45.0
        >>> wrap_latitude(-135.0, radians=False)
        45.0

    """
    # Put phi in range -pi <= phi < pi.
    phi = wrap_longitude(phi, radians=radians)

    if radians:
        half_range = pi
    else:
        half_range = 180

    if abs(phi) <= half_range / 2:
        result = phi
    else:
        result = phi - copysign(half_range, phi)
    return result


def _wrap_longitude_array(lam: FloatArray, radians: bool = False) -> FloatArray:
    """
    ``wrap_longitude`` applied elementwise to an array.
    """
    lam = np.asarray(lam, dtype=np.float64)
    half_range = pi if radians else 180
    wrapped = lam % (2 * half_range)
    wrapped = np.where(wrapped >= half_range, wrapped - 2 * half_range, wrapped)
    out_of_range = (lam < -half_range) | (lam >= half_range)
    return np.asarray(np.where(out_of_range, wrapped, lam), dtype=np.float64)


def _wrap_latitude_array(phi: FloatArray, radians: bool = False) -> FloatArray:
    """
    ``wrap_latitude`` applied elementwise to an array.
    """
    phi = _wrap_longitude_array(phi, radians=radians)
    half_range = pi if radians else 180
    reflected = phi - np.copysign(half_range, phi)
    return np.asarray(
        np.where(np.abs(phi) <= half_range / 2, phi, reflected), dtype=np.float64
    )


def _auth_lat_coefficients(
    n: float, inverse: bool
) -> tuple[float, float, float, float, float, float]:
    """
    The power-series coefficients of sin(2 phi), ..., sin(12 phi) used by
    ``auth_lat`` for third flattening `n`: equation A19 (forward) or A20
    (inverse) of https://doi.org/10.48550/arXiv.2212.05818.
    """
    if not inverse:
        return (
            n
            * (
                -4 / 3
                + n
                * (
                    -4 / 45
                    + n
                    * (
                        88 / 315
                        + n
                        * (538 / 4725 + n * (20824 / 467775 + n * (-44732 / 2837835)))
                    )
                )
            ),
            n
            * (
                n
                * (
                    34 / 45
                    + n
                    * (
                        8 / 105
                        + n
                        * (
                            -2482 / 14175
                            + n * (-37192 / 467775 + n * (-12467764 / 212837625))
                        )
                    )
                )
            ),
            n
            * (
                n
                * (
                    n
                    * (
                        -1532 / 2835
                        + n
                        * (
                            -898 / 14175
                            + n * (54968 / 467775 + n * 100320856 / 1915538625)
                        )
                    )
                )
            ),
            n
            * (
                n
                * (
                    n
                    * (
                        n
                        * (
                            6007 / 14175
                            + n * (24496 / 467775 + n * (-5884124 / 70945875))
                        )
                    )
                )
            ),
            n * (n * (n * (n * (n * (-23356 / 66825 + n * (-839792 / 19348875)))))),
            n * (n * (n * (n * (n * (n * 570284222 / 1915538625))))),
        )
    return (
        n
        * (
            4 / 3
            + n
            * (
                4 / 45
                + n
                * (
                    -16 / 35
                    + n
                    * (-2582 / 14175 + n * (60136 / 467775 + n * 28112932 / 212837625))
                )
            )
        ),
        n
        * (
            n
            * (
                46 / 45
                + n
                * (
                    152 / 945
                    + n
                    * (
                        -11966 / 14175
                        + n * (-21016 / 51975 + n * 251310128 / 638512875)
                    )
                )
            )
        ),
        n
        * (
            n
            * (
                n
                * (
                    3044 / 2835
                    + n
                    * (3802 / 14175 + n * (-94388 / 66825 + n * (-8797648 / 10945935)))
                )
            )
        ),
        n
        * (
            n
            * (
                n
                * (
                    n
                    * (
                        6059 / 4725
                        + n * (41072 / 93555 + n * (-1472637812 / 638512875))
                    )
                )
            )
        ),
        n * (n * (n * (n * (n * (768272 / 467775 + n * 455935736 / 638512875))))),
        n * (n * (n * (n * (n * (n * 4210684958 / 1915538625))))),
    )


def auth_lat(
    phi: float, e: float, inverse: bool = False, radians: bool = False
) -> float:
    """
    Given a point of geographic latitude `phi` on an ellipse of
    eccentricity `e`, return the authalic latitude of the point.
    If `inverse` =True, then compute its inverse approximately.

    EXAMPLES::

        >>> beta = auth_lat(pi/3, 0.08181919104281579, radians=True)
        >>> print(my_round(beta, 15))
        1.045256493205824

        >>> print(my_round(auth_lat(beta, 0.08181919104281579, radians=True, inverse=True), 15))
        1.047197551196598

        >>> print(my_round(pi/3, 15))
        1.047197551196598

    NOTES:

    For small flattenings f (f < 1/150), when calculating
    authalic from common latitude, power series approximation
    (from https://doi.org/10.48550/arXiv.2212.05818) gives
    more accurate results than direct formula (for
    double-precission accuracy). For the inverse, again power
    series approximation is used, which is standard in
    cartography for small flattenings. The one used in this
    case is from https://doi.org/10.48550/arXiv.2212.05818
    """
    if e == 0:
        return phi
    # Compute flattening f and third flattening n from eccentricity e.
    f = 1 - sqrt(1 - e**2)
    n = (1 - sqrt(1 - e**2)) / (1 + sqrt(1 - e**2))

    if not inverse:
        # Compute authalic latitude from common latitude phi.
        # For large flattenings (f > 1/150) use direct formula,
        # for small flattenings (f <= 1/150) use power series.
        if abs(f) > 1 / 150:
            # Use direct formula for large flattenings.
            if not radians:
                # Convert to radians to do calculations below.
                phi = phi * pi / 180
            # Compute authalic latitude from latitude phi.
            q = ((1 - e**2) * sin(phi)) / (1 - (e * sin(phi)) ** 2) - (1 - e**2) / (
                2.0 * e
            ) * log((1 - e * sin(phi)) / (1 + e * sin(phi)))
            qp = 1 - (1 - e**2) / (2.0 * e) * log((1.0 - e) / (1.0 + e))
            ratio = q / qp
            # Avoid rounding errors.
            if abs(ratio) > 1:
                # Make abs(ratio) = 1
                ratio = copysign(1, ratio)
            result = asin(ratio)
            if not radians:
                result = result * 180 / pi
            return result
        else:
            # Use power series approximation for small flattenings (f <= 1/150).
            # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A19)
            if not radians:
                phi = phi * pi / 180

            c2, c4, c6, c8, c10, c12 = _auth_lat_coefficients(n, inverse=False)
            authalic_lat = phi + (
                c2 * sin(2 * phi)
                + c4 * sin(4 * phi)
                + c6 * sin(6 * phi)
                + c8 * sin(8 * phi)
                + c10 * sin(10 * phi)
                + c12 * sin(12 * phi)
            )

            if not radians:
                authalic_lat = authalic_lat * 180 / pi

            return authalic_lat
    else:
        # Compute common latitude from authalic latitude phi.
        # Power series expansion taken from https://doi.org/10.48550/arXiv.2212.05818 (Equation A20)
        if not radians:
            phi = phi * pi / 180

        c2, c4, c6, c8, c10, c12 = _auth_lat_coefficients(n, inverse=True)
        common_lat = phi + (
            c2 * sin(2 * phi)
            + c4 * sin(4 * phi)
            + c6 * sin(6 * phi)
            + c8 * sin(8 * phi)
            + c10 * sin(10 * phi)
            + c12 * sin(12 * phi)
        )

        if not radians:
            common_lat = common_lat * 180 / pi

        return common_lat


def _auth_lat_array(
    phi: FloatArray, e: float, inverse: bool = False, radians: bool = False
) -> FloatArray:
    """
    ``auth_lat`` applied elementwise to an array. Every expression mirrors
    the scalar function's form so the two agree to the last bit wherever
    the platform's array and scalar transcendental kernels do.
    """
    phi = np.asarray(phi, dtype=np.float64)
    if e == 0:
        return phi.copy()
    f = 1 - sqrt(1 - e**2)
    n = (1 - sqrt(1 - e**2)) / (1 + sqrt(1 - e**2))
    if not radians:
        phi = phi * pi / 180
    if not inverse and abs(f) > 1 / 150:
        s = np.sin(phi)
        q = ((1 - e**2) * s) / (1 - np.float_power(e * s, 2)) - (1 - e**2) / (
            2.0 * e
        ) * np.log((1 - e * s) / (1 + e * s))
        qp = 1 - (1 - e**2) / (2.0 * e) * log((1.0 - e) / (1.0 + e))
        ratio = q / qp
        ratio = np.where(np.abs(ratio) > 1, np.copysign(1, ratio), ratio)
        result = np.arcsin(ratio)
    else:
        c2, c4, c6, c8, c10, c12 = _auth_lat_coefficients(n, inverse)
        result = phi + (
            c2 * np.sin(2 * phi)
            + c4 * np.sin(4 * phi)
            + c6 * np.sin(6 * phi)
            + c8 * np.sin(8 * phi)
            + c10 * np.sin(10 * phi)
            + c12 * np.sin(12 * phi)
        )
    if not radians:
        result = result * 180 / pi
    return np.asarray(result, dtype=np.float64)


def auth_rad(a: float, e: float, inverse: bool = False) -> float:
    """
    Return the radius of the authalic sphere of the ellipsoid with major
    radius `a` and eccentricity `e`.
    If `inverse` = True, then return the major radius of the ellipsoid
    with authalic radius `a` and eccentricity `e`.

    EXAMPLES::

        >>> auth_rad(1, 0)
        1
        >>> for i in range(2, 11):
        ...     e = 1.0/i**2
        ...     print(my_round((e, auth_rad(1, 1.0/i**2)), 15))
        (0.25, 0.989393259670095)
        (0.111111111111111, 0.997935147429943)
        (0.0625, 0.999348236455825)
        (0.04, 0.99973321235361)
        (0.027777777777778, 0.99987137105188)
        (0.020408163265306, 0.999930576285614)
        (0.015625, 0.999959307080847)
        (0.012345679012346, 0.999974596271211)
        (0.01, 0.999983332861089)

    """
    if e == 0:
        return a
    k = sqrt(0.5 * (1 - (1 - e**2) / (2 * e) * log((1 - e) / (1 + e))))
    if not inverse:
        # The expression below is undefined when e=0 (sphere),
        # but its limit as e tends to 0 is a, as expected.
        return a * k
    else:
        # Then a is the authalic radius and output major radius of ellipsoid.
        return a / k
