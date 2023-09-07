"""
This Python 3.11 module implements the rHEALPix map projection.

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
from numpy import pi, sign, array, identity, dot, deg2rad, rad2deg
from typing import Callable

# Import my modules.
from rhealpixdggs.pj_healpix import (
    healpix_sphere,
    healpix_sphere_inverse,
    healpix_ellipsoid,
    healpix_ellipsoid_inverse,
)
from rhealpixdggs.utils import my_round, auth_rad

# Matrix for anticlockwise rotation by pi/2:
ROTATE1 = array([[0, -1], [1, 0]])
# Matrix for anticlockwise rotation by pi:
ROTATE2 = dot(ROTATE1, ROTATE1)
# Matrix for anticlockwise rotation by 3*pi/2.
ROTATE3 = dot(ROTATE2, ROTATE1)
# Dictionary of all powers of ROTATE1 and its inverse (ROTATE3).
ROTATE = {
    0: identity(2, int),
    1: ROTATE1,
    2: ROTATE2,
    3: ROTATE3,
    -1: ROTATE3,
    -2: ROTATE2,
    -3: ROTATE1,
}


def combine_triangles(
    x: float,
    y: float,
    north_square: int = 0,
    south_square: int = 0,
    inverse: bool = False,
) -> tuple[float, float]:
    """
    Rearrange point `(x, y)` in the HEALPix projection by
    combining the polar triangles into two polar squares.
    Put the north polar square in position `north_square` and
    the south polar square in position `south_square`.
    If `inverse` = True, uncombine the polar triangles.

    INPUT:

    - `x, y` - Coordinates in the HEALPix projection of the unit sphere.
    - `north_square, south_square` - Integers between 0 and 3 indicating
      the positions of the north_square polar square and south_square polar
      square respectively.
      See rhealpix_sphere() docstring for a diagram.
    - `inverse` - (Optional; default = False) Boolean. If False, then compute
      forward function. If True, then compute inverse function.

    EXAMPLES::

        >>> u, v = -pi/4, pi/3
        >>> x, y = combine_triangles(u, v)
        >>> print(my_round((x, y), 15))
        (-1.832595714594046, 1.570796326794896)
        >>> print(my_round(combine_triangles(x, y, inverse=True), 15))
        (-0.785398163397448, 1.047197551196598)
        >>> print(my_round((u, v), 15))
        (-0.785398163397448, 1.047197551196598)

    """
    # Ensure north_square and south_square lie in {0, 1, 2, 3}.
    north_square = north_square % 4
    south_square = south_square % 4
    c, region = triangle(
        x, y, north_square=north_square, south_square=south_square, inverse=inverse
    )
    if region == "equatorial":
        # (x,y) remains fixed
        return x, y
    xy = array((x, y))
    tc = array((-3 * pi / 4 + c * pi / 2, sign(y) * pi / 2))
    if not inverse:
        # Forward function.
        # Rotate (x, y) about tc and then translate it to
        # the tip u of the polar triangle it will be assembled upon.
        if region == "north_polar":
            u = array((-3 * pi / 4 + north_square * pi / 2, pi / 2))
            x, y = dot(ROTATE[c - north_square], xy - tc) + u
        elif region == "south_polar":
            u = array((-3 * pi / 4 + south_square * pi / 2, -pi / 2))
            x, y = dot(ROTATE[-(c - south_square)], xy - tc) + u
    else:
        # Inverse function.
        # Unrotate (x, y) about u and then translate it to tc.
        if region == "north_polar":
            u = array((-3 * pi / 4 + north_square * pi / 2, pi / 2))
            x, y = dot(ROTATE[-(c - north_square)], xy - u) + tc
        elif region == "south_polar":
            u = array((-3 * pi / 4 + south_square * pi / 2, -pi / 2))
            x, y = dot(ROTATE[c - south_square], xy - u) + tc
    return x, y


def triangle(
    x: float,
    y: float,
    north_square: int = 0,
    south_square: int = 0,
    inverse: bool = False,
) -> tuple[int, str]:
    """
    Return the number of the polar triangle and region that `(x, y)` lies in.
    If `inverse` = False, then assume `(x,y)` lies in the image of the HEALPix
    projection of the unit sphere.
    If `inverse` = True, then assume `(x,y)` lies in the image of the
    `(north_square, south_square)`-rHEALPix projection of the unit sphere.

    INPUT:

    - `x, y` - Coordinates in the HEALPix or rHEALPix (if `inverse` = True)
      projection of the unit sphere.
    - `north_square, south_square` - Integers between 0 and 3 indicating the
      positions of the north_square pole square and south_square pole square
      respectively.
      See rhealpix_sphere() docstring for a diagram.
    - `inverse` - (Optional; default = False) Boolean. If False, then compute
      forward function. If True, then compute inverse function.

    OUTPUT:

    The pair (triangle_number, region).
    Here region equals 'north_polar' (polar), 'south_polar' (polar), or
    'equatorial', indicating where `(x, y)` lies.
    If region = 'equatorial', then triangle_number = None.
    Suppose now that region != 'equatorial'.
    If `inverse` = False, then triangle_number is the number (0, 1, 2, or 3) of
    the HEALPix polar triangle Z that `(x, y)` lies in.
    If `inverse` = True, then triangle_number is the number (0, 1, 2, or 3) of
    the HEALPix polar triangle that `(x, y)` will get moved into.

    EXAMPLES::

        >>> triangle(-pi/4, pi/4 + 0.1)
        (1, 'north_polar')
        >>> triangle(-3*pi/4 + 0.1, pi/2, inverse=True)
        (1, 'north_polar')

    NOTES:

    In the HEALPix projection, the polar triangles are labeled 0--3 from
    east to west like this::

            *       *       *       *
          * 0 *   * 1 *   * 2 *   * 3 *
        *-------*-------*-------*-------*
        |       |       |       |       |
        |       |       |       |       |
        |       |       |       |       |
        *-------*-------*-------*-------*
          * 0 *   * 1 *   * 2 *   * 3 *
            *       *       *       *

    In the rHEALPix projection these polar triangles get rearranged
    into a square with the triangles numbered `north_square` and `south_square`
    remaining fixed.
    For example, if `north_square` = 1 and `south_square` = 3,
    then the triangles get rearranged this way::

        North polar square:     *-------*
                                | * 3 * |
                                | 0 * 2 |
                                | * 1 * |
                            ----*-------*----

        South polar square: ----*-------*----
                                | * 3 * |
                                | 2 * 0 |
                                | * 1 * |
                                *-------*

    """
    if not inverse:
        # Forward function.
        # Find the region (x, y) lies in.
        if y > pi / 4:
            region = "north_polar"
        elif y < -pi / 4:
            region = "south_polar"
        else:
            region = "equatorial"
        # Find the triangle number of (x, y) in the image of the
        # HEALPix projection.
        if region == "equatorial":
            triangle_number = None
        else:
            if x < -pi / 2:
                triangle_number = 0
            elif x >= -pi / 2 and x < 0:
                triangle_number = 1
            elif x >= 0 and x < pi / 2:
                triangle_number = 2
            else:
                triangle_number = 3
    else:
        # Inverse function.
        # Find the region (x, y) lies in.
        if y > pi / 4:
            region = "north_polar"
        elif y < -pi / 4:
            region = "south_polar"
        else:
            region = "equatorial"
        # Find HEALPix polar triangle number that (x, y) moves to when
        # the rHEALPix polar square is disassembled.
        eps = 1e-15  # Fuzz to avoid some rounding errors.
        if region == "equatorial":
            triangle_number = None
        elif region == "north_polar":
            L1 = x - (-3 * pi / 4 + (north_square - 1) * pi / 2)
            L2 = -x + (-3 * pi / 4 + (north_square + 1) * pi / 2)
            if y < L1 - eps and y >= L2 - eps:
                triangle_number = (north_square + 1) % 4
            elif y >= L1 - eps and y > L2 + eps:
                triangle_number = (north_square + 2) % 4
            elif y > L1 + eps and y <= L2 + eps:
                triangle_number = (north_square + 3) % 4
            else:
                triangle_number = north_square
        else:
            # region == 'south_polar':
            L1 = x - (-3 * pi / 4 + (south_square + 1) * pi / 2)
            L2 = -x + (-3 * pi / 4 + (south_square - 1) * pi / 2)
            if y <= L1 + eps and y > L2 + eps:
                triangle_number = (south_square + 1) % 4
            elif y < L1 - eps and y <= L2 + eps:
                triangle_number = (south_square + 2) % 4
            elif y >= L1 - eps and y < L2 - eps:
                triangle_number = (south_square + 3) % 4
            else:
                triangle_number = south_square
    return triangle_number, region


def rhealpix_sphere(
    lam: float,
    phi: float,
    north_square: int = 0,
    south_square: int = 0,
    region: str = "none",
) -> tuple[float, float]:
    """
    Compute the signature functions of the rHEALPix map projection of
    the unit sphere.
    The north polar square is put in position `north_square`, and the
    south polar square is put in position `south_square`.

    INPUT:

    - `lam, phi` -Geographic longitude-latitude coordinates in radians.
      Assume -pi <= `lam` < pi and -pi/2 <= `phi` <= pi/2.
    - `north_square, south_square` - (Optional; defaults = 0, 0) Integers
      between 0 and 3 indicating positions of north polar and
      south polar squares, respectively.

    EXAMPLES::

        >>> print(my_round(rhealpix_sphere(0, pi/4), 15))
        (-1.619978633413937, 2.307012183573304)

    NOTE:

    The polar squares are labeled 0, 1, 2, 3 from east to west like this::

        east         west
        *---*---*---*---*
        | 0 | 1 | 2 | 3 |
        *---*---*---*---*
        |   |   |   |   |
        *---*---*---*---*
        | 0 | 1 | 2 | 3 |
        *---*---*---*---*
    """
    x, y = healpix_sphere(lam, phi)
    return combine_triangles(x, y, north_square=north_square, south_square=south_square)


def rhealpix_sphere_inverse(
    x: float,
    y: float,
    north_square: int = 0,
    south_square: int = 0,
    region: str = "none",
) -> tuple[float, float]:
    """
    Compute the inverse of rhealpix_sphere().

    EXAMPLES::

        >>> p = (0, pi/4)
        >>> q = rhealpix_sphere(*p)
        >>> print(my_round(rhealpix_sphere_inverse(*q), 15))
        (0.0, 0.785398163397448)
        >>> print(my_round(p, 15))
        (0, 0.785398163397448)

    """
    # Throw error if input coordinates are out of bounds.
    if not in_rhealpix_image(
        x, y, south_square=south_square, north_square=north_square
    ):
        print("Error (rsi): input coordinates (%.20f,%.20f) are out of bounds" % (x, y))
        return
    if region != "equatorial":
        x, y = combine_triangles(
            x, y, north_square=north_square, south_square=south_square, inverse=True
        )
    return healpix_sphere_inverse(x, y)


def rhealpix_ellipsoid(
    lam: float,
    phi: float,
    e: float = 0,
    north_square: int = 0,
    south_square: int = 0,
    region: str = "none",
) -> tuple[float, float]:
    """
    Compute the signature functions of the rHEALPix map
    projection of an oblate ellipsoid with eccentricity `e` whose
    authalic sphere is the unit sphere.
    The north polar square is put in position `north_square`,
    and the south polar square is put in position `south_square`.
    Works when `e` = 0 (spherical case) too.

    INPUT:

    - `lam, phi` - Geographic longitude-latitude coordinates in radian.
      Assume -pi <= `lam` < pi and -pi/2 <= `phi` <= pi/2.
    - `e` - Eccentricity of the ellipsoid.
    - `north_square, south_square` - (Optional; defaults = 0, 0) Integers
      between 0 and 3 indicating positions of north polar and
      south polar squares, respectively.
      See rhealpix_sphere() docstring for a diagram.

    EXAMPLES::

        >>> from numpy import arcsin
        >>> print(my_round(rhealpix_ellipsoid(0, arcsin(2.0/3)), 15))
        (0, 0.785398163397448)

    """
    # Ensure north_square and south_square lie in {0, 1,2, 3}.
    x, y = healpix_ellipsoid(lam, phi, e)
    if region != "equatorial":
        x, y = combine_triangles(
            x, y, north_square=north_square, south_square=south_square
        )
    return x, y


def rhealpix_ellipsoid_inverse(
    x: float,
    y: float,
    e: float = 0,
    north_square: int = 0,
    south_square: int = 0,
    region: str = "none",
) -> tuple[float, float]:
    """
    Compute the inverse of rhealpix_ellipsoid().

    EXAMPLES::

        >>> p = (0, pi/4)
        >>> q = rhealpix_ellipsoid(*p)
        >>> print(my_round(rhealpix_ellipsoid_inverse(*q), 15))
        (0.0, 0.785398163397448)
        >>> print(my_round(p, 15))
        (0, 0.785398163397448)

    """
    # Throw error if input coordinates are out of bounds.
    if not in_rhealpix_image(
        x, y, south_square=south_square, north_square=north_square
    ):
        print("Error (rei): input coordinates (%.20f,%.20f) are out of bounds" % (x, y))
        return

    if region != "equatorial":
        x, y = combine_triangles(
            x, y, north_square=north_square, south_square=south_square, inverse=True
        )

    return healpix_ellipsoid_inverse(x, y, e=e)


def in_rhealpix_image(
    x: float, y: float, north_square: int = 0, south_square: int = 0
) -> bool:
    """
    Return True if and only if the point `(x, y)` lies in the image of
    the rHEALPix projection of the unit sphere.

    EXAMPLES::

        >>> eps = 0     # Test boundary points.
        >>> north_square, south_square = 0, 0
        >>> rhp = [
        ... (-pi - eps, pi/4 + eps),
        ... (-pi + north_square*pi/2 - eps, pi/4 + eps),
        ... (-pi + north_square*pi/2 - eps, 3*pi/4 + eps),
        ... (-pi + (north_square + 1)*pi/2 + eps, 3*pi/4 + eps),
        ... (-pi + (north_square + 1)*pi/2 + eps, pi/4 + eps),
        ... (pi + eps, pi/4 + eps),
        ... (pi + eps,-pi/4 - eps),
        ... (-pi + (south_square + 1)*pi/2 + eps,-pi/4 - eps),
        ... (-pi + (south_square + 1)*pi/2 + eps,-3*pi/4 - eps),
        ... (-pi + south_square*pi/2 - eps,-3*pi/4 - eps),
        ... (-pi + south_square*pi/2 -eps,-pi/4 - eps),
        ... (-pi - eps,-pi/4 - eps)
        ... ]
        >>> for p in rhp:
        ...     if not in_rhealpix_image(*p):
        ...             print('Fail')
        ...
        >>> print(in_rhealpix_image(0, 0))
        True
        >>> print(in_rhealpix_image(0, pi/4 + 0.1))
        False

    """
    # matplotlib is a third-party module.
    from matplotlib.path import Path

    # Fuzz to slightly expand rHEALPix image so that
    # points on the boundary count as lying in the image.
    eps = 1e-15
    vertices = [
        (-pi - eps, pi / 4 + eps),
        (-pi + north_square * pi / 2 - eps, pi / 4 + eps),
        (-pi + north_square * pi / 2 - eps, 3 * pi / 4 + eps),
        (-pi + (north_square + 1) * pi / 2 + eps, 3 * pi / 4 + eps),
        (-pi + (north_square + 1) * pi / 2 + eps, pi / 4 + eps),
        (pi + eps, pi / 4 + eps),
        (pi + eps, -pi / 4 - eps),
        (-pi + (south_square + 1) * pi / 2 + eps, -pi / 4 - eps),
        (-pi + (south_square + 1) * pi / 2 + eps, -3 * pi / 4 - eps),
        (-pi + south_square * pi / 2 - eps, -3 * pi / 4 - eps),
        (-pi + south_square * pi / 2 - eps, -pi / 4 - eps),
        (-pi - eps, -pi / 4 - eps),
    ]
    poly = Path(vertices)
    return bool(poly.contains_point([x, y]))


def rhealpix_vertices(
    north_square: int = 0, south_square: int = 0
) -> list[tuple[float, float]]:
    """
    Return a list of the planar vertices of the rHEALPix projection of
    the unit sphere.
    """
    vertices = [
        (pi, pi / 4),
        (-pi + (north_square + 1) * pi / 2, pi / 4),
        (-pi + (north_square + 1) * pi / 2, 3 * pi / 4),
        (-pi + north_square * pi / 2, 3 * pi / 4),
        (-pi + north_square * pi / 2, pi / 4),
        (-pi, pi / 4),
        (-pi, -pi / 4),
        (-pi + south_square * pi / 2, -pi / 4),
        (-pi + south_square * pi / 2, -3 * pi / 4),
        (-pi + (south_square + 1) * pi / 2, -3 * pi / 4),
        (-pi + (south_square + 1) * pi / 2, -pi / 4),
        (pi, -pi / 4),
    ]
    # Delete unnecessary non-vertices.
    if north_square == 3:
        vertices.remove((pi, pi / 4))
        vertices.remove((pi, pi / 4))
    elif north_square == 0:
        vertices.remove((-pi, pi / 4))
        vertices.remove((-pi, pi / 4))
    if south_square == 3:
        vertices.remove((pi, -pi / 4))
        vertices.remove((pi, -pi / 4))
    elif south_square == 0:
        vertices.remove((-pi, -pi / 4))
        vertices.remove((-pi, -pi / 4))

    return vertices


def rhealpix(
    a: float = 1,
    e: float = 0,
    north_square: int = 0,
    south_square: int = 0,
    region: str = "none",
) -> Callable[[float, float, bool, bool], tuple[float, float]]:
    """
    Return a function object that wraps the rHEALPix projection and its inverse
    of an ellipsoid with major radius `a` and eccentricity `e`.

    EXAMPLES::

        >>> f = rhealpix(a=2, e=0, north_square=1, south_square=2)
        >>> print(my_round(f(0, pi/3, radians=True), 15))
        (-0.574951359778215, 2.145747686573111)
        >>> p = (0, 60)
        >>> q = f(*p, radians=False)
        >>> print(my_round(q, 15))
        (-0.574951359778215, 2.145747686573111)
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
            return tuple(
                R_A
                * array(
                    rhealpix_ellipsoid(
                        lam,
                        phi,
                        e=e,
                        north_square=north_square,
                        south_square=south_square,
                        region=region,
                    )
                )
            )
        else:
            # Scale down to R_A = 1.
            x, y = array((u, v)) / R_A
            lam, phi = array(
                rhealpix_ellipsoid_inverse(
                    x,
                    y,
                    e=e,
                    north_square=north_square,
                    south_square=south_square,
                    region=region,
                )
            )

            if not radians:
                # Convert to degrees.
                lam, phi = rad2deg([lam, phi])
            return lam, phi

    return f
