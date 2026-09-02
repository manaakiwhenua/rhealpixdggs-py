"""
This Python 3.11 module implements the rHEALPix discrete global grid system.

- Alexander Raichev (AR), 2012-11-12: Initial version based upon grids.py.

NOTES:

All lengths are measured in meters and all angles are measured in radians
unless indicated otherwise.

By 'ellipsoid' throughout, I mean an ellipsoid of revolution and *not* a general (triaxial) ellipsoid.

Points lying on the plane are given in rectangular (horizontal, vertical) coordinates, and points lying on the ellipsoid are given in geodetic (longitude, latitude) coordinates unless indicated otherwise.

DGGS abbreviates 'discrete global grid system'.

Except when manipulating positive integers, I avoid the modulo function '%'
and insted write everything in terms of 'floor()'.
This is because Python interprets the sign of '%' differently than
Java or C, and I don't want to confuse people who are translating this code
to those languages.

EXAMPLES:

Create the (1, 2)-rHEALPix DGGS with N_side = 3 that is based on the WGS84 ellipsoid. Use degrees instead of the default radians for angular measurements ::

    >>> from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
    >>> E = WGS84_ELLIPSOID
    >>> rdggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=2, N_side=3)
    >>> print(rdggs)
    rHEALPix DGGS:
        N_side = 3
        north_square = 1
        south_square = 2
        max_areal_resolution = 1
        max_resolution = 15
        ellipsoid:
            R_A = 6371007.180918476
            a = 6378137.0
            b = 6356752.314245179
            e = 0.08181919084262149
            f = 0.0033528106647474805
            lat_0 = 0
            lon_0 = 0
            radians = False
            sphere = False

Pick a (longitude-latitude) point on the ellipsoid and find the resolution 1 cell that contains it ::

    >>> p = (0, 45)
    >>> c = rdggs.cell_from_point(1, p, plane=False); print(c)
    N8

Find the ellipsoidal (edge) neighbors of this cell ::

    >>> for (direction, cell) in sorted(c.neighbors(plane=False).items()):
    ...     print(direction, cell)
    east N5
    south_east Q0
    south_west P2
    west N7

Find the planar (edge) neighbors of this cell ::

    >>> for (direction, cell) in sorted(c.neighbors('plane').items()):
    ...     print(direction, cell)
    down P2
    left N7
    right Q0
    up N5

Find all the resolution 1 cells intersecting the longitude-latitude aligned ellipsoidal quadrangle with given northwest and southeast corners ::

    >>> nw = (0, 45)
    >>> se = (90, 0)
    >>> cells = rdggs.cells_from_region(1, nw, se, plane=False)
    >>> for row in cells:
    ...     print([str(cell) for cell in row])
    ['N8', 'N5', 'N2']
    ['Q0', 'Q1', 'Q2', 'R0']
    ['Q3', 'Q4', 'Q5', 'R3']

Compute the ellipsoidal nuclei of these cells ::

    >>> expected_results = [
    ...    [
    ...        (1.90833280887811e-14, 58.52801748206219),
    ...        (45.00000000000002, 58.52801748206219),
    ...        (89.99999999999997, 58.52801748206219)
    ...    ], [
    ...        (14.999999999999998, 26.490118751439734),
    ...        (45.0, 26.490118751439734),
    ...        (74.99999999999999, 26.490118751439734),
    ...        (105.00000000000001, 26.490118751439734)
    ...    ], [
    ...        (14.999999999999998, 0),
    ...        (45.0, 0),
    ...        (74.99999999999999, 0),
    ...        (105.00000000000001, 0)
    ...    ]]
    >>> for i, row in enumerate(cells):
    ...     for j, cell in enumerate(row):
    ...         print(cell, assert_allclose(cell.nucleus(plane=False), expected_results[i][j], rtol=1e-15, atol=0) == None)
    N8 True
    N5 True
    N2 True
    Q0 True
    Q1 True
    Q2 True
    R0 True
    Q3 True
    Q4 True
    Q5 True
    R3 True

Create a (0, 0)-rHEALPix DGGS with N_side = 3 based on the WGS84 ellipsoid.
Use degrees instead of the default radians for angular measurements and
orient the DGGS so that the planar origin (0, 0) is on Auckland, New Zealand ::

    >>> p = (174, -37)  # Approximate Auckland lon-lat coordinates
    >>> from rhealpixdggs.ellipsoids import *
    >>> E = Ellipsoid(a=WGS84_A, f=WGS84_F, radians=False, lon_0=p[0], lat_0=p[1])
    >>> rdggs = RHEALPixDGGS(E, N_side=3, north_square=0, south_square=0)
    >>> print(rdggs)
    rHEALPix DGGS:
        N_side = 3
        north_square = 0
        south_square = 0
        max_areal_resolution = 1
        max_resolution = 15
        ellipsoid:
            R_A = 6371007.180918476
            a = 6378137.0
            b = 6356752.314245179
            e = 0.08181919084262149
            f = 0.0033528106647474805
            lat_0 = -37
            lon_0 = 174
            radians = False
            sphere = False

    >>> print(rdggs.cell_from_point(1, p, plane=False))
    Q3

"""

# *****************************************************************************
#       Copyright (C) 2012 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http: //www.gnu.org/licenses/
# *****************************************************************************
from itertools import pairwise, product
from math import asin, copysign, floor
from random import randint
from typing import Literal, overload

from numpy import array, base_repr, ceil, log, pi

# assert_allclose is doctest-only: the doctests use it from the module globals.
from numpy.testing import assert_allclose  # noqa: F401

import rhealpixdggs.pj_rhealpix as pjr
import rhealpixdggs.projection_wrapper as pw
from rhealpixdggs.cell import CELLS0, Cell
from rhealpixdggs.ellipsoids import (
    UNIT_SPHERE,
    UNIT_SPHERE_RADIANS,
    WGS84_ELLIPSOID,
    WGS84_ELLIPSOID_RADIANS,
)

# my_round is doctest-only: the doctests use it from the module globals.
from rhealpixdggs.utils import auth_lat, my_round  # noqa: F401


class RHEALPixDGGS:
    """
    Represents an rHEALPix DGGS on a given ellipsoid.

    INSTANCE ATTRIBUTES:

    - `ellipsoid` - The underlying ellipsoid (Ellipsoid instance).
    - `N_side` - An integer of size at least 2.
      Each planar cell has N_side x N_side child cells.
    - `(north_square, south_square)` - Integers between 0 and 3 indicating
      the positions of north polar and south polar squares, respectively,
      of the rHEALPix projection used.
    - `max_areal_resolution` - An area measured in square meters that
      upper bounds the area of the smallest ellipsoidal grid cells.
    - `max_resolution` - A nonnegative integer that is the maximum grid
      resolution needed to have ellipsoidal cells of area at most
      `max_areal_resolution`.
    - `child_order` - A dictionary of the ordering (Morton order) of child
      cells of a cell in terms of the row-column coordinates in the matrix
      of child cells.  Child cell are numbered 0 to `N_side**2 -1` from
      left to right and top to bottom.
    - `ul_vertex` - A dictionary with key-value pairs (c, (x, y)), where
      c is an element of `CELLS0` and (x, y) is the upper left corner point
      of the resolution 0 planar cell c.
    - `atomic_neighbors` - A dictionary with key-value pairs
      (n, {'up': a, 'down': b, 'left': c, 'right': d}),
      where n, a, b, c, and d are elements of `CELLS0` or
      {0, 1, ..., `N_side**2 -1`}.
      Describes the planar (edge) neighbors of cell0 letter / child cell number
      n.

    NOTE:

    Several RHEALPixDGGS methods have the keyword argument 'plane'.
    Setting it to True indicates that all input and output points and cells are
    interpreted as lying in the planar DGGS.
    Setting it to False indicates that they are interpreted as lying in
    the ellipsoidal DGGS.
    """

    def __init__(
        self,
        ellipsoid=WGS84_ELLIPSOID,
        N_side=3,
        north_square=0,
        south_square=0,
        max_areal_resolution=1,  # square metres
    ):
        self.N_side = N_side
        self.north_square = north_square % 4  # = 0, 1, 2, or 3.
        self.south_square = south_square % 4  # = 0, 1, 2, or 3.
        self.max_areal_resolution = max_areal_resolution
        # Find the maximum grid resolution needed to have ellipsoidal
        # cells of area at most max_areal_resolution.
        self.max_resolution = int(
            ceil(
                log(ellipsoid.R_A**2 * (2 * pi / 3) / max_areal_resolution)
                / (2 * log(N_side))
            )
        )

        self.ellipsoid = ellipsoid
        self._projection_cache: dict[str, pw.Projection] = {}

        # Dictionary of the ordering (Morton order) of child cells of a cell
        # in terms of the row-column coordinates in the matrix of child cells.
        # Child cell are numbered 0 to N_side**2 -1 from left to right and top
        # to bottom.
        # Here's a diagram of the ordering and coordinates for N_side=3:
        #
        # 0 | 0 1 2
        # 1 | 3 4 5
        # 2 | 6 7 8
        #   --------
        #     0 1 2
        #
        child_order = {}
        for row, col in product(list(range(N_side)), repeat=2):
            order = row * N_side + col
            # Handy to have both coordinates and order as dictionary keys.
            child_order[(row, col)] = order
            child_order[order] = (row, col)
        self.child_order = child_order

        # Store the upper left vertices of the resolution 0 cells
        # in the rHEALPix grid hierarchy for this ellipsoid.
        # The default layout by CELLS0 index is
        #
        # 0
        # 1 2 3 4
        # 5.
        #
        ul_vertex = {  # Location for radius = 1
            CELLS0[0]: (-pi + self.north_square * pi / 2, 3 * pi / 4),
            CELLS0[1]: (-pi, pi / 4),
            CELLS0[2]: (-pi / 2, pi / 4),
            CELLS0[3]: (0, pi / 4),
            CELLS0[4]: (pi / 2, pi / 4),
            CELLS0[5]: (-pi + self.south_square * pi / 2, -pi / 4),
        }
        # Scale up ul_vertex by authalic radius of ellipsoid.
        self.ul_vertex = {}
        for k in list(ul_vertex.keys()):
            self.ul_vertex[k] = tuple(self.ellipsoid.R_A * array(ul_vertex[k]))

        # Initialize atomic neighbor relationships among cells.
        # Dictionary of up, right, down, and left neighbors of
        # resolution 0 cells and their subcells 0--(N_side**2 -1),
        # aka the atomic neighbors.
        # Based on the layouts
        #
        #   0
        #   1 2 3 4   (but folded into a cube) and
        #   5
        #
        #   0 1 2
        #   3 4 5
        #   6 7 8   (example for N_side=3).
        #
        an = {}
        # Neighbors of CELLS0[1], ..., CELLS0[4]
        an[CELLS0[1]] = {
            "left": CELLS0[4],
            "right": CELLS0[2],
            "down": CELLS0[5],
            "up": CELLS0[0],
        }
        an[CELLS0[2]] = {
            "left": CELLS0[1],
            "right": CELLS0[3],
            "down": CELLS0[5],
            "up": CELLS0[0],
        }
        an[CELLS0[3]] = {
            "left": CELLS0[2],
            "right": CELLS0[4],
            "down": CELLS0[5],
            "up": CELLS0[0],
        }
        an[CELLS0[4]] = {
            "left": CELLS0[3],
            "right": CELLS0[1],
            "down": CELLS0[5],
            "up": CELLS0[0],
        }
        # Neighbors of CELLS0[0] and CELLS0[5] depend on
        # volues of north_square and south_square, respectively.
        nn = self.north_square
        an[CELLS0[0]] = {
            "down": CELLS0[(nn + 0) % 4 + 1],
            "right": CELLS0[(nn + 1) % 4 + 1],
            "up": CELLS0[(nn + 2) % 4 + 1],
            "left": CELLS0[(nn + 3) % 4 + 1],
        }
        ss = self.south_square
        an[CELLS0[5]] = {
            "up": CELLS0[(ss + 0) % 4 + 1],
            "right": CELLS0[(ss + 1) % 4 + 1],
            "down": CELLS0[(ss + 2) % 4 + 1],
            "left": CELLS0[(ss + 3) % 4 + 1],
        }

        N = self.N_side
        # Neighbors of 0, 1, ..., N**2 - 1.
        for i in range(N**2):
            an[i] = {
                "left": i - 1,
                "right": i + 1,
                "up": (i - N) % N**2,
                "down": (i + N) % N**2,
            }
        # Adjust left and right edge cases.
        for i in range(0, N**2, N):
            an[i]["left"] = an[i]["left"] + N
        for i in range(N - 1, N**2, N):
            an[i]["right"] = an[i]["right"] - N
        self.atomic_neighbors = an

    def __str__(self):
        result = ["rHEALPix DGGS:"]
        result.append(f"    N_side = {self.N_side}")
        result.append(f"    north_square = {self.north_square}")
        result.append(f"    south_square = {self.south_square}")
        result.append(f"    max_areal_resolution = {self.max_areal_resolution}")
        result.append(f"    max_resolution = {self.max_resolution}")
        result.append("    ellipsoid:")
        for k, v in sorted(self.ellipsoid.__dict__.items()):
            if k == "phi_0":
                continue
            result.append(" " * 8 + k + " = " + str(v))
        return "\n".join(result)

    def __eq__(self, other):
        return (
            other is not None
            and self.ellipsoid == other.ellipsoid
            and self.N_side == other.N_side
            and self.north_square == other.north_square
            and self.south_square == other.south_square
            and self.max_resolution == other.max_resolution
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def healpix(self, u: float, v: float, inverse: bool = False) -> tuple[float, float]:
        """
        Return the HEALPix projection of point `(u, v)` (or its inverse if
        `inverse` = True) appropriate to this rHEALPix DGGS.

        EXAMPLES::

            >>> rdggs = UNIT_003_RADIANS
            >>> print(tuple(x.tolist() for x in my_round(rdggs.healpix(-pi, pi/2), 14)))
            (-2.35619449019234, 1.5707963267949)

        NOTE:

        Uses ``pj_healpix`` instead of the PROJ.4 version of HEALPix.
        """
        if "healpix" not in self._projection_cache:
            self._projection_cache["healpix"] = pw.Projection(
                ellipsoid=self.ellipsoid, proj="healpix"
            )
        return self._projection_cache["healpix"](u, v, inverse=inverse)

    def rhealpix(
        self, u: float, v: float, inverse: bool = False, region: str = "none"
    ) -> tuple[float, float]:
        """
        Return the rHEALPix projection of the point `(u, v)` (or its inverse if
        `inverse` = True) appropriate to this rHEALPix DGGS.

        EXAMPLES::

            >>> rdggs = UNIT_003_RADIANS
            >>> print(tuple(x.tolist() for x in my_round(rdggs.rhealpix(0, pi/3), 14)))
            (-1.858272006684, 2.06871881030324)

        NOTE:

        Uses ``pj_rhealpix`` instead of the PROJ.4 version of rHEALPix.
        """
        if region not in self._projection_cache:
            self._projection_cache[region] = pw.Projection(
                ellipsoid=self.ellipsoid,
                proj="rhealpix",
                north_square=self.north_square,
                south_square=self.south_square,
                region=region,
            )
        return self._projection_cache[region](u, v, inverse=inverse)

    def combine_triangles(
        self, u: float, v: float, inverse: bool = False, region: str = "none"
    ) -> tuple[float, float]:
        """
        Return the combine_triangles() transformation of the point `(u, v)`
        (or its inverse if `inverse` = True) appropriate to the underlying
        ellipsoid.
        It maps the HEALPix projection to the rHEALPix projection.

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> p = (0, 0)
            >>> q = (-pi/4, pi/2)
            >>> print(tuple(x.tolist() for x in rdggs.combine_triangles(*p)))
            (0.0, 0.0)
            >>> print(tuple(x.tolist() for x in my_round(rdggs.combine_triangles(*q), 14)))
            (-2.35619449019234, 1.5707963267949)
        """
        R_A = self.ellipsoid.R_A
        ns = self.north_square
        ss = self.south_square
        # Scale down.
        u, v = array((u, v)) / R_A
        # Combine triangles.
        if region != "equatorial":
            u, v = pjr.combine_triangles(
                u, v, inverse=inverse, north_square=ns, south_square=ss
            )
        # Scale up.
        return tuple(R_A * array((u, v)))

    def triangle(
        self, x: float, y: float, inverse: bool = True
    ) -> tuple[int | None, str]:
        """
        If `inverse` = False, then assume `(x,y)` lies in the image of the
        HEALPix projection that comes with this DGGS, and
        return the number of the HEALPix polar triangle (0, 1, 2, 3, or None)
        and the region ('north_polar', 'south_polar', or 'equatorial') that
        `(x, y)` lies in.
        If `inverse` = True, then assume `(x, y)` lies in the image of
        the rHEALPix projection that comes with this DGGS, map `(x, y)`
        to its HEALPix image (x', y'), and return the number of the HEALPix
        polar triangle and the region that (x', y') lies in.
        If `(x, y)` lies in the equatorial region, then the triangle number
        returned is None.

        OUTPUT:

        The pair (triangle_number, region).

        NOTES:

        This is a wrapper for pjr.triangle().

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> c = rdggs.cell(['N', 7])
            >>> print(rdggs.triangle(*c.nucleus(), inverse=True))
            (0, 'north_polar')

            >>> c = rdggs.cell(['N', 3])
            >>> print(rdggs.triangle(*c.nucleus(), inverse=True))
            (3, 'north_polar')

            >>> c = rdggs.cell(['P', 3])
            >>> print(rdggs.triangle(*c.nucleus(), inverse=True))
            (None, 'equatorial')

            >>> c = rdggs.cell(['S', 5, 2])
            >>> print(rdggs.triangle(*c.nucleus(), inverse=True))
            (1, 'south_polar')

        """
        R_A = self.ellipsoid.R_A
        ns = self.north_square
        ss = self.south_square
        # Scale down.
        x, y = array((x, y)) / R_A
        # Get triangle.
        return pjr.triangle(x, y, inverse=inverse, north_square=ns, south_square=ss)

    def xyz(
        self, u: float, v: float, lonlat: bool = False
    ) -> tuple[float, float, float]:
        """
        Given a point `(u, v)` in the planar image of the rHEALPix projection,
        project it back to the ellipsoid and return its 3D rectangular
        coordinates.
        If `lonlat` = True, then assume `(u, v)` is a longitude-latitude
        point.

        EXAMPLES::

            >>> rdggs = UNIT_003_RADIANS
            >>> print(tuple(x.tolist() for x in my_round(rdggs.xyz(0, pi/4, lonlat=True), 14)))
            (0.70710678118655, 0.0, 0.70710678118655)
        """
        if lonlat:
            lam, phi = u, v
        else:
            lam, phi = self.rhealpix(u, v, inverse=True)
        return self.ellipsoid.xyz(lam, phi)

    def xyz_cube(
        self, u: float, v: float, lonlat: bool = False
    ) -> tuple[float, float, float]:
        """
        Given a point `(u, v)` in the planar version of this rHEALPix DGGS,
        fold the rHEALPix image into a cube centered at the origin,
        and return the resulting point's 3D rectangular coordinates.
        If `lonlat` = True, then assume `(u, v)` is a longitude-latitude
        point.

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> print(tuple(x.tolist() for x in my_round(rdggs.xyz_cube(0, 0), 14)))
            (0.78539816339745, 0.0, -0.78539816339745)

        """
        if lonlat:
            x, y = self.rhealpix(u, v)
        else:
            x, y = u, v
        w = self.cell_width(0)
        north = self.north_square
        south = self.south_square
        # Shift rHEALPix projection (with (x, y) in it) so that cell O
        # has downleft corner (0, 0).
        x, y = array((x, y)) + array((2 * w, w / 2))
        q: tuple[float, float, float]
        # Fold projection.
        if y < 0:
            # S
            x += -south * w
            if south == 0:
                q = (x, 0, y)
            elif south == 1:
                q = (y + w, 0, -x)
            elif south == 2:
                q = (w - x, 0, -y - w)
            else:
                q = (-y, 0, x - w)
        elif y > w:
            # N
            x += -north * w
            if north == 0:
                q = (x, w, -y + w)
            elif north == 1:
                q = (-y + 2 * w, w, -x)
            elif north == 2:
                q = (-x + w, w, y - 2 * w)
            else:
                q = (y - w, w, x - w)
        elif x < w:
            # O
            q = (x, y, 0)
        elif (x >= w) and (x < 2 * w):
            # P
            x += -w
            q = (w, y, -x)
        elif (x >= 2 * w) and (x < 3 * w):
            # Q
            x += -2 * w
            q = (w - x, y, -w)
        else:
            # R
            x += -3 * w
            q = (0, y, x - w)
        # Translate the cube's center to (0, 0).
        centered = array(q) + (w / 2) * array((-1, -1, 1))
        return tuple(centered)

    def cell(self, suid=None, level_order_index=None, post_order_index=None):
        """
        Return a cell (Cell instance) of this DGGS either from its ID or
        from its resolution and index.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> c = rdggs.cell(('N', 4, 5))
            >>> print(isinstance(c, Cell))
            True
            >>> print(c)
            N45

        """
        return Cell(self, suid, level_order_index, post_order_index)

    def grid(self, resolution: int):
        """
        Generator function for all the cells at resolution `resolution`.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> grid0 = rdggs.grid(0)
            >>> print([str(x) for x in grid0])
            ['N', 'O', 'P', 'Q', 'R', 'S']

        """
        suid = [CELLS0[0]] + [0 for i in range(resolution)]
        c = self.cell(suid)
        yield c
        cs = c.successor(resolution)
        while cs:
            yield cs
            cs = cs.successor(resolution)

    def num_cells(
        self, res_1: int, res_2: int | None = None, subcells: bool = False
    ) -> int:
        """
        Return the number of cells of resolutions `res_1` to `res_2`
        (inclusive).
        Assume `res_1 <= res_2`.
        If `subcells` = True, then return the number of subcells at resolutions
        `res_1` to `res_2` (inclusive) of a cell at resolution `res_1`.
        If `res_2=None` and `subcells=False, then return the number of
        cells at resolution `res_1`.
        If `res_2=None` and `subcells` = True, then return the number of
        subcells from resolution `res_1` to resolution `self.max_resolution`.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> rdggs.num_cells(0)
            6
            >>> rdggs.num_cells(0, 1)
            60
            >>> rdggs.num_cells(0, subcells=True)
            231627523606480
            >>> rdggs.num_cells(0, 1, subcells=True)
            10
            >>> rdggs.num_cells(5, 6, subcells=True)
            10

        """
        k = self.N_side**2
        if subcells:
            if (res_2 is None) or (res_2 < res_1):
                res_2 = self.max_resolution
            num = int((k ** (res_2 - res_1 + 1) - 1) / (k - 1))
        else:
            if (res_2 is None) or (res_2 < res_1):
                res_2 = res_1
            num = int(6 * (k ** (res_2 + 1) - k**res_1) / (k - 1))
        return num

    @overload
    def cell_width(self, resolution: int, plane: Literal[True] = ...) -> float: ...

    @overload
    def cell_width(self, resolution: int, plane: Literal[False]) -> None: ...

    def cell_width(self, resolution: int, plane: bool = True) -> float | None:
        """
        Return the width of a planar cell at the given resolution.
        If `plane` = False, then return None,
        because the ellipsoidal cells don't have constant width.

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> print(rdggs.cell_width(0) == pi/2)
            True
            >>> print(rdggs.cell_width(1) == pi/6)
            True

        """
        if plane:
            return self.ellipsoid.R_A * (pi / 2) * self.N_side ** (-resolution)
        return None

    def cell_area(self, resolution: int, plane: bool = True) -> float:
        """
        Return the area of a planar or ellipsoidal cell at the given
        resolution.

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> a = rdggs.cell_area(1)
            >>> print(a == (pi/6)**2)
            True
            >>> print(rdggs.cell_area(1, plane=False) == 8/(3*pi)*a)
            True

        """
        w = self.cell_width(resolution)
        if plane:
            return w**2
        else:
            return 8 / (3 * pi) * w**2

    def area_error_budget(self) -> dict[int, dict]:
        """
        Return an analytical error budget for cell area equality testing at
        each resolution from 0 to ``max_resolution``.

        rHEALPix is theoretically equal-area: all cells at a given resolution
        have identical area by construction.  Floating-point arithmetic
        introduces a small representational error when that area is computed;
        this method exposes it so callers can apply a principled tolerance when
        testing whether two computed cell areas should be considered equal.

        The tolerance is derived analytically by tracing the rounding errors
        through ``cell_area(r, plane=False) = 8/(3π) × (R_A × π/2 × N_side⁻ʳ)²``:

        * ``N_side⁻ʳ``         — one division of an exact integer, ≤ ε rel. error
        * ``R_A × (π/2)``      — two multiplications, ≤ 2ε
        * ``× N_side⁻ʳ``       — one more multiplication, ≤ 3ε cumulative
        * squaring (``w²``)    — doubles rel. error, ≤ 7ε
        * ``8/(3π) × w²``      — two multiplications and one division, ≤ 10ε

        The relative tolerance is therefore bounded by 10 × machine epsilon
        regardless of resolution.

        Returns a ``dict`` keyed by resolution.  Each value is a ``dict`` with:

        ``cell_area_m2``
            Theoretical equal area per cell in square metres.
        ``abs_tolerance``
            Absolute area tolerance (m²) for equal-area comparisons.
        ``rel_tolerance``
            Relative tolerance (dimensionless); constant across resolutions.

        EXAMPLES::

            >>> rdggs = WGS84_003
            >>> budget = rdggs.area_error_budget()
            >>> all(r in budget for r in range(rdggs.max_resolution + 1))
            True
            >>> budget[0]['cell_area_m2'] > budget[1]['cell_area_m2']
            True
            >>> budget[0]['rel_tolerance'] == budget[1]['rel_tolerance']
            True

        """
        import sys

        # Conservative upper bound: 10 rounding steps in the area formula.
        # See docstring for the derivation.
        rel_tol = 10 * sys.float_info.epsilon
        budget = {}
        for r in range(self.max_resolution + 1):
            area = self.cell_area(r, plane=False)
            budget[r] = {
                "cell_area_m2": area,
                "abs_tolerance": area * rel_tol,
                "rel_tolerance": rel_tol,
            }
        return budget

    def interval(self, a, b):
        """
        Generator function for all the resolution
        `max(a.resolution, b.resolution)` cells between cell
        `a` and cell `b` (inclusive and with respect to the
        postorder ordering on cells).
        Note that `a` and `b` don't have to lie at the same resolution.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> a = rdggs.cell(('N', 1))
            >>> b = rdggs.cell(('N',))
            >>> print([str(c) for c in list(rdggs.interval(a, b))])
            ['N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8']

        """
        # Choose the starting cell, which might not be A.
        resolution = max(a.resolution, b.resolution)
        if a.resolution < resolution:
            cell = a.successor(resolution)
        else:
            cell = Cell(self, a.suid[: resolution + 1])
        while cell <= b:
            yield cell
            cell = cell.successor(resolution)

    def cell_from_point(
        self, resolution: int, p: tuple[float, float], plane: bool = True
    ) -> Cell | None:
        """
        Return the resolution `resolution` cell that contains the point `p`.
        If `plane` = True, then `p` and the output cell lie in the
        planar DGGS.
        Otherwise, `p` and the output cell lie in the ellipsoidal DGGS.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> p = (0, 0)
            >>> c = rdggs.cell_from_point(1, p)
            >>> print(c)
            Q3
            >>> rdggs = RHEALPixDGGS(N_side=15)
            >>> p = (80, -20)
            >>> c = rdggs.cell_from_point(1, p, plane=False)
            >>> print(c)
            (Q, 178)

        """
        # Get the rectangular coordinates of p.
        if plane:
            x, y = p
        else:
            x, y = self.rhealpix(*p)

        # Determine the resolution 0 cell c0 that (x, y) lies in,
        # since resolution 0 cells are anamolous.
        ns = self.north_square
        ss = self.south_square
        R = self.ellipsoid.R_A
        if (
            y > R * pi / 4
            and y < R * 3 * pi / 4
            and x > R * (-pi + ns * (pi / 2))
            and x < R * (-pi / 2 + ns * (pi / 2))
        ):
            s0 = CELLS0[0]
        elif (
            y > -R * 3 * pi / 4
            and y < -R * pi / 4
            and x > R * (-pi + ss * (pi / 2))
            and x < R * (-pi / 2 + ss * (pi / 2))
        ):
            s0 = CELLS0[5]
        elif y >= -R * pi / 4 and y <= R * pi / 4 and x >= -R * pi and x < -R * pi / 2:
            s0 = CELLS0[1]
        elif y >= -R * pi / 4 and y <= R * pi / 4 and x >= -R * pi / 2 and x < 0:
            s0 = CELLS0[2]
        elif y >= -R * pi / 4 and y <= R * pi / 4 and x >= 0 and x < R * pi / 2:
            s0 = CELLS0[3]
        elif y >= -R * pi / 4 and y <= R * pi / 4 and x >= R * pi / 2 and x < R * pi:
            s0 = CELLS0[4]
        else:
            # (x, y) doesn't lie in the DGGS.
            return None

        suid = [s0]
        if resolution == 0:
            # Done.
            return Cell(self, suid)

        # Compute the horizontal and vertical distances between (x, y) and
        # the ul_vertex of c0 as fractions of the width of c0.
        w = self.cell_width(0)
        dx = abs(x - self.ul_vertex[suid[0]][0]) / w
        dy = abs(y - self.ul_vertex[suid[0]][1]) / w
        if dx == 1:
            # This case is analytically impossible
            # but, i guess, numerically possible because of rounding errors.
            # Border case. Take a smidgen off dx.
            dx -= 0.5 * self.cell_width(self.max_resolution) / w
        if dy == 1:
            # Border case. Take a smidgen off dy.
            dy -= 0.5 * self.cell_width(self.max_resolution) / w

        N = self.N_side
        # Compute the base N expansions of dx and dy and truncate them
        # at index resolution to get the row and column SUIDs of
        # the resolution resolution cell c containing (x,y).
        suid_row = base_repr(int(float(str(dy * N**resolution))), N)
        suid_col = base_repr(int(float(str(dx * N**resolution))), N)
        # Using int(float(str(.))) instead of the straightforward int(.),
        # because the latter gave me rounding errors.
        # Prefix with the appropriate amount of zeros.
        suid_row = "0" * (resolution - len(suid_row)) + suid_row
        suid_col = "0" * (resolution - len(suid_col)) + suid_col

        # Use the column and row SUIDs of c to get the SUID of c.
        for i in range(resolution):
            suid.append(self.child_order[(int(suid_row[i], N), int(suid_col[i], N))])
        return Cell(self, suid)

    def cell_from_region(
        self, ul: tuple[float, float], dr: tuple[float, float], plane: bool = True
    ) -> Cell | None:
        """
        Return the smallest planar or ellipsoidal cell wholly containing
        the region bounded by the axis-aligned rectangle with upper left
        and lower right vertices given by the the points `ul` and `dr`,
        respectively.
        If such as cell does not exist, then return None.
        If `plane` = True, then `ul` and `dr` and the returned cell
        lie in the planar DGGS.
        Otherwise, `ul` and `dr` and the returned cell lie in the ellipsoidal
        DGGS.

        To specify an ellipsoidal cap region, set `ul` = (-pi, pi/2) and
        `dr` = (-pi, phi) for a northern cap from latitudes pi/2 to phi, or
        set `ul` = (-pi, phi) and `dr` = (-pi, -pi/2) for a southern cap from
        latitudes phi to -pi/2.
        (As usual, if `self.ellipsoid.radians` = False,
        then use degrees instead of radians when specifying ul and dr.)

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> p = (0, pi/12)
            >>> q = (pi/6 - 1e-6, 0)
            >>> c = rdggs.cell_from_region(p, q)
            >>> print(c)
            Q3

        """
        if not plane:
            # Compute planar ul and dr as follows.
            # Get all four vertices of the ellipsoidal cap or quadrangle.
            PI = self.ellipsoid.pi()
            if ul == (-PI, PI / 2) or dr == (-PI, -PI / 2):
                # Cap.
                if dr[1] != -PI / 2:
                    phi = dr[1]
                else:
                    phi = ul[1]
                vertices = [
                    (-3 * PI / 4, phi),
                    (-PI / 4, phi),
                    (PI / 4, phi),
                    (3 * PI / 4, phi),
                ]
            else:
                # Quadrangle.
                vertices = [ul, (ul[0], dr[1]), dr, (dr[0], ul[1])]
            # Project the vertices onto the plane.
            vertices = [self.rhealpix(*p) for p in vertices]
            # Find the upper left and lower right vertices of the
            # planar bounding rectangle.
            ul = (min([p[0] for p in vertices]), max([p[1] for p in vertices]))
            dr = (max([p[0] for p in vertices]), min([p[1] for p in vertices]))

        # Find the resolution max_resolution cells containing ul and dr.
        resolution = self.max_resolution
        ul_cell = self.cell_from_point(resolution, ul)
        dr_cell = self.cell_from_point(resolution, dr)
        if ul_cell is None or dr_cell is None:
            return None
        ul_suid = ul_cell.suid
        dr_suid = dr_cell.suid

        # Find the longest common prefix of ul_suid and dr_suid.
        least = resolution + 1  # Default if the suids agree everywhere
        for i in range(resolution + 1):
            if ul_suid[i] != dr_suid[i]:
                least = i
                break
        if least == 0:
            # No one cell contains R.
            return None
        else:
            return self.cell(ul_suid[:least])

    def cell_latitudes(
        self,
        resolution: int,
        phi_min: float,
        phi_max: float,
        nucleus: bool = True,
        plane: bool = True,
    ) -> list[float]:
        """
        Return a list of every latitude phi whose parallel intersects
        a resolution `resolution` cell nucleus and satisfies
        `phi_min` < phi < `phi_max`.
        If `plane` = True, then use rHEALPix y-coordinates for `phi_min`,
        `phi_max`, and the result. Return the list in increasing order.
        If `nucleus` = False, then return a list of every latitude phi whose
        parallel intersects the north or south boundary of a resolution
        `resolution` cell and that satisfies `phi_min` < phi < `phi_max`.

        NOTE:

        By convention, the pole latitudes pi/2 and -pi/2 (or their
        corresponding rHEALPix y-coordinates) will be excluded.

        There are 2*self.N_side**resolution - 1 nuclei
        latitudes between the poles if self.N_side is odd and
        2*self.N_side**resolution if self.N_side is even.
        Consequently, there are 2*self.N_side**resolution
        boundary latitudes between the poles if self.N_side is odd and
        2*self.N_side**resolution - 1 boundary latitudes if self.N_side is
        even.

        EXAMPLES::

            >>> rdggs = WGS84_003_RADIANS
            >>> for phi in rdggs.cell_latitudes(1, -pi/2, pi/2, plane=False):
            ...     print(my_round(phi, 14))
            -1.02150660972679
            -0.46233979145691
            0.0
            0.46233979145691
            1.02150660972679

            >>> for phi in rdggs.cell_latitudes(1, -pi/2, pi/2, nucleus=False, plane=False):
            ...     print(my_round(phi, 14))
            -1.29894395947616
            -0.73195363195267
            -0.22506566919844
            0.22506566919844
            0.73195363195267
            1.29894395947616

        """
        if phi_min > phi_max:
            return []
        # Work in the plane first, because that's easier.
        R = self.ellipsoid.R_A
        # Set y_min and y_max.
        if plane:
            y_min = phi_min
            y_max = phi_max
        else:
            y_min = self.healpix(0, phi_min)[1]
            y_max = self.healpix(0, phi_max)[1]
        w = self.cell_width(resolution)
        # Set first y, which is the only step that depends on the
        # nucleus keyword.
        if self.N_side % 2 == 1:
            if nucleus:
                y = -R * pi / 2 + w
            else:
                y = -R * pi / 2 + w / 2
        else:
            if resolution == 0:
                # Anomalous.
                if nucleus:
                    y = 0
                else:
                    y = -R * pi / 4
            else:
                if nucleus:
                    y = -R * pi / 2 + w / 2
                else:
                    y = -R * pi / 2 + w
        # Start y above y_min.
        if y <= y_min:
            dy = y_min - y
            y = max(y + int(ceil(dy / w)) * w, y + w)
        # Collect the ys.
        result = []
        while y < y_max:
            result.append(y)
            y += w
        # Convert to latitudes if desired.
        if not plane:
            result = [self.healpix(R * pi / 4, y, inverse=True)[1] for y in result]
        return result

    def cells_from_meridian(
        self, resolution: int, lam: float, phi_min: float, phi_max: float
    ) -> list[Cell]:
        """
        Return a list of the resolution `resolution` cells that intersect
        the meridian segment of longitude `lam` whose least latitude is
        `phi_min` and whose greatest latitude is `phi_max`.
        Sort the cells from north to south and west to east in case
        two cells with the same nucleus latitude intersect the
        meridian.

        EXAMPLES::

            >>> rdggs = WGS84_003_RADIANS
            >>> cells = rdggs.cells_from_meridian(1, 0.1, -pi/2, pi/2)
            >>> print([str(cell) for cell in cells])
            ['N4', 'N2', 'N1', 'Q0', 'Q3', 'Q6', 'S8', 'S7', 'S4']

        """
        if phi_min > phi_max:
            return []
        start = self.cell_from_point(resolution, (lam, phi_max), plane=False)
        end = self.cell_from_point(resolution, (lam, phi_min), plane=False)
        # Points on the ellipsoid always land in a cell.
        assert start is not None and end is not None
        if start == end:
            return [start]
        # Get latitudes of cell nuclei that lie ibetween start and end.
        phis = self.cell_latitudes(resolution, phi_min, phi_max, True, plane=False)
        if not phis:
            return [start, end]
        # Will have at least three cells in the final list.
        # Collect all the cells that contain the points
        # (lam, phi) for phi in phis.
        result = []
        for phi in reversed(phis):
            c = self.cell_from_point(resolution, (lam, phi), plane=False)
            assert c is not None  # points on the ellipsoid always land in a cell
            new_cells = [c]
            if c.ellipsoidal_shape in ["dart", "skew_quad"]:
                # Either the east or the west neighbor of c
                # might also intersect the meridian.
                # So include the neighbor too.
                west = c.neighbor("west", plane=False)
                east = c.neighbor("east", plane=False)
                if west.intersects_meridian(lam):
                    new_cells = [west, c]
                elif east.intersects_meridian(lam):
                    new_cells = [c, east]
            result.extend(new_cells)
        # Add start and end if they weren't added in the for loop.
        if start not in result[0:2]:
            result.insert(0, start)
        if end not in result[-1:-3]:
            result.append(end)
        return result

    def cells_from_parallel(
        self, resolution: int, phi: float, lam_min: float, lam_max: float
    ) -> list[Cell]:
        """
        Return a list of the resolution `resolution` cells that intersect
        the parallel segment of latitude `phi` whose least longitude is
        `lam_min` and whose greatest longitude is `lam_max`.
        Sort the list from west to east.

        EXAMPLES::

            >>> rdggs = WGS84_003_RADIANS
            >>> cells = rdggs.cells_from_parallel(1, pi/3, -pi, pi)
            >>> print([str(cell) for cell in cells])
            ['N6', 'N7', 'N8', 'N5', 'N2', 'N1', 'N0', 'N3']

        """
        if lam_min > lam_max:
            return []
        start = self.cell_from_point(resolution, (lam_min, phi), plane=False)
        end = self.cell_from_point(resolution, (lam_max, phi), plane=False)
        # Points on the ellipsoid always land in a cell.
        assert start is not None and end is not None
        PI = self.ellipsoid.pi()
        if start == end:
            if start.ellipsoidal_shape == "cap" or lam_max - lam_min < PI / 2:
                return [start]
            else:
                # Need to wrap all the way around globe.
                end = start.neighbor("west", plane=False)
        result = []
        current = start
        while current != end:
            result.append(current)
            current = current.neighbor("east", plane=False)
        result.append(end)
        return result

    def cells_from_line(
        self,
        resolution: int,
        lstart: tuple[float, float],
        lend: tuple[float, float],
        plane: bool = True,
        wrap_antimeridian: bool = False,
    ) -> list[Cell]:
        """
        Return the ordered list of resolution `resolution` cells that the
        line segment from `lstart` to `lend` passes through.

        The segment is straight in the given coordinate space: planar
        coordinates if `plane` = True, longitude-latitude coordinates
        otherwise (so it is a plate carree straight line, not a geodesic;
        to trace a geodesic, densify it into short segments first). In
        particular, a longitude-latitude segment spanning more than half
        a turn of longitude does not, by default, wrap around the
        antimeridian: a segment from longitude 179 to longitude -179 runs
        the long way around, through longitude 0 -- the literal planar
        reading of the coordinates. Splitting inputs at the antimeridian
        (as GeoJSON's RFC 7946 prescribes for data producers) avoids the
        ambiguity entirely; alternatively, set `wrap_antimeridian` = True
        to trace such a segment the short way, across the antimeridian,
        which traces exactly the same cells as splitting it there. The
        flag is meaningless (and ignored) for `plane` = True.

        The sequence is exact, computed by a sweep over the segment's
        planar image: every parameter value where the projected segment
        crosses a cell edge is found (the projection of a longitude-
        latitude segment is piecewise smooth, with its pieces' boundaries
        -- region and polar-triangle changes and the longitude wrap --
        known in closed form, and each planar coordinate at most singly
        non-monotone per piece, so every crossing is bracketed and solved
        to machine precision), and each inter-crossing interval's midpoint
        is located with `cell_from_point`, which is exact for every cell
        shape, including polar cap cells. Cells the segment meets in a
        single point only (passing exactly through a cell corner) are not
        included.

        If either endpoint lies outside the grid (`cell_from_point`
        returns None there), return []. If some middle stretch of a
        *planar* segment leaves the grid's cross-shaped image, the cells
        on both sides of the gap are still returned, in order.

        EXAMPLES::

            >>> rdggs = WGS84_003
            >>> cells = rdggs.cells_from_line(3, (-89.669615, 86.549596), (-134, 86), False)
            >>> print([str(cell) for cell in cells])
            ['N448', 'N447']

        """
        if wrap_antimeridian and not plane:
            # Take the short way in longitude: shift the end longitude by
            # a full turn so the segment crosses the antimeridian. The
            # projection wraps longitudes, so the out-of-range value
            # locates to the same cells throughout.
            full = 2 * pi if self.ellipsoid.radians else 360.0
            dlon = lend[0] - lstart[0]
            if abs(dlon) > full / 2:
                lend = (lend[0] - copysign(full, dlon), lend[1])

        start = self.cell_from_point(resolution, lstart, plane)
        end = self.cell_from_point(resolution, lend, plane)
        if start is None or end is None:
            return []
        if start == end:
            return [start]

        R = self.ellipsoid.R_A
        w = self.cell_width(resolution)

        def point_at(t):
            return (
                lstart[0] + t * (lend[0] - lstart[0]),
                lstart[1] + t * (lend[1] - lstart[1]),
            )

        if plane:
            q = point_at
        else:
            proj = self.rhealpix

            def q(t):
                return proj(*point_at(t))

        # Piece boundaries: parameter values where the planar image of
        # the segment may kink or jump. For a planar segment there are
        # none; for a longitude-latitude segment they are the crossings
        # of the equatorial/polar region boundaries, the polar-triangle
        # edges and face columns (longitude multiples of a quarter turn,
        # relative to the ellipsoid's lon_0), and the longitude/latitude
        # wraps -- all exact, since longitude and latitude are linear in
        # the parameter.
        breakpoints = {0.0, 1.0}

        def add_linear_crossings(f0, f1, targets):
            for target in targets:
                if f1 != f0:
                    t = (target - f0) / (f1 - f0)
                    if 0.0 < t < 1.0:
                        breakpoints.add(t)

        if not plane:
            ell = self.ellipsoid
            half = pi if ell.radians else 180.0
            quarter = half / 2
            # Geodetic latitude of the region boundary (the authalic
            # latitude arcsin(2/3)).
            phi_reg = auth_lat(asin(2.0 / 3), ell.e, inverse=True, radians=True)
            if not ell.radians:
                phi_reg = phi_reg * 180 / pi
            add_linear_crossings(
                lstart[1],
                lend[1],
                [
                    ell.lat_0 + phi_reg,
                    ell.lat_0 - phi_reg,
                    ell.lat_0 + quarter,
                    ell.lat_0 - quarter,
                ],
            )
            lo = min(lstart[0], lend[0]) - ell.lon_0
            hi = max(lstart[0], lend[0]) - ell.lon_0
            k0, k1 = floor(lo / quarter), int(ceil(hi / quarter))
            add_linear_crossings(
                lstart[0],
                lend[0],
                [ell.lon_0 + k * quarter for k in range(k0, k1 + 1)],
            )

        # All planar cell edges of this resolution lie on one global
        # lattice: face origins are exact multiples of the face width,
        # which is an exact multiple of the cell width.
        x_anchor = -pi * R
        y_anchor = -3 * pi * R / 4

        def lattice_lines_between(a, b, anchor):
            lo, hi = (a, b) if a <= b else (b, a)
            i0 = floor((lo - anchor) / w) + 1
            i1 = floor((hi - anchor) / w)
            return [anchor + i * w for i in range(i0, i1 + 1)]

        def root(f, ta, tb, fa):
            # Bisection to machine precision on a bracketed sign change.
            for _ in range(120):
                tm = 0.5 * (ta + tb)
                if tm == ta or tm == tb:
                    return tm
                fm = float(f(tm))
                if fm == 0:
                    return tm
                if (fa < 0) != (fm < 0):
                    tb = tm
                else:
                    ta, fa = tm, fm
            return 0.5 * (ta + tb)

        def monotone_runs(axis, ta, tb):
            # Split [ta, tb] into runs on which q(t)[axis] is monotone:
            # scan for direction flips, then pin each extremum by ternary
            # search. Within one smooth piece each planar coordinate has
            # the form c + l(t)*sigma(t) with l linear and sigma monotone,
            # so it has at most one interior extremum, which a 64-point
            # scan brackets comfortably.
            N = 64
            ts = [ta + (tb - ta) * i / N for i in range(N + 1)]
            vs = [float(q(t)[axis]) for t in ts]
            cuts = [ta]
            direction = 0
            for i in range(1, N + 1):
                d = (vs[i] > vs[i - 1]) - (vs[i] < vs[i - 1])
                if d == 0:
                    continue
                if direction == 0:
                    direction = d
                elif d != direction:
                    lo, hi = ts[max(i - 2, 0)], ts[i]
                    for _ in range(200):
                        if hi - lo < 1e-14:
                            break
                        m1 = lo + (hi - lo) / 3
                        m2 = hi - (hi - lo) / 3
                        if direction * float(q(m1)[axis] - q(m2)[axis]) < 0:
                            lo = m1
                        else:
                            hi = m2
                    cuts.append(0.5 * (lo + hi))
                    direction = d
            cuts.append(tb)
            return list(pairwise(cuts))

        crossings = set()
        for pa, pb in zip(sorted(breakpoints), sorted(breakpoints)[1:]):
            if pb - pa < 1e-14:
                continue
            # Evaluate strictly inside the piece, clear of its kinks.
            eps = (pb - pa) * 1e-12
            ta, tb = pa + eps, pb - eps
            for axis in (0, 1):
                anchor = x_anchor if axis == 0 else y_anchor
                for ra, rb in monotone_runs(axis, ta, tb):
                    va = float(q(ra)[axis])
                    vb = float(q(rb)[axis])
                    for line in lattice_lines_between(va, vb, anchor):

                        def f(t, line=line, axis=axis):
                            return q(t)[axis] - line

                        crossings.add(root(f, ra, rb, va - line))
            if pb < 1.0:
                # The piece boundary itself may be a cell change (a face
                # jump, or a kink lying exactly on a cell edge).
                crossings.add(pb)

        ts = [0.0] + sorted(crossings) + [1.0]
        line_cells = [start]
        for a, b in pairwise(ts):
            cell = self.cell_from_point(resolution, point_at(0.5 * (a + b)), plane)
            if cell is not None and cell != line_cells[-1]:
                line_cells.append(cell)
        if line_cells[-1] != end:
            line_cells.append(end)
        return line_cells

    def cell_boundaries(
        self, cells, n: int = 2, plane: bool = True
    ) -> dict[Cell, list]:
        """
        Return a dictionary mapping each cell in `cells` to its boundary
        points -- each value agreeing with that cell's own
        ``boundary(n=n, plane=plane)`` in point count, order, and
        coordinates (up to floating point) -- while computing the
        projection of every shared boundary point only once.

        Adjacent cells share their common edge's points (and cells
        meeting at a corner share that corner point), so computing each
        cell's boundary independently projects every interior edge of a
        map of cells twice. This method projects each distinct planar
        boundary point once and reuses it, roughly halving the projection
        work for contiguous sets of cells; as a corollary, two adjacent
        cells' copies of their shared points are identical floats rather
        than two independently computed (and potentially last-digit
        different) values, which helps downstream consumers that dissolve
        or snap cell geometries.

        Points are shared only between cells of the same region
        ('equatorial', 'north_polar', 'south_polar'): the inverse
        projection takes a per-cell region hint that legitimately
        disambiguates points lying exactly on a region boundary, so
        edges along those parallels are computed per region, exactly as
        ``boundary()`` computes them.

        `cells` may mix resolutions; sharing happens per resolution.
        For `plane` = True there is no projection work to share and this
        is simply a convenience over calling ``boundary()`` per cell.

        EXAMPLES::

            >>> rdggs = WGS84_003
            >>> cells = list(rdggs.cell(('P', 0)).subcells())
            >>> boundaries = rdggs.cell_boundaries(cells, n=3, plane=False)
            >>> all(len(b) == 4 * 3 - 4 for b in boundaries.values())
            True

        """
        if plane:
            return {cell: cell.boundary(n=n, plane=True) for cell in cells}
        n = max(n, 2)
        R = self.ellipsoid.R_A
        x_anchor = -pi * R
        y_anchor = -3 * pi * R / 4
        cache = {}
        result = {}
        for cell in cells:
            # The same planar points, in the same order, as
            # cell.boundary(n=n, plane=False) computes: the planar
            # boundary, reordered to start at the northwest vertex.
            planar = cell.boundary(n=n, plane=True)
            v = cell.vertices(plane=True)
            nw = cell.nw_vertex(plane=True)
            i = (n - 1) * v.index(nw)
            planar = planar[i:] + planar[:i]
            region = cell.region()
            # All of this cell's boundary points lie on the fine lattice
            # of pitch w/(n - 1) anchored at the planar image's corner,
            # shared with every same-resolution neighbor's points, so an
            # integer lattice key identifies coincident points robustly.
            pitch = cell.width(plane=True) / (n - 1)
            points = []
            for p in planar:
                key = (
                    cell.resolution,
                    region,
                    round((p[0] - x_anchor) / pitch),
                    round((p[1] - y_anchor) / pitch),
                )
                if key not in cache:
                    cache[key] = self.rhealpix(p[0], p[1], inverse=True, region=region)
                points.append(cache[key])
            result[cell] = points
        return result

    def cells_from_region(
        self,
        resolution: int,
        ul: tuple[float, float],
        dr: tuple[float, float],
        plane: bool = True,
    ) -> list[list[Cell]]:
        """
        If `plane` = True, then return a list of lists of resolution
        `resolution` cells that cover the axis-aligned rectangle whose
        upper left and lower right vertices are the points `ul` and `dr`,
        respectively.
        In the output, sort each sublist of cells from left to right (in the
        planar DGGS) and sort the sublists from top to bottom.

        If `plane` = False, then return a list of lists of resolution
        `resolution`
        cells that cover the longitude-latitude aligned ellipsoidal quadrangle
        whose northwest and southeast vertices are the points `ul` and `dr`,
        respectively. Defunct quads with `ul` = (stuff, pi/2) or
        `dr` = (stuff, -pi/2) also work
        (and rely on the fact that the north and south pole can both
        be specified by infinitely many longitudes).

        To specify an ellipsoidal cap region, set `ul` = (-pi, pi/2) and
        `dr` = (-pi, phi) for a northern cap from latitudes pi/2 to phi, or
        set `ul` = (-pi, phi) and `dr` = (-pi, -pi/2) for a southern cap from
        latitudes phi to -pi/2. (As usual, if `self.ellipsoid.radians` = False,
        then use degrees instead of radians when specifying ul and dr.)

        In the output, sort each sublist of cells from west to east (in the
        ellipsoidal DGGS) and sort the sublists from north to south.

        Return the empty list if if `ul[0] > dr[0]` or `ul[1] < dr[1]`.

        NOTE:

        If `plane` = True, then the resulting list is a matrix, that is,
        each sublist has the same length.  This is not necessarily so if
        `plane` = False; see the examples below.

        EXAMPLES::

            >>> rdggs = WGS84_003_RADIANS
            >>> R_A = rdggs.ellipsoid.R_A
            >>> ul = R_A*array((-0.1, pi/4))
            >>> dr = R_A*array((0.1, -pi/4))  # Rectangle
            >>> M = rdggs.cells_from_region(1, ul, dr)
            >>> for row in M:
            ...     print([str(cell) for cell in row])
            ['P2', 'Q0']
            ['P5', 'Q3']
            ['P8', 'Q6']

            >>> ul = (0, pi/3)
            >>> dr = (pi/2, 0)  # Quad
            >>> M = rdggs.cells_from_region(1, ul, dr, plane=False)
            >>> for row in M:
            ...     print([str(cell) for cell in row])
            ['N2', 'N1', 'N0']
            ['Q0', 'Q1', 'Q2', 'R0']
            ['Q3', 'Q4', 'Q5', 'R3']

            >>> ul = (0, -pi/6)
            >>> dr = (pi/2, -pi/2)  # Defunct quad / lune segment
            >>> M = rdggs.cells_from_region(1, ul, dr, plane=False)
            >>> for row in M:
            ...     print([str(cell) for cell in row])
            ['Q6', 'Q7', 'Q8', 'R6']
            ['S8', 'S7', 'S6']
            ['S4']

            >>> ul = (-pi, -pi/5)
            >>> dr = (-pi, -pi/2)  # Cap
            >>> M = rdggs.cells_from_region(1, ul, dr, plane=False)
            >>> for row in M:
            ...     print([str(cell) for cell in row])
            ['O6', 'O7', 'O8', 'P6', 'P7', 'P8', 'Q6', 'Q7', 'Q8', 'R6', 'R7', 'R8']
            ['S0', 'S1', 'S2', 'S5', 'S8', 'S7', 'S6', 'S3']
            ['S4']

        """
        if ul[0] > dr[0] or ul[1] < dr[1]:
            return []
        if plane:
            # Rectangle region.
            # Get the four corner cells.
            ur_cell = self.cell_from_point(resolution, (dr[0], ul[1]), plane)
            dl_cell = self.cell_from_point(resolution, (ul[0], dr[1]), plane)
            ul_cell = self.cell_from_point(resolution, ul, plane)
            dr_cell = self.cell_from_point(resolution, dr, plane)
            if ul_cell is None or ur_cell is None or dl_cell is None or dr_cell is None:
                return []
            if ul_cell == dr_cell:
                return [[ul_cell]]
            # Starting from ul, collect cells from left to right and
            # then top to bottom, ending at dr.
            result = []
            row_start = ul_cell
            row_end = ur_cell
            while True:
                row = []
                current = row_start
                while current != row_end:
                    row.append(current)
                    current = current.neighbor("right", plane)
                row.append(current)
                result.append(row)
                if current == dr_cell:
                    # Done.
                    break
                # Update row start and end cells to their down neighbors,
                # and collect another row of cells.
                row_start = row_start.neighbor("down", plane)
                row_end = row_end.neighbor("down", plane)
                current = row_start
            return result
        # Ellipsoid: quad or cap region.
        phi_min = dr[1]
        phi_max = ul[1]
        # Get the cell nucleus latitudes between phi_min and phi_max.
        phis = self.cell_latitudes(resolution, phi_min, phi_max, True, plane)
        # Collect the cells along the parallels of phis.
        PI = self.ellipsoid.pi()
        if (ul == (-PI, PI / 2) and dr[0] == -PI) or (
            dr == (-PI, -PI / 2) and ul[0] == -PI
        ):
            # Cap.
            lam_min = -PI
            lam_max = PI
        else:
            # Quad.
            lam_min = ul[0]
            lam_max = dr[0]
        result = []
        for phi in reversed(phis):
            cells = self.cells_from_parallel(resolution, phi, lam_min, lam_max)
            result.append(cells)
        # Add the cells along parallels phi_min and phi_max if necessary.
        ul_cell = self.cell_from_point(resolution, ul, plane)
        if not result or result[0][0] != ul_cell:
            # Add cells along phi_max parallel.
            cells = self.cells_from_parallel(resolution, phi_max, lam_min, lam_max)
            result.insert(0, cells)
            ul_cell = self.cell_from_point(resolution, ul, plane)
        dl_cell = self.cell_from_point(resolution, (ul[0], dr[1]), plane)
        if not result or result[-1][0] != dl_cell:
            # Add cells along phi_min parallel.
            cells = self.cells_from_parallel(resolution, phi_min, lam_min, lam_max)
            result.append(cells)
        return result

    def random_point(self, plane: bool = True) -> tuple[float, float]:
        """
        Return a point in this DGGS sampled uniformly at
        random from the plane or from the ellipsoid.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> print(E.random_point()) # doctest: +SKIP
            (-1.0999574573422948, 0.21029104897701129)

        """
        # Pick a random resolution 0 cell.
        n = randint(0, 5)
        c = self.cell([CELLS0[n]])
        # Pick a random point in that cell.
        return c.random_point(plane=plane)

    def random_cell(self, resolution: int | None = None) -> Cell:
        """
        Return a cell of the given resolution chosen uniformly at random
        from all cells at that resolution.
        If `resolution=None`, then the cell resolution is first chosen
        uniformly at random from [0,..,self.max_resolution].

        EXAMPLES::

            >>> print(RHEALPixDGGS().random_cell()) # doctest: +SKIP
            S480586367780080

        """
        if resolution == None:
            resolution = randint(0, self.max_resolution)
        suid: list[str | int] = []
        suid.append(CELLS0[randint(0, 5)])
        for i in range(1, resolution + 1):
            suid.append(randint(0, self.N_side**2 - 1))
        return Cell(self, suid)

    def minimal_cover(
        self, resolution: int, points: list[tuple[float, float]], plane: bool = True
    ) -> list[Cell]:
        """
        Find the minimal set of resolution `resolution` cells that covers
        the list of points `points`.
        If `plane` = True, then assume `points` is a list of x-y
        coordinates in the planar DGGS.
        If `plane` = False, then assume `points` is a list of
        longitude-latitude coordinates in the ellipsoidal DGGS.
        This method will be made redundant by standard GIS rasterization
        tools that implement the rHEALPix projection.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> c1 = rdggs.cell(['N', 0, 2, 1])
            >>> c2 = rdggs.cell(['P', 7, 3, 3])
            >>> points = [c.nucleus() for c in [c1, c2]]
            >>> for r in range(5):
            ...     cover = sorted(rdggs.minimal_cover(r, points))
            ...     print([str(c) for c in cover])
            ['N', 'P']
            ['N0', 'P7']
            ['N02', 'P73']
            ['N021', 'P733']
            ['N0214', 'P7334']

        """
        cover: dict[str, Cell] = {}  # Use a dictionary to ignore repeated cells.
        for p in points:
            c = self.cell_from_point(resolution, p, plane=plane)
            if c is not None:
                # nuc = c.nucleus(plane=plane)
                cover[str(c)] = c  # (c, nuc[0], nuc[1])
        return list(cover.values())
        # Sort cells by nuclei y-coordinate and then by x-coordinate.
        # cover.sort(key=lambda x: (x[2], -x[1]), reverse=True)
        # return [t[0] for t in cover]


# Some common rHEALPix DGGSs.
WGS84_002 = RHEALPixDGGS(
    ellipsoid=WGS84_ELLIPSOID, north_square=0, south_square=0, N_side=2
)
WGS84_003 = RHEALPixDGGS(
    ellipsoid=WGS84_ELLIPSOID, north_square=0, south_square=0, N_side=3
)
WGS84_003_RADIANS = RHEALPixDGGS(
    ellipsoid=WGS84_ELLIPSOID_RADIANS, north_square=0, south_square=0, N_side=3
)
UNIT_003 = RHEALPixDGGS(ellipsoid=UNIT_SPHERE, north_square=0, south_square=0, N_side=3)
UNIT_003_RADIANS = RHEALPixDGGS(
    ellipsoid=UNIT_SPHERE_RADIANS, north_square=0, south_square=0, N_side=3
)
