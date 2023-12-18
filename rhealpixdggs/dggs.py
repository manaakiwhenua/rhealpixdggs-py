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
            R_A = 6374581.467096525
            a = 6378137.0
            b = 6356752.314140356
            e = 0.0578063088401125
            f = 0.003352810681182319
            lat_0 = 0
            lon_0 = 0
            radians = False
            sphere = False

NOTES::  .. Issue #1 was ..
        ellipsoid:
            R_A = 6374581.4671 *
            a = 6378137.0
            b = 6356752.314140356 *
            e = 0.0578063088401
            f = 0.003352810681182319

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
    ...        (0.0, 58.47067782962736),
    ...        (45.000000000000036, 58.47067782962734),
    ...        (89.99999999999996, 58.47067782962736)
    ...    ], [
    ...        (14.999999999999998, 26.438744923100096),
    ...        (45.0, 26.438744923100096),
    ...        (74.99999999999999, 26.438744923100096),
    ...        (105.00000000000001, 26.438744923100096)
    ...    ], [
    ...        (14.999999999999998, 3.560649871414923e-15),
    ...        (45.0, 3.560649871414923e-15),
    ...        (74.99999999999999, 3.560649871414923e-15),
    ...        (105.00000000000001, 3.560649871414923e-15)
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
            R_A = 6374581.467096525
            a = 6378137.0
            b = 6356752.314140356
            e = 0.0578063088401125
            f = 0.003352810681182319
            lat_0 = -37
            lon_0 = 174
            radians = False
            sphere = False

NOTES::  .. Issue #1 was ..
        ellipsoid:
            R_A = 6374581.4671 *
            a = 6378137.0
            b = 6356752.314140356
            e = 0.0578063088401 *
            f = 0.003352810681182319

    >>> print(rdggs.cell_from_point(1, p, plane=False))
    Q3

"""
# *****************************************************************************
#       Copyright (C) 2012 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http: //www.gnu.org/licenses/
# *****************************************************************************
# Import third-party modules.
from numpy import array, base_repr, ceil, log, pi, sign
from scipy import integrate

# Import standard modules.
from itertools import product
from random import uniform, randint
from colorsys import hsv_to_rgb

# Import my modules.
import rhealpixdggs.pj_rhealpix as pjr
import rhealpixdggs.projection_wrapper as pw
from rhealpixdggs.ellipsoids import (
    WGS84_ELLIPSOID,
    WGS84_ELLIPSOID_RADIANS,
    UNIT_SPHERE,
    UNIT_SPHERE_RADIANS,
)
from rhealpixdggs.utils import my_round
from numpy.testing import assert_allclose


class RHEALPixDGGS(object):
    """
    Represents an rHEALPix DGGS on a given ellipsoid.

    CLASS ATTRIBUTES:

    - `cells0` - A list of the resolution 0 cell IDs (strings).

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
      c is an element of `cells0` and (x, y) is the upper left corner point
      of the resolution 0 planar cell c.
    - `atomic_neighbors` - A dictionary with key-value pairs
      (n, {'up': a, 'down': b, 'left': c, 'right': d}),
      where n, a, b, c, and d are elements of `cells0` or
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

    # Level 0 cell IDs, which are anomalous.
    cells0 = ["N", "O", "P", "Q", "R", "S"]

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
        # The default layout by cells0 index is
        #
        # 0
        # 1 2 3 4
        # 5.
        #
        cells0 = RHEALPixDGGS.cells0
        ul_vertex = {  # Location for radius = 1
            cells0[0]: (-pi + self.north_square * pi / 2, 3 * pi / 4),
            cells0[1]: (-pi, pi / 4),
            cells0[2]: (-pi / 2, pi / 4),
            cells0[3]: (0, pi / 4),
            cells0[4]: (pi / 2, pi / 4),
            cells0[5]: (-pi + self.south_square * pi / 2, -pi / 4),
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
        # Neighbors of cells0[1], ..., cells0[4]
        an[cells0[1]] = {
            "left": cells0[4],
            "right": cells0[2],
            "down": cells0[5],
            "up": cells0[0],
        }
        an[cells0[2]] = {
            "left": cells0[1],
            "right": cells0[3],
            "down": cells0[5],
            "up": cells0[0],
        }
        an[cells0[3]] = {
            "left": cells0[2],
            "right": cells0[4],
            "down": cells0[5],
            "up": cells0[0],
        }
        an[cells0[4]] = {
            "left": cells0[3],
            "right": cells0[1],
            "down": cells0[5],
            "up": cells0[0],
        }
        # Neighbors of cells0[0] and cells0[5] depend on
        # volues of north_square and south_square, respectively.
        nn = self.north_square
        an[cells0[0]] = {
            "down": cells0[(nn + 0) % 4 + 1],
            "right": cells0[(nn + 1) % 4 + 1],
            "up": cells0[(nn + 2) % 4 + 1],
            "left": cells0[(nn + 3) % 4 + 1],
        }
        ss = self.south_square
        an[cells0[5]] = {
            "up": cells0[(ss + 0) % 4 + 1],
            "right": cells0[(ss + 1) % 4 + 1],
            "down": cells0[(ss + 2) % 4 + 1],
            "left": cells0[(ss + 3) % 4 + 1],
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
        result.append("    N_side = %s" % self.N_side)
        result.append("    north_square = %s" % self.north_square)
        result.append("    south_square = %s" % self.south_square)
        result.append("    max_areal_resolution = %s" % self.max_areal_resolution)
        result.append("    max_resolution = %s" % self.max_resolution)
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

    def healpix(self, u, v, inverse=False):
        """
        Return the HEALPix projection of point `(u, v)` (or its inverse if
        `inverse` = True) appropriate to this rHEALPix DGGS.

        EXAMPLES::

            >>> rdggs = UNIT_003_RADIANS
            >>> print(my_round(rdggs.healpix(-pi, pi/2), 14))
            (-2.35619449019234, 1.5707963267949)

        NOTES:: Issue #1 was ..
            (-2.35619449019234, 1.5707963267949001) *

        NOTE:

        Uses ``pj_healpix`` instead of the PROJ.4 version of HEALPix.
        """
        f = pw.Projection(ellipsoid=self.ellipsoid, proj="healpix")
        return f(u, v, inverse=inverse)

    def rhealpix(self, u, v, inverse=False, region="none"):
        """
        Return the rHEALPix projection of the point `(u, v)` (or its inverse if
        `inverse` = True) appropriate to this rHEALPix DGGS.

        EXAMPLES::

            >>> rdggs = UNIT_003_RADIANS
            >>> print(my_round(rdggs.rhealpix(0, pi/3), 14))
            (-1.858272006684, 2.06871881030324)

        NOTES:: Issue #1 was ..
            (-1.8582720066839999, 2.0687188103032401)

        NOTE:

        Uses ``pj_rhealpix`` instead of the PROJ.4 version of rHEALPix.
        """
        f = pw.Projection(
            ellipsoid=self.ellipsoid,
            proj="rhealpix",
            north_square=self.north_square,
            south_square=self.south_square,
            region=region,
        )
        return f(u, v, inverse=inverse, region=region)

    def combine_triangles(self, u, v, inverse=False, region="none"):
        """
        Return the combine_triangles() transformation of the point `(u, v)`
        (or its inverse if `inverse` = True) appropriate to the underlying
        ellipsoid.
        It maps the HEALPix projection to the rHEALPix projection.

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> p = (0, 0)
            >>> q = (-pi/4, pi/2)
            >>> print(rdggs.combine_triangles(*p))
            (0.0, 0.0)
            >>> print(my_round(rdggs.combine_triangles(*q), 14))
            (-2.35619449019234, 1.5707963267949)

        NOTES:: Issue #1 was ..
            (-2.35619449019234, 1.5707963267949001)

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

    def triangle(self, x, y, inverse=True):
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

    def xyz(self, u, v, lonlat=False):
        """
        Given a point `(u, v)` in the planar image of the rHEALPix projection,
        project it back to the ellipsoid and return its 3D rectangular
        coordinates.
        If `lonlat` = True, then assume `(u, v)` is a longitude-latitude
        point.

        EXAMPLES::

            >>> rdggs = UNIT_003_RADIANS
            >>> print(my_round(rdggs.xyz(0, pi/4, lonlat=True), 14))
            (0.70710678118655, 0.0, 0.70710678118655)

        NOTES:: Issue #1 was ..
            (0.70710678118655002, 0.0, 0.70710678118655002)

        """
        if lonlat:
            lam, phi = u, v
        else:
            lam, phi = self.rhealpix(u, v, inverse=True)
        return self.ellipsoid.xyz(lam, phi)

    def xyz_cube(self, u, v, lonlat=False):
        """
        Given a point `(u, v)` in the planar version of this rHEALPix DGGS,
        fold the rHEALPix image into a cube centered at the origin,
        and return the resulting point's 3D rectangular coordinates.
        If `lonlat` = True, then assume `(u, v)` is a longitude-latitude
        point.

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> print(my_round(rdggs.xyz_cube(0, 0), 14))
            (0.78539816339745, 0.0, -0.78539816339745)

        NOTES:: Issue #1 was ..
            (0.78539816339745006, 0.0, -0.78539816339745006)

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
        q = array(q) + (w / 2) * array((-1, -1, 1))
        return tuple(q)

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

    def grid(self, resolution):
        """
        Generator function for all the cells at resolution `resolution`.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> grid0 = rdggs.grid(0)
            >>> print([str(x) for x in grid0])
            ['N', 'O', 'P', 'Q', 'R', 'S']

        """
        suid = [RHEALPixDGGS.cells0[0]] + [0 for i in range(resolution)]
        c = self.cell(suid)
        yield c
        cs = c.successor(resolution)
        while cs:
            yield cs
            cs = cs.successor(resolution)

    def num_cells(self, res_1, res_2=None, subcells=False):
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

    def cell_width(self, resolution, plane=True):
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

    def cell_area(self, resolution, plane=True):
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

    def cell_from_point(self, resolution, p, plane=True):
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
            s0 = RHEALPixDGGS.cells0[0]
        elif (
            y > -R * 3 * pi / 4
            and y < -R * pi / 4
            and x > R * (-pi + ss * (pi / 2))
            and x < R * (-pi / 2 + ss * (pi / 2))
        ):
            s0 = RHEALPixDGGS.cells0[5]
        elif y >= -R * pi / 4 and y <= R * pi / 4 and x >= -R * pi and x < -R * pi / 2:
            s0 = RHEALPixDGGS.cells0[1]
        elif y >= -R * pi / 4 and y <= R * pi / 4 and x >= -R * pi / 2 and x < 0:
            s0 = RHEALPixDGGS.cells0[2]
        elif y >= -R * pi / 4 and y <= R * pi / 4 and x >= 0 and x < R * pi / 2:
            s0 = RHEALPixDGGS.cells0[3]
        elif y >= -R * pi / 4 and y <= R * pi / 4 and x >= R * pi / 2 and x < R * pi:
            s0 = RHEALPixDGGS.cells0[4]
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
            suid.append(self.child_order[(int(suid_row[i]), int(suid_col[i]))])
        return Cell(self, suid)

    def cell_from_region(self, ul, dr, plane=True):
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

    def cell_latitudes(self, resolution, phi_min, phi_max, nucleus=True, plane=True):
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
            -1.02050584399985
            -0.46144314900303
            -0.0
            0.46144314900303
            1.02050584399985
            1.5707963267949

        NOTES:: .. Issue @1 was ..
            -1.020505844
            -0.461443149003
            -0
            0.461443149003
            1.020505844
            1.57079632679

            >>> for phi in rdggs.cell_latitudes(1, -pi/2, pi/2, nucleus=False, plane=False):
            ...     print(my_round(phi, 14))
            -1.29836248988994
            -0.73083668811327
            -0.22457715619516
            0.22457715619516
            0.73083668811327
            1.29836248988994

        NOTES:: .. Issue @1 was ..
            -1.29836248989
            -0.730836688113
            -0.224577156195
            0.224577156195
            0.730836688113
            1.29836248989

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

    def cells_from_meridian(self, resolution, lam, phi_min, phi_max):
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
            new_cells = [c]
            if c.ellipsoidal_shape() in ["dart", "skew_quad"]:
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

    def cells_from_parallel(self, resolution, phi, lam_min, lam_max):
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
        PI = self.ellipsoid.pi()
        if start == end:
            if start.ellipsoidal_shape() == "cap" or lam_max - lam_min < PI / 2:
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

    def cells_from_region(self, resolution, ul, dr, plane=True):
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
            ur = self.cell_from_point(resolution, (dr[0], ul[1]), plane)
            dl = self.cell_from_point(resolution, (ul[0], dr[1]), plane)
            ul = self.cell_from_point(resolution, ul, plane)
            dr = self.cell_from_point(resolution, dr, plane)
            if ur is None or dl is None:
                return []
            if ul == dr:
                return [[ul]]
            # Starting from ul, collect cells from left to right and
            # then top to bottom, ending at dr.
            result = []
            row_start = ul
            row_end = ur
            while True:
                row = []
                current = row_start
                while current != row_end:
                    row.append(current)
                    current = current.neighbor("right", plane)
                row.append(current)
                result.append(row)
                if current == dr:
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

    def random_point(self, plane=True):
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
        c = self.cell([RHEALPixDGGS.cells0[n]])
        # Pick a random point in that cell.
        return c.random_point(plane=plane)

    def random_cell(self, resolution=None):
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
        suid = []
        suid.append(RHEALPixDGGS.cells0[randint(0, 5)])
        for i in range(1, resolution + 1):
            suid.append(randint(0, self.N_side**2 - 1))
        return Cell(self, suid)

    def minimal_cover(self, resolution, points, plane=True):
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
        cover = dict()  # Use a dictionary to ignore repeated cells.
        for p in points:
            c = self.cell_from_point(resolution, p, plane=plane)
            # nuc = c.nucleus(plane=plane)
            cover[str(c)] = c  # (c, nuc[0], nuc[1])
        cover = list(cover.values())
        return cover
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

from functools import total_ordering


@total_ordering
class Cell(object):
    """
    Represents a cell of the planar or ellipsoidal rHEALPix grid hierarchies.
    Cell identifiers are of the form (p_0, p_1,...,p_l), where p_0 is one of
    the characters 'A', 'B', 'C', 'D', 'E', 'F' and p_i for i > 0 is one of
    the integers 0, 1,..., N_side**2 - 1, where N_side is the instance
    attribute from RHEALPixDGGS (the number of children cells along a cell's
    side).

    INSTANCE ATTRIBUTES:

    - `rdggs` - The DGGS that the cell comes from.
    - `ellipsoid` - The underlying ellipsoid of the DGGS.
    - `N_side` - The N_side attribute of the DGGS
    - `suid` - The cell's ID (tuple).  SUID = spatially unique identifier.
      ('id' is a reserved word in Python)
    - `resolution` - The cell's resolution (nonnegative integer).

    NOTE:

    Several Cell methods have the keyword argument 'plane'.
    Setting it to True indicates that all input and output points and cells are
    to be interpreted as lying in the planar DGGS.
    Setting it to False indicates that they are to be interpreted as lying in
    the ellipsoidal DGGS.
    """

    @staticmethod
    def suid_from_index(rdggs, index, order="resolution"):
        """
        Return the suid of a cell from its index.
        The index is according to the cell ordering `order`,
        which can be 'resolution' (default) or 'post'.
        See the `index()` docstring for more details on orderings.
        For internal use.
        """
        from math import log  # Allows for different bases.

        if order == "post":
            # Compute suid from post order index one character at a time.
            suid = []
            p = index

            def num(k):
                return rdggs.num_cells(res_1=k, subcells=True)

            # Consider the tree T of all cells.
            # The indices of the cells in the six subtrees rooted at the
            # resolution 0 cells lie in the intervals
            # [0, num(0)), [num(0), 2*num(0)),..., [5*num(0), 6*num(0)),
            # respectively.
            # So computing p // num(0) gives us the first character of
            # the suid of our cell c.
            # Setting p = p % num(0) gives us the post order index
            # of c relative to the subtree rooted at cell suid[0].
            # The relative indices of the cells in the nine subtree rooted at
            # the nine children of suid[0] lie in the intervals
            # [0, num(1)), [num(1), 2*num(1)),...,
            # [(N_side**2 - 1)*num(1), N_side**2*num(1)),
            # respectively.
            # So computing p // num(1) gives us suid[1].
            # Repeating this procedure until p == num(i) - 1 for some i
            # (which will happen when i = max_resolution at the latest)
            # gives us all the characters of suid.
            for i in range(rdggs.max_resolution + 1):
                n = num(i)
                q, r = divmod(p, n)
                suid.append(q)
                p = r
                if p == n - 1:
                    break
            # Replace first digit with appropriate letter.
            suid[0] = RHEALPixDGGS.cells0[suid[0]]
            suid = tuple(suid)
        else:
            b = rdggs.N_side**2

            # Compute suid from level order index.
            def ind(k):
                """
                Return the level order index of the first cell at
                resolution k.
                """
                return int(6 * ((b**k - 1) / (b - 1)))

            # The cells at resolution L have indices in the interval
            # [ind(L), ind(L + 1)).
            k = int(log((b - 1) * (index / 6.0) + 1, b))
            # k = L  or L + 1. Find out which one.
            remainder = index - ind(k)
            if remainder >= 0:
                # Then k = L
                L = k
            else:
                # Then k = L + 1
                L = k - 1
                remainder = index - ind(L)
            # Now compute cell suid from remainder.
            suid = base_repr(remainder, b)
            # If necessary, prepend with zeros to get a length L + 1 string.
            suid = "0" * (L + 1 - len(suid)) + suid
            suid = [int(s) for s in suid]
            # Replace first digit with appropriate letter.
            suid[0] = RHEALPixDGGS.cells0[suid[0]]
            suid = tuple(suid)
        return suid

    def __init__(
        self, rdggs=WGS84_003, suid=None, level_order_index=None, post_order_index=None
    ):
        """
        Create a cell either from its suid or from its level order or
        post order index.

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> print(rdggs.max_resolution)
            1
            >>> c = Cell(rdggs, ['N', 3])
            >>> print(c)
            N3
            >>> loi = 6 + 6*9 - 1
            >>> c = Cell(rdggs, level_order_index=loi)
            >>> print(c)
            S8
            >>> c = Cell(rdggs, post_order_index=3)
            >>> print(c)
            N3

        """
        self.rdggs = rdggs
        self.ellipsoid = rdggs.ellipsoid
        self.N_side = rdggs.N_side
        self.suid = ()  # Spatially unique identifier of self.
        self.resolution = None  # Level of self in grid hierarchy.
        if suid is not None:
            # A little error checking.
            assert isinstance(suid, list) or isinstance(suid, tuple), (
                "Cell suid must be a list or tuple. Got %s." % suid
            )
            assert suid[0] in RHEALPixDGGS.cells0, "suid[0] must lie in %s. Got %s." % (
                RHEALPixDGGS.cells0,
                suid[0],
            )
            digits = set(range(self.N_side**2))
            for x in suid[1:]:
                assert x in digits, "Digits of suid must lie in %s" % digits
            assert (len(suid) > 0) and (
                len(suid) <= rdggs.max_resolution + 1
            ), "Need 0 < len(suid) <= %s. Got %s." % (rdggs.max_resolution + 1, suid)

            self.suid = [suid[0]] + [int(n) for n in suid[1:]]
            self.suid = tuple(self.suid)
        elif level_order_index is not None:
            self.suid = Cell.suid_from_index(
                self.rdggs, level_order_index, order="resolution"
            )
        elif post_order_index is not None:
            self.suid = Cell.suid_from_index(self.rdggs, post_order_index, order="post")
        self.resolution = len(self.suid) - 1

    def __bool__(self):
        return bool(self.suid)

    def __str__(self):
        if (self.rdggs.N_side) ** 2 < 10:
            s0 = self.suid[0]
            s = "".join(str(n) for n in self.suid[1:])
            return s0 + s
        else:
            # Comma separate suid entries.
            return "(" + self.suid[0] + str(self.suid)[4:]

    def __eq__(self, other):
        return (
            (other is not None)
            and (self.rdggs == other.rdggs)
            and (self.suid == other.suid)
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __le__(self, other):
        """
        The (strictly) less-than relation on cells.
        Derived from the post order traversal of the tree T of all cells
        defined in the `index()` docstring.
        Return True if (`self.suid < other.suid` and
        `self.suid` is not a prefix of `other.suid`) or
        `self` is a subcell of `other`.
        Here < is the lexicographic order.
        Returns False otherwise.
        """
        s = ",".join([str(x) for x in self.suid])
        t = ",".join([str(x) for x in other.suid])
        if (s <= t and not t.startswith(s)) or s.startswith(t):
            return True
        else:
            return False

    def index(self, order="resolution"):
        """
        Return the index of `self` when it's ordered according to `order`.
        Here `order` can be 'resolution' (default) or 'post'.
        Indices start at 0.
        The empty cell has index None.

        The ordering comes from the way of traversing the tree T
        of all cells defined as follows.
        The root of T is a non-cell place holder.
        The children of the root are the cells A < B < ... < F.
        The children of a cell in T with suid s are s0 < s1 < ... < sn,
        where n = self.N_side**2.

        The level order index of a nonempty cell is its position
        (starting from 0) in the level order traversal of T starting at cell A.

        The post order index of a nonempty cell is its position
        (starting from 0) in the post order traversal of T.

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> c = Cell(rdggs, ['N', 2])
            >>> print(c.index(order='resolution'))
            8
            >>> print(c.index(order='post'))
            2

        """
        L = self.resolution
        if not self.suid:
            return None
        s = list(self.suid)
        s[0] = RHEALPixDGGS.cells0.index(s[0])
        if order == "post":

            def num(k):
                return self.rdggs.num_cells(res_1=k, subcells=True)

            result = sum(s[i] * num(i) for i in range(L + 1)) + num(L) - 1
        else:
            # Calculate level order index.
            b = self.N_side**2
            n = len(s)
            result = self.rdggs.num_cells(res_1=0, res_2=L - 1) + sum(
                [s[n - 1 - i] * b**i for i in range(n)]
            )
        return result

    def suid_rowcol(self):
        """
        Return the pair of row- and column-suids of `self`, each as tuples.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> c = Cell(rdggs, ['N', 7, 3])
            >>> rsuid, csuid = c.suid_rowcol()
            >>> print(rsuid == ('N', 2, 1))
            True
            >>> print(csuid == ('N', 1, 0))
            True

        """
        suid_row = [self.suid[0]]
        suid_col = [self.suid[0]]
        for n in self.suid[1:]:
            row, col = self.rdggs.child_order[n]
            suid_row.append(row)
            suid_col.append(col)
        return tuple(suid_row), tuple(suid_col)

    def width(self, plane=True):
        """
        Return the width of this cell.
        If `plane` = False, then return None, because ellipsoidal cells
        don't have a fixed width.

        EXAMPLES::

            >>> c = Cell(UNIT_003, ('N', 8))
            >>> print(c)
            N8
            >>> c.width() == pi/2*3**(-1)
            True

        """
        return self.rdggs.cell_width(self.resolution, plane=plane)

    def area(self, plane=True):
        """
        Return the area of this cell.
        """
        return self.rdggs.cell_area(self.resolution, plane=plane)

    def successor(self, resolution=None):
        """
        Return the least resolution `resolution` cell greater than `self`.
        Note: `self` need not be a resolution `resolution` cell.

        EXAMPLES::

            >>> c = Cell(RHEALPixDGGS(), ('N', 8, 2))
            >>> print(c.successor())
            N83
            >>> print(c.successor(0))
            O
            >>> print(c.successor(1))
            O0
            >>> print(c.successor(3))
            N830

        """
        suid = list(self.suid)
        if resolution is None:
            resolution = self.resolution
        if resolution < self.resolution:
            # Truncate suid at resolution resolution and return its successor.
            return Cell(self.rdggs, suid[: resolution + 1]).successor()
        elif resolution > self.resolution:
            # Find the resolution self.resolution successor of suid
            # and pad it with zeros.
            suid = list(self.successor().suid) + [
                0 for i in range(resolution - self.resolution)
            ]
            return Cell(self.rdggs, suid)

        # Can now assume resolution = self.resolution.
        # First, find the greatest index i such that suid[i] != M.
        M = self.N_side**2 - 1
        greatest = 0
        for i in reversed(list(range(1, resolution + 1))):
            if suid[i] != M:
                greatest = i
                break

        # Second, increment suid[greatest] and append all zeros
        # if possible.
        if greatest == 0:
            # suid[greatest] is a letter.
            if suid[0] == RHEALPixDGGS.cells0[-1]:
                # End of the line. No successor.
                return None
            else:
                i = RHEALPixDGGS.cells0.index(suid[0])
                suid = [RHEALPixDGGS.cells0[i + 1]] + [0 for j in range(resolution)]
        else:
            # suid[greatest] is a number.
            suid = (
                suid[0:greatest]
                + [suid[greatest] + 1]
                + [0 for j in range(resolution - greatest)]
            )
        return Cell(self.rdggs, suid)

    def predecessor(self, resolution=None):
        """
        Return the greatest resolution `resolution` cell less than `self`.
        Note: `self` need not be a resolution `resolution` cell.

        EXAMPLES::

            >>> c = Cell(RHEALPixDGGS(), ('N', 0, 8))
            >>> print(c.predecessor())
            N07
            >>> print(c.predecessor(0))
            None
            >>> print(c.predecessor(1))
            None
            >>> print(c.predecessor(3))
            N088

        """
        M = self.N_side**2 - 1
        suid = list(self.suid)
        if resolution is None:
            resolution = self.resolution
        if resolution < self.resolution:
            # Return predecessor of suid[:resolution + 1]
            return Cell(self.rdggs, suid[: resolution + 1]).predecessor()
        elif resolution > self.resolution:
            # Return suid padded with Ms.
            suid = suid + [M for i in range(resolution - self.resolution)]
            return Cell(self.rdggs, suid)

        # Can now assume resolution = self.resolution.
        # Find the predecessor of suid.
        # First, find the greatest index i such that suid[i] != 0.
        greatest = 0
        for i in reversed(list(range(1, resolution + 1))):
            if suid[i] != 0:
                greatest = i
                break

        # Second, decrement suid[greatest] and append all Ms
        # if possible.
        if greatest == 0:
            # suid[greatest] is a letter.
            i = RHEALPixDGGS.cells0.index(suid[greatest])
            if i == 0:
                # End of the line. No predecessor.
                return None
            else:
                suid = [RHEALPixDGGS.cells0[i - 1]] + [M for i in range(resolution)]
        else:
            # nome[greatest] is a number > 0.
            suid = (
                suid[0:greatest]
                + [suid[greatest] - 1]
                + [M for i in range(resolution - greatest)]
            )
        return Cell(self.rdggs, suid)

    def subcell(self, other):
        """
        Subcell (subset) relation on cells.

        EXAMPLES::

            >>> a = Cell(RHEALPixDGGS(), ('N', 1))
            >>> b = Cell(RHEALPixDGGS(), ['N'])
            >>> print(a.subcell(b))
            True
            >>> print(b.subcell(a))
            False

        """
        s = ",".join([str(x) for x in self.suid])
        t = ",".join([str(x) for x in other.suid])
        return s.startswith(t)

    def subcells(self, resolution=None):
        """
        Generator function for the set of all resolution `resolution` subcells
        of this cell.
        If `resolution=None`, then return a generator function for the children
        of this cell.

        EXAMPLES::

            >>> c = Cell(RHEALPixDGGS(), ['N'])
            >>> print([str(cell) for cell in c.subcells()])
            ['N0', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8']

        """
        L = self.resolution
        if resolution is None:
            resolution = L + 1
        if resolution < L:
            return  # Stop iteration
        if resolution == L:
            yield self
            return
        N = self.N_side
        for t in product(list(range(N**2)), repeat=resolution - L):
            yield Cell(self.rdggs, list(self.suid) + list(t))

    def ul_vertex(self, plane=True):
        """
        If `plane` = True, then return the upper left vertex of this
        planar cell.
        If `plane` = False, then return the projection onto the ellipsoid
        of the planar upper left vertex.
        Note that for polar cells, this projection is not necessarily the
        northwest vertex.
        For the latter vertex use nw_vertex().

        WARNING: The upper left vertex of a cell might not lie in the cell,
        because not all cells contain their boundary.

        EXAMPLES::

            >>> c = Cell(UNIT_003, ['N', 0])
            >>> print(c.ul_vertex() == (-pi, 3*pi/4))
            True

        """
        # Call this cell c.
        # Find the location of the resolution 0 cell c0 containing c.
        x0, y0 = self.rdggs.ul_vertex[self.suid[0]]
        resolution = self.resolution

        # The column and row SUIDs of c give the the horizontal and vertical
        # distances, respectively, between the ul_vertex of c0 and
        # the ul_vertex of c as fractions of the width of c0.
        suid_row, suid_col = self.suid_rowcol()
        N = self.N_side
        dx = sum(
            N ** (resolution - i) * suid_col[i] for i in range(1, resolution + 1)
        ) * N ** (-resolution)
        dy = sum(
            N ** (resolution - i) * suid_row[i] for i in range(1, resolution + 1)
        ) * N ** (-resolution)
        # Multiplied by N**(-resolution) at the end instead of at each term to
        # avoid rounding errors.

        # Use dx, dy, and the width of c0 to calculate the ul_vertex of c.
        E = self.rdggs
        x = x0 + E.cell_width(0) * dx
        y = y0 - E.cell_width(0) * dy

        # Project onto ellipsoid if necessary.
        if not plane:
            x, y = self.rdggs.rhealpix(x, y, inverse=True)
        return x, y

    def nw_vertex(self, plane=True):
        """
        If `plane` = False, then return the northwest vertex of this
        ellipsoidal cell.
        If `plane` = True, then return the projection onto the plane
        of the ellipsoidal northwest vertex.
        On quad cells and cap cells, this function returns the same output
        as ul_vertex().
        On skew quad cells and dart cells, this function returns output
        different from ul_vertex().

        WARNING: The northwest vertex of a cell might not lie in the cell,
        because not all cells contain their boundary.

        EXAMPLES::
            >>> rdggs = RHEALPixDGGS()
            >>> c = rdggs.cell(['P', 5, 7]) # Quad cell.
            >>> print(my_round(c.ul_vertex(plane=True), 14))
            (-2225148.7007489, -556287.1751872245)
            >>> print(my_round(c.nw_vertex(plane=True), 14))
            (-2225148.7007489, -556287.1751872245)

        NOTES:: .. Issue #1 was ..
            (-2225148.7007489, -556287.17518722452) * for both above tests

            >>> c = rdggs.cell(['S', 4])  # Cap cell.
            >>> print(my_round(c.ul_vertex(plane=True), 14))
            (-16688615.255616743, -8344307.627808371)
            >>> print(my_round(c.nw_vertex(plane=True), 14))
            (-16688615.255616743, -8344307.627808371)

        NOTES:: .. Issue #1 was ..
            (-16688615.255616743, -8344307.6278083706) * for both above tests

            >>> c = rdggs.cell(['N', 4, 3]) # Skew quad cell.
            >>> print(my_round(c.ul_vertex(plane=True), 14))
            (-16688615.255616743, 10569456.32855727)
            >>> print(my_round(c.nw_vertex(plane=True), 14))
            (-15576040.905242294, 10569456.32855727)

            >>> c = rdggs.cell(['S', 4, 3])  # Skew quad cell.
            >>> print(my_round(c.ul_vertex(plane=True), 14))
            (-16688615.255616743, -9456881.97818282)
            >>> print(my_round(c.nw_vertex(plane=True), 14))
            (-16688615.255616743, -10569456.32855727)

        NOTES:: .. Issue #1 was ..
            (-16688615.255616743, -9456881.9781828206) *

            >>> c = rdggs.cell(['N', 6, 2])  # Dart cell.
            >>> print(my_round(c.ul_vertex(plane=True), 14))
            (-17801189.605991192, 8344307.627808372)
            >>> print(my_round(c.nw_vertex(plane=True), 14))
            (-16688615.255616743, 8344307.627808372)

        NOTES:: .. Issue #1 was ..
            (-16688615.255616743, 8344307.6278083716) * for both above tests

            >>> c = rdggs.cell(['S', 6, 2])  # Dart cell.
            >>> print(my_round(c.ul_vertex(plane=True), 14))
            (-17801189.605991192, -11682030.678931719)
            >>> print(my_round(c.nw_vertex(plane=True), 14))
            (-16688615.255616743, -12794605.029306168)

        """
        v = self.vertices(plane=True)  # Planar vertices.
        shape = self.ellipsoidal_shape()
        if shape == "quad" or shape == "cap":
            # Northwest vertex is the upper left vertex.
            result = v[0]
        elif shape == "skew_quad":
            # Could project cell to ellipsoid and then sort vertices
            # by latitude and then by longitude, but handling the possible
            # rounding errors gets fiddly.
            # Instead determine northwest vertex based on what HEALPix triangle
            # the cell lies in.
            rdggs = self.rdggs
            triangle, region = rdggs.triangle(*self.nucleus(plane=True))
            if region == "north_polar":
                ns = rdggs.north_square
                i = (triangle - ns) % 4
                result = v[-i]
            else:
                ss = rdggs.south_square
                i = (triangle - ss) % 4
                result = v[i]
        else:
            # shape == 'dart':
            # Map cell to ellipsoid and get the polewards vertex.
            ev = [self.rdggs.rhealpix(*vv, inverse=True) for vv in v]
            i = max((abs(ev[j][1]), j) for j in range(4))[1]
            if self.region() == "north_polar":
                # Northwest vertex is the polewards vertex.
                result = v[i]
            else:
                # Northwest vertex is one step clockwise
                # from the polewards vertex.
                result = v[(i + 1) % 4]
        if not plane:
            result = self.rdggs.rhealpix(*result, inverse=True)
        return result

    def nucleus(self, plane=True):
        """
        Return the nucleus and vertices of this planar or ellipsoidal cell
        in the order (nucleus, upper left corner, lower left corner,
        lower right corner, upper right corner) with reference to the
        planar cell.
        The output for ellipsoidal cells is the projection onto the ellipsoid
        of the output for planar cells.  In particular, while the
        nucleus of a planar cell is its centroid, the nucleus
        of an ellipsoidal cell is not its centroid.
        To compute the centroid of a cell, use centroid() below.

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> c = rdggs.cell(['N'])
            >>> print(my_round(c.nucleus(), 14))
            (-2.35619449019234, 1.5707963267949)

        NOTES:: .. Issue #1 was ..
            (-2.35619449019234, 1.5707963267949001) *

        """
        ul = self.ul_vertex(plane=True)
        w = self.width()
        result = (ul[0] + w / 2, ul[1] - w / 2)
        if not plane:
            # Project to ellipsoid.
            result = self.rdggs.rhealpix(*result, inverse=True)
        return result

    def vertices(self, plane=True, trim_dart=False):
        """
        If `plane` = True, then assume this cell is planar and return
        its four vertices in the order (upper left corner, upper right corner,
        lower right corner, lower left corner).
        If `plane` = False, then assume this cell is ellipsoidal and return
        the projection of the planar vertices in the order
        (northwest, northeast, southeast, southwest).
        If `plane` = False, this cell is a dart cell, and
        `trim_dart` = True, then remove the one non-vertex point from
        the output. (Dart cells only have three vertices.)

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> c = rdggs.cell(['N'])
            >>> for p in c.vertices():
            ...     print(my_round(p, 14))
            (-3.14159265358979, 2.35619449019234)
            (-1.5707963267949, 2.35619449019234)
            (-1.5707963267949, 0.78539816339745)
            (-3.14159265358979, 0.78539816339745)

            >>> rdggs = WGS84_003
            >>> c = rdggs.cell(['N', 0])
            >>> for p in c.vertices(plane=False):
            ...     print(my_round(p, 13))
            (89.9999999999999, 74.3906909487906)
            (120.0, 41.8738577422094)
            (90.0, 41.8738577422094)
            (60.0, 41.8738577422094)
            >>> for p in c.vertices(plane=False, trim_dart=True):
            ...     print(my_round(p, 13))
            (89.9999999999999, 74.3906909487906)
            (120.0, 41.8738577422094)
            (60.0, 41.8738577422094)

            >>> c = rdggs.cell(['S', 0])
            >>> for p in c.vertices(plane=False):
            ...     print(my_round(p, 14))
            (149.99999999999997, -41.87385774220941)
            (-180.0, -41.87385774220941)
            (-150.0, -41.87385774220941)
            (-180.0, -74.39069094879065)
            >>> for p in c.vertices(plane=False, trim_dart=True):
            ...     print(my_round(p, 14))
            (149.99999999999997, -41.87385774220941)
            (-150.0, -41.87385774220941)
            (-180.0, -74.39069094879065)
        """
        ul = self.ul_vertex(plane=True)
        w = self.width()
        ur = (ul[0] + w, ul[1])
        dr = (ul[0] + w, ul[1] - w)
        dl = (ul[0], ul[1] - w)
        result = [ul, ur, dr, dl]
        if not plane:
            # Reorder result so that it starts with the northwest vertex.
            # Clockwise ordering is preserved when mapping from plane to
            # ellipsoid.
            nw = self.nw_vertex(plane=True)
            i = result.index(nw)
            result = result[i:] + result[:i]
            # Project to ellipsoid.
            region = self.region()
            result = [
                self.rdggs.rhealpix(*p, inverse=True, region=region) for p in result
            ]
            if trim_dart and self.ellipsoidal_shape() == "dart":
                # Remove non-vertex point.
                if self.region() == "north_polar":
                    result.pop(2)
                else:
                    result.pop(1)
        return result

    def xy_range(self):
        """
        Return the x- and y-coordinate extremes of the planar version of
        this cell in the format ((x_min, x_max), (y_min, y_max)).

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> c = rdggs.cell(['N'])
            >>> c.xy_range() == ((-pi, -pi/2), (pi/4, 3*pi/4))
            True

        """
        ul = self.ul_vertex(plane=True)
        w = self.width()
        x_min = ul[0]
        x_max = x_min + w
        y_max = ul[1]
        y_min = y_max - w
        return (x_min, x_max), (y_min, y_max)

    def boundary(self, n=2, plane=True, interior=False):
        """
        Return a list of `4*n - 4` boundary points of this cell,
        `n` on each edge, where `n` >= 2.
        List the points in clockwise order starting from the cell's upper left
        corner if `plane` = True, or from the cell's northwest corner
        if `plane` = False.

        If `n` = 2, then the output is the same as vertices().
        If `interior` = True, then push the boundary points slighly into the
        interior of the cell, which is convenient for some graphics methods.

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> c = rdggs.cell(['N', 6])
            >>> c.boundary(n=2, plane=True) == c.vertices(plane=True)
            True
            >>> for p in c.boundary(n=3, plane=True):
            ...     print(my_round(p, 14))
            (-3.14159265358979, 1.30899693899575)
            (-2.87979326579064, 1.30899693899575)
            (-2.61799387799149, 1.30899693899575)
            (-2.61799387799149, 1.0471975511966)
            (-2.61799387799149, 0.78539816339745)
            (-2.87979326579064, 0.78539816339745)
            (-3.14159265358979, 0.78539816339745)
            (-3.14159265358979, 1.0471975511966)

        NOTES::  .. Issue #1 was ..
            (-3.14159265358979, 1.3089969389957501) *
            (-2.87979326579064, 1.3089969389957501) *
            (-2.61799387799149, 1.3089969389957501) *
            (-2.61799387799149, 1.0471975511966001) *
            (-2.61799387799149, 0.78539816339745006) *
            (-2.87979326579064, 0.78539816339745006) *
            (-3.14159265358979, 0.78539816339745006) *
            (-3.14159265358979, 1.0471975511966001) *


            >>> for p in c.boundary(n=3, plane=False):
            ...     print(my_round(p, 14))
            (-180.0, 74.35752898700072)
            (-157.50000000000003, 58.41366190347208)
            (-150.0, 41.8103148957786)
            (-165.00000000000003, 41.8103148957786)
            (-180.0, 41.8103148957786)
            (165.0, 41.8103148957786)
            (149.99999999999997, 41.8103148957786)
            (157.49999999999997, 58.41366190347208)

        NOTES::  .. Issue #1 was ..
            (-180.0, 74.35752898700072)
            (-157.50000000000003, 58.413661903472082) *
            (-150.0, 41.810314895778603) *
            (-165.00000000000003, 41.810314895778603) *
            (-180.0, 41.810314895778603) *
            (165.0, 41.810314895778603) *
            (149.99999999999997, 41.810314895778603) *
            (157.49999999999997, 58.413661903472082) *

        """
        ul = self.ul_vertex(plane=True)
        w = self.width(plane=True)
        if n < 2:
            n == 2
        if interior:
            eps = w / 10000  # A smidgen.
        else:
            eps = 0
        delta = (w - 2 * eps) / (n - 1)
        point = (ul[0] + eps, ul[1] - eps)
        result = [point]
        for direction in [(1, 0), (0, -1), (-1, 0), (0, 1)]:
            for j in range(1, n):
                temp = array(point) + j * delta * array(direction)
                result.append(tuple(temp))
            point = result[-1]
        # Remove the last point because it's the first point.
        result.pop()
        if not plane:
            # Reorder result so that it starts with the northwest vertex.
            # Clockwise ordering is preserved when mapping from plane to
            # ellipsoid.
            v = self.vertices(plane=True)
            nw = self.nw_vertex(plane=True)
            i = v.index(nw)  # Index of northwest vertex in planar vertex list
            i = (n - 1) * i  # Index of northwest vertex in result.
            result = result[i:] + result[:i]
            # Project to ellipsoid.
            region = self.region()
            result = [
                self.rdggs.rhealpix(*p, inverse=True, region=region) for p in result
            ]
        return result

    def interior(self, n=2, plane=True, flatten=False):
        """
        Return an `n` x `n` matrix of interior points of this cell.
        If the cell is planar, space the interior points on a regular
        square grid.
        List the points in standard, row-major matrix order.
        If the cell is ellipsoidal, project the matrix of points to the
        ellipsoid (longitude-latitude points).
        If `flatten` = True, then flatten the matrix into a one dimensional
        array of pairs.

        EXAMPLES::

            >>> rdggs = UNIT_003
            >>> c = rdggs.cell(['N'])
            >>> for p in c.interior(n=2, plane=False, flatten=True):
            ...     print(my_round(p, 13))
            (90.0, 41.8103801453539)
            (-180.0, 41.8103801453539)
            (-0.0, 41.8103801453539)
            (-90.0, 41.8103801453539)

            >>> all([c.contains(p) for p in c.interior(n=5, plane=True, flatten=True)])
            True

        """
        ul = self.ul_vertex(plane=True)
        w = self.width(plane=True)
        eps = 1e-6
        delta = (w - 2 * eps) / (n - 1)

        def g(x, y):
            if plane:
                return (x, y)
            else:
                return self.rdggs.rhealpix(x, y, inverse=True)

        if flatten:
            result = [
                g(ul[0] + eps + delta * j, ul[1] - eps - delta * i)
                for j in range(n)
                for i in range(n)
            ]
        else:
            result = [
                [g(ul[0] + eps + delta * j, ul[1] - eps - delta * i) for j in range(n)]
                for i in range(n)
            ]
        return result

    def contains(self, p, plane=True):
        """
        Return True if this cell contains point `p`, and return False
        otherwise.

        EXAMPLES::

            >>> rdggs = WGS84_003_RADIANS
            >>> p = (pi/4, 0)
            >>> c = rdggs.cell_from_point(2, p, plane=False)
            >>> print(c)
            Q44
            >>> print(c.contains(p, plane=False))
            True

        """
        # Calling cell_from_point() is the simplest (but maybe not the
        # fastest) way to do this check, because given a planar cell,
        # deciding which of its edges it contains involves several cases,
        # because the rHEALPix map projection does not contain all of its
        # edges.
        return self.rdggs.cell_from_point(self.resolution, p, plane=plane) == self

    def intersects_meridian(self, lam):
        """
        Return True if this ellipsoidal cell's boundary intersects the
        meridian of longitude `lam`, and return False otherwise.

        EXAMPLES::

            >>> rdggs = WGS84_003_RADIANS
            >>> c = rdggs.cell(['N', 6])
            >>> print(c.intersects_meridian(-pi))
            True
            >>> print(c.intersects_meridian(-pi/2))
            False

        """
        if self.ellipsoidal_shape() == "cap":
            return True
        # Not a cap cell.
        vertices = self.vertices(plane=False)
        lon_min = min([v[0] for v in vertices])
        lon_max = max([v[0] for v in vertices])
        PI = self.ellipsoid.pi()
        if abs(lon_min - lon_max) > PI:
            # Exceptional case of a dart cell with nucleus at longitude -pi.
            # The cell straddles the -pi/pi boundary.
            lon_min = -lon_max
            return lon_max <= lam or lam <= lon_min
        else:
            # Typical case.
            return lon_min <= lam and lam <= lon_max

    def intersects_parallel(self, phi):
        """
        Return True if this cell's boundary intersects the parallel of latitude
        `phi`, and return False otherwise.
        """
        # Cell's boundary intersects parallel iff its extreme latitudes lie on
        # opposite sides of parallel.
        vertices = self.vertices(plane=False)
        lat_min = min([v[1] for v in vertices])
        lat_max = max([v[1] for v in vertices])
        if self.ellipsoidal_shape() == "cap":
            if self.region() == "north_polar":
                return phi >= lat_min
            else:
                return phi <= lat_max
        else:
            return lat_min <= phi and lat_max >= phi

    def overlaps(self, other_cell):
        """
        Determines whether two DGGS cells overlap.
        Where cells are of different resolution, they will have different suid lengths. The zip function truncates the longer
        to be the same length as the shorter, producing two lists for comparison. If these lists are equal, the cells overlap.
        :param cell_one: the first DGGS cell
        :param cell_two: the second DGGS cell
        :return: True if overlaps
        """
        assert self.suid is not tuple()  # cell cannot be empty
        for i, j in zip(self.suid, other_cell.suid):
            if i != j:
                return False
        return True

    def region_overlaps(self, region: list):
        """
        Determine whether a cell overlaps with any cell in a list of cells
        :param cell: a DGGS cell
        :param region: a list of DGGS cells
        :return: True if any overlapping cells
        """
        for component_cell in region:
            if self.overlaps(component_cell):
                return True
        return False

    def region(self):
        """
        Return the region of this cell: 'equatorial', 'north_polar', or
        'south_polar'.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> print(Cell(rdggs, ['P', 2]).region())
            equatorial
            >>> print(Cell(rdggs, ['N', 2]).region())
            north_polar

        """
        if self.suid[0] == RHEALPixDGGS.cells0[0]:
            return "north_polar"
        elif self.suid[0] == RHEALPixDGGS.cells0[5]:
            return "south_polar"
        else:
            return "equatorial"

    def ellipsoidal_shape(self):
        """
        Return the shape of this cell ('quad', 'cap', 'dart', or
        'skew_quad') when viewed on the ellipsoid.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> print(Cell(rdggs, ['P', 2]).ellipsoidal_shape())
            quad
            >>> print(Cell(rdggs, ['N', 2]).ellipsoidal_shape())
            dart

        """
        suid = self.suid
        if suid[0] in RHEALPixDGGS.cells0[1:5]:
            return "quad"
        N = self.N_side
        # Cap check.
        cap = True
        if N % 2 != 1:
            cap = False
        if cap:
            for n in suid[1:]:
                if n != (N**2 - 1) // 2:
                    cap = False
                    break
        if cap:
            return "cap"
        # Dart check 1.
        dart = True
        S = set([i * (N + 1) for i in range(N)])
        for n in suid[1:]:
            if n not in S:
                dart = False
                break
        if dart:
            return "dart"
        # Dark check 2.
        dart = True
        S = set([(i + 1) * (N - 1) for i in range(N)])
        for n in suid[1:]:
            if n not in S:
                dart = False
                break
        if dart:
            return "dart"
        # Must be a skew quad then.
        return "skew_quad"

    def centroid(self, plane=True):
        """
        Return the centroid of this planar or ellipsoidal cell.

        EXAMPLES::

            >>> rdggs = RHEALPixDGGS()
            >>> c = Cell(rdggs, ['P', 0, 2])
            >>> centroid = c.centroid()
            >>> nucleus = c.nucleus()
            >>> print(centroid == nucleus)
            True

        """
        if plane:
            # Then this cell's centroid is its nucleus.
            return self.nucleus(plane=True)

        # This cell is ellipsoidal.
        # So we have to do some work.
        nucleus = self.nucleus(plane=False)
        vertices = self.vertices(plane=False)
        shape = self.ellipsoidal_shape()
        if shape == "cap":
            return nucleus
        if shape == "quad":
            lam_bar = nucleus[0]
            phi_bar = sum([v[1] for v in vertices]) / 4
            return lam_bar, phi_bar
        planar_vertices = self.vertices(plane=True)
        x1 = min([v[0] for v in planar_vertices])
        x2 = max([v[0] for v in planar_vertices])
        y1 = min([v[1] for v in planar_vertices])
        y2 = max([v[1] for v in planar_vertices])
        area = (x2 - x1) ** 2

        def lam(x, y):
            return self.rdggs.rhealpix(x, y, inverse=True)[0]

        def phi(x, y):
            return self.rdggs.rhealpix(x, y, inverse=True)[1]

        if shape == "dart":
            lam_bar = nucleus[0]
            phi_bar = (1 / area) * integrate.dblquad(
                phi, y1, y2, lambda x: x1, lambda x: x2
            )[0]
            return lam_bar, phi_bar
        # Now shape == 'skew_quad'.
        # phi_bar formula same as dart case.
        phi_bar = (1 / area) * integrate.dblquad(
            phi, y1, y2, lambda x: x1, lambda x: x2
        )[0]
        # lam_bar formula changes.
        # Option 1 (clean, possibly slow):
        # Compute lam_bar by numerical integration.
        lam_bar = (1 / area) * integrate.dblquad(
            lam, y1, y2, lambda x: x1, lambda x: x2
        )[0]
        # Option 2 (messy, possibly fast):
        # Evaluate the integral symbolically and then plug in values.
        # w = x2 - x1 # Cell width.
        # R_A = self.rdggs.ellipsoid.R_A
        # hx0, hy0 = self.rdggs.healpix(*nucleus)
        # # x and y extremes of the HEALPix projection of this cell's interior:
        # hx1 = hx0 - w/2
        # hx2 = hx0 + w/2
        # # Without loss of generality, force HEALPix y coordinates into
        # # the northern hemisphere:
        # hy1 = abs(hy0) - w/2
        # hy2 = abs(hy0) + w/2
        # # Compute xc.
        # cap_number = floor(2*hx0/(pi*R_A) + 2)
        # if cap_number >= 4:
        #     # Rounding error.
        #     cap_number = 3
        # xc = -3*pi/4 + (pi/2)*cap_number
        # integral = lambda x, y: (pi/8)*x*(2*R_A*xc - x)*\
        #            log(1 - 2*y/(pi*R_A)) + xc*x*y
        # lam_bar = (1/area)*\
        #           (integral(hx2, hy2) - integral(hx1, hy2) -\
        #            integral(hx2, hy1) + integral(hx1, hy1))
        # if not self.rdggs.ellipsoid.radians:
        #     # Convert to degrees.
        #     lam_bar = rad2deg(lam_bar)
        return lam_bar, phi_bar

    def rotate_entry(self, x, quarter_turns):
        """
        Let N = self.N_side and rotate the N x N matrix of subcell numbers ::

            0        1          ... N - 1
            N        N+1        ... 2*N - 1
            ...
            (N-1)*N  (N-1)*N+1  ... N**2-1

        anticlockwise by `quarter_turns` quarter turns to obtain a
        new table with entries f(0), f(1), ..., f(N**2 - 1) read from
        left to right and top to bottom.
        Given entry number `x` in the original matrix, return `f(x)`.
        Used in rotate().

        INPUT:

        - `x` - A letter from RHEALPixDGGS.cells0 or one of the integers
          0, 1, ..., N**2 - 1.
        - `quarter_turns` - 0, 1, 2, or 3.

        EXAMPLES::

            >>> c = Cell(RHEALPixDGGS(), ['P', 2])
            >>> print([c.rotate_entry(0, t) for t in range(4)])
            [0, 2, 8, 6]

        NOTES:

        Operates on letters from RHEALPixDGGS.cells0 too.
        They stay fixed under f.
        Only depends on `self` through `self.N_side`.
        """
        N = self.N_side
        # Original matrix of subcell numbers as drawn in the docstring.
        A = self.rdggs.child_order
        # Function (written as a dictionary) describing action of rotating A
        # one quarter turn anticlockwise.
        f = dict()
        for i in range(N):
            for j in range(N):
                n = A[(i, j)]
                f[n] = A[(j, N - 1 - i)]
        # Level 0 cell names stay the same.
        for c in self.rdggs.cells0:
            f[c] = c

        quarter_turns = quarter_turns % 4
        if quarter_turns == 1:
            return f[x]
        elif quarter_turns == 2:
            return f[f[x]]
        elif quarter_turns == 3:
            return f[f[f[x]]]
        else:
            return x

    def rotate(self, quarter_turns):
        """
        Return the cell that is the result of rotating this cell's
        resolution 0 supercell by `quarter_turns` quarter turns anticlockwise.
        Used in neighbor().

        EXAMPLES::

            >>> c = Cell(RHEALPixDGGS(), ['N', 0])
            >>> print([str(c.rotate(t)) for t in range(4)])
            ['N0', 'N2', 'N8', 'N6']

        """
        suid = [self.rotate_entry(x, quarter_turns) for x in self.suid]
        return Cell(self.rdggs, suid)

    def neighbor(self, direction, plane=True):
        """
        Return this cell's (edge) neighbor in the given direction.
        If `plane` = True, then the direction is one of the strings
        'up', 'right', 'down', 'left', which indicates the desired neighbor
        relative to x-y coordinates
        in the following planar neighbor diagram, (drawn for self.N_side = 3)
        where `self` is the middle cell ::

                            up
                          *-----*
                          |     |
                          |     |
                          |     |
                    *-----*-----*-----*
                    |     | 012 |     |
               left |     | 345 |     | right
                    |     | 678 |     |
                    *-----*-----*-----*
                          |     |
                          |     |
                          |     |
                          *-----*
                           down

        If `plane` = False, then the direction is relative to
        longitude-latitude coordinates and is one of the strings
        'west', 'east', 'north', 'south' for a quad or skew quad cell;
        'west', 'east', 'southwest', 'southeast' for a northern dart cell;
        'west', 'east', 'northwest', 'northeast' for a southern dart cell;
        'south_0', 'south_1', 'south_2', 'south_3' for a northern cap cell;
        'north_0', 'north_1', 'north_2', 'north_3' for a southern cap cell;
        For a cap cell, neighbor directions are numbered in increasing
        longitude, so that the longitude of the (nucleus of) north_0 is less
        than the longitude of north_1 is less than the longitude of north_2 is
        less than the longitude of north_3, and the longitude of the south_0 is
        less than the longitude of south_1, etc.

        The tricky part in the planar scenario is that the neighbor
        relationships of the six resolution 0 cells is determined by the
        positions of those cells on the surface of a cube,
        one cell on each face, and not on a plane.
        So sometimes rotating cells is needed to compute neighbors.

        Return None if the given direction is invalid for this cell.

        EXAMPLES::

            >>> c = Cell(RHEALPixDGGS(), ['N', 0])
            >>> print(c.neighbor('down'))
            N3

        """
        if plane:
            if direction not in {"left", "right", "down", "up"}:
                return None
            an = self.rdggs.atomic_neighbors
            # First, compute the neighbor of self naively, that is,
            # without considering rotations.
            self_suid = self.suid
            neighbor_suid = []
            N = self.N_side
            up_border = set(range(N))
            down_border = set([(N - 1) * N + i for i in range(N)])
            left_border = set([i * N for i in range(N)])
            right_border = set([(i + 1) * N - 1 for i in range(N)])
            border = {
                "left": left_border,
                "right": right_border,
                "up": up_border,
                "down": down_border,
            }
            crossed_all_borders = False
            # Scan from the back to the front of suid.
            for i in reversed(list(range(len(self_suid)))):
                n = self_suid[i]
                if crossed_all_borders:
                    neighbor_suid.append(n)
                else:
                    neighbor_suid.append(an[n][direction])
                    if n not in border[direction]:
                        crossed_all_borders = True
            neighbor_suid.reverse()
            neighbor = Cell(self.rdggs, neighbor_suid)

            # Second, rotate the neighbor if necessary.
            # If self is a polar cell and neighbor is not, or vice versa,
            # then rotate neighbor accordingly.
            self0 = self_suid[0]
            neighbor0 = neighbor_suid[0]
            cells0 = RHEALPixDGGS.cells0
            if (
                (self0 == cells0[5] and neighbor0 == an[self0]["left"])
                or (self0 == an[cells0[5]]["right"] and neighbor0 == cells0[5])
                or (self0 == cells0[0] and neighbor0 == an[self0]["right"])
                or (self0 == an[cells0[0]]["left"] and neighbor0 == cells0[0])
            ):
                neighbor = neighbor.rotate(1)
            elif (
                (self0 == cells0[5] and neighbor0 == an[self0]["down"])
                or (self0 == an[cells0[5]]["down"] and neighbor0 == cells0[5])
                or (self0 == cells0[0] and neighbor0 == an[self0]["up"])
                or (self0 == an[cells0[0]]["up"] and neighbor0 == cells0[0])
            ):
                neighbor = neighbor.rotate(2)
            elif (
                (self0 == cells0[5] and neighbor0 == an[self0]["right"])
                or (self0 == an[cells0[5]]["left"] and neighbor0 == cells0[5])
                or (self0 == cells0[0] and neighbor0 == an[self0]["left"])
                or (self0 == an[cells0[0]]["right"] and neighbor0 == cells0[0])
            ):
                neighbor = neighbor.rotate(3)
        else:
            # Ellipsoid.
            # Call neighbors() to do all the work.
            neighbors = self.neighbors(plane=False)
            try:
                neighbor = neighbors[direction]
            except KeyError:
                # Invalid direction given.
                neighbor = None
        return neighbor

    def neighbors(self, plane=True):
        """
        Return this cell's planar or ellipsoidal (edge) neighbors
        as a dictionary whose keys are the directions of the neighbors.
        See neighbor() for a list of valid directions.

        EXAMPLES::

            >>> c = Cell(RHEALPixDGGS(), ['N', 0])
            >>> for k, v in sorted(c.neighbors().items()):
            ...     print(k, v)
            ...
            down N3
            left R0
            right N1
            up Q2

        """
        plane_neighbors = dict()
        for d in ["left", "right", "down", "up"]:
            plane_neighbors[d] = self.neighbor(d, "plane")
        if plane:
            return plane_neighbors
        # Ellipsoid case.
        result = dict()
        shape = self.ellipsoidal_shape()
        if shape == "quad":
            result["north"] = plane_neighbors["up"]
            result["south"] = plane_neighbors["down"]
            result["west"] = plane_neighbors["left"]
            result["east"] = plane_neighbors["right"]
        elif shape == "cap":
            # Sort neighbors by nuclei longitudes.
            nuc_cell = []
            for cell in list(plane_neighbors.values()):
                nucleus = cell.nucleus(plane=False)
                nuc_cell.append((nucleus[0], nucleus[1], cell))
            nuc_cell.sort()
            if self.region() == "north_polar":
                result["south_0"] = nuc_cell[0][2]
                result["south_1"] = nuc_cell[1][2]
                result["south_2"] = nuc_cell[2][2]
                result["south_3"] = nuc_cell[3][2]
            else:
                result["north_0"] = nuc_cell[0][2]
                result["north_1"] = nuc_cell[1][2]
                result["north_2"] = nuc_cell[2][2]
                result["north_3"] = nuc_cell[3][2]
        elif shape == "skew_quad":
            # To avoid east-west longitude wrapping, move prime meridian
            # so that nucleus of this cell is at longitude 0.
            old_lon_0 = self.rdggs.ellipsoid.lon_0
            self.rdggs.ellipsoid.lon_0 = -self.nucleus(plane=False)[0]
            # Get lon-lat coordinates of neighbor centroids.
            nuc_cell = []
            for cell in list(plane_neighbors.values()):
                nucleus = cell.nucleus(plane=False)
                nuc_cell.append((nucleus[0], nucleus[1], cell))
            # Max latitude cell is north neighbor:
            north = max(nuc_cell, key=lambda x: x[1])
            result["north"] = north[2]
            nuc_cell.remove(north)
            # Min latitude cell is south neighbor:
            south = min(nuc_cell, key=lambda x: x[1])
            result["south"] = south[2]
            nuc_cell.remove(south)
            # Max longitude cell is east neighbor
            # (because i moved the prime meridian):
            result["east"] = max(nuc_cell, key=lambda x: x[0])[2]
            # Min longitude cell is west neighbor
            # (because i moved the prime meridian and removed cap cells):
            result["west"] = min(nuc_cell, key=lambda x: x[0])[2]
            # Return prime meridian to its original position.
            self.rdggs.ellipsoid.lon_0 = old_lon_0
        else:
            # Dart cell.
            # To avoid east-west longitude wrapping, move prime meridian
            # so that nucleus of this cell is at longitude 0.
            old_lon_0 = self.rdggs.ellipsoid.lon_0
            self.rdggs.ellipsoid.lon_0 = -self.nucleus(plane=False)[0]
            nuc_cell = []
            for cell in list(plane_neighbors.values()):
                nucleus = cell.nucleus(plane=False)
                nuc_cell.append((nucleus[0], nucleus[1], cell))
            # Sort cells by longitude. Works because moved prime meridian.
            nuc_cell.sort()
            if self.region() == "north_polar":
                result["west"] = nuc_cell[0][2]
                result["south_west"] = nuc_cell[1][2]
                result["south_east"] = nuc_cell[2][2]
                result["east"] = nuc_cell[3][2]
            else:
                result["west"] = nuc_cell[0][2]
                result["north_west"] = nuc_cell[1][2]
                result["north_east"] = nuc_cell[2][2]
                result["east"] = nuc_cell[3][2]
            # Return prime meridian to its original position.
            self.rdggs.ellipsoid.lon_0 = old_lon_0
        return result

    def random_point(self, plane=True):
        """
        Return a random point in this cell.
        If `plane` = True, then choose the point from
        the planar cell.
        Otherwise, choose the point from the ellipsoidal cell.

        EXAMPLES::

            >>> c = Cell(RHEALPixDGGS(), ['N', 0])
            >>> print(c.random_point(plane=False))  # doctest: +SKIP
            (1.4840291937583836, 0.90042819146088571)

        """
        vertices = self.vertices(plane=plane)
        u_min = min([p[0] for p in vertices])
        u_max = max([p[0] for p in vertices])
        v_min = min([p[1] for p in vertices])
        v_max = max([p[1] for p in vertices])
        if plane:
            return uniform(u_min, u_max), uniform(v_min, v_max)
        else:
            if self.ellipsoidal_shape() == "cap":
                # Need to adjust extremes.
                PI = self.ellipsoid.pi()
                u_max = PI
                if v_min > 0:
                    v_max = PI / 2
                else:
                    v_min = -PI / 2
            # Sample longitude and latitude within extremes, but reject if
            # they don't lie in the cell.
            # Rejection can happen for polar cells, because they are not
            # rectangular.
            while True:
                lam, phi = self.ellipsoid.random_point(u_min, u_max, v_min, v_max)
                if self.contains((lam, phi), plane=False):
                    # Success
                    return lam, phi

    def color(self, saturation=0.5):
        """
        Return a unique RGB color tuple for this cell.
        Inessential graphics method.
        """
        suid = self.suid
        N = self.rdggs.N_side
        hue_resolution0 = dict(
            [(v, k / 6.0) for (k, v) in enumerate(RHEALPixDGGS.cells0)]
        )
        hue = hue_resolution0[suid[0]]
        n = len(suid)
        if n > 1:
            hue += sum([suid[i] * N ** (-2 * i) for i in range(1, n)]) / 6.0
            # hue += sum([suid[i + 1]*N**(resolution - 1 - i)
            #             for i in range(resolution)])/\
            #        float(6*N**(2*resolution))
        return hsv_to_rgb(hue, saturation, 1)


class RhealPolygon(object):
    """ """

    def __init__(self, rdggs=WGS84_003, suid_list=None):
        self.rdggs = rdggs
        self.ellipsoid = rdggs.ellipsoid
        self.N_side = rdggs.N_side
        self.suid = ()  # Spatially unique identifier of self.
        self.suid_list = suid_list
