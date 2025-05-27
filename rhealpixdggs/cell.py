# from rhealpixdggs.dggs import WGS84_003

from numpy import array, base_repr, pi  # pi is just for the doctests
from scipy import integrate
from itertools import product
from random import uniform
from colorsys import hsv_to_rgb
from functools import total_ordering

# Level 0 cell IDs, which are anomalous.
CELLS0 = ["N", "O", "P", "Q", "R", "S"]


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
            suid[0] = CELLS0[suid[0]]
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
            suid[0] = CELLS0[suid[0]]
            suid = tuple(suid)
        return suid

    def __init__(
        self, rdggs=None, suid=None, level_order_index=None, post_order_index=None
    ):
        """
        Create a cell either from its suid or from its level order or
        post order index.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import UNIT_003
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
            assert (len(suid) > 0) and (
                len(suid) <= rdggs.max_resolution + 1
            ), "Need 0 < len(suid) <= %s. Got %s." % (rdggs.max_resolution + 1, suid)
            assert suid[0] in CELLS0, "suid[0] must lie in %s. Got %s." % (
                CELLS0,
                suid[0],
            )
            digits = set(range(self.N_side**2))
            for x in suid[1:]:
                assert x in digits, "Digits of suid must lie in %s" % digits

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

            >>> from rhealpixdggs.dggs import UNIT_003
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
        s[0] = CELLS0.index(s[0])
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

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
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

            >>> from rhealpixdggs.dggs import UNIT_003
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

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
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
            if suid[0] == CELLS0[-1]:
                # End of the line. No successor.
                return None
            else:
                i = CELLS0.index(suid[0])
                suid = [CELLS0[i + 1]] + [0 for j in range(resolution)]
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

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
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
            i = CELLS0.index(suid[greatest])
            if i == 0:
                # End of the line. No predecessor.
                return None
            else:
                suid = [CELLS0[i - 1]] + [M for i in range(resolution)]
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

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
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

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
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

            >>> from rhealpixdggs.dggs import UNIT_003
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
            >>> from rhealpixdggs.utils import my_round
            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> c = rdggs.cell(['P', 5, 7]) # Quad cell.
            >>> print(tuple(x.tolist() for x in my_round(c.ul_vertex(plane=True), 14)))
            (-2223901.0395045895, -555975.2598761478)
            >>> print(tuple(x.tolist() for x in my_round(c.nw_vertex(plane=True), 14)))
            (-2223901.0395045895, -555975.2598761478)

            >>> c = rdggs.cell(['S', 4])  # Cap cell.
            >>> print(tuple(x.tolist() for x in my_round(c.ul_vertex(plane=True), 14)))
            (-16679257.796284417, -8339628.898142208)
            >>> print(tuple(x.tolist() for x in my_round(c.nw_vertex(plane=True), 14)))
            (-16679257.796284417, -8339628.898142208)

            >>> c = rdggs.cell(['N', 4, 3]) # Skew quad cell.
            >>> print(tuple(x.tolist() for x in my_round(c.ul_vertex(plane=True), 14)))
            (-16679257.796284417, 10563529.937646799)
            >>> print(tuple(x.tolist() for x in my_round(c.nw_vertex(plane=True), 14)))
            (-15567307.276532121, 10563529.937646799)

            >>> c = rdggs.cell(['S', 4, 3])  # Skew quad cell.
            >>> print(tuple(x.tolist() for x in my_round(c.ul_vertex(plane=True), 14)))
            (-16679257.796284417, -9451579.417894501)
            >>> print(tuple(x.tolist() for x in my_round(c.nw_vertex(plane=True), 14)))
            (-16679257.796284417, -10563529.937646795)

            >>> c = rdggs.cell(['N', 6, 2])  # Dart cell.
            >>> print(tuple(x.tolist() for x in my_round(c.ul_vertex(plane=True), 14)))
            (-17791208.316036712, 8339628.898142208)
            >>> print(tuple(x.tolist() for x in my_round(c.nw_vertex(plane=True), 14)))
            (-16679257.796284417, 8339628.898142208)

            >>> c = rdggs.cell(['S', 6, 2])  # Dart cell.
            >>> print(tuple(x.tolist() for x in my_round(c.ul_vertex(plane=True), 14)))
            (-17791208.316036712, -11675480.457399093)
            >>> print(tuple(x.tolist() for x in my_round(c.nw_vertex(plane=True), 14)))
            (-16679257.796284417, -12787430.977151388)

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

            >>> from rhealpixdggs.utils import my_round
            >>> from rhealpixdggs.dggs import UNIT_003
            >>> rdggs = UNIT_003
            >>> c = rdggs.cell(['N'])
            >>> print(tuple(x.tolist() for x in my_round(c.nucleus(), 14)))
            (-2.35619449019234, 1.5707963267949)

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

            >>> from rhealpixdggs.utils import my_round
            >>> from rhealpixdggs.dggs import UNIT_003, WGS84_003
            >>> rdggs = UNIT_003
            >>> c = rdggs.cell(['N'])
            >>> for p in c.vertices():
            ...     print(tuple(x.tolist() for x in my_round(p, 14)))
            (-3.14159265358979, 2.35619449019234)
            (-1.5707963267949, 2.35619449019234)
            (-1.5707963267949, 0.78539816339745)
            (-3.14159265358979, 0.78539816339745)

            >>> rdggs = WGS84_003
            >>> c = rdggs.cell(['N', 0])
            >>> for p in c.vertices(plane=False):
            ...     print(tuple(x.tolist() for x in my_round(p, 14)))
            (89.99999999999993, 74.42400670199599)
            (119.99999999999999, 41.93785391016014)
            (90.0, 41.93785391016014)
            (59.99999999999999, 41.93785391016014)
            >>> for p in c.vertices(plane=False, trim_dart=True):
            ...     print(tuple(x.tolist() for x in my_round(p, 14)))
            (89.99999999999993, 74.42400670199599)
            (119.99999999999999, 41.93785391016014)
            (59.99999999999999, 41.93785391016014)

            >>> c = rdggs.cell(['S', 0])
            >>> for p in c.vertices(plane=False):
            ...     print(tuple(x.tolist() for x in my_round(p, 13)))
            (150.0, -41.9378539101601)
            (-180.0, -41.9378539101601)
            (-150.0, -41.9378539101601)
            (-179.9999999999999, -74.424006701996)
            >>> for p in c.vertices(plane=False, trim_dart=True):
            ...     print(tuple(x.tolist() for x in my_round(p, 13)))
            (150.0, -41.9378539101601)
            (-150.0, -41.9378539101601)
            (-179.9999999999999, -74.424006701996)

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

            >>> from rhealpixdggs.dggs import UNIT_003
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

            >>> from rhealpixdggs.utils import my_round
            >>> from rhealpixdggs.dggs import UNIT_003
            >>> rdggs = UNIT_003
            >>> c = rdggs.cell(['N', 6])
            >>> c.boundary(n=2, plane=True) == c.vertices(plane=True)
            True
            >>> for p in c.boundary(n=3, plane=True):
            ...     print(tuple(x.tolist() for x in my_round(p, 14)))
            (-3.14159265358979, 1.30899693899575)
            (-2.87979326579064, 1.30899693899575)
            (-2.61799387799149, 1.30899693899575)
            (-2.61799387799149, 1.0471975511966)
            (-2.61799387799149, 0.78539816339745)
            (-2.87979326579064, 0.78539816339745)
            (-3.14159265358979, 0.78539816339745)
            (-3.14159265358979, 1.0471975511966)

            >>> for p in c.boundary(n=3, plane=False):
            ...     print(tuple(x.tolist() for x in my_round(p, 14)))
            (-180.0, 74.35752898700072)
            (-157.50000000000003, 58.41366190347208)
            (-150.0, 41.8103148957786)
            (-165.00000000000003, 41.8103148957786)
            (-180.0, 41.8103148957786)
            (165.0, 41.8103148957786)
            (149.99999999999997, 41.8103148957786)
            (157.49999999999997, 58.41366190347208)

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

            >>> from rhealpixdggs.utils import my_round
            >>> from rhealpixdggs.dggs import UNIT_003
            >>> rdggs = UNIT_003
            >>> c = rdggs.cell(['N'])
            >>> for p in c.interior(n=2, plane=False, flatten=True):
            ...     print(tuple(x.tolist() for x in my_round(p, 14)))
            (90.0, 41.8103801453539)
            (-180.0, 41.8103801453539)
            (-1e-14, 41.8103801453539)
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

            >>> from rhealpixdggs.dggs import WGS84_003_RADIANS
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

            >>> from rhealpixdggs.dggs import WGS84_003_RADIANS
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

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> print(Cell(rdggs, ['P', 2]).region())
            equatorial
            >>> print(Cell(rdggs, ['N', 2]).region())
            north_polar

        """
        if self.suid[0] == CELLS0[0]:
            return "north_polar"
        elif self.suid[0] == CELLS0[5]:
            return "south_polar"
        else:
            return "equatorial"

    def ellipsoidal_shape(self):
        """
        Return the shape of this cell ('quad', 'cap', 'dart', or
        'skew_quad') when viewed on the ellipsoid.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> print(Cell(rdggs, ['P', 2]).ellipsoidal_shape())
            quad
            >>> print(Cell(rdggs, ['N', 2]).ellipsoidal_shape())
            dart

        """
        suid = self.suid
        if suid[0] in CELLS0[1:5]:
            return "quad"
        # Cap check.
        if suid == tuple(CELLS0[0]) or suid == tuple(CELLS0[5]):
            return "cap"
        N = self.N_side
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

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
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

        - `x` - A letter from CELLS0 or one of the integers
          0, 1, ..., N**2 - 1.
        - `quarter_turns` - 0, 1, 2, or 3.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> c = Cell(RHEALPixDGGS(), ['P', 2])
            >>> print([c.rotate_entry(0, t) for t in range(4)])
            [0, 2, 8, 6]

        NOTES:

        Operates on letters from CELLS0 too.
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
        for c in CELLS0:
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

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
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

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
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
            if (
                (self0 == CELLS0[5] and neighbor0 == an[self0]["left"])
                or (self0 == an[CELLS0[5]]["right"] and neighbor0 == CELLS0[5])
                or (self0 == CELLS0[0] and neighbor0 == an[self0]["right"])
                or (self0 == an[CELLS0[0]]["left"] and neighbor0 == CELLS0[0])
            ):
                neighbor = neighbor.rotate(1)
            elif (
                (self0 == CELLS0[5] and neighbor0 == an[self0]["down"])
                or (self0 == an[CELLS0[5]]["down"] and neighbor0 == CELLS0[5])
                or (self0 == CELLS0[0] and neighbor0 == an[self0]["up"])
                or (self0 == an[CELLS0[0]]["up"] and neighbor0 == CELLS0[0])
            ):
                neighbor = neighbor.rotate(2)
            elif (
                (self0 == CELLS0[5] and neighbor0 == an[self0]["right"])
                or (self0 == an[CELLS0[5]]["left"] and neighbor0 == CELLS0[5])
                or (self0 == CELLS0[0] and neighbor0 == an[self0]["left"])
                or (self0 == an[CELLS0[0]]["right"] and neighbor0 == CELLS0[0])
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

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
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

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
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
        hue_resolution0 = dict([(v, k / 6.0) for (k, v) in enumerate(CELLS0)])
        hue = hue_resolution0[suid[0]]
        n = len(suid)
        if n > 1:
            hue += sum([suid[i] * N ** (-2 * i) for i in range(1, n)]) / 6.0
            # hue += sum([suid[i + 1]*N**(resolution - 1 - i)
            #             for i in range(resolution)])/\
            #        float(6*N**(2*resolution))
        return hsv_to_rgb(hue, saturation, 1)
