# from rhealpixdggs.dggs import WGS84_003

import enum
from collections.abc import Iterator
from colorsys import hsv_to_rgb
from functools import cache, cached_property, total_ordering
from itertools import product
from math import fsum, hypot
from random import uniform
from typing import TYPE_CHECKING, ClassVar, Literal, cast, overload

if TYPE_CHECKING:
    from rhealpixdggs.dggs import RHEALPixDGGS
    from rhealpixdggs.zoneset import ZoneSet

# pi is doctest-only: the doctests use it from the module globals.
import numpy as np
import pyproj
from numpy import base_repr, pi  # noqa: F401
from scipy import integrate, optimize

from rhealpixdggs.utils import (
    FloatArray,
    _normalise_antimeridian_rings,
    wrap_longitude,
)

# Level 0 cell IDs, which are anomalous.
CELLS0 = ["N", "O", "P", "Q", "R", "S"]

# Points per axis of the Gauss-Legendre product rule that integrates the
# centroid of a dart or skew quad cell over its planar square (issue #120).
# The integrand is smooth within one polar triangle, so the rule converges
# exponentially: 20 points agree with 40 to about 1e-14 degrees, the
# projection's own floating-point floor.
_CENTROID_QUADRATURE_ORDER = 20


@cache
def _gauss_legendre_unit(n: int) -> tuple[np.ndarray, np.ndarray]:
    """
    The `n`-point Gauss-Legendre nodes and weights on [0, 1]; the weights
    sum to 1, so a weighted sum of function values is the mean over [0, 1].
    """
    nodes, weights = np.polynomial.legendre.leggauss(n)
    return (nodes + 1) / 2, weights / 2


@cache
def _geod(a: float, f: float) -> pyproj.Geod:
    """
    The geodesic calculator for the ellipsoid of semi-major axis `a` and
    flattening `f` (a sphere when `f` is 0), shared by every cell of every
    grid on that ellipsoid.
    """
    return pyproj.Geod(a=a, f=f)


class RelativePosition(enum.Enum):
    """
    The relative position of two cells projected onto one direction, as
    ``Cell.relative_position`` reports it: the enumeration of OGC Topic 21
    v2.0 Table 54, whose names follow OGC 16-071r3 (OWL-Time). Each cell
    projects to a closed interval from `begin` to `end` along the
    direction, and the thirteen primary members are the interval
    relations of Allen (1983), here with `self` the cell asked and
    `other` the argument:

    =============  =================================================
    member         condition
    =============  =================================================
    BEFORE         self.end < other.begin
    AFTER          self.begin > other.end
    MEETS          self.end == other.begin
    MET_BY         self.begin == other.end
    OVERLAPS       self.begin < other.begin < self.end < other.end
    OVERLAPPED_BY  other.begin < self.begin < other.end < self.end
    STARTS         self.begin == other.begin and self.end < other.end
    STARTED_BY     self.begin == other.begin and self.end > other.end
    DURING         other.begin < self.begin and self.end < other.end
    CONTAINS       self.begin < other.begin and other.end < self.end
    FINISHES       self.end == other.end and self.begin > other.begin
    FINISHED_BY    self.end == other.end and self.begin < other.begin
    EQUALS         self.begin == other.begin and self.end == other.end
    =============  =================================================

    (Table 54 prints FINISHED_BY with its `begin` inequality reversed and
    OVERLAPPED_BY without its `begin` clause; these are the standard
    relations.) The table's two further members are groupings: IN is
    STARTS, DURING or FINISHES, and DISJOINT is BEFORE or AFTER.
    ``relative_position`` never returns them; ``relate_position`` accepts
    them, as do the ``is_in`` and ``is_disjoint`` properties.

    Read along the planar x axis the words describe east-west
    arrangement in the equatorial belt: BEFORE is "entirely west of",
    MEETS "shares an edge on the east side", DURING "lies within the
    other's east-west span". Along planar y they describe south to
    north.
    """

    BEFORE = "Before"
    AFTER = "After"
    MEETS = "Meets"
    MET_BY = "MetBy"
    OVERLAPS = "Overlaps"
    OVERLAPPED_BY = "OverlappedBy"
    STARTS = "Starts"
    STARTED_BY = "StartedBy"
    DURING = "During"
    CONTAINS = "Contains"
    FINISHES = "Finishes"
    FINISHED_BY = "FinishedBy"
    EQUALS = "Equals"
    IN = "In"
    DISJOINT = "Disjoint"

    @property
    def is_in(self) -> bool:
        """True for STARTS, DURING and FINISHES, the members grouped as IN."""
        return self in (
            RelativePosition.STARTS,
            RelativePosition.DURING,
            RelativePosition.FINISHES,
        )

    @property
    def is_disjoint(self) -> bool:
        """True for BEFORE and AFTER, the members grouped as DISJOINT."""
        return self in (RelativePosition.BEFORE, RelativePosition.AFTER)


@total_ordering
class Cell:
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
    def suid_from_index(
        rdggs: "RHEALPixDGGS", index: int, order: str = "resolution"
    ) -> tuple[str | int, ...]:
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
            digits: list[int] = []
            p = index

            def num(k: int) -> int:
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
                digits.append(q)
                p = r
                if p == n - 1:
                    break
        else:
            b = rdggs.N_side**2

            # Compute suid from level order index.
            def ind(k: int) -> int:
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
            digits_str = base_repr(remainder, b)
            # If necessary, prepend with zeros to get a length L + 1 string.
            digits_str = "0" * (L + 1 - len(digits_str)) + digits_str
            digits = [int(s) for s in digits_str]
        # Replace the first digit with the appropriate letter.
        return (CELLS0[digits[0]], *digits[1:])

    def __init__(
        self,
        rdggs: "RHEALPixDGGS",
        suid: list[str | int] | tuple[str | int, ...] | None = None,
        level_order_index: int | None = None,
        post_order_index: int | None = None,
    ) -> None:
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
        # Spatially unique identifier of self: a face letter followed by
        # int digits, or empty for the (rarely used) empty cell.
        self.suid: tuple[str | int, ...] = ()
        # Level of self in grid hierarchy; None only for the empty cell.
        self.resolution: int | None = None
        if suid is not None:
            # A little error checking.
            if not isinstance(suid, (list, tuple)):
                raise TypeError(f"Cell suid must be a list or tuple. Got {suid}.")
            if not (0 < len(suid) <= rdggs.max_resolution + 1):
                raise ValueError(
                    f"Need 0 < len(suid) <= {rdggs.max_resolution + 1}. Got {suid}."
                )
            if suid[0] not in CELLS0:
                raise ValueError(f"suid[0] must lie in {CELLS0}. Got {suid[0]}.")
            digits = set(range(self.N_side**2))
            for x in suid[1:]:
                if x not in digits:
                    raise ValueError(f"Digits of suid must lie in {digits}. Got {x}.")

            self.suid = tuple([suid[0]] + [int(n) for n in suid[1:]])
        elif level_order_index is not None:
            self.suid = Cell.suid_from_index(
                self.rdggs, level_order_index, order="resolution"
            )
        elif post_order_index is not None:
            self.suid = Cell.suid_from_index(self.rdggs, post_order_index, order="post")
        self.resolution = len(self.suid) - 1

    def __bool__(self) -> bool:
        return bool(self.suid)

    def __str__(self) -> str:
        """
        The cell's index string, ``rdggs.format_index(suid)``, for a grid
        with index strings (``N_side`` 2 or 3); for any other grid a display
        form such as ``'(N, 10, 15)'``, which is not an identifier (issue
        #146).
        """
        if self.rdggs.has_index_strings:
            return self.rdggs.format_index(self.suid)
        return "(" + str(self.suid[0]) + str(self.suid)[4:]

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Cell)
            and (self.rdggs == other.rdggs)
            and (self.suid == other.suid)
        )

    def __ne__(self, other: object) -> bool:
        return not self.__eq__(other)

    def __hash__(self) -> int:
        # Consistent with __eq__ (equal cells have equal suids), and
        # makes cells usable as dictionary keys and set members. A
        # cell's suid never changes after construction.
        return hash(self.suid)

    def __le__(self, other: "Cell") -> bool:
        """
        The (strictly) less-than relation on cells.
        Derived from the post order traversal of the tree T of all cells
        defined in the `index()` docstring.
        Return True if (`self.suid < other.suid` and
        `self.suid` is not a prefix of `other.suid`) or
        `self` is a subcell of `other`.
        Here < is the lexicographic order on the suid tuple itself (face
        letter, then digits compared as the integers they are -- not on
        any string rendering of it, which for `N_side >= 4` would compare
        multi-character digits like 10-15 character-by-character rather
        than numerically).
        Returns False otherwise.
        """
        s = self.suid
        t = other.suid
        t_starts_with_s = t[: len(s)] == s
        s_starts_with_t = s[: len(t)] == t
        return (s <= t and not t_starts_with_s) or s_starts_with_t

    def index(self, order: str = "resolution") -> int | None:
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
        if not self.suid:
            return None
        L = self.resolution
        assert L is not None  # non-empty cells always have a resolution
        # suid[0] is the face letter; the rest are int digits already.
        s = [CELLS0.index(str(self.suid[0]))] + [int(d) for d in self.suid[1:]]
        if order == "post":

            def num(k: int) -> int:
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

    def suid_rowcol(
        self,
    ) -> tuple[tuple[str | int, ...], tuple[str | int, ...]]:
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
        suid_row: list[str | int] = [self.suid[0]]
        suid_col: list[str | int] = [self.suid[0]]
        for n in self.suid[1:]:
            row, col = cast("tuple[int, int]", self.rdggs.child_order[int(n)])
            suid_row.append(row)
            suid_col.append(col)
        return tuple(suid_row), tuple(suid_col)

    def ring(self) -> int:
        """
        Return the isolatitude ring of this cell's nucleus. Rings are
        numbered 0 at the north pole through ``n + 2 * q - 1`` at the south
        pole, where ``n = N_side ** resolution`` cells span a base cell's
        side and ``q = ceil(n / 2)`` rings fill each polar cap. Every cell on
        a ring, in any base cell, has the same nucleus latitude. Ring numbers
        are only comparable within one resolution. See :doc:`isolatitude`.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> suids = [['N', 4], ['N', 0], ['P', 4], ['S', 8]]
            >>> [Cell(WGS84_003, s).ring() for s in suids]
            [0, 1, 3, 5]

        """
        N = self.rdggs.N_side
        n = N ** (len(self.suid) - 1)
        q = -(-n // 2)
        row = col = 0
        for digit in self.suid[1:]:
            r, c = cast("tuple[int, int]", self.rdggs.child_order[int(digit)])
            row = row * N + r
            col = col * N + c
        if self.suid[0] not in (CELLS0[0], CELLS0[5]):
            return q + row
        m = max(abs(2 * row - (n - 1)), abs(2 * col - (n - 1))) // 2
        if self.suid[0] == CELLS0[0]:
            return m
        return q + n + (q - 1 - m)

    @overload
    def width(self, plane: Literal[True] = ...) -> float: ...

    @overload
    def width(self, plane: Literal[False]) -> None: ...

    def width(self, plane: bool = True) -> float | None:
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
        assert self.resolution is not None
        return self.rdggs.cell_width(self.resolution, plane=plane)

    def area(self, plane: bool = True) -> float:
        """
        Return the area of this cell.
        """
        assert self.resolution is not None
        return self.rdggs.cell_area(self.resolution, plane=plane)

    def successor(self, resolution: int | None = None) -> "Cell | None":
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
        assert self.resolution is not None
        suid = list(self.suid)
        if resolution is None:
            resolution = self.resolution
        if resolution < self.resolution:
            # Truncate suid at resolution resolution and return its successor.
            return Cell(self.rdggs, suid[: resolution + 1]).successor()
        elif resolution > self.resolution:
            # Find the resolution self.resolution successor of suid
            # and pad it with zeros.
            base = self.successor()
            assert base is not None
            suid = list(base.suid) + [0 for i in range(resolution - self.resolution)]
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
                i = CELLS0.index(str(suid[0]))
                suid = [CELLS0[i + 1], *(0 for j in range(resolution))]
        else:
            # suid[greatest] is a number.
            suid = (
                suid[0:greatest]
                + [int(suid[greatest]) + 1]
                + [0 for j in range(resolution - greatest)]
            )
        return Cell(self.rdggs, suid)

    def predecessor(self, resolution: int | None = None) -> "Cell | None":
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
        assert self.resolution is not None
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
            i = CELLS0.index(str(suid[greatest]))
            if i == 0:
                # End of the line. No predecessor.
                return None
            else:
                suid = [CELLS0[i - 1], *(M for i in range(resolution))]
        else:
            # nome[greatest] is a number > 0.
            suid = (
                suid[0:greatest]
                + [int(suid[greatest]) - 1]
                + [M for i in range(resolution - greatest)]
            )
        return Cell(self.rdggs, suid)

    def subcell(self, other: "Cell") -> bool:
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
        # Compare the suid tuples directly (not a string rendering of
        # them): for N_side >= 4, some digits are multi-character (e.g.
        # 10-15), and comparing stringified, comma-joined suids with
        # startswith() can misfire, e.g. treating ('N', 15) as a subcell
        # of ('N', 1), since "N,15" starts with "N,1" as a string even
        # though these are same-resolution siblings, not ancestor and
        # descendant.
        return self.suid[: len(other.suid)] == other.suid

    def subcells(self, resolution: int | None = None) -> Iterator["Cell"]:
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
        assert L is not None  # the empty cell has no subcells
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

    def ul_vertex(self, plane: bool = True) -> tuple[float, float]:
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
        x0, y0 = self.rdggs.ul_vertex[str(self.suid[0])]
        resolution = self.resolution
        assert resolution is not None  # the empty cell has no vertices

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

    def nw_vertex(self, plane: bool = True) -> tuple[float, float]:
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
        shape = self.ellipsoidal_shape
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
            assert triangle is not None
            if region == "north_polar":
                ns = rdggs.north_square
                i = (triangle - ns) % 4
                result = v[-i]
            else:
                ss = rdggs.south_square
                i = (triangle - ss) % 4
                result = v[i]
        else:
            # shape == 'dart': find the polewards vertex in the plane.
            # Parallels in a polar square are concentric squares about its
            # centre, the pole, so the vertex nearest the pole is the one
            # with the smallest Chebyshev distance to that centre. It is
            # unique: only the dart's two base vertices, farthest from the
            # pole, are equidistant.
            cx, cy = self.rdggs.cell([self.suid[0]]).nucleus(plane=True)
            i = min(range(4), key=lambda j: max(abs(v[j][0] - cx), abs(v[j][1] - cy)))
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

    def nucleus(self, plane: bool = True) -> tuple[float, float]:
        """
        Return the nucleus of this cell: the centre of its planar square,
        or with `plane` = False the projection of that point onto the
        ellipsoid. The nucleus is the grid's indexing point: it lies
        inside the cell by construction, ``cell_from_point`` maps it back
        to the cell, and the nuclei of a resolution lie on rings of
        constant latitude, one per ``ring()`` number (see
        :doc:`isolatitude`). It is not the cell's published position:
        while the nucleus of a planar cell is its centroid, the nucleus of
        an ellipsoidal cell is not, and the direct (representative)
        position of a cell in the OGC Topic 21 sense is `centroid()`.

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

    def vertices(
        self, plane: bool = True, trim_dart: bool = False
    ) -> list[tuple[float, float]]:
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
        A vertex on the antimeridian takes the sign that keeps the cell's
        longitude span under half a turn, as ``boundary()`` describes.

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
            result = self._normalise_ring(result)
            if trim_dart and self.ellipsoidal_shape == "dart":
                # Remove non-vertex point.
                if self.region() == "north_polar":
                    result.pop(2)
                else:
                    result.pop(1)
        return result

    def xy_range(self) -> tuple[tuple[float, float], tuple[float, float]]:
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

    def boundary(
        self, n: int = 2, plane: bool = True, interior: bool = False
    ) -> list[tuple[float, float]]:
        """
        Return a list of `4*n - 4` boundary points of this cell,
        `n` on each edge, where `n` >= 2.
        List the points in clockwise order starting from the cell's upper left
        corner if `plane` = True, or from the cell's northwest corner
        if `plane` = False.

        If `n` = 2, then the output is the same as vertices().
        If `interior` = True, then push the boundary points slighly into the
        interior of the cell, which is convenient for some graphics methods.
        When `plane` = False, a point on the antimeridian takes the sign
        that keeps the ring's longitude span under half a turn: +180 for a
        cell just west of the antimeridian, -180 for one just east. Cells
        that straddle it, and cap cells, keep both signs.

        When `plane` = False, the cost scales with `n` because each point
        requires an inverse projection call, except on quad cells. A quad
        cell lies entirely in the equatorial region, where the inverse
        projection is separable: longitude depends only on `x` and latitude
        only on `y`. Its east and west edges are meridians sharing the same
        `n - 2` interior latitudes, and its north and south edges are
        parallels sharing the same `n - 2` interior longitudes. Only the four
        corners and the interior points of the west and north edges are
        projected (`2*n` calls instead of `4*n - 4`); the east and south
        edges reuse those values. Every coordinate returned is one the
        projection computed, and adjacent quad cells get bit-identical
        shared points, except that an edge on the antimeridian is +180 in
        the western cell's ring and -180 in the eastern's (see
        ``boundary()``).

        For quad and cap cells with `n` = 2 and `interior` = False the result
        is ``vertices(plane=False)``. Cap cells with `n` > 2 take the general
        per-point path: a cap's boundary is a single parallel and could be
        computed directly too, but there are only two cap cells per
        resolution, so it isn't worth a third code path.

        EXAMPLES::

            >>> from rhealpixdggs.utils import my_round
            >>> from rhealpixdggs.dggs import UNIT_003
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
        n = max(n, 2)
        # Quad and cap cells have straight or rotationally-symmetric edges on
        # the ellipsoid, so at n=2 extra boundary points would add no accuracy.
        # Fall back to vertices() and avoid the per-point projection cost
        # entirely.
        if (
            not plane
            and n == 2
            and not interior
            and self.ellipsoidal_shape in ("quad", "cap")
        ):
            return self.vertices(plane=False)
        ul = self.ul_vertex(plane=True)
        w = self.width(plane=True)
        if interior:
            eps = w / 10000  # A smidgen.
        else:
            eps = 0
        if not plane and self.ellipsoidal_shape == "quad":
            return self._quad_boundary(n, eps)
        delta = (w - 2 * eps) / (n - 1)
        x, y = float(ul[0]) + eps, float(ul[1]) - eps
        result = [(x, y)]
        for dx, dy in [(1, 0), (0, -1), (-1, 0), (0, 1)]:
            for j in range(1, n):
                result.append((x + j * delta * dx, y + j * delta * dy))
            x, y = result[-1]
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
            xs = np.array([p[0] for p in result])
            ys = np.array([p[1] for p in result])
            lons, lats = self.rdggs.rhealpix(xs, ys, inverse=True, region=self.region())
            return self._normalise_ring(list(zip(lons, lats)))
        return result

    def _quad_boundary(self, n: int, eps: float) -> list[tuple[float, float]]:
        """
        Return ``boundary(n, plane=False)`` for this quad cell, its planar
        square shrunk inward by `eps`, projecting only the north and west
        edges in one array call. See ``boundary()``.
        """
        ul = self.ul_vertex(plane=True)
        w = self.width(plane=True)
        x_west, y_north = ul[0] + eps, ul[1] - eps
        x_east, y_south = ul[0] + w - eps, ul[1] - w + eps
        delta = (w - 2 * eps) / (n - 1)
        region = self.region()

        # One batch: the north edge west to east (n points, ending exactly at
        # x_east), then the west edge below the north-west corner (n - 1
        # points, ending exactly at y_south).
        xs = np.concatenate(
            [x_west + delta * np.arange(n - 1), [x_east], np.full(n - 1, x_west)]
        )
        ys = np.concatenate(
            [np.full(n, y_north), y_north - delta * np.arange(1, n - 1), [y_south]]
        )
        all_lons, all_lats = self.rdggs.rhealpix(xs, ys, inverse=True, region=region)
        lons = list(all_lons[1 : n - 1])
        lats = list(all_lats[n:-1])
        nw = (all_lons[0], all_lats[0])
        ne = (all_lons[n - 1], all_lats[n - 1])
        sw = (all_lons[-1], all_lats[-1])
        # Longitude depends only on x and latitude only on y in the
        # equatorial region, so the south-east corner needs no projection.
        se = (ne[0], sw[1])
        return self._normalise_ring(
            [nw]
            + [(lon, nw[1]) for lon in lons]
            + [ne]
            + [(ne[0], lat) for lat in lats]
            + [se]
            + [(lon, se[1]) for lon in reversed(lons)]
            + [sw]
            + [(sw[0], lat) for lat in reversed(lats)]
        )

    def _normalise_ring(
        self, ring: list[tuple[float, float]]
    ) -> list[tuple[float, float]]:
        """
        The ring with any point on the antimeridian given the sign that
        keeps the ring's longitude span under half a turn. See
        ``utils._normalise_antimeridian_rings``.
        """
        lons = _normalise_antimeridian_rings(
            np.array([p[0] for p in ring]), radians=self.rdggs.ellipsoid.radians
        )
        return [(lon, p[1]) for lon, p in zip(lons, ring)]

    def interior(
        self, n: int = 2, plane: bool = True, flatten: bool = False
    ) -> list[tuple[float, float]] | list[list[tuple[float, float]]]:
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
        if plane:
            if flatten:
                return [
                    (ul[0] + eps + delta * j, ul[1] - eps - delta * i)
                    for j in range(n)
                    for i in range(n)
                ]
            return [
                [(ul[0] + eps + delta * j, ul[1] - eps - delta * i) for j in range(n)]
                for i in range(n)
            ]
        # Project the whole grid in one call; row i, column j is (lons[i, j],
        # lats[i, j]). The flattened order is column-major, as it always was.
        xs, ys = np.meshgrid(
            ul[0] + eps + delta * np.arange(n), ul[1] - eps - delta * np.arange(n)
        )
        lons, lats = self.rdggs.rhealpix(xs.ravel(), ys.ravel(), inverse=True)
        lons, lats = lons.reshape(n, n), lats.reshape(n, n)
        if flatten:
            return [(lons[i, j], lats[i, j]) for j in range(n) for i in range(n)]
        return [[(lons[i, j], lats[i, j]) for j in range(n)] for i in range(n)]

    def contains(self, p: tuple[float, float], plane: bool = True) -> bool:
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
        assert self.resolution is not None
        return self.rdggs.cell_from_point(self.resolution, p, plane=plane) == self

    def intersects_meridian(self, lam: float) -> bool:
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
        if self.ellipsoidal_shape == "cap":
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
            return lon_min <= lam <= lon_max

    def intersects_parallel(self, phi: float) -> bool:
        """
        Return True if this cell's boundary intersects the parallel of latitude
        `phi`, and return False otherwise.
        """
        # Cell's boundary intersects parallel iff its extreme latitudes lie on
        # opposite sides of parallel.
        vertices = self.vertices(plane=False)
        lat_min = min([v[1] for v in vertices])
        lat_max = max([v[1] for v in vertices])
        if self.ellipsoidal_shape == "cap":
            if self.region() == "north_polar":
                return phi >= lat_min
            else:
                return phi <= lat_max
        else:
            return lat_min <= phi and lat_max >= phi

    def _nests(self, other: "Cell") -> bool:
        """
        Return True if one of `self` and `other` is the other or an
        ancestor of it, i.e. one suid is a prefix of the other's. Assumes
        both cells belong to the same grid and neither is empty; the
        public predicates check that first.
        """
        n = min(len(self.suid), len(other.suid))
        return self.suid[:n] == other.suid[:n]

    def region(self) -> str:
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

    @cached_property
    def ellipsoidal_shape(self) -> str:
        """
        Return the shape of this cell ('quad', 'cap', 'dart', or
        'skew_quad') when viewed on the ellipsoid.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> print(Cell(rdggs, ['P', 2]).ellipsoidal_shape)
            quad
            >>> print(Cell(rdggs, ['N', 2]).ellipsoidal_shape)
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
        S = {i * (N + 1) for i in range(N)}
        for n in suid[1:]:
            if n not in S:
                dart = False
                break
        if dart:
            return "dart"
        # Dark check 2.
        dart = True
        S = {(i + 1) * (N - 1) for i in range(N)}
        for n in suid[1:]:
            if n not in S:
                dart = False
                break
        if dart:
            return "dart"
        # Must be a skew quad then.
        return "skew_quad"

    def centroid(self, plane: bool = True) -> tuple[float, float]:
        """
        Return the centroid of this planar or ellipsoidal cell: the
        cell's direct position, in the sense of OGC Topic 21 v2.0
        requirement 27 (the DirectPosition of an equal-area cell is its
        centroid, the geodesic centre of surface area). In the plane it is
        the centre of the square, the nucleus. On the ellipsoid it is the
        area-weighted mean of longitude and latitude over the cell,
        integrated in the equal-area plane, where equal planar area is
        equal surface area: it lies inside the cell but, except for cap
        cells and on the equator, not at the nucleus. `nucleus()` remains
        the indexing point.

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
        shape = self.ellipsoidal_shape
        if shape == "cap":
            return nucleus
        planar_vertices = self.vertices(plane=True)
        x1 = min([v[0] for v in planar_vertices])
        x2 = max([v[0] for v in planar_vertices])
        y1 = min([v[1] for v in planar_vertices])
        y2 = max([v[1] for v in planar_vertices])

        if shape == "quad":
            # A quad cell is symmetric about its nucleus meridian, and its
            # meridians are equally spaced in planar x, so the mean
            # longitude is the nucleus longitude. Latitude is independent
            # of planar x on a quad cell, so the area-weighted mean
            # latitude reduces to a single integral over y. Note it is
            # NOT the midpoint of the two edge latitudes: latitude is a
            # nonlinear function of planar y, so the mean sits closer to
            # the equator than the midpoint (by up to ~0.6 degrees for
            # resolution 1 cells).
            lam_bar = float(nucleus[0])
            # Integrate along the nucleus meridian's planar x, which is
            # safely interior to the cell (any x in the cell would do,
            # since latitude doesn't depend on it). Fixed-order
            # Gauss-Legendre quadrature is effectively exact here -- the
            # integrand is smooth, and n=20 agrees with adaptive
            # quadrature to machine precision -- while avoiding adaptive
            # quadrature's error estimation, whose requested tolerances
            # collide with the projection stack's floating-point noise
            # floor and trigger spurious IntegrationWarnings.
            x_mid = self.nucleus(plane=True)[0]

            def phi_of_y(ys: np.ndarray) -> np.ndarray:
                return self.rdggs.rhealpix(np.full_like(ys, x_mid), ys, inverse=True)[1]

            phi_bar = float(
                (1 / (y2 - y1)) * integrate.fixed_quad(phi_of_y, y1, y2, n=20)[0]
            )
            return lam_bar, phi_bar
        # Dart or skew quad: the mean latitude (and, for a skew quad, the
        # mean longitude) is an area-weighted integral over the planar
        # square, evaluated by a fixed Gauss-Legendre product rule on one
        # batch of projected points.
        xs, ys, weights = self._centroid_quadrature(x1, x2, y1, y2)
        lons, lats = self.rdggs.rhealpix(xs, ys, inverse=True, region=self.region())
        # fsum is exactly rounded, so the result is the same on every
        # platform regardless of summation order.
        phi_bar = fsum(weights * lats)
        if shape == "dart":
            # A dart is symmetric about its nucleus meridian.
            return float(nucleus[0]), phi_bar
        return fsum(weights * lons), phi_bar

    def _centroid_quadrature(
        self, x1: float, x2: float, y1: float, y2: float
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Planar sample points and weights (summing to 1) whose weighted sum
        of a function's values is its mean over this dart or skew quad
        cell's planar square [x1, x2] x [y1, y2].

        A skew quad lies within one polar triangle, where the inverse
        projection is smooth, so a plain product rule over the square
        serves. A dart straddles a diagonal of its polar square, along
        which the projection has a kink; its square is split along that
        diagonal (which passes through the dart's nucleus and so runs
        corner to corner of the square) into two triangles, each mapped
        from the unit square by the Duffy transform (u, v) -> (u, u*v),
        whose Jacobian u folds into the weights.
        """
        t, w = _gauss_legendre_unit(_CENTROID_QUADRATURE_ORDER)
        u, v = np.meshgrid(t, t, indexing="ij")
        if self.ellipsoidal_shape == "skew_quad":
            s, r = u.ravel(), v.ravel()
            weights = np.outer(w, w).ravel()
        else:
            # Which diagonal: the nucleus's offset from the polar square's
            # centre has equal-sign components on the rising (slope +1)
            # diagonal and opposite-sign components on the falling one.
            cx, cy = self.rdggs.cell([self.suid[0]]).nucleus(plane=True)
            nx, ny = self.nucleus(plane=True)
            rising = (nx - cx) * (ny - cy) > 0
            # In unit-square coordinates (s, r) the rising diagonal is
            # r = s; the triangle below it is the Duffy image (u, u*v) and
            # the one above is its transpose (u*v, u).
            a, b = u.ravel(), (u * v).ravel()
            s = np.concatenate([a, b])
            r = np.concatenate([b, a])
            if not rising:
                # The falling diagonal r = 1 - s is the mirror image.
                s = 1 - s
            triangle_weights = (np.outer(w, w) * t[:, None]).ravel()
            weights = np.concatenate([triangle_weights, triangle_weights])
        xs = x1 + (x2 - x1) * s
        ys = y1 + (y2 - y1) * r
        return xs, ys, weights

    def rotate_entry(self, x: str | int, quarter_turns: int) -> str | int:
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
        f: dict[str | int, str | int] = {}
        for i in range(N):
            for j in range(N):
                n = cast(int, A[(i, j)])
                f[n] = cast(int, A[(j, N - 1 - i)])
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

    def rotate(self, quarter_turns: int) -> "Cell":
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

    def neighbor(self, direction: str, plane: bool = True) -> "Cell | None":
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
            down_border = {(N - 1) * N + i for i in range(N)}
            left_border = {i * N for i in range(N)}
            right_border = {(i + 1) * N - 1 for i in range(N)}
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

    def neighbors(self, plane: bool = True) -> "dict[str, Cell]":
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
        plane_neighbors: dict[str, Cell] = {}
        for d in ["left", "right", "down", "up"]:
            neighbor = self.neighbor(d, plane=True)
            assert neighbor is not None
            plane_neighbors[d] = neighbor
        if plane:
            return plane_neighbors
        # Ellipsoid case.
        result = {}
        shape = self.ellipsoidal_shape
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
            # Get lon-lat coordinates of neighbor centroids, paired with
            # longitude relative to this cell's nucleus (see the helper's
            # docstring for why "relative" instead of raw longitude).
            nuc_cell = self._neighbor_nuclei_by_relative_longitude(plane_neighbors)
            # Max latitude cell is north neighbor:
            north = max(nuc_cell, key=lambda x: x[1])
            result["north"] = north[2]
            nuc_cell.remove(north)
            # Min latitude cell is south neighbor:
            south = min(nuc_cell, key=lambda x: x[1])
            result["south"] = south[2]
            nuc_cell.remove(south)
            # Max relative longitude cell is east neighbor:
            result["east"] = max(nuc_cell, key=lambda x: x[0])[2]
            # Min relative longitude cell is west neighbor:
            result["west"] = min(nuc_cell, key=lambda x: x[0])[2]
        else:
            # Dart cell.
            nuc_cell = self._neighbor_nuclei_by_relative_longitude(plane_neighbors)
            # Sort cells by longitude relative to this cell's nucleus.
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
        return result

    def _neighbor_nuclei_by_relative_longitude(
        self, plane_neighbors: "dict[str, Cell]"
    ) -> "list[tuple[float, float, Cell]]":
        """
        Return a list of `(relative_longitude, latitude, cell)` triples,
        one per cell in `plane_neighbors`, where `relative_longitude` is
        each neighbor's nucleus longitude minus this cell's own nucleus
        longitude, wrapped into `(-pi, pi]` (or `(-180, 180]` in degrees
        mode). Used by `neighbors()` for the dart and skew_quad cases to
        compare/sort neighbors by longitude without an east-west
        antimeridian wrap-around artefact -- e.g. a neighbor at -179
        degrees is 2 degrees *east* of one at 179 degrees, not far to the
        west, and comparing raw longitudes would get that backwards.
        Computed relative to this cell's own nucleus rather than any
        global reference, so the result is independent of the ellipsoid's
        `lon_0` and touches no shared state.
        """
        self_lon = self.nucleus(plane=False)[0]
        radians = self.rdggs.ellipsoid.radians
        nuc_cell = []
        for cell in plane_neighbors.values():
            nucleus = cell.nucleus(plane=False)
            rel_lon = wrap_longitude(nucleus[0] - self_lon, radians=radians)
            nuc_cell.append((rel_lon, nucleus[1], cell))
        return nuc_cell

    # Diagonal (corner-touching only) directions, each a (row, column)
    # direction pair, keyed to match neighbor()'s own plane direction names.
    _DIAGONAL_DIRECTIONS: ClassVar[dict[str, tuple[str, str]]] = {
        "up_left": ("up", "left"),
        "up_right": ("up", "right"),
        "down_left": ("down", "left"),
        "down_right": ("down", "right"),
    }

    def diagonal_neighbor(self, direction: str) -> "Cell | None":
        """
        Return this cell's diagonal (corner-touching only, not sharing an
        edge) planar neighbor in the given `direction`, one of
        'up_left', 'up_right', 'down_left', 'down_right'.

        Return `None` if this cell sits at a genuine cube corner in that
        direction (i.e. all of its ancestors, down to and including its
        resolution 0 cell, are positioned exactly in that corner). Cube
        corners are 3-valent -- exactly 3 cells meet there, not 4 -- so
        there is no diagonal 4th cell distinct from the two already
        reachable via `neighbor()`.

        Unlike `neighbor()`, this has no `plane` parameter: cell adjacency
        is a property of the grid's topology, independent of whether it's
        expressed in planar or ellipsoidal coordinates, so this is always
        computed in the (topologically primary) planar grid.

        Note this is not simply `self.neighbor(d1, plane=True).neighbor(d2,
        plane=True)` (nor the same composed in the other order): composing
        two independently rotation-corrected edge steps does not reliably
        commute near a face boundary, since crossing one face in the first
        step can change what "the other direction" even means for the
        second step. This computes both directions' digit-carrying in a
        single pass, matching how `neighbor()` itself only ever applies one
        rotation correction, and applies that same correction once, for
        whichever single direction (if any) actually crosses a face.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> c = Cell(RHEALPixDGGS(), ['N', 4, 0])
            >>> print(c.diagonal_neighbor('up_left'))
            N08

        """
        row_dir, col_dir = Cell._DIAGONAL_DIRECTIONS[direction]
        an = self.rdggs.atomic_neighbors
        child_order = self.rdggs.child_order
        N = self.N_side
        self_suid = list(self.suid)
        neighbor_suid = list(self_suid)
        row_active = True
        col_active = True
        # Scan the numeric digits from finest to coarsest, carrying each
        # dimension (row, column) independently into the parent for as
        # long as it keeps landing on that dimension's border -- exactly
        # neighbor()'s own carrying rule, just tracked for two dimensions
        # at once instead of one.
        for i in reversed(range(1, len(self_suid))):
            row, col = cast("tuple[int, int]", child_order[int(self_suid[i])])
            if row_active:
                if (row == 0 and row_dir == "up") or (
                    row == N - 1 and row_dir == "down"
                ):
                    row = N - 1 if row_dir == "up" else 0
                else:
                    row = row - 1 if row_dir == "up" else row + 1
                    row_active = False
            if col_active:
                if (col == 0 and col_dir == "left") or (
                    col == N - 1 and col_dir == "right"
                ):
                    col = N - 1 if col_dir == "left" else 0
                else:
                    col = col - 1 if col_dir == "left" else col + 1
                    col_active = False
            neighbor_suid[i] = cast(int, child_order[(row, col)])
            if not row_active and not col_active:
                # Both dimensions resolved locally: everything coarser,
                # including the resolution 0 cell, is unchanged.
                return Cell(self.rdggs, neighbor_suid)

        # Every numeric digit was on the border in at least one dimension,
        # so the resolution 0 cell itself must change in whichever
        # dimension(s) are still active.
        if row_active and col_active:
            return None
        face = self_suid[0]
        single_dir = row_dir if row_active else col_dir
        neighbor_suid[0] = an[face][single_dir]
        neighbor = Cell(self.rdggs, neighbor_suid)

        # Apply the same rotation correction neighbor() applies for a
        # single step in single_dir, since crossing between a polar
        # (N or S) and non-polar resolution 0 cell needs it regardless of
        # whether that crossing happens alone or as half of a diagonal step.
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
        return neighbor

    def _check_comparable(self, other: "Cell", verb: str) -> None:
        """
        Raise `ValueError` if `self` and `other` aren't cells of the same
        `RHEALPixDGGS`, or if either is the empty cell. Shared precondition
        check for the topological predicates below, whose result would
        otherwise be meaningless (comparing SUIDs from two differently
        configured grids) or undefined (an empty cell has no spatial
        extent to compare).
        """
        if self.rdggs != other.rdggs:
            raise ValueError(
                f"Cannot test {verb} between cells of different RHEALPixDGGS "
                "instances."
            )
        if not self.suid or not other.suid:
            raise ValueError(f"Cannot test {verb} for an empty cell.")

    def equals(self, other: "Cell") -> bool:
        """
        DE-9IM `equals` predicate: return True if this cell and `other`
        are the same cell, and False otherwise. For two cells of one grid
        this is `self == other`; unlike `==`, which is False across grids,
        it raises `ValueError` for cells of different grids or for the
        empty cell, as every predicate here does.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> Cell(rdggs, ['N', 0]).equals(Cell(rdggs, ['N', 0]))
            True
            >>> Cell(rdggs, ['N', 0]).equals(Cell(rdggs, ['N', 1]))
            False

        """
        self._check_comparable(other, "equality")
        return self.suid == other.suid

    def contains_cell(self, other: "Cell") -> bool:
        """
        DE-9IM-style `contains` predicate for a pair of cells: return True
        if `other` is this cell or a descendant of it, and False
        otherwise.

        Named `contains_cell` rather than `contains` to avoid confusion
        with the pre-existing `contains()` method, which tests whether
        this cell contains a *point*, not another cell.

        Note this coincides with `covers()`: because cells form a strict
        hierarchical partition, a descendant cell's boundary that touches
        its ancestor's boundary is still entirely contained in the
        ancestor's closed region, so there's no DE-9IM-style distinction
        here between "contains" (usually excludes boundary-touching) and
        "covers" (usually includes it) the way there can be for general
        geometries.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> Cell(rdggs, ['N']).contains_cell(Cell(rdggs, ['N', 0]))
            True
            >>> Cell(rdggs, ['N', 0]).contains_cell(Cell(rdggs, ['N']))
            False

        """
        self._check_comparable(other, "containment")
        return other.suid[: len(self.suid)] == self.suid

    def within(self, other: "Cell") -> bool:
        """
        DE-9IM `within` predicate: return True if this cell is `other` or
        a descendant of it, and False otherwise. The converse of
        `contains_cell()`; see that method's docstring for why this also
        coincides with `covered_by()` for cells.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> Cell(rdggs, ['N', 0]).within(Cell(rdggs, ['N']))
            True

        """
        return other.contains_cell(self)

    def covers(self, other: "Cell") -> bool:
        """
        DE-9IM `covers` predicate. See `contains_cell()`'s docstring for
        why this is the same relation as `contains_cell()` for cells.
        """
        return self.contains_cell(other)

    def covered_by(self, other: "Cell") -> bool:
        """
        DE-9IM `coveredBy` predicate. See `within()`'s docstring for why
        this is the same relation as `within()` for cells.
        """
        return other.contains_cell(self)

    def touches(self, other: "Cell") -> bool:
        """
        DE-9IM `touches` predicate: return True if this cell and `other`
        share at least one boundary point but neither contains the other
        (equivalently, no ancestor/descendant relationship and no shared
        interior), and False otherwise. `self` and `other` may be of
        different resolutions.

        Two cells of the *same* resolution touch exactly when one is an
        edge or diagonal (corner-only) neighbor of the other -- see
        `neighbor()` and `diagonal_neighbor()`.

        For cells of *different* resolutions (informally, "cousins" --
        neither an ancestor of the other, but nested in siblings that are
        themselves edge/diagonal neighbors, possibly several levels up),
        this finds their two ancestors at the shallower of the two
        resolutions and checks those for edge/diagonal adjacency. If
        they're edge-adjacent, the deeper cell touches the shallower one
        exactly when every one of the deeper cell's digits below that
        ancestor lies in the row or column bordering the shared edge (so
        the deeper cell never strays from that edge); if they're
        diagonally adjacent, the same must hold for the single digit
        value that is that shared corner, since a corner is a single
        point rather than a whole edge. Since the shallower cell exposes
        its entire edge or corner at that resolution, nothing else needs
        checking on its side.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> Cell(rdggs, ['N', 0]).touches(Cell(rdggs, ['N', 1]))
            True
            >>> Cell(rdggs, ['O']).touches(Cell(rdggs, ['P', 3]))
            True
            >>> Cell(rdggs, ['N', 0]).touches(Cell(rdggs, ['N', 8]))
            False

        """
        self._check_comparable(other, "touches")
        return self._touch_dimension(other) is not None

    def _touch_dimension(self, other: "Cell") -> int | None:
        """
        The dimension of the boundary this cell and `other` share when
        they touch: 1 for a segment of an edge, 0 for a single corner
        point, and None when they do not touch (they nest or are
        disjoint). Assumes both cells belong to the same grid and neither
        is empty; `touches()` documents the algorithm.
        """
        if self._nests(other):
            # One is an ancestor of (or the same cell as) the other: their
            # closed regions share interior points, so this isn't touches.
            return None
        assert self.resolution is not None and other.resolution is not None
        r = min(self.resolution, other.resolution)
        shallow, deep = (self, other) if self.resolution == r else (other, self)
        deep_ancestor = Cell(self.rdggs, deep.suid[: r + 1])
        # deep_ancestor != shallow is guaranteed here: if they were equal,
        # the cells would nest and we would have returned above.
        tail = deep.suid[r + 1 :]
        child_order = self.rdggs.child_order
        N = self.N_side
        row_edge = {"up": 0, "down": N - 1}
        col_edge = {"left": 0, "right": N - 1}
        for direction in ("up", "down", "left", "right"):
            if deep_ancestor.neighbor(direction, plane=True) == shallow:
                if direction in row_edge:
                    on_edge = all(
                        cast("tuple[int, int]", child_order[int(d)])[0]
                        == row_edge[direction]
                        for d in tail
                    )
                else:
                    on_edge = all(
                        cast("tuple[int, int]", child_order[int(d)])[1]
                        == col_edge[direction]
                        for d in tail
                    )
                return 1 if on_edge else None
        for direction in ("up_left", "up_right", "down_left", "down_right"):
            if deep_ancestor.diagonal_neighbor(direction) == shallow:
                target_row = 0 if direction.startswith("up") else N - 1
                target_col = 0 if direction.endswith("left") else N - 1
                at_corner = all(
                    child_order[int(d)] == (target_row, target_col) for d in tail
                )
                return 0 if at_corner else None
        return None

    def disjoint(self, other: "Cell") -> bool:
        """
        DE-9IM `disjoint` predicate: return True if this cell and `other`
        share no point at all (no shared interior and no shared
        boundary), and False otherwise.

        Two cells of one hierarchical grid always either nest (one
        contains the other), touch along their boundary only, or are
        disjoint; see `overlaps()` for why partial overlap is impossible.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> Cell(rdggs, ['N', 0]).disjoint(Cell(rdggs, ['S', 0]))
            True
            >>> Cell(rdggs, ['N', 0]).disjoint(Cell(rdggs, ['N', 1]))
            False

        """
        self._check_comparable(other, "disjoint")
        return not (self._nests(other) or self.touches(other))

    def overlaps(self, other: "Cell") -> bool:
        """
        DE-9IM `overlaps` predicate: return True if the interiors of this
        cell and `other` intersect and neither contains the other. Two
        cells of one hierarchical grid always either nest, touch along
        their boundary only, or are disjoint, so this is False for every
        pair of cells; it exists because OGC Topic 21 mandates the
        predicate under this name, and it raises `ValueError` for cells
        of different grids or the empty cell like the other predicates.

        Before 0.9.0 this method answered containment in either direction
        (one suid a prefix of the other's). For that relation use
        ``a.contains_cell(b) or a.within(b)``.

        .. versionchanged:: 0.9.0
           DE-9IM semantics; formerly containment in either direction.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> Cell(rdggs, ['N']).overlaps(Cell(rdggs, ['N', 0]))
            False
            >>> Cell(rdggs, ['N']).contains_cell(Cell(rdggs, ['N', 0]))
            True

        """
        self._check_comparable(other, "overlaps")
        return False

    def region_overlaps(self, region: "list[Cell]") -> bool:
        """
        DE-9IM `overlaps` predicate between this cell and the union of the
        cells in `region`: return True if they share interior points but
        neither covers the other, and False otherwise (including for an
        empty list). That happens exactly when `region` holds descendants
        of this cell that do not tile it completely, and also at least one
        cell outside it. A region containing this cell or an ancestor, or
        whose cells inside this cell tile it, covers the cell (`within`);
        a region entirely inside this cell is covered by it
        (`contains`); a region with no descendant of this cell shares no
        interior with it. Raises `ValueError` if any cell belongs to a
        different grid or is empty.

        With `a` the cell P0 and `region` given as index strings:

        ==================================  ======  ====================
        `region`                            result  why
        ==================================  ======  ====================
        [S8, P03]                           True    partial fill, and S8
        [P0]                                False   covers `a` exactly
        [P03]                               False   `a` covers it
        the nine children of P0, and S8     False   the children tile `a`
        [P1, S8]                            False   no shared interior
        ==================================  ======  ====================

        .. versionchanged:: 0.9.0
           DE-9IM semantics; formerly True if any cell of `region`
           contained, or was contained by, this cell.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> a = Cell(rdggs, ['P', 0])
            >>> a.region_overlaps([Cell(rdggs, ['S', 8]), Cell(rdggs, ['P', 0, 3])])
            True
            >>> a.region_overlaps([Cell(rdggs, ['P', 0, 3])])
            False
            >>> a.region_overlaps(list(a.subcells()) + [Cell(rdggs, ['S', 8])])
            False

        """
        for cell in region:
            self._check_comparable(cell, "overlaps")
        depth = len(self.suid)
        if any(len(cell.suid) <= depth and self._nests(cell) for cell in region):
            return False  # A region cell is this cell or an ancestor.
        inside = [cell for cell in region if self._nests(cell)]
        if not inside or len(inside) == len(region):
            return False  # No shared interior, or the region is inside.
        # Does the region tile this cell? Drop cells nested in another
        # region cell, then compare areas exactly in units of the finest
        # cell present: a cell at resolution r covers N**(2*(rmax - r))
        # such units.
        suids = sorted({cell.suid for cell in inside}, key=len)
        disjoint: list[tuple[str | int, ...]] = []
        for suid in suids:
            if not any(suid[: len(kept)] == kept for kept in disjoint):
                disjoint.append(suid)
        N2 = self.N_side**2
        deepest = max(len(suid) for suid in disjoint)
        covered = sum(N2 ** (deepest - len(suid)) for suid in disjoint)
        return covered != N2 ** (deepest - depth)

    def intersects(self, other: "Cell") -> bool:
        """
        DE-9IM `intersects` predicate: return True if this cell and
        `other` share at least one point, which for cells means they nest
        or touch, and False otherwise. The negation of `disjoint()`.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> Cell(rdggs, ['N', 0]).intersects(Cell(rdggs, ['N', 1]))
            True
            >>> Cell(rdggs, ['N']).intersects(Cell(rdggs, ['N', 0]))
            True
            >>> Cell(rdggs, ['N', 0]).intersects(Cell(rdggs, ['S', 0]))
            False

        """
        self._check_comparable(other, "intersects")
        return self._nests(other) or self.touches(other)

    def crosses(self, other: "Cell") -> bool:
        """
        DE-9IM `crosses` predicate: return True if the interiors of this
        cell and `other` intersect in something of lower dimension than
        at least one of them, which needs inputs of different dimension
        (or two lines). Two cells are regions of the same dimension, so
        this is False for every pair; it exists because OGC Topic 21
        mandates the predicate, and it validates its inputs like the
        other predicates.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> Cell(rdggs, ['N']).crosses(Cell(rdggs, ['N', 0]))
            False

        """
        self._check_comparable(other, "crosses")
        return False

    def distance(self, other: "Cell", plane: bool = True, n: int = 8) -> float:
        """
        Return the distance between this cell and `other`: the infimum
        of the distance between a point of one and a point of the other.
        It is 0 when the cells nest or touch, and otherwise the least
        distance between their boundaries. With `plane` = True that is
        the Euclidean distance in the rHEALPix plane; with `plane` =
        False it is the geodesic (shortest surface) distance on the
        ellipsoid, from ``pyproj.Geod``. Both are in the ellipsoid's
        length units: metres on the WGS84 grids, and on the unit sphere
        the central angle in radians.

        The plane is one unfolding of the cube, so `plane` = True can put
        cells that are near each other on the ellipsoid far apart: cells
        close to the seam at x = +-pi R from opposite sides, or a polar
        cell and an equatorial face its polar square is not attached to
        in the layout. (Cells that touch across such a seam still give 0,
        since touching is decided from the grid's topology.) Use `plane` =
        False for a surface distance.

        On the ellipsoid each cell's boundary is sampled at `n` points
        per edge, the nearest pair of samples is found, and the result is
        refined by minimising the geodesic distance along the nearest
        pair of edges, which are smooth curves, so the answer is accurate
        to the optimiser's tolerance (millimetres on WGS84) rather than
        to the sample spacing. This is the ``distance`` operation of OGC
        Topic 21 v2.0 Table 53; its ``projectTo`` argument selects a
        dimension, and here the distance is taken in the two surface
        dimensions.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import UNIT_003
            >>> from rhealpixdggs.utils import my_round
            >>> rdggs = UNIT_003
            >>> my_round(Cell(rdggs, ['P', 0]).distance(Cell(rdggs, ['P', 2])), 9)
            0.523598776
            >>> Cell(rdggs, ['P', 0]).distance(Cell(rdggs, ['P', 1]))
            0.0
            >>> my_round(Cell(rdggs, ['N']).distance(Cell(rdggs, ['S']), plane=False), 9)
            1.459455312

        """
        self._check_comparable(other, "distance")
        if self._nests(other) or self.touches(other):
            return 0.0
        if plane:
            (ax0, ax1), (ay0, ay1) = self.xy_range()
            (bx0, bx1), (by0, by1) = other.xy_range()
            dx = max(0.0, ax0 - bx1, bx0 - ax1)
            dy = max(0.0, ay0 - by1, by0 - ay1)
            return float(hypot(dx, dy))
        return self._geodesic_distance(other, max(n, 2))

    def _boundary_samples(self, n: int) -> tuple[FloatArray, FloatArray, np.ndarray]:
        """
        `n` points along each of the four planar edges of this cell,
        projected to the ellipsoid: their longitudes, their latitudes and
        the index (0 to 3, clockwise from the upper-left corner) of the
        edge each lies on.
        """
        corners = np.array(self.vertices(plane=True))
        ends = np.roll(corners, -1, axis=0)
        t = np.linspace(0.0, 1.0, n)
        x = (corners[:, None, 0] + t * (ends[:, None, 0] - corners[:, None, 0])).ravel()
        y = (corners[:, None, 1] + t * (ends[:, None, 1] - corners[:, None, 1])).ravel()
        lon, lat = self.rdggs.rhealpix(x, y, inverse=True, region=self.region())
        return lon, lat, np.repeat(np.arange(4), n)

    def _geodesic_distance(self, other: "Cell", n: int) -> float:
        """
        The least geodesic distance between the boundaries of this cell
        and `other`, which do not intersect: a nearest pair among `n`
        samples per edge, then a bounded minimisation along every pair of
        edges whose samples could hide the true minimum.
        """
        E = self.ellipsoid
        geod = _geod(E.a, E.f)
        radians = E.radians
        lon_a, lat_a, edge_a = self._boundary_samples(n)
        lon_b, lat_b, edge_b = other._boundary_samples(n)

        def between(
            u0: FloatArray, v0: FloatArray, u1: FloatArray, v1: FloatArray
        ) -> FloatArray:
            return np.asarray(geod.inv(u0, v0, u1, v1, radians=radians)[2], dtype=float)

        m, k = len(lon_a), len(lon_b)
        dist = between(
            np.repeat(lon_a, k),
            np.repeat(lat_a, k),
            np.tile(lon_b, m),
            np.tile(lat_b, m),
        ).reshape(m, k)
        # Any boundary point lies within half a sample step of a sample, so
        # the sampled minimum overestimates the true one by at most the
        # largest step on either boundary; only edge pairs within that
        # margin of the best sample can hold the true minimum.
        spacing = max(
            float(between(lon_a[:-1], lat_a[:-1], lon_a[1:], lat_a[1:]).max()),
            float(between(lon_b[:-1], lat_b[:-1], lon_b[1:], lat_b[1:]).max()),
        )
        best = float(dist.min())
        corners_a = np.array(self.vertices(plane=True))
        corners_b = np.array(other.vertices(plane=True))
        region_a, region_b = self.region(), other.region()

        def objective(
            t: np.ndarray,
            a0: np.ndarray,
            a1: np.ndarray,
            b0: np.ndarray,
            b1: np.ndarray,
        ) -> float:
            # Geodesic distance between the points at fractions t[0] and
            # t[1] along the planar edges a0-a1 of this cell and b0-b1 of
            # the other.
            ax, ay = a0 + t[0] * (a1 - a0)
            bx, by = b0 + t[1] * (b1 - b0)
            pa = self.rdggs.rhealpix(
                float(ax), float(ay), inverse=True, region=region_a
            )
            pb = other.rdggs.rhealpix(
                float(bx), float(by), inverse=True, region=region_b
            )
            return float(geod.inv(pa[0], pa[1], pb[0], pb[1], radians=radians)[2])

        result = best
        for i in range(4):
            a0, a1 = corners_a[i], corners_a[(i + 1) % 4]
            for j in range(4):
                block = dist[np.ix_(edge_a == i, edge_b == j)]
                if block.min() > best + 2 * spacing:
                    continue
                b0, b1 = corners_b[j], corners_b[(j + 1) % 4]
                row, col = np.unravel_index(block.argmin(), block.shape)
                start = np.array([row, col], dtype=float) / (n - 1)
                fit = optimize.minimize(
                    objective,
                    start,
                    args=(a0, a1, b0, b1),
                    method="L-BFGS-B",
                    bounds=[(0.0, 1.0), (0.0, 1.0)],
                )
                result = min(result, float(fit.fun))
        return result

    def within_distance(
        self, other: "Cell", dist: float, plane: bool = True, n: int = 8
    ) -> bool:
        """
        Return True if `distance(other, plane, n)` is less than `dist`,
        and False otherwise: the ``withinDistance`` operation of OGC Topic
        21 v2.0 Table 53, with its strict inequality. See `distance()` for
        the definition and units.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import UNIT_003
            >>> rdggs = UNIT_003
            >>> Cell(rdggs, ['P', 0]).within_distance(Cell(rdggs, ['P', 2]), 0.6)
            True
            >>> Cell(rdggs, ['P', 0]).within_distance(Cell(rdggs, ['P', 2]), 0.5)
            False

        """
        return self.distance(other, plane=plane, n=n) < dist

    def relative_position(
        self, other: "Cell", direction: tuple[float, float] = (1.0, 0.0)
    ) -> RelativePosition:
        """
        Return the position of this cell relative to `other` along
        `direction`, as a `RelativePosition`: project both planar cells
        onto the line through the origin with that direction, giving each
        a closed interval from its smallest to its largest coordinate,
        and classify the two intervals. This is the ``relativePosition``
        operation of OGC Topic 21 v2.0 Table 53; `direction` is the
        spatial part of its ``projectTo`` vector, which for a
        two-dimensional grid can only have those two components nonzero.

        The default `direction` (1, 0) is the planar x axis, west to east
        in the equatorial belt; (0, 1) is the planar y axis, south to
        north there. Polar cells' orientation depends on
        ``north_square`` and ``south_square``, and the plane is one
        unfolding of the cube, so cells that meet across the seam at
        x = +-pi R are far apart along x. Along either axis two cells of
        one hierarchy only ever give EQUALS, the nesting relations,
        MEETS/MET_BY or BEFORE/AFTER, because their edges line up;
        OVERLAPS needs a diagonal direction. Raise `ValueError` for a zero
        vector.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> p0, p1 = Cell(rdggs, ['P', 0]), Cell(rdggs, ['P', 1])
            >>> p0.relative_position(p1).name
            'MEETS'
            >>> p1.relative_position(Cell(rdggs, ['P'])).name
            'DURING'
            >>> p0.relative_position(p1, direction=(1, 1)).name
            'OVERLAPS'
            >>> Cell(rdggs, ['P', 3]).relative_position(p0, direction=(0, 1)).name
            'MEETS'

        """
        self._check_comparable(other, "relative position")
        dx, dy = float(direction[0]), float(direction[1])
        norm = hypot(dx, dy)
        if norm == 0:
            raise ValueError("direction must be a nonzero vector.")

        def interval(cell: "Cell") -> tuple[float, float]:
            values = [x * dx + y * dy for x, y in cell.vertices(plane=True)]
            return min(values), max(values)

        a0, a1 = interval(self)
        b0, b1 = interval(other)
        tol = 1e-9 * min(self.width(), other.width()) * norm

        def same(u: float, v: float) -> bool:
            return abs(u - v) <= tol

        RP = RelativePosition
        if a1 < b0 - tol:
            return RP.BEFORE
        if a0 > b1 + tol:
            return RP.AFTER
        if same(a1, b0):
            return RP.MEETS
        if same(a0, b1):
            return RP.MET_BY
        begin, end = same(a0, b0), same(a1, b1)
        if begin and end:
            return RP.EQUALS
        if begin:
            return RP.STARTS if a1 < b1 else RP.STARTED_BY
        if end:
            return RP.FINISHES if a0 > b0 else RP.FINISHED_BY
        if a0 > b0 and a1 < b1:
            return RP.DURING
        if a0 < b0 and a1 > b1:
            return RP.CONTAINS
        return RP.OVERLAPS if a0 < b0 else RP.OVERLAPPED_BY

    def relate_position(
        self,
        other: "Cell",
        relate: RelativePosition,
        direction: tuple[float, float] = (1.0, 0.0),
    ) -> bool:
        """
        Return True if `relative_position(other, direction)` is `relate`,
        or belongs to it when `relate` is one of the groupings
        ``RelativePosition.IN`` (STARTS, DURING or FINISHES) and
        ``RelativePosition.DISJOINT`` (BEFORE or AFTER). The
        ``relatePosition`` operation of OGC Topic 21 v2.0 Table 53.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> p1 = Cell(rdggs, ['P', 1])
            >>> p1.relate_position(Cell(rdggs, ['P']), RelativePosition.IN)
            True
            >>> p1.relate_position(Cell(rdggs, ['P', 0]), RelativePosition.DISJOINT)
            False

        """
        position = self.relative_position(other, direction)
        if relate is RelativePosition.IN:
            return position.is_in
        if relate is RelativePosition.DISJOINT:
            return position.is_disjoint
        return position is relate

    def _on_boundary_of(self, ancestor: "Cell") -> bool:
        """
        True if this cell, a proper descendant of `ancestor`, touches the
        ancestor's boundary: along some side, every digit below the
        ancestor keeps the extreme row or column.
        """
        child_order = self.rdggs.child_order
        N = self.N_side
        rows_cols = [
            cast("tuple[int, int]", child_order[int(d)])
            for d in self.suid[len(ancestor.suid) :]
        ]
        return any(
            all(rc[axis] == extreme for rc in rows_cols)
            for axis in (0, 1)
            for extreme in (0, N - 1)
        )

    def _de9im(self, other: "Cell") -> str:
        """
        The DE-9IM matrix of this cell against `other` as nine
        characters, row by row: this cell's interior, boundary and
        exterior against `other`'s, each entry the dimension of the
        intersection (0, 1 or 2) or F for empty. Assumes both cells belong
        to the same grid and neither is empty.
        """
        if self.suid == other.suid:
            return "2FFF1FFF2"
        if self._nests(other):
            if len(self.suid) < len(other.suid):
                shared = other._on_boundary_of(self)
                return "212F11FF2" if shared else "212FF1FF2"
            shared = self._on_boundary_of(other)
            return "2FF11F212" if shared else "2FF1FF212"
        touch = self._touch_dimension(other)
        if touch is None:
            return "FF2FF1212"
        return "FF2F11212" if touch == 1 else "FF2F01212"

    def relate(self, other: "Cell", matrix: str) -> bool:
        """
        DE-9IM `relate` predicate: return True if the DE-9IM matrix of
        this cell against `other` matches the nine-character pattern
        `matrix`, read row by row as this cell's interior, boundary and
        exterior against `other`'s. Pattern characters are ``T`` (a
        non-empty intersection of any dimension), ``F`` (empty), ``*``
        (anything) and ``0``, ``1``, ``2`` (that dimension exactly). The
        ``relate`` operation of OGC Topic 21 v2.0 Table 53.

        Between two cells of one hierarchy only six matrices occur: equal
        cells ``2FFF1FFF2``; an ancestor against a descendant
        ``212F11FF2`` when the descendant touches the ancestor's boundary
        and ``212FF1FF2`` when it is interior, and their transposes for a
        descendant against its ancestor; cells touching along an edge
        ``FF2F11212``, at a corner ``FF2F01212``; and disjoint cells
        ``FF2FF1212``.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> p0, p1 = Cell(rdggs, ['P', 0]), Cell(rdggs, ['P', 1])
            >>> p0.relate(p1, 'FF*F1****')  # touching along an edge
            True
            >>> p0.relate(p1, 'T********')  # interiors intersect
            False
            >>> Cell(rdggs, ['P']).relate(p0, 'T*****FF*')  # contains
            True

        """
        self._check_comparable(other, "relate")
        if len(matrix) != 9 or any(c not in "TF*012" for c in matrix):
            raise ValueError(
                "matrix must be nine characters from 'T', 'F', '*', '0', '1', "
                f"'2', not {matrix!r}."
            )
        for want, have in zip(matrix, self._de9im(other)):
            if want == "*":
                continue
            if want == "T":
                if have == "F":
                    return False
            elif want != have:
                return False
        return True

    def parent_of(self, other: "Cell", inherit_id: bool = False) -> bool:
        """
        Return True if this cell is the parent of `other`, one level up:
        ``parentOf`` of OGC Topic 21 v2.0 Table 53. `inherit_id` is
        accepted for the spec's signature; in rHEALPix's nested hierarchy
        a child always inherits its parent's identifier prefix, so it
        changes nothing.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> Cell(rdggs, ['P', 4]).parent_of(Cell(rdggs, ['P', 4, 0]))
            True
            >>> Cell(rdggs, ['P', 4]).parent_of(Cell(rdggs, ['P', 4, 0, 0]))
            False

        """
        self._check_comparable(other, "parenthood")
        return len(other.suid) == len(self.suid) + 1 and other.suid[:-1] == self.suid

    def child_of(self, other: "Cell", inherit_id: bool = False) -> bool:
        """
        Return True if this cell is a child of `other`, one level down:
        ``childOf`` of Table 53, the converse of `parent_of()`.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> Cell(rdggs, ['P', 4, 0]).child_of(Cell(rdggs, ['P', 4]))
            True

        """
        return other.parent_of(self, inherit_id)

    def sibling_of(self, other: "Cell", inherit_id: bool = False) -> bool:
        """
        Return True if `other` is a sibling of this cell, itself included:
        ``siblingOf`` of Table 53. With `inherit_id` True siblings share a
        parent; with `inherit_id` False (the default) inheritance is
        ignored and siblings are the cells at the same resolution that
        share an edge or a corner, so the neighbours across a parent's
        edge count too, as the spec's example ``40.siblingOf(31)`` shows.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import RHEALPixDGGS
            >>> rdggs = RHEALPixDGGS()
            >>> p40, p32 = Cell(rdggs, ['P', 4, 0]), Cell(rdggs, ['P', 3, 2])
            >>> p40.sibling_of(Cell(rdggs, ['P', 4, 1]))
            True
            >>> p40.sibling_of(p32), p40.sibling_of(p32, inherit_id=True)
            (True, False)

        """
        self._check_comparable(other, "siblinghood")
        if len(self.suid) != len(other.suid):
            return False
        if self.suid == other.suid:
            return True
        if inherit_id:
            return len(self.suid) > 1 and self.suid[:-1] == other.suid[:-1]
        return self._touch_dimension(other) is not None

    def union(
        self,
        other: "Cell | ZoneSet",
        min_res: int | None = None,
        max_res: int | None = None,
    ) -> "ZoneSet":
        """The ``union`` of Table 53: see ``ZoneSet.union``."""
        from rhealpixdggs.zoneset import ZoneSet

        return ZoneSet(self.rdggs, [self]).union(other, min_res, max_res)

    def intersection(
        self,
        other: "Cell | ZoneSet",
        min_res: int | None = None,
        max_res: int | None = None,
    ) -> "ZoneSet":
        """The ``intersection`` of Table 53: see ``ZoneSet.intersection``."""
        from rhealpixdggs.zoneset import ZoneSet

        return ZoneSet(self.rdggs, [self]).intersection(other, min_res, max_res)

    def difference(
        self,
        other: "Cell | ZoneSet",
        min_res: int | None = None,
        max_res: int | None = None,
    ) -> "ZoneSet":
        """The ``difference`` of Table 53: see ``ZoneSet.difference``."""
        from rhealpixdggs.zoneset import ZoneSet

        return ZoneSet(self.rdggs, [self]).difference(other, min_res, max_res)

    def sym_difference(
        self,
        other: "Cell | ZoneSet",
        min_res: int | None = None,
        max_res: int | None = None,
    ) -> "ZoneSet":
        """The ``symDifference`` of Table 53: see ``ZoneSet.sym_difference``."""
        from rhealpixdggs.zoneset import ZoneSet

        return ZoneSet(self.rdggs, [self]).sym_difference(other, min_res, max_res)

    def buffer(self, dist: float, plane: bool = True, n: int = 8) -> "ZoneSet":
        """The ``buffer`` of Table 53: see ``ZoneSet.buffer``."""
        from rhealpixdggs.zoneset import ZoneSet

        return ZoneSet(self.rdggs, [self]).buffer(dist, plane, n)

    def parent(self, levels: int = 1, inherit_id: bool = False) -> "ZoneSet":
        """The ``parent`` of Table 53: see ``ZoneSet.parent``."""
        from rhealpixdggs.zoneset import ZoneSet

        return ZoneSet(self.rdggs, [self]).parent(levels, inherit_id)

    def child(self, levels: int = 1, inherit_id: bool = False) -> "ZoneSet":
        """The ``child`` of Table 53: see ``ZoneSet.child``."""
        from rhealpixdggs.zoneset import ZoneSet

        return ZoneSet(self.rdggs, [self]).child(levels, inherit_id)

    def sibling(self, levels: int = 1, inherit_id: bool = False) -> "ZoneSet":
        """The ``sibling`` of Table 53: see ``ZoneSet.sibling``."""
        from rhealpixdggs.zoneset import ZoneSet

        return ZoneSet(self.rdggs, [self]).sibling(levels, inherit_id)

    def random_point(self, plane: bool = True) -> tuple[float, float]:
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
            if self.ellipsoidal_shape == "cap":
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

    def color(self, saturation: float = 0.5) -> tuple[float, float, float]:
        """
        Return a unique RGB color tuple for this cell.
        Inessential graphics method.
        """
        suid = self.suid
        N = self.rdggs.N_side
        hue_resolution0 = {v: k / 6.0 for (k, v) in enumerate(CELLS0)}
        hue = hue_resolution0[str(suid[0])]
        n = len(suid)
        if n > 1:
            hue += sum([suid[i] * N ** (-2 * i) for i in range(1, n)]) / 6.0
            # hue += sum([suid[i + 1]*N**(resolution - 1 - i)
            #             for i in range(resolution)])/\
            #        float(6*N**(2*resolution))
        return hsv_to_rgb(hue, saturation, 1)
