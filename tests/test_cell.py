# Import third-party modules.
# Import standard modules
import unittest
from itertools import product
from math import pi

from scipy.spatial.distance import euclidean, norm

# Import my modules.
from rhealpixdggs.cell import CELLS0, Cell
from rhealpixdggs.dggs import WGS84_003, WGS84_003_RADIANS, RHEALPixDGGS
from rhealpixdggs.ellipsoids import (
    WGS84_ASPHERE_RADIANS,
    WGS84_ELLIPSOID,
    WGS84_ELLIPSOID_RADIANS,
    Ellipsoid,
)

# Level 0 cell names
N = CELLS0[0]
O = CELLS0[1]
P = CELLS0[2]
Q = CELLS0[3]
R = CELLS0[4]
S = CELLS0[5]

# Some DGGSs used below.
WGS84_123 = RHEALPixDGGS(
    ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=2, N_side=3
)
WGS84_123_RADIANS = RHEALPixDGGS(
    ellipsoid=WGS84_ELLIPSOID_RADIANS, north_square=1, south_square=2, N_side=3
)
WGS84_122 = RHEALPixDGGS(
    ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=2, N_side=2
)


# Relative error function.
def rel_err(get, expect):
    a = euclidean(get, expect)
    b = norm(expect)
    if b == 0:
        return a
    else:
        return a / b


class SCENZGridCELLTestCase(unittest.TestCase):
    def test_suid_from_index(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            # Test order='level'.
            # The suid with index 0 should be N.
            suid = (N,)
            get = Cell.suid_from_index(rdggs, 0, order="level")
            self.assertEqual(get, suid)

            # Should be correct on cell P1.
            suid = (P, 1)
            i = 25
            get = Cell.suid_from_index(rdggs, i, order="level")
            self.assertEqual(get, suid)

            # Test order='post'.
            # The suid with index 0 should be N00...0.
            suid = [N] + [0 for i in range(rdggs.max_resolution)]
            suid = tuple(suid)
            get = Cell.suid_from_index(rdggs, 0, order="post")
            self.assertEqual(get, suid)

            # Should be correct on cell P1.
            suid = (P, 1)

            def num(k):
                return WGS84_123.num_cells(res_1=k, subcells=True)

            i = 2 * num(0) + 1 * num(1) + num(1) - 1
            get = Cell.suid_from_index(rdggs, i, order="post")
            self.assertEqual(get, suid)

    def test_Cell_init(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            # Should set the rdggs, suid, and resolution attributes correctly.
            suid = (S, 0, 1, 2, 6, 3, 1)
            resolution = 6
            C = Cell(rdggs, suid)
            self.assertEqual(C.rdggs, rdggs)
            self.assertEqual(C.suid, suid)
            self.assertEqual(C.resolution, resolution)

            # Should not create invalid cells.
            suid = (P, rdggs.N_side**2)
            self.assertRaises(ValueError, Cell, rdggs, suid)

            # Should create cell P1.
            expect = (P, 1)
            i = 25
            get = Cell(rdggs, level_order_index=i).suid
            self.assertEqual(get, expect)

            expect = (P, 1)

            def num(k):
                # rdggs is stable across this loop iteration, and num() is
                # only called within it.
                return rdggs.num_cells(res_1=k, subcells=True)  # noqa: B023

            i = 2 * num(0) + 1 * num(1) + num(1) - 1
            get = Cell(rdggs, post_order_index=i).suid
            self.assertEqual(get, expect)

    def test_Cell_init_invalid_suid_raises(self):
        # Regression test for issue #54: these were all bare `assert`s,
        # which python -O silently compiles out, and which raised
        # AssertionError regardless of what actually went wrong. Each
        # should now raise a specific, conventional exception.
        rdggs = WGS84_123
        # Wrong type entirely.
        with self.assertRaises(TypeError):
            Cell(rdggs, suid="N0")
        with self.assertRaises(TypeError):
            Cell(rdggs, suid={"N", 0})
        # Right type (tuple/list), wrong length.
        with self.assertRaises(ValueError):
            Cell(rdggs, suid=())
        with self.assertRaises(ValueError):
            Cell(rdggs, suid=(N,) * (rdggs.max_resolution + 2))
        # suid[0] not a valid level-0 cell name.
        with self.assertRaises(ValueError):
            Cell(rdggs, suid=("X", 0))
        # A later digit out of range for this DGGS's N_side.
        with self.assertRaises(ValueError):
            Cell(rdggs, suid=(P, rdggs.N_side**2))

    def test_suid_rowcol(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            # Should work for resolution 0 cells.
            for suid in CELLS0:
                self.assertEqual(suid, Cell(rdggs, (suid,)).suid_rowcol()[0][0])
                self.assertEqual(suid, Cell(rdggs, (suid,)).suid_rowcol()[1][0])

            # Should work on an arbitrary cell.
            n = rdggs.N_side
            c = Cell(rdggs, (N, 1, 3, n**2 - 1, 2 * n))
            row, col = c.suid_rowcol()
            for i in range(1, c.resolution + 1):
                self.assertEqual(c.suid[i], n * row[i] + col[i])

    def test_eq(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            A = Cell(rdggs, (O, 1, 2, 3))
            AA = Cell(rdggs, (O, 1, 2, 3))
            B = Cell(rdggs, (N, 2, 8, 7, 5))
            self.assertTrue(A == AA)
            self.assertFalse(A == B)

    def test_ne(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            A = Cell(rdggs, (O, 1, 2, 3))
            AA = Cell(rdggs, (O, 1, 2, 3))
            B = Cell(rdggs, (N, 2, 8, 7, 5))
            self.assertFalse(A != AA)
            self.assertTrue(A != B)

    def test_lt(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            a = Cell(rdggs, (N, 7, 6, 8, 1))
            b = Cell(rdggs, (O, 1, 2, 3))
            c = Cell(rdggs, (O, 1, 2))
            self.assertTrue(a < b)
            self.assertTrue(a < c)
            self.assertFalse(a < a)  # noqa: PLR0124 -- tests < irreflexivity
            self.assertFalse(b < a)
            self.assertFalse(c < a)

    def test_le(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            a = Cell(rdggs, (N, 7, 6, 8, 1))
            b = Cell(rdggs, (O, 1, 2, 3))
            bb = Cell(rdggs, (O, 1, 2, 3))
            c = Cell(rdggs, (O, 1, 2))
            self.assertTrue(b <= bb)
            self.assertTrue(bb <= b)
            self.assertTrue(a <= b)
            self.assertTrue(a <= c)
            self.assertFalse(b <= a)
            self.assertFalse(c <= a)

    def test_le_N_side_ge_4(self):
        # Regression test for issue #71: __le__ used to lexicographically
        # compare a string rendering of the suid (comma-joined), which
        # breaks for N_side >= 4 once some digits become multi-character
        # (e.g. 10-15 for N_side=4): "9" > "1" as strings, so same-
        # resolution siblings (N, 9) and (N, 10) came out mis-ordered
        # despite 9 < 10 numerically. Compare the suid tuple itself
        # instead.
        rdggs = RHEALPixDGGS(N_side=4)
        c9 = Cell(rdggs, (N, 9))
        c10 = Cell(rdggs, (N, 10))
        self.assertTrue(c9 <= c10)
        self.assertFalse(c10 <= c9)
        self.assertTrue(c9 < c10)
        self.assertFalse(c10 < c9)

    def test_gt(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            a = Cell(rdggs, (N, 7, 6, 8, 1))
            b = Cell(rdggs, (O, 1, 2, 3))
            c = Cell(rdggs, (O, 1, 2))
            self.assertFalse(a > b)
            self.assertFalse(a > c)
            self.assertFalse(a > a)  # noqa: PLR0124 -- tests > irreflexivity
            self.assertTrue(b > a)
            self.assertTrue(c > a)

    def test_ge(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            a = Cell(rdggs, (N, 7, 6, 8, 1))
            b = Cell(rdggs, (O, 1, 2, 3))
            bb = Cell(rdggs, (O, 1, 2, 3))
            c = Cell(rdggs, (O, 1, 2))
            self.assertTrue(b >= bb)
            self.assertTrue(bb >= a)
            self.assertFalse(a >= b)
            self.assertFalse(a >= c)
            self.assertTrue(b >= a)
            self.assertTrue(c >= a)

    def test_index(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            # Empty cell should have index None
            empty = rdggs.cell()
            self.assertEqual(empty.index(order="level"), None)
            self.assertEqual(empty.index(order="post"), None)

            # Test order='level'.
            # A cell should have index one greater than its predecessor
            # at the same resolution.
            c = rdggs.cell((N, 2, 7, 4, 8))
            b = c.predecessor()
            self.assertEqual(c.index(order="level"), b.index(order="level") + 1)

            # It should invert suid_from_index().
            # A fixed index; originally drawn as
            # randint(0, rdggs.num_cells(0, rdggs.max_resolution) - 1).
            k = 3616048
            get = Cell(rdggs, Cell.suid_from_index(rdggs, k, order="level")).index(
                order="level"
            )
            self.assertEqual(get, k)

            # Test order='post'.
            # The cell N00...0 should have index 0.
            suid = [N] + [0] * rdggs.max_resolution
            a = Cell(rdggs, suid)
            self.assertEqual(a.index(order="post"), 0)

            # A cell should have index one greater than its last child.
            cc = list(c.subcells())[8]
            self.assertEqual(c.index(order="post"), cc.index(order="post") + 1)

            # It should invert suid_from_index().
            k = 3616048  # randint(0, M - 1)
            get = Cell(rdggs, Cell.suid_from_index(rdggs, k, order="post")).index(
                order="post"
            )
            self.assertEqual(get, k)

    def test_contains(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            # A cell should contain its nucleus, but not my test point p.
            # Assume that nucleus() and vertices() work.
            for suid in [[N, 3, 1], [P, 5, 7, 5, 1, 3], [S, 0]]:
                c = Cell(rdggs, suid)
                for plane in [True, False]:
                    nucleus = c.nucleus(plane=plane)
                    vertices = c.vertices(plane=plane)
                    p = (max([v[0] for v in vertices]) + 1, vertices[3][1])
                    self.assertTrue(c.contains(nucleus, plane=plane))
                    self.assertFalse(c.contains(p, plane=plane))

    def test_intersects_meridian(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            PI = rdggs.ellipsoid.pi()
            c = rdggs.cell([S, 4])
            self.assertTrue(c.intersects_meridian(-PI))
            self.assertTrue(c.intersects_meridian(-PI / 2))
            self.assertTrue(c.intersects_meridian(0))
            self.assertTrue(c.intersects_meridian(PI / 2))

            c = rdggs.cell([N, 6])
            self.assertTrue(c.intersects_meridian(-PI))
            self.assertTrue(c.intersects_meridian(-PI + 0.1))
            self.assertTrue(c.intersects_meridian(PI - 0.1))
            self.assertFalse(c.intersects_meridian(-PI / 2))
            self.assertFalse(c.intersects_meridian(0))
            self.assertFalse(c.intersects_meridian(PI / 2))

            c = rdggs.cell([P])
            self.assertFalse(c.intersects_meridian(-PI))
            self.assertTrue(c.intersects_meridian(-PI / 2))
            self.assertTrue(c.intersects_meridian(-PI / 2 + 0.1))
            self.assertTrue(c.intersects_meridian(0))
            self.assertFalse(c.intersects_meridian(0.1))
            self.assertFalse(c.intersects_meridian(PI / 2))

    def test_intersects_parallel(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            PI = rdggs.ellipsoid.pi()
            c = rdggs.cell([N])
            self.assertFalse(c.intersects_parallel(0.1))
            self.assertTrue(c.intersects_parallel(PI / 3))
            self.assertTrue(c.intersects_parallel(PI / 2))

            c = rdggs.cell([S])
            self.assertFalse(c.intersects_parallel(0.1))
            self.assertTrue(c.intersects_parallel(-PI / 3))
            self.assertTrue(c.intersects_parallel(-PI / 2))

            c = rdggs.cell([P])
            self.assertTrue(c.intersects_parallel(0))
            self.assertTrue(c.intersects_parallel(PI / 5))
            self.assertTrue(c.intersects_parallel(-PI / 5))
            self.assertFalse(c.intersects_parallel(PI / 3))
            self.assertFalse(c.intersects_parallel(-PI / 3))

    def test_successor(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            A = Cell(rdggs, (N, 8, 1))
            As0 = Cell(rdggs, (O,))
            As1 = Cell(rdggs, (O, 0))
            As2 = Cell(rdggs, (N, 8, 2))
            As3 = Cell(rdggs, (N, 8, 2, 0))
            self.assertTrue(A.successor(0) == As0)
            self.assertTrue(A.successor(1) == As1)
            self.assertTrue(A.successor(2) == As2)
            self.assertTrue(A.successor(3) == As3)
            A = Cell(rdggs, (S, 8, 8))
            self.assertTrue(A.successor(1) == None)

    def test_predecessor(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            A = Cell(rdggs, (O, 0, 1))
            Ap0 = Cell(rdggs, (N,))
            Ap1 = Cell(rdggs, (N, 8))
            Ap2 = Cell(rdggs, (O, 0, 0))
            Ap3 = Cell(rdggs, (O, 0, 1, 8))
            self.assertTrue(A.predecessor(0) == Ap0)
            self.assertTrue(A.predecessor(1) == Ap1)
            self.assertTrue(A.predecessor(2) == Ap2)
            self.assertTrue(A.predecessor(3) == Ap3)
            A = Cell(rdggs, (N, 0, 0))
            self.assertTrue(A.predecessor(1) == None)

    def test_subcell(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            C = Cell(rdggs, (S, 1, 2, 0))
            D = Cell(rdggs, (S, 1, 2, 0, 5))
            self.assertTrue(D.subcell(C))
            self.assertFalse(C.subcell(D))

    def test_subcell_N_side_ge_4(self):
        # Regression test for issue #71: subcell() used to compare a
        # string rendering of the suid (comma-joined) with startswith(),
        # which breaks for N_side >= 4 once some digits become
        # multi-character (e.g. 10-15 for N_side=4): "N,15" starts with
        # "N,1" as a string, incorrectly making (N, 15) look like a
        # subcell of (N, 1) even though they're same-resolution siblings,
        # not ancestor and descendant.
        rdggs = RHEALPixDGGS(N_side=4)
        a = Cell(rdggs, (N, 15))
        b = Cell(rdggs, (N, 1))
        self.assertFalse(a.subcell(b))
        self.assertFalse(b.subcell(a))
        # A genuine ancestor/descendant pair should still work.
        parent = Cell(rdggs, (N, 1))
        child = Cell(rdggs, (N, 1, 15))
        self.assertTrue(child.subcell(parent))
        self.assertFalse(parent.subcell(child))

    def test_subcells(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            l = 6
            C = Cell(rdggs, (S, 1, 0, 5, 7, 7, 3))
            s = list(C.subcells(l + 1))
            t = [Cell(rdggs, list(C.suid) + [i]) for i in range(9)]
            for i in range(9):
                self.assertTrue(s[i] == t[i])

    def test_ul_vertex(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            # Should work for resolution 0 cells.
            for X in rdggs.grid(0):
                get = X.ul_vertex()
                expect = rdggs.ul_vertex[X.suid[0]]
                self.assertEqual(get, expect)

            for X in rdggs.grid(0):
                get = X.ul_vertex(plane=False)
                p = rdggs.ul_vertex[X.suid[0]]
                expect = rdggs.rhealpix(*p, inverse=True)
                self.assertEqual(get, expect)

            # Should work on children cells.
            a = Cell(rdggs, (S, 2, 3))
            l = a.resolution
            x, y = a.ul_vertex()
            w = rdggs.cell_width(l + 1)
            error = 1e-10  # Error tolerance.
            for i, j in product(list(range(3)), repeat=2):
                b = Cell(rdggs, list(a.suid) + [i + 3 * j])
                xx, yy = b.ul_vertex()
                xp, yp = (x + i * w, y - j * w)
                self.assertTrue(rel_err([xx, yy], [xp, yp]) < error)

            a = Cell(rdggs, (S, 2, 3))
            l = a.resolution
            x, y = a.ul_vertex()
            w = rdggs.cell_width(l + 1)
            error = rdggs.ellipsoid.R_A * 1e-15  # Error tolerance.
            for i, j in product(list(range(3)), repeat=2):
                b = Cell(rdggs, list(a.suid) + [i + 3 * j])
                xx, yy = b.ul_vertex(plane=False)
                xp, yp = rdggs.rhealpix(x + i * w, y - j * w, inverse=True)
                self.assertTrue(rel_err([xx, yy], [xp, yp]) < error)

    def test_boundary(self):
        # Regression test for the `n < 2` guard (was `n == 2`, a
        # comparison instead of an assignment, so it never actually
        # clamped n). Before the fix, boundary(n=1, plane=True) crashed
        # with ZeroDivisionError and boundary(n=0, plane=True) silently
        # returned an empty list instead of being clamped to n=2. See
        # issue #48.
        rdggs = RHEALPixDGGS()
        a = Cell(rdggs, [N, 2])
        self.assertEqual(a.ellipsoidal_shape, "dart")
        expect = a.boundary(n=2, plane=True)
        for n in (0, 1):
            self.assertEqual(a.boundary(n=n, plane=True), expect)

    def test_boundary_quad_cap_n_contract(self):
        # Regression test for the boundary(plane=False) short-circuit that
        # silently returned only 4 points for quad/cap cells regardless of
        # n, instead of the documented 4*n - 4. See issue #49.
        rdggs = RHEALPixDGGS()
        quad = Cell(rdggs, [P, 2])
        cap = Cell(rdggs, [N])
        for c in (quad, cap):
            for n in (2, 3, 5):
                b = c.boundary(n=n, plane=False)
                self.assertEqual(len(b), max(4 * n - 4, 4))
            # n=2 (the short-circuited case) must still agree with vertices().
            self.assertEqual(c.boundary(n=2, plane=False), c.vertices(plane=False))

    def test_nucleus(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            # Nuclei of children should be in correct position
            # relative to parent cell in rHEALPix projection.
            a = Cell(rdggs, (S, 7, 4, 1, 2, 1))  # Arbitrary cell.
            w = a.width()
            (x, y) = a.ul_vertex()
            error = 1e-10
            for row, col in product(list(range(3)), repeat=2):
                # Child cell in (row, column) position relative to a:
                b = Cell(rdggs, list(a.suid) + [3 * row + col])
                get = b.nucleus()
                expect = (x + w / 6 + (w / 3) * col, y - w / 6 - (w / 3) * row)
                self.assertTrue(rel_err(get, expect) < error)

    def test_rotate_entry(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            X = Cell(rdggs, (N,))
            # Should return correct values.
            s = []  # Function values.
            t = []  # Correct values.
            for q in range(4):
                s.append([X.rotate_entry(x, q) for x in range(9)])
            t.append([0, 1, 2, 3, 4, 5, 6, 7, 8])
            t.append([2, 5, 8, 1, 4, 7, 0, 3, 6])
            t.append([8, 7, 6, 5, 4, 3, 2, 1, 0])
            t.append([6, 3, 0, 7, 4, 1, 8, 5, 2])
            for q in range(4):
                self.assertEqual(s[q], t[q])

    def test_neighbor(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            # Plane test.
            c = Cell(rdggs, (N, 0))
            get = c.neighbor("left")
            expect = Cell(rdggs, (O, 0))
            self.assertEqual(get, expect)

            # Ellipsoid test.
            c = Cell(rdggs, (O, 0))
            get = c.neighbor("east", plane=False)
            expect = Cell(rdggs, (O, 1))
            self.assertEqual(get, expect)

            c = Cell(rdggs, (N, 6, 4))
            get = c.neighbor("south_west", plane=False)
            expect = Cell(rdggs, (N, 6, 3))
            self.assertEqual(get, expect)

            c = Cell(rdggs, (P, 2, 2))
            get = c.neighbor("north", plane=False)
            expect = Cell(rdggs, (N, 8, 8))
            self.assertEqual(get, expect)

    def test_neighbors(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            # Plane test.
            # Should work on resolution 1 cells with suids of the form s0
            get = []
            for s in CELLS0:
                nb = Cell(rdggs, [s] + [0]).neighbors()
                for d in ["up", "left", "down", "right"]:
                    get.append(nb[d])
            expect = [
                (Q, 2),
                (R, 0),
                (N, 3),
                (N, 1),
                (N, 6),
                (R, 2),
                (O, 3),
                (O, 1),
                (N, 8),
                (O, 2),
                (P, 3),
                (P, 1),
                (N, 2),
                (P, 2),
                (Q, 3),
                (Q, 1),
                (N, 0),
                (Q, 2),
                (R, 3),
                (R, 1),
                (O, 6),
                (R, 8),
                (S, 3),
                (S, 1),
            ]
            expect = [Cell(rdggs, s) for s in expect]
            for i in range(len(expect)):
                self.assertEqual(get[i], expect[i])

            # Ellipsoid test.
            # Quad.
            c = Cell(rdggs, (O, 0))
            get = c.neighbors(plane=False)
            expect = {}
            expect["north"] = Cell(rdggs, (N, 6))
            expect["south"] = Cell(rdggs, (O, 3))
            expect["west"] = Cell(rdggs, (R, 2))
            expect["east"] = Cell(rdggs, (O, 1))
            for k in list(get.keys()):
                self.assertEqual(get[k], expect[k])
            # Cap.
            c = Cell(rdggs, (S, 4))
            get = c.neighbors(plane=False)
            expect = {}
            expect["north_0"] = Cell(rdggs, (S, 1))
            expect["north_1"] = Cell(rdggs, (S, 5))
            expect["north_2"] = Cell(rdggs, (S, 7))
            expect["north_3"] = Cell(rdggs, (S, 3))
            for k in list(get.keys()):
                self.assertEqual(get[k], expect[k])
            # Dart.
            c = Cell(rdggs, (N, 6))
            get = c.neighbors(plane=False)
            expect = {}
            expect["west"] = Cell(rdggs, (N, 3))
            expect["east"] = Cell(rdggs, (N, 7))
            expect["south_west"] = Cell(rdggs, (R, 2))
            expect["south_east"] = Cell(rdggs, (O, 0))
            for k in list(get.keys()):
                self.assertEqual(get[k], expect[k])
            # Skew quad.
            c = Cell(rdggs, (N, 3))
            get = c.neighbors(plane=False)
            expect = {}
            expect["north"] = Cell(rdggs, (N, 4))
            expect["south"] = Cell(rdggs, (R, 1))
            expect["west"] = Cell(rdggs, (N, 0))
            expect["east"] = Cell(rdggs, (N, 6))
            for k in list(get.keys()):
                self.assertEqual(get[k], expect[k])

    def test_neighbors_does_not_mutate_ellipsoid(self):
        # Regression test for issue #53: neighbors() used to temporarily
        # reassign self.rdggs.ellipsoid.lon_0 (a singleton shared by every
        # Cell built from this rdggs) and restore it afterwards. Confirm
        # it's never touched at all, for both the dart and skew_quad
        # branches where the mutation used to happen.
        ellipsoid = Ellipsoid(a=WGS84_ELLIPSOID.a, f=WGS84_ELLIPSOID.f, lon_0=-131.25)
        rdggs = RHEALPixDGGS(ellipsoid=ellipsoid, N_side=3)
        before = ellipsoid.lon_0
        for suid in [(N, 6), (N, 3)]:  # dart, skew_quad
            c = Cell(rdggs, suid)
            self.assertIn(c.ellipsoidal_shape, ("dart", "skew_quad"))
            c.neighbors(plane=False)
            self.assertEqual(ellipsoid.lon_0, before)

    def test_neighbors_independent_of_lon_0(self):
        # Regression test for issue #53. neighbors() labels ("east",
        # "west", "north", "south", "south_east", etc.) are documented as
        # compass directions, which are physically meaningful and must
        # not depend on where the ellipsoid's prime meridian (lon_0) is
        # drawn. The previous try/finally-free mutation-based
        # implementation got this wrong specifically when lon_0 was near
        # +-180 degrees (confirmed by testing before this fix): the same
        # physical cell's "east" and "west" neighbors could come out
        # swapped depending on lon_0 alone, with the SUID grid and all
        # other inputs unchanged.
        dart_suid = (N, 6)
        skew_quad_suid = (N, 3)
        for suid in [dart_suid, skew_quad_suid]:
            results = []
            for lon_0 in (0, 47.9, -131.25, 179.5, -179.5, 180, -180):
                ellipsoid = Ellipsoid(
                    a=WGS84_ELLIPSOID.a, f=WGS84_ELLIPSOID.f, lon_0=lon_0
                )
                rdggs = RHEALPixDGGS(ellipsoid=ellipsoid, N_side=3)
                c = Cell(rdggs, suid)
                result = {k: v.suid for k, v in c.neighbors(plane=False).items()}
                results.append(result)
            for result in results[1:]:
                self.assertEqual(result, results[0])

    def test_neighbors_thread_safety(self):
        # neighbors() no longer mutates any shared state, so concurrent
        # calls against cells sharing the same rdggs/ellipsoid should be
        # safe and should all agree with the single-threaded result.
        # (Python's GIL means this isn't a guaranteed way to catch a
        # reintroduced race, but it's a real exercise of concurrent
        # access, and combined with test_neighbors_does_not_mutate_
        # ellipsoid -- which shows there's no shared state to race on in
        # the first place -- it's a reasonable belt-and-braces check.)
        import concurrent.futures

        rdggs = WGS84_003
        suids = [(N, 6), (N, 3), (O, 0), (S, 4), (Q, 2, 5)]
        cells = [Cell(rdggs, suid) for suid in suids]
        expected = [c.neighbors(plane=False) for c in cells]

        def compute(i):
            return cells[i].neighbors(plane=False)

        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            futures = [pool.submit(compute, i % len(cells)) for i in range(200)]
            for i, future in enumerate(futures):
                result = future.result()
                self.assertEqual(result, expected[i % len(cells)])

    def test_diagonal_neighbor(self):
        rdggs = WGS84_003
        # Interior case: no face crossing.
        c = rdggs.cell((N, 4, 0))
        self.assertEqual(c.diagonal_neighbor("up_left"), rdggs.cell((N, 0, 8)))

        # Genuine cube corner: exactly 3 faces meet there (self plus its
        # up- and left-neighbors), so there's no distinct 4th diagonal
        # cell. Confirmed against real vertex geometry: N, Q, and R all
        # share this exact corner point, with no other cell reaching it.
        self.assertIsNone(rdggs.cell((N,)).diagonal_neighbor("up_left"))

        # Face-crossing case needing a rotation correction. This is the
        # case that plain composition of neighbor() calls gets wrong in
        # about a third of tested cases (self.neighbor('up').neighbor(
        # 'left') disagrees with self.neighbor('left').neighbor('up')
        # here): confirmed against real vertex geometry that ('N', 1, 0),
        # its up-neighbor, its left-neighbor, and this diagonal candidate
        # all meet at a single shared point.
        c = rdggs.cell((N, 1, 0))
        self.assertEqual(c.diagonal_neighbor("up_left"), rdggs.cell((Q, 2, 0)))

        # Invalid direction.
        with self.assertRaises(KeyError):
            c.diagonal_neighbor("north")

    def test_overlaps(self):
        rdggs = WGS84_003
        a = rdggs.cell((P, 0))
        descendant = rdggs.cell((P, 0, 3))
        sibling = rdggs.cell((P, 1))
        self.assertTrue(a.overlaps(a))
        self.assertTrue(a.overlaps(descendant))
        self.assertTrue(descendant.overlaps(a))
        self.assertFalse(a.overlaps(sibling))
        # Regression test for issue #54: this was a bare `assert`, which
        # python -O silently compiles out, and which raised
        # AssertionError rather than a conventional exception.
        empty = rdggs.cell()
        with self.assertRaises(ValueError):
            empty.overlaps(a)

    def test_equals(self):
        rdggs = WGS84_003
        a = rdggs.cell((P, 0))
        self.assertTrue(a.equals(rdggs.cell((P, 0))))
        self.assertFalse(a.equals(rdggs.cell((P, 1))))
        # Different RHEALPixDGGS instances: never equal, even with the
        # same suid, since __eq__ (which equals() delegates to) also
        # compares rdggs.
        other_rdggs = RHEALPixDGGS(N_side=4)  # A genuinely different RHEALPixDGGS.
        self.assertFalse(a.equals(Cell(other_rdggs, (P, 0))))

    def test_contains_cell_and_within(self):
        rdggs = WGS84_003
        parent = rdggs.cell((P, 0))
        child = rdggs.cell((P, 0, 3))
        sibling = rdggs.cell((P, 1))

        self.assertTrue(parent.contains_cell(child))
        self.assertTrue(parent.contains_cell(parent))  # A cell contains itself.
        self.assertFalse(child.contains_cell(parent))
        self.assertFalse(parent.contains_cell(sibling))

        self.assertTrue(child.within(parent))
        self.assertTrue(parent.within(parent))
        self.assertFalse(parent.within(child))
        self.assertFalse(sibling.within(parent))

        # covers()/covered_by() are the same relations for cells; see
        # contains_cell()'s docstring for why.
        self.assertTrue(parent.covers(child))
        self.assertTrue(child.covered_by(parent))

        empty = rdggs.cell()
        with self.assertRaises(ValueError):
            empty.contains_cell(parent)
        other_rdggs = RHEALPixDGGS(N_side=4)  # A genuinely different RHEALPixDGGS.
        with self.assertRaises(ValueError):
            parent.contains_cell(Cell(other_rdggs, (P, 0)))

    def test_touches(self):
        rdggs = WGS84_003

        # Same-resolution edge-adjacent siblings.
        self.assertTrue(rdggs.cell((P, 0)).touches(rdggs.cell((P, 1))))
        # Same-resolution corner-adjacent (diagonal-only) siblings.
        self.assertTrue(rdggs.cell((P, 0)).touches(rdggs.cell((P, 4))))
        # Same-resolution siblings that are neither: too far apart.
        self.assertFalse(rdggs.cell((P, 0)).touches(rdggs.cell((P, 8))))
        # Ancestor/descendant: touching is impossible (interiors overlap).
        self.assertFalse(rdggs.cell((P,)).touches(rdggs.cell((P, 0))))
        # The same cell: not touches (that's equals).
        self.assertFalse(rdggs.cell((P, 0)).touches(rdggs.cell((P, 0))))

        # Cross-resolution "cousins": O (a whole resolution 0 face) sits
        # immediately left of P, so O's entire shared edge touches any
        # P-descendant that never strays from P's own left-hand column.
        self.assertTrue(rdggs.cell((O,)).touches(rdggs.cell((P, 3))))
        self.assertTrue(rdggs.cell((O,)).touches(rdggs.cell((P, 6))))
        # P0 isn't in that left-hand column, so no touch.
        self.assertFalse(rdggs.cell((O,)).touches(rdggs.cell((P, 1))))
        # Symmetric regardless of argument order.
        self.assertTrue(rdggs.cell((P, 3)).touches(rdggs.cell((O,))))

        # Cross-resolution cousin that starts on the shared edge (O2 is
        # in O's right-hand column, bordering P) but then strays from it
        # one level deeper (O2's own child 0 is back in its left column):
        # no longer touches.
        self.assertTrue(rdggs.cell((O, 2)).touches(rdggs.cell((P,))))
        self.assertFalse(rdggs.cell((O, 2, 0)).touches(rdggs.cell((P,))))

        empty = rdggs.cell()
        with self.assertRaises(ValueError):
            empty.touches(rdggs.cell((P, 0)))
        other_rdggs = RHEALPixDGGS(N_side=4)  # A genuinely different RHEALPixDGGS.
        with self.assertRaises(ValueError):
            rdggs.cell((P, 0)).touches(Cell(other_rdggs, (P, 1)))

    def test_disjoint(self):
        rdggs = WGS84_003
        # Opposite poles: nowhere near each other.
        self.assertTrue(rdggs.cell((N, 0)).disjoint(rdggs.cell((S, 0))))
        # Edge-adjacent siblings: touching, so not disjoint.
        self.assertFalse(rdggs.cell((P, 0)).disjoint(rdggs.cell((P, 1))))
        # Ancestor/descendant: sharing interior, so not disjoint.
        self.assertFalse(rdggs.cell((P,)).disjoint(rdggs.cell((P, 0))))
        # Same-face siblings far enough apart to neither touch nor nest.
        self.assertTrue(rdggs.cell((P, 0)).disjoint(rdggs.cell((P, 8))))

    def test_area(self):
        rdggs = WGS84_003
        for resolution in (0, 1, 3):
            c = rdggs.cell((P,) + (0,) * resolution)
            # Planar cells are squares of side width().
            self.assertEqual(c.area(plane=True), c.width() ** 2)
            # Must agree with the DGGS-level formula it delegates to.
            self.assertEqual(
                c.area(plane=False), rdggs.cell_area(resolution, plane=False)
            )
            # Independent check: the grid is an equal-area partition of the
            # ellipsoid, so the 6 * N_side**(2r) ellipsoidal cells at any
            # resolution must sum to the surface area of the ellipsoid's
            # authalic sphere, 4*pi*R_A**2.
            total = 6 * rdggs.N_side ** (2 * resolution) * c.area(plane=False)
            sphere = 4 * pi * rdggs.ellipsoid.R_A**2
            self.assertAlmostEqual(total / sphere, 1, places=12)

    def test_color(self):
        rdggs = WGS84_003
        cells = list(rdggs.grid(1))
        colors = [c.color() for c in cells]
        # Deterministic, in-range RGB.
        for c, rgb in zip(cells, colors):
            self.assertEqual(len(rgb), 3)
            for component in rgb:
                self.assertGreaterEqual(component, 0)
                self.assertLessEqual(component, 1)
            self.assertEqual(c.color(), rgb)
        # Documented as "a unique RGB color tuple for this cell": distinct
        # cells at the same resolution get distinct colors.
        self.assertEqual(len(set(colors)), len(cells))
        # The saturation parameter is honored.
        self.assertNotEqual(cells[7].color(saturation=0.2), cells[7].color(0.9))

    def test_region_overlaps(self):
        rdggs = WGS84_003
        a = rdggs.cell((P, 0))
        descendant = rdggs.cell((P, 0, 3))
        sibling = rdggs.cell((P, 1))
        far = rdggs.cell((S, 8))
        self.assertTrue(a.region_overlaps([far, descendant]))
        self.assertTrue(a.region_overlaps([a]))
        self.assertFalse(a.region_overlaps([sibling, far]))
        self.assertFalse(a.region_overlaps([]))
        empty = rdggs.cell()
        with self.assertRaises(ValueError):
            empty.region_overlaps([a])

    def test_region(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            c = rdggs.cell((P, 0))
            self.assertEqual(c.region(), "equatorial")
            c = rdggs.cell((S, 7, 8))
            self.assertEqual(c.region(), "south_polar")
            c = rdggs.cell((N, 7, 8))
            self.assertEqual(c.region(), "north_polar")

    def test_ellipsoidal_shape(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            for s in CELLS0[1:5]:
                X = rdggs.cell([s])
                self.assertEqual(X.ellipsoidal_shape, "quad")
            for s in [CELLS0[0], CELLS0[5]]:
                for t in [[], [4], [4, 4]]:
                    u = [s] + t
                    X = rdggs.cell(u)
                    self.assertEqual(X.ellipsoidal_shape, "cap")
                for t in [[4, 0, 8], [0], [4, 4, 4, 4, 6], [2, 4]]:
                    u = [s] + t
                    X = rdggs.cell(u)
                    self.assertEqual(X.ellipsoidal_shape, "dart")
                for t in [[4, 1, 8], [1], [4, 6, 5], [2, 3]]:
                    u = [s] + t
                    X = rdggs.cell(u)
                    self.assertEqual(X.ellipsoidal_shape, "skew_quad")
        rdggs = WGS84_122
        cell_n = rdggs.cell((CELLS0[0],))
        cell_s = rdggs.cell((CELLS0[5],))
        self.assertEqual(cell_n.ellipsoidal_shape, "cap")
        self.assertEqual(cell_s.ellipsoidal_shape, "cap")

    @staticmethod
    def monte_carlo_mean_lon_lat(cell, n=4000, seed=20260811):
        """
        Estimate the mean longitude-latitude of ellipsoidal cell `cell` --
        i.e. its centroid, by the definition centroid() implements -- by
        sampling `n` points uniformly at random from the planar cell and
        projecting them onto the ellipsoid. The projection is equal-area,
        so planar-uniform samples are uniform on the ellipsoidal cell, and
        their sample mean estimates the centroid with standard error
        (sample standard deviation)/sqrt(n), independently of the
        integration centroid() itself performs. The seed is fixed so the
        estimate (and hence the test outcome) is deterministic.

        Longitudes are handled relative to the cell's nucleus meridian and
        wrapped, so a cell straddling the +-180 degree antimeridian
        doesn't produce a meaningless raw average. Returns
        (mean_lon, mean_lat, standard_error_lon, standard_error_lat).
        """
        from random import Random
        from statistics import fmean, stdev

        from rhealpixdggs.utils import wrap_longitude

        rng = Random(seed)
        rdggs = cell.rdggs
        nucleus_lon = cell.nucleus(plane=False)[0]
        ul = cell.ul_vertex()
        w = cell.width()
        lons, lats = [], []
        for _ in range(n):
            x = rng.uniform(ul[0], ul[0] + w)
            y = rng.uniform(ul[1] - w, ul[1])
            lon, lat = rdggs.rhealpix(x, y, inverse=True)
            lons.append(wrap_longitude(lon - nucleus_lon, radians=False))
            lats.append(lat)
        mean_lon = wrap_longitude(fmean(lons) + nucleus_lon, radians=False)
        return (
            mean_lon,
            fmean(lats),
            stdev(lons) / n**0.5,
            stdev(lats) / n**0.5,
        )

    def test_centroid(self):
        rdggs = WGS84_003
        # The centroid of a planar cell is its nucleus, whatever the
        # cell's ellipsoidal shape.
        for suid in [(Q, 7), (S, 2, 2), (N,), (N, 6), (N, 7)]:
            X = rdggs.cell(suid)
            self.assertEqual(X.centroid(plane=True), X.nucleus(plane=True))

        # The centroid of an ellipsoidal cap cell is its nucleus (the
        # pole), by symmetry.
        for suid in [(N, 4), (S,)]:
            X = rdggs.cell(suid)
            self.assertEqual(X.centroid(plane=False), X.nucleus(plane=False))

        # For the ellipsoidal shapes whose centroid is computed by
        # numerical integration -- dart and skew_quad -- check both
        # coordinates against an independent fixed-seed Monte Carlo
        # estimate, within 6 standard errors (a bound the estimate has
        # essentially no chance of missing unless the integration itself
        # is wrong). The dart cell at longitude -180 also exercises the
        # antimeridian-straddling case.
        from rhealpixdggs.utils import wrap_longitude

        for suid in [(N, 6), (S, 2, 4), (N, 7), (N, 7, 3)]:
            X = rdggs.cell(suid)
            self.assertIn(X.ellipsoidal_shape, ("dart", "skew_quad"))
            lon, lat = X.centroid(plane=False)
            mc_lon, mc_lat, se_lon, se_lat = self.monte_carlo_mean_lon_lat(X)
            lon_gap = abs(wrap_longitude(lon - mc_lon, radians=False))
            self.assertLess(lon_gap, 6 * se_lon, msg=str(suid))
            self.assertLess(abs(lat - mc_lat), 6 * se_lat, msg=str(suid))

    def test_centroid_quad(self):
        # Regression test for issue #75: the centroid latitude of an
        # ellipsoidal quad cell is the area-weighted mean latitude over
        # the cell, NOT the midpoint of its two edge latitudes -- latitude
        # is a nonlinear function of planar y, so the two differ, by up to
        # ~0.63 degrees for resolution 1 quads (the midpoint is what
        # centroid() used to return, following an erratum in the founding
        # paper's summary table that contradicts the paper's own integral
        # definition of the centroid).
        from scipy import integrate

        from rhealpixdggs.utils import wrap_longitude

        rdggs = WGS84_003
        for suid in [(Q, 7), (O, 0), (P, 3, 1), (Q, 4)]:
            X = rdggs.cell(suid)
            self.assertEqual(X.ellipsoidal_shape, "quad")
            lon, lat = X.centroid(plane=False)

            # The centroid longitude is the nucleus longitude (meridians
            # are equally spaced in planar x, so mean == midpoint there).
            self.assertEqual(lon, X.nucleus(plane=False)[0])

            # Deterministic ground truth for the latitude: adaptive
            # quadrature of the mean-latitude integral, evaluated along
            # the cell's left edge x (latitude is independent of x on a
            # quad, and adaptive quadrature at a different abscissa is an
            # implementation-independent path from centroid()'s own
            # fixed-order rule at the nucleus meridian).
            pv = X.vertices(plane=True)
            x1 = min(v[0] for v in pv)
            y1 = min(v[1] for v in pv)
            y2 = max(v[1] for v in pv)

            def phi(y):
                # rdggs and x1 are stable across this loop iteration, and
                # phi() is only integrated within it.
                return rdggs.rhealpix(x1, y, inverse=True)[1]  # noqa: B023

            expected_lat = (
                integrate.quad(phi, y1, y2, epsabs=1e-5, epsrel=1e-10)[0]
            ) / (y2 - y1)
            self.assertAlmostEqual(lat, expected_lat, places=9, msg=str(suid))

            # And explicitly: NOT the edge-latitude midpoint, except for
            # cells symmetric about the equator, where the two coincide.
            midpoint = (phi(y1) + phi(y2)) / 2
            if abs(midpoint) > 1e-12:
                self.assertNotAlmostEqual(lat, midpoint, places=3, msg=str(suid))

        # Fully independent statistical check for the deepest-affected
        # case: a fixed-seed Monte Carlo mean with enough samples that its
        # 6-standard-error bound (~0.35 degrees) is tighter than the
        # ~0.63 degree error this test guards against.
        X = rdggs.cell((Q, 7))
        lon, lat = X.centroid(plane=False)
        mc_lon, mc_lat, se_lon, se_lat = self.monte_carlo_mean_lon_lat(X, n=20000)
        self.assertLess(abs(wrap_longitude(lon - mc_lon, radians=False)), 6 * se_lon)
        self.assertLess(abs(lat - mc_lat), 6 * se_lat)

    def test_random_point(self):
        # Output should lie in the cell at least.
        for E in [WGS84_ASPHERE_RADIANS, WGS84_ELLIPSOID]:
            rdggs = RHEALPixDGGS(E)
            for plane in [True, False]:
                c = Cell(rdggs, [N, 8, 7])
                p = c.random_point(plane=plane)
                self.assertTrue(
                    c.contains(p, plane=plane)
                )  # ------------------------------------------------------------------------------
