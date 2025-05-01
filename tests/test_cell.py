# Import third-party modules.
from scipy.spatial.distance import euclidean, norm

# Import standard modules
import unittest
from itertools import product

# Import my modules.
from rhealpixdggs.cell import Cell, CELLS0

from rhealpixdggs.dggs import RHEALPixDGGS, WGS84_003, WGS84_003_RADIANS
from rhealpixdggs.ellipsoids import (
    WGS84_ASPHERE_RADIANS,
    WGS84_ELLIPSOID,
    WGS84_ELLIPSOID_RADIANS,
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
            self.assertRaises(AssertionError, Cell, rdggs, suid)

            # Should create cell P1.
            expect = (P, 1)
            i = 25
            get = Cell(rdggs, level_order_index=i).suid
            self.assertEqual(get, expect)

            expect = (P, 1)

            def num(k):
                return rdggs.num_cells(res_1=k, subcells=True)

            i = 2 * num(0) + 1 * num(1) + num(1) - 1
            get = Cell(rdggs, post_order_index=i).suid
            self.assertEqual(get, expect)

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
            self.assertFalse(a < a)
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

    def test_gt(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            a = Cell(rdggs, (N, 7, 6, 8, 1))
            b = Cell(rdggs, (O, 1, 2, 3))
            c = Cell(rdggs, (O, 1, 2))
            self.assertFalse(a > b)
            self.assertFalse(a > c)
            self.assertFalse(a > a)
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
            M = rdggs.num_cells(0, rdggs.max_resolution)
            k = 3616048  # randint(0, M - 1)
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
                w = c.width()
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

    def test_nucleus(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            # Nuclei of children should be in correct position
            # relative to parent cell in rHEALPix projection.
            a = Cell(rdggs, (S, 7, 4, 1, 2, 1))  # Arbitrary cell.
            w = a.width()
            (x, y) = a.ul_vertex()
            error = 1e-10
            for row, col in product(list(range(3)), repeat=2):
                s = str(row * 3 + col)
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
            expect = dict()
            expect["north"] = Cell(rdggs, (N, 6))
            expect["south"] = Cell(rdggs, (O, 3))
            expect["west"] = Cell(rdggs, (R, 2))
            expect["east"] = Cell(rdggs, (O, 1))
            for k in list(get.keys()):
                self.assertEqual(get[k], expect[k])
            # Cap.
            c = Cell(rdggs, (S, 4))
            get = c.neighbors(plane=False)
            expect = dict()
            expect["north_0"] = Cell(rdggs, (S, 1))
            expect["north_1"] = Cell(rdggs, (S, 5))
            expect["north_2"] = Cell(rdggs, (S, 7))
            expect["north_3"] = Cell(rdggs, (S, 3))
            for k in list(get.keys()):
                self.assertEqual(get[k], expect[k])
            # Dart.
            c = Cell(rdggs, (N, 6))
            get = c.neighbors(plane=False)
            expect = dict()
            expect["west"] = Cell(rdggs, (N, 3))
            expect["east"] = Cell(rdggs, (N, 7))
            expect["south_west"] = Cell(rdggs, (R, 2))
            expect["south_east"] = Cell(rdggs, (O, 0))
            for k in list(get.keys()):
                self.assertEqual(get[k], expect[k])
            # Skew quad.
            c = Cell(rdggs, (N, 3))
            get = c.neighbors(plane=False)
            expect = dict()
            expect["north"] = Cell(rdggs, (N, 4))
            expect["south"] = Cell(rdggs, (R, 1))
            expect["west"] = Cell(rdggs, (N, 0))
            expect["east"] = Cell(rdggs, (N, 6))
            for k in list(get.keys()):
                self.assertEqual(get[k], expect[k])

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
                self.assertEqual(X.ellipsoidal_shape(), "quad")
            for s in [CELLS0[0], CELLS0[5]]:
                for t in [[], [4], [4, 4]]:
                    u = [s] + t
                    X = rdggs.cell(u)
                    self.assertEqual(X.ellipsoidal_shape(), "cap")
                for t in [[4, 0, 8], [0], [4, 4, 4, 4, 6], [2, 4]]:
                    u = [s] + t
                    X = rdggs.cell(u)
                    self.assertEqual(X.ellipsoidal_shape(), "dart")
                for t in [[4, 1, 8], [1], [4, 6, 5], [2, 3]]:
                    u = [s] + t
                    X = rdggs.cell(u)
                    self.assertEqual(X.ellipsoidal_shape(), "skew_quad")
        rdggs = WGS84_122
        cell_n = rdggs.cell((CELLS0[0],))
        cell_s = rdggs.cell((CELLS0[5],))
        self.assertEqual(cell_n.ellipsoidal_shape(), "cap")
        self.assertEqual(cell_s.ellipsoidal_shape(), "cap")

    def test_centroid(self):
        # Warning: This test is slow.
        # Uncomment below if you want to test it and wait.
        pass
        # print
        # print 'Testing centroid() method now. Takes about 2 minutes.'
        #
        # # For non-cap ellipsoidal cells, test centroid() against a Monte Carlo
        # # approximation of the centroid.
        # def monte_carlo_centroid(cell):
        #     nv = cell.nucleus_and_vertices()
        #     rdggs = cell.rdggs
        #     lam_nucleus = rdggs.rhealpix(*nv[0], inverse=True)[0]
        #     vertices = nv[1:]
        #     x1, x2 = vertices[0][0], vertices[3][0]
        #     y1, y2 = vertices[1][1], vertices[0][1]
        #     N = 10000
        #     sample_points = []
        #     for i in range(N):
        #         x, y = uniform(x1, x2), uniform(y1, y2)
        #         lam, phi = rdggs.rhealpix(x, y, inverse=True)
        #         sample_points.append(array((lam, phi)))
        #     lam_bar, phi_bar = sum(sample_points)/N
        #     sample_var = sum([array((
        #                              (p[0] - lam_bar)**2,
        #                              (p[1] - phi_bar)**2))
        #                       for p in sample_points])/(N - 1)
        #     lam_bar_err = sqrt(sample_var[0]/N) # Approximately
        #     phi_bar_err = sqrt(sample_var[1]/N) # Approximately
        #     PI = cell.rdggs.ellipsoid.pi()
        #     if cell.ellipsoidal_shape() == 'dart':
        #         lam_bar = lam_nucleus
        #     return lam_bar, phi_bar, abs(lam_bar_err), abs(phi_bar_err)
        #
        # for rdggs in [WGS84_003, WGS84_003_RADIANS]:
        #     # The centroid of a planar cell is its nucleus.
        #     for suid in [(Q, 7), (S, 2, 2)]:
        #         X = rdggs.cell(suid)
        #         centroid = X.centroid()
        #         nucleus = X.nucleus_and_vertices()[0]
        #         self.assertEqual(centroid, nucleus)
        #
        #     # The centroid of a ellipsoidal cap cell is also its nucleus.
        #     for suid in [(N, 4, 4, 4), [S]]:
        #         X = rdggs.cell(suid)
        #         centroid = X.centroid(plane=False)
        #         nucleus = X.nucleus_and_vertices(plane=False)[0]
        #         self.assertEqual(centroid, nucleus)
        #
        #     for suid in [
        #       (Q, 7), # quad
        #       (O, 5, 8), # quad
        #       (N, 6), # dart
        #       (N, 6, 2), # dart
        #       #(N, 6, 2, 4), # dart
        #       (N, 7), # skew quad
        #       (N, 7, 3), # skew quad
        #       #(N, 7, 3, 5), # skew quad
        #       (S, 2), # dart
        #       (S, 2, 4), # dart
        #       #(S, 2, 4, 4), # dart
        #       (S, 1), # skew quad
        #       (S, 1, 1), # skew quad
        #       #(S, 1, 1, 1), # skew quad
        #       ]:
        #         X = rdggs.cell(suid)
        #         lam_bar, phi_bar = X.centroid(plane=False)
        #         lam_bar_approx, phi_bar_approx, lam_bar_err, phi_bar_err = \
        #         monte_carlo_centroid(X)
        #         # print "Testing centroid(plane=False) for %s cell %s..."\
        #         #  % (X.ellipsoidal_shape(), X)
        #         # print 'lam:', lam_bar, lam_bar_approx, lam_bar_err
        #         # print 'phi:', phi_bar, phi_bar_approx, phi_bar_err
        #         self.assertTrue(euclidean(lam_bar, lam_bar_approx) <\
        #                                                   10*lam_bar_err)
        #         self.assertTrue(euclidean(phi_bar, phi_bar_approx) <\
        #                                                   10*phi_bar_err)

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
