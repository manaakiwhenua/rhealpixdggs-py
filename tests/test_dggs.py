"""
This Python 3.11 code tests the ``dggs_rhealpix`` module.
Beware, while these tests cover most functions and scenarios, they don't cover them all.
Keep adding tests!

- Alexander Raichev (AR), 2013-01-26: Initial version based on previous tests.

"""

# *****************************************************************************
#       Copyright (C) 2011 Alexander Raichev <tortoise.said@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import third-party modules.
from numpy import array, pi, sin, cos, tan, arctan, arcsin
from scipy.spatial.distance import euclidean, norm

# Import standard modules
import unittest
from itertools import product
from random import randint, uniform

# Import my modules.
import rhealpixdggs.dggs as dggs
from rhealpixdggs.cell import Cell, CELLS0
from rhealpixdggs.dggs import RHEALPixDGGS, WGS84_003, WGS84_003_RADIANS
from rhealpixdggs.ellipsoids import (
    Ellipsoid,
    WGS84_ASPHERE_RADIANS,
    WGS84_ELLIPSOID,
    WGS84_ELLIPSOID_RADIANS,
)


# Relative error function.
def rel_err(get, expect):
    a = euclidean(get, expect)
    b = norm(expect)
    if b == 0:
        return a
    else:
        return a / b


# Level 0 cell names
N = CELLS0[0]
O = CELLS0[1]
P = CELLS0[2]
Q = CELLS0[3]
R = CELLS0[4]
S = CELLS0[5]

# Some DGGSs used below.
WGS84_123 = dggs.RHEALPixDGGS(
    ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=2, N_side=3
)
WGS84_123_RADIANS = dggs.RHEALPixDGGS(
    ellipsoid=WGS84_ELLIPSOID_RADIANS, north_square=1, south_square=2, N_side=3
)
WGS84_122 = dggs.RHEALPixDGGS(
            ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=2, N_side=2
        )


class SCENZGridTestCase(unittest.TestCase):
    # ------------------------------------------------------------------------------
    # Test Cell methods
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

    # Test RHEALPixDGGS methods.
    def test_RHEALPixDGGS_init(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            # Atomic neighbors of integers 0,..., 8 should be consistent.
            # For example, the down neighbor of the up neighbor of k
            # should be k.
            an = rdggs.atomic_neighbors
            n = rdggs.N_side**2
            for k in range(n):
                self.assertEqual(an[an[k]["up"]]["down"], k)
                self.assertEqual(an[an[k]["right"]]["left"], k)

    def test_grid(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            # Should output the correct number of cells.
            resolution = 3
            cells = list(rdggs.grid(resolution))
            n = rdggs.N_side**2
            self.assertEqual(len(cells), 6 * n**resolution)

            # Pick a random cell of index k in cells and check that cells[k+1]
            # is the successor cell.
            k = randint(0, 6 * n**resolution - 2)
            a = cells[k]
            b = cells[k + 1]
            c = a.successor()
            self.assertEqual(str(b), str(c))

    def test_interval(self):
        for rdggs in [WGS84_123, WGS84_123_RADIANS]:
            # Should produce the correct number of cells
            # and each cell should have level order index one greater
            # than the previous.
            count = 0
            A = rdggs.cell((N, 1))
            B = rdggs.cell((S, 2, 8))
            start_index = rdggs.cell((N, 2, 0)).index(order="level")
            end_index = rdggs.cell((S, 2, 8)).index(order="level")
            correct_count = end_index - start_index + 1
            count = 0
            old_index = start_index - 1
            for X in rdggs.interval(A, B):
                new_index = X.index(order="level")
                self.assertEqual(new_index, old_index + 1)
                old_index = new_index
                count += 1
            self.assertEqual(count, correct_count)

            # Should produce the correct number of cells
            # and each cell should have level order index one greater
            # than the previous.
            count = 0
            A = rdggs.cell((N, 2, 0))
            B = rdggs.cell((S, 2))
            start_index = rdggs.cell((N, 2, 0)).index(order="level")
            end_index = rdggs.cell((S, 2, 8)).index(order="level")
            correct_count = end_index - start_index + 1
            count = 0
            old_index = start_index - 1
            for X in rdggs.interval(A, B):
                new_index = X.index(order="level")
                self.assertEqual(new_index, old_index + 1)
                old_index = new_index
                count += 1
            self.assertEqual(count, correct_count)

    def test_cell_from_point(self):
        # The nucleus of a cell should yield the cell.
        for plane in [True, False]:
            rdggs = WGS84_123
            for x in rdggs.grid(1):
                c = x.nucleus(plane=plane)
                y = rdggs.cell_from_point(x.resolution, c, plane=plane)
                self.assertEqual(x.suid, y.suid)

        # Test an edge case.
        rdggs = WGS84_123
        R_A = rdggs.ellipsoid.R_A
        p = R_A * array((-0.1, -pi / 4))
        get = rdggs.cell_from_point(1, p)
        expect = rdggs.cell((P, 8))
        self.assertEqual(get, expect)

        # Test with origin translation.
        p = (0.1, 0.7)
        E = Ellipsoid(lon_0=p[0], lat_0=p[1])
        rdggs = RHEALPixDGGS(E)
        get = rdggs.cell_from_point(1, p, plane=False)
        expect = rdggs.cell_from_point(1, (0, 0))
        self.assertEqual(get, expect)

        # Test with origin translation and degrees mode.
        p = (-50, 20)
        E = Ellipsoid(lon_0=p[0], lat_0=p[1], radians=False)
        rdggs = RHEALPixDGGS(E)
        get = rdggs.cell_from_point(1, p, plane=False)
        expect = rdggs.cell_from_point(1, (0, 0))
        self.assertEqual(get, expect)

        # Test extreme value (point outside the planar DGGS)
        p = (11500249, 56898969)
        E = Ellipsoid(lon_0=0, lat_0=0, radians=False)
        rdggs = RHEALPixDGGS(E)
        get = rdggs.cell_from_point(0, p, plane=True)
        self.assertIsNone(get)

        # Assert inequality of non-coincident points
        p1 = (0, 0)
        p2 = (p1[0]+45,p1[1]+360)
        E = Ellipsoid(lon_0=0, lat_0=0, radians=False)
        rdggs = RHEALPixDGGS(E)
        c1 = rdggs.cell_from_point(1, p1, plane=False)
        c2 = rdggs.cell_from_point(1, p2, plane=False)
        self.assertNotEqual(c1, c2)

        # Assert equality of coincident points
        # (Given a comparison point with >90 latitude, >180 longitude; i.e. assert that lat/lon wrap)
        p1 = (0, 0)
        p2 = (p1[0]+360,p1[1]+360)
        E = Ellipsoid(lon_0=0, lat_0=0, radians=False)
        rdggs = RHEALPixDGGS(E)
        c1 = rdggs.cell_from_point(1, p1, plane=False)
        c2 = rdggs.cell_from_point(1, p2, plane=False)
        self.assertEqual(c1, c2)

    def test_cell_from_region(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            # For any planar cell X with nucleus c and width w,
            # a square centered at c with width w/2 should yield X.
            X = rdggs.cell((S, 4, 7))
            c = X.nucleus()
            w = X.width()
            ul = (c[0] - w / 4, c[1] + w / 4)
            dr = (c[0] + w / 4, c[1] - w / 4)
            Y = rdggs.cell_from_region(ul, dr)
            self.assertEqual(X.suid, Y.suid)

            # The nuclei of two cells at the same resolution that have
            # the same parent cell should yield the parent cell.
            A = rdggs.cell((O, 3, 0, 1))
            B = rdggs.cell((O, 3, 0, 5))
            expect = rdggs.cell((O, 3, 0))
            Ac = A.nucleus()
            Bc = B.nucleus()
            get = rdggs.cell_from_region(Ac, Bc)
            self.assertEqual(get, expect)

            # Should work on defunct cases.
            # Rectangle with only two vertices.
            A = rdggs.cell((O, 3, 0, 1))
            B = rdggs.cell((O, 3, 0, 2))
            expect = rdggs.cell((O, 3, 0))
            Ac = A.nucleus()
            Bc = B.nucleus()
            get = rdggs.cell_from_region(Ac, Bc)
            self.assertEqual(get, expect)
            # Rectangle with only one vertex.
            A = rdggs.cell((O, 3, 0, 1))
            B = rdggs.cell((O, 3, 0, 1))
            expect = [O, 3, 0, 1] + [4 for i in range(rdggs.max_resolution - 3)]
            expect = rdggs.cell(expect)
            Ac = A.nucleus()
            Bc = B.nucleus()
            get = rdggs.cell_from_region(Ac, Bc)
            self.assertEqual(get, expect)

            # Test ellipsoidal regions.
            # An ellipsoidal cap contained in S4 but not in S44.
            PI = rdggs.ellipsoid.pi()
            X = rdggs.cell((S, 4))
            vertices = X.vertices(plane=False)
            phi = vertices[0][1] - PI / 100
            ul = (-PI, phi)
            dr = (-PI, -PI / 2)
            get = rdggs.cell_from_region(ul, dr, plane=False)
            expect = rdggs.cell((S, 4))
            self.assertEqual(get, expect)

            # An ellipsoidal quadrangle sharing two vertices with an
            # ellipsoidal cell.
            X = rdggs.cell((N, 7))  # A skew quad cell
            vertices = X.vertices(plane=False)
            fuzz = 1e-3  # To avoid rounding errors
            ul = vertices[3][0] + fuzz, vertices[0][1] - fuzz
            dr = vertices[2][0] - fuzz, vertices[2][1] + fuzz
            get = rdggs.cell_from_region(ul, dr, plane=False)
            expect = X
            self.assertEqual(get, expect)

    def test_cell_latitudes(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            # Test nuclei latitudes.
            R_A = rdggs.ellipsoid.R_A
            resolution = 1
            get = rdggs.cell_latitudes(
                resolution, -R_A * pi / 2, R_A * pi / 2, True, plane=True
            )
            expect = R_A * array(
                [
                    -pi / 3,
                    -pi / 6,
                    0,
                    pi / 6,
                    pi / 3,
                ]
            )
            for i in range(len(expect)):
                self.assertAlmostEqual(get[i], expect[i])

            get = rdggs.cell_latitudes(
                resolution, -R_A * pi / 3, R_A * pi / 3, True, plane=True
            )
            expect = R_A * array([-pi / 6, 0, pi / 6])
            for i in range(len(expect)):
                self.assertAlmostEqual(get[i], expect[i])

            # Test boundary latitudes:
            resolution = 1
            get = rdggs.cell_latitudes(
                resolution, -R_A * pi / 2, R_A * pi / 2, False, plane=True
            )
            expect = R_A * array(
                [-5 * pi / 12, -pi / 4, -pi / 12, pi / 12, pi / 4, 5 * pi / 12]
            )
            for i in range(len(get)):
                self.assertAlmostEqual(get[i], expect[i])

            get = rdggs.cell_latitudes(
                resolution, -R_A * pi / 4, R_A * 5 * pi / 12, False, plane=True
            )
            expect = R_A * array([-pi / 4, -pi / 12, pi / 12, pi / 4])
            for i in range(len(expect)):
                self.assertAlmostEqual(get[i], expect[i])

    def test_cells_from_meridian(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            PI = rdggs.ellipsoid.pi()
            lam = -PI / 2 + 0.2
            phi_min = -PI / 2 + 0.1
            phi_max = PI / 2 - 0.1
            get = rdggs.cells_from_meridian(1, lam, phi_min, phi_max)
            expect = [
                rdggs.cell([N, 4]),
                rdggs.cell([N, 8]),
                rdggs.cell([N, 5]),
                rdggs.cell([P, 0]),
                rdggs.cell([P, 3]),
                rdggs.cell([P, 6]),
                rdggs.cell([S, 2]),
                rdggs.cell([S, 5]),
                rdggs.cell([S, 4]),
            ]
            self.assertEqual(get, expect)

    def test_cells_from_parallel(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            PI = rdggs.ellipsoid.pi()
            phi = PI / 3
            lam_min = -PI
            lam_max = PI
            end_points = [(lam_min, phi), (lam_max, phi)]
            get = rdggs.cells_from_parallel(1, phi, lam_min, lam_max)
            expect = [
                rdggs.cell([N, 6]),
                rdggs.cell([N, 7]),
                rdggs.cell([N, 8]),
                rdggs.cell([N, 5]),
                rdggs.cell([N, 2]),
                rdggs.cell([N, 1]),
                rdggs.cell([N, 0]),
                rdggs.cell([N, 3]),
            ]
            self.assertEqual(get, expect)

            phi = PI / 3
            lam_min = -PI
            lam_max = -PI + 0.1
            end_points = [(lam_min, phi), (lam_max, phi)]
            get = rdggs.cells_from_parallel(1, phi, lam_min, lam_max)
            expect = [rdggs.cell([N, 6])]
            self.assertEqual(get, expect)

            phi = PI / 3
            lam_min = -PI
            lam_max = 0
            end_points = [(lam_min, phi), (lam_max, phi)]
            get = rdggs.cells_from_parallel(1, phi, lam_min, lam_max)
            expect = [
                rdggs.cell([N, 6]),
                rdggs.cell([N, 7]),
                rdggs.cell([N, 8]),
                rdggs.cell([N, 5]),
                rdggs.cell([N, 2]),
            ]
            self.assertEqual(get, expect)

    def test_cells_from_region(self):
        for rdggs in [WGS84_003, WGS84_003_RADIANS]:
            PI = rdggs.ellipsoid.pi()
            # Rectangle.
            R_A = rdggs.ellipsoid.R_A
            ul = R_A * array((-0.1, pi / 4))
            dr = R_A * array((0.1, -pi / 4))
            get = rdggs.cells_from_region(1, ul, dr, plane=True)
            expect = [
                [rdggs.cell((P, 2)), rdggs.cell((Q, 0))],
                [rdggs.cell((P, 5)), rdggs.cell((Q, 3))],
                [rdggs.cell((P, 8)), rdggs.cell((Q, 6))],
            ]
            self.assertEqual(get, expect)

            # Quad.
            ul = (0, PI / 3)
            dr = (PI / 2, PI / 6)
            get = rdggs.cells_from_region(1, ul, dr, plane=False)
            expect = [
                [rdggs.cell((N, 2)), rdggs.cell((N, 1)), rdggs.cell((N, 0))],
                [
                    rdggs.cell((Q, 0)),
                    rdggs.cell((Q, 1)),
                    rdggs.cell((Q, 2)),
                    rdggs.cell((R, 0)),
                ],
            ]

            # Defunct quad.
            ul = (0, -PI / 6)
            dr = (PI / 2, -PI / 2)
            get = rdggs.cells_from_region(1, ul, dr, plane=False)
            expect = [
                [
                    rdggs.cell((Q, 6)),
                    rdggs.cell((Q, 7)),
                    rdggs.cell((Q, 8)),
                    rdggs.cell((R, 6)),
                ],
                [rdggs.cell((S, 8)), rdggs.cell((S, 7)), rdggs.cell((S, 6))],
                [rdggs.cell((S, 4))],
            ]

            # Cap.
            ul = (-PI, PI / 2)
            dr = (-PI, PI / 7)
            get = rdggs.cells_from_region(1, ul, dr, plane=False)
            expect = [
                [rdggs.cell((N, 4))],
                [
                    rdggs.cell((N, 6)),
                    rdggs.cell((N, 7)),
                    rdggs.cell((N, 8)),
                    rdggs.cell((N, 5)),
                    rdggs.cell((N, 2)),
                    rdggs.cell((N, 1)),
                    rdggs.cell((N, 0)),
                    rdggs.cell((N, 3)),
                ],
                [
                    rdggs.cell((O, 0)),
                    rdggs.cell((O, 1)),
                    rdggs.cell((O, 2)),
                    rdggs.cell((P, 0)),
                    rdggs.cell((P, 1)),
                    rdggs.cell((P, 2)),
                    rdggs.cell((Q, 0)),
                    rdggs.cell((Q, 1)),
                    rdggs.cell((Q, 2)),
                    rdggs.cell((R, 0)),
                    rdggs.cell((R, 1)),
                    rdggs.cell((R, 2)),
                ],
            ]
            self.assertEqual(get, expect)

            # Cap.
            ul = (-PI, -PI / 3)
            dr = (-PI, -PI / 2)
            get = rdggs.cells_from_region(1, ul, dr, plane=False)
            expect = [
                [
                    rdggs.cell((S, 0)),
                    rdggs.cell((S, 1)),
                    rdggs.cell((S, 2)),
                    rdggs.cell((S, 5)),
                    rdggs.cell((S, 8)),
                    rdggs.cell((S, 7)),
                    rdggs.cell((S, 6)),
                    rdggs.cell((S, 3)),
                ],
                [rdggs.cell((S, 4))],
            ]
            self.assertEqual(get, expect)


# ------------------------------------------------------------------------------
#     # Test CellFamily methods
#     def test_CellFamily_init(self):
#         X = WGS84_123.cell((S, 0))
#         Y = WGS84_123.cell((P, 3, 3))
#         Z = WGS84_123.cell('S')
#         F = CellFamily([X, Y, X, Z, Y])
#         cells = [Y, Z]
#         min_resolution = 0
#         max_resolution = 2
#         self.assertEqual(len(F.cells), len(cells))
#         for i in range(len(F.cells)):
#             self.assertEqual(F.cells[i].suid, cells[i].suid)
#         self.assertEqual(F.min_resolution, min_resolution)
#         self.assertEqual(F.max_resolution, max_resolution)
#
#     def test_CellFamily_contains(self):
#         f = CellFamily(list(WGS84_123.grid(1)))
#         c = WGS84_123.cell((O, 5))
#         d = WGS84_123.cell((O, 5, 5))
#         self.assertTrue(c in f)
#         self.assertTrue(d not in f)
#
#     def test_CellFamily_delitem(self):
#         f = CellFamily(list(WGS84_123.grid(1)))
#         c = WGS84_123.cell((N, 0))
#         del(f[0])
#         self.assertTrue(c not in f)
#
#     def test_CellFamily_add(self):
#         f = CellFamily(list(WGS84_123.grid(0)))
#         del f[2]    # Remove cell C.
#
#         # Should add cell C7 to f in correct spot.
#         c = WGS84_123.cell((P, 7))
#         f.add(c)
#         self.assertTrue(str(f[2]), str(c))
#
#         # Should not add cell B7 to f since it's a subcell of a cell in f.
#         d = WGS84_123.cell((O, 7))
#         f.add(d)
#         self.assertTrue(d not in f)
#
#     def test_CellFamily_minimize(self):
#         # The minimizeion of the empty cell family should be itself.
#         f = CellFamily()
#         fc = f.minimize()
#         self.assertTrue(f is fc)
#         self.assertEqual(fc.min_resolution, None)
#         self.assertEqual(fc.max_resolution, None)
#
#         # Minimize.
#         A00 = WGS84_123.cell((S, 0, 0))
#         B = WGS84_123.cell('O')
#         C00 = WGS84_123.cell((P, 0, 0))
#         C02 = WGS84_123.cell((P, 0, 2))
#         f = CellFamily([A00, C00, C02] + list(B.subcells(3)))
#         fc = f.minimize()
#         g = CellFamily([A00, B, C00, C02])
#         self.assertTrue(fc == g)
#         self.assertEqual(fc.min_resolution, 0)
#         self.assertEqual(fc.max_resolution, 2)
#
#         # This family can't be minimizeed.
#         cells = list(B.subcells(1))
#         cells.pop(7)
#         f = CellFamily(cells)
#         fc = f.minimize()
#         self.assertTrue(f == fc)
#         self.assertEqual(fc.min_resolution, 1)
#         self.assertEqual(fc.max_resolution, 1)
#
#     def test_CellFamily_union(self):
#         # Union with the empty cell family.
#         f = CellFamily()
#         g = CellFamily(list(WGS84_123.grid(0)))
#         self.assertTrue(f.union(g) == g)
#
#         # Union.
#         p17 = WGS84_123.cell((P, 1, 7))
#         p23 = WGS84_123.cell((P, 2, 3))
#         p345 = WGS84_123.cell((P, 3, 4, 5))
#         p346 = WGS84_123.cell((P, 3, 4, 6))
#         p381 = WGS84_123.cell((P, 3, 8, 1))
#         p2 = WGS84_123.cell((P, 2))
#         p3 = WGS84_123.cell((P, 3))
#         f1 = CellFamily([p17, p23, p345, p346, p381])
#         f2 = CellFamily([p2, p3])
#         u = f1.union(f2)
#         uu = CellFamily([p17, p2, p3])
#         self.assertTrue(u == uu)
#
#         # Union with filtering.
#         u = f1.union(f2, filter_resolution=1)
#         uu = CellFamily([p2, p3])
#         self.assertTrue(u == uu)
#
#         u = f1.union(f2, filter_resolution=0)
#         uu = CellFamily()
#         self.assertTrue(u == uu)
#
#     def test_CellFamily_union_all(self):
#         # Union.
#         p17 = WGS84_123.cell((P, 1, 7))
#         p23 = WGS84_123.cell((P, 2, 3))
#         p345 = WGS84_123.cell((P, 3, 4, 5))
#         p346 = WGS84_123.cell((P, 3, 4, 6))
#         p381 = WGS84_123.cell((P, 3, 8, 1))
#         p2 = WGS84_123.cell((P, 2))
#         p3 = WGS84_123.cell((P, 3))
#         p34 = WGS84_123.cell((P, 3, 4))
#         p38 = WGS84_123.cell((P, 3, 8))
#         f1 = CellFamily([p17, p23, p345, p346, p381])
#         f2 = CellFamily([p2, p3])
#         f3 = CellFamily([p2, p34, p38])
#         u = f1.union_all([f2, f3])
#         uu = CellFamily([p17, p2, p3])
#         self.assertTrue(u == uu)
#
#         # Union with filtering.
#         u = f1.union_all([f2, f3], filter_resolution=1)
#         uu = CellFamily([p2, p3])
#         self.assertTrue(u == uu)
#
#         u = f1.union_all([f2, f3], filter_resolution=0)
#         uu = CellFamily()
#         self.assertTrue(u == uu)
#
#     def test_CellFamily_intersect(self):
#         # Intersect with an empty cell family.
#         f = CellFamily()
#         g = CellFamily(list(WGS84_123.grid(0)))
#         self.assertTrue(f.intersect(g) == f)
#
#         # Intersect.
#         p17 = WGS84_123.cell((P, 1, 7))
#         p23 = WGS84_123.cell((P, 2, 3))
#         p345 = WGS84_123.cell((P, 3, 4, 5))
#         p346 = WGS84_123.cell((P, 3, 4, 6))
#         p381 = WGS84_123.cell((P, 3, 8, 1))
#         p2 = WGS84_123.cell((P, 2))
#         p3 = WGS84_123.cell((P, 3))
#         f1 = CellFamily([p17, p23, p345, p346, p381])
#         f2 = CellFamily([p2, p3])
#         v = f1.intersect(f2)
#         vv = CellFamily([p23, p345, p346, p381])
#         self.assertTrue(v == vv)
#
#         # Intersect with filtering.
#         v = f1.intersect(f2, filter_resolution=2)
#         vv = CellFamily([p23])
#         self.assertTrue(v == vv)
#
#         v = f1.intersect(f2, filter_resolution=1)
#         vv = CellFamily()
#         self.assertTrue(v == vv)
#
#     def test_CellFamily_intersect_all(self):
#         # Intersect.
#          p17 = WGS84_123.cell((P, 1, 7))
#          p23 = WGS84_123.cell((P, 2, 3))
#          p345 = WGS84_123.cell((P, 3, 4, 5))
#          p346 = WGS84_123.cell((P, 3, 4, 6))
#          p381 = WGS84_123.cell((P, 3, 8, 1))
#          p2 = WGS84_123.cell((P, 2))
#          p3 = WGS84_123.cell((P, 3))
#          p34 = WGS84_123.cell((P, 3, 4))
#          p38 = WGS84_123.cell((P, 3, 8))
#          f1 = CellFamily([p17, p23, p345, p346, p381])
#          f2 = CellFamily([p2, p3])
#          f3 = CellFamily([p2, p34, p38])
#          v = f1.intersect_all([f2, f3])
#          vv = CellFamily([p23, p345, p346, p381])
#          self.assertTrue(v == vv)
#
#          # Intersect with filtering.
#          v = f1.intersect_all([f2, f3], filter_resolution=2)
#          vv = CellFamily([p23])
#          self.assertTrue(v == vv)
#
#          v = f1.intersect_all([f2, f3], filter_resolution=1)
#          vv = CellFamily()
#          self.assertTrue(v == vv)
#
# #------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
