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
from numpy import array, pi

# Import standard modules
import unittest
from random import randint  # , uniform

# Import my modules.
# import rhealpixdggs.dggs as dggs
from rhealpixdggs.cell import CELLS0
from rhealpixdggs.dggs import RHEALPixDGGS, WGS84_003, WGS84_003_RADIANS
from rhealpixdggs.ellipsoids import (
    Ellipsoid,
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


class SCENZGridRHEALPixDGGSTestCase(unittest.TestCase):
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
        p2 = (p1[0] + 45, p1[1] + 360)
        E = Ellipsoid(lon_0=0, lat_0=0, radians=False)
        rdggs = RHEALPixDGGS(E)
        c1 = rdggs.cell_from_point(1, p1, plane=False)
        c2 = rdggs.cell_from_point(1, p2, plane=False)
        self.assertNotEqual(c1, c2)

        # Assert equality of coincident points
        # (Given a comparison point with >90 latitude, >180 longitude; i.e. assert that lat/lon wrap)
        p1 = (0, 0)
        p2 = (p1[0] + 360, p1[1] + 360)
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
