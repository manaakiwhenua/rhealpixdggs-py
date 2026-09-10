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

import itertools
import unittest
from random import randint  # , uniform

import numpy as np
from numpy import array, pi

# import rhealpixdggs.dggs as dggs
from rhealpixdggs.cell import CELLS0
from rhealpixdggs.dggs import WGS84_003, WGS84_003_RADIANS, RHEALPixDGGS
from rhealpixdggs.ellipsoids import (
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

    def test_area_error_budget(self):
        import sys

        rdggs = WGS84_003
        budget = rdggs.area_error_budget()
        # One entry per resolution, none missing, none extra.
        self.assertEqual(sorted(budget.keys()), list(range(rdggs.max_resolution + 1)))
        rel_tol = 10 * sys.float_info.epsilon
        previous_area = None
        for r in range(rdggs.max_resolution + 1):
            entry = budget[r]
            # The budgeted area is exactly the ellipsoidal cell area at
            # that resolution.
            self.assertEqual(entry["cell_area_m2"], rdggs.cell_area(r, plane=False))
            # Documented invariants: the relative tolerance is 10 machine
            # epsilons at every resolution, and the absolute tolerance is
            # the area scaled by it.
            self.assertEqual(entry["rel_tolerance"], rel_tol)
            self.assertEqual(entry["abs_tolerance"], entry["cell_area_m2"] * rel_tol)
            # Areas strictly decrease with resolution.
            if previous_area is not None:
                self.assertLess(entry["cell_area_m2"], previous_area)
            previous_area = entry["cell_area_m2"]

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

    def test_interval_N_side_ge_4(self):
        # Regression test for issue #71: interval() walks cells via
        # cell.successor() and terminates on `cell <= b`, and __le__ used
        # to compare a string rendering of the suid, which breaks for
        # N_side >= 4 (some digits, e.g. 10-15, are multi-character).
        # Confirm interval() still produces the correct, fully-ordered
        # sequence of cells in that case.
        rdggs = RHEALPixDGGS(N_side=4)
        A = rdggs.cell((N, 9))
        B = rdggs.cell((N, 12))
        start_index = A.index(order="level")
        end_index = B.index(order="level")
        correct_count = end_index - start_index + 1
        count = 0
        old_index = start_index - 1
        for X in rdggs.interval(A, B):
            new_index = X.index(order="level")
            self.assertEqual(new_index, old_index + 1)
            old_index = new_index
            count += 1
        self.assertEqual(count, correct_count)
        self.assertEqual(old_index, end_index)

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

        # Test with the central meridian rotated: a point on the meridian
        # lon_0 lands where the same latitude on the prime meridian lands
        # in the standard grid. (Only longitude recentring exists; see
        # test_ellipsoid_rejects_nonzero_lat_0.)
        for lon_0, lat, radians in [(0.1, 0.7, True), (-50, 20, False)]:
            E = Ellipsoid(lon_0=lon_0, radians=radians)
            rdggs = RHEALPixDGGS(E)
            standard = RHEALPixDGGS(Ellipsoid(radians=radians))
            get = rdggs.cell_from_point(1, (lon_0, lat), plane=False)
            expect = standard.cell_from_point(1, (0, lat), plane=False)
            self.assertEqual(get, expect)
            # And the rotated grid's planar origin is (lon_0, 0).
            get = rdggs.cell_from_point(1, (lon_0, 0), plane=False)
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

    def test_cells_in_box(self):
        # cells_in_box enumerates candidate cells with array arithmetic; it
        # must contain every cell cells_from_region returns for the same box,
        # in all three regions, across the wrap meridian of a rotated grid,
        # and for boxes touching the poles or longitude 180.
        from random import Random

        rotated = RHEALPixDGGS(
            Ellipsoid(lon_0=129, radians=False), north_square=0, south_square=0
        )
        rng = Random(20260908)
        for rdggs in (WGS84_003, rotated):
            for _ in range(60):
                lon1 = rng.uniform(-180, 170)
                lon2 = min(180, lon1 + rng.uniform(0.5, 60))
                lat1 = rng.uniform(-90, 85)
                lat2 = min(90, lat1 + rng.uniform(0.5, 40))
                res = rng.randint(1, 3)
                rows = rdggs.cells_from_region(
                    res, (lon1, lat2), (lon2, lat1), plane=False
                )
                want = {str(c) for row in rows for c in row}
                have = set(
                    rdggs.cells_in_box(res, (lon1, lat2), (lon2, lat1), plane=False)
                )
                self.assertTrue(
                    want <= have, (lon1, lat2, lon2, lat1, res, want - have)
                )
        # Whole-globe box: every cell.
        self.assertEqual(
            len(WGS84_003.cells_in_box(2, (-180, 90), (180, -90), plane=False)), 6 * 81
        )
        # Planar box inside one face, and one reaching outside the image.
        R = WGS84_003.ellipsoid.R_A
        inside = set(WGS84_003.cells_in_box(1, (0.1 * R, 0.5 * R), (0.6 * R, 0.1 * R)))
        rows = WGS84_003.cells_from_region(1, (0.1 * R, 0.5 * R), (0.6 * R, 0.1 * R))
        self.assertTrue({str(c) for row in rows for c in row} <= inside)
        # Reaching above the band, where cells_from_region finds no corner
        # cell: the cells of the in-image part are still enumerated, with
        # the one-cell margin reaching into face P but not into a polar
        # square.
        beyond = set(WGS84_003.cells_in_box(1, (0.1 * R, 3 * R), (0.6 * R, 0.1 * R)))
        self.assertTrue(inside <= beyond)
        self.assertEqual({index[0] for index in beyond}, {"P", "Q"})

    def test_cells_from_line(self):
        rdggs = WGS84_003

        # Both endpoints in the same cell.
        cells = rdggs.cells_from_line(1, (10, 10), (11, 11), plane=False)
        self.assertEqual([str(c) for c in cells], ["Q3"])

        # An endpoint outside the (planar) grid: no trace.
        R = rdggs.ellipsoid.R_A
        outside = (0.9 * pi * R, 0.6 * pi * R)
        inside = (0.1 * R, 0.1 * R)
        self.assertEqual(rdggs.cells_from_line(1, inside, outside, plane=True), [])

        # A planar segment crossing a void of the cross-shaped image:
        # cells on both sides of the gap, in order, and nothing invented
        # in between. From the N square rightward across the void into
        # the R face's top row.
        n_square_point = (-0.6 * pi * R, 1.5 * R)  # in N (the polar square)
        r_face_point = (0.75 * pi * R, 0.7 * R)  # in R (equatorial band)
        cells = rdggs.cells_from_line(0, n_square_point, r_face_point, plane=True)
        self.assertEqual([str(c) for c in cells], ["N", "R"])

        # Antimeridian: longitude-latitude segments are straight in
        # coordinate space and by default don't wrap, so a segment from
        # longitude 179 to -179 runs the long way around through
        # longitude 0.
        cells = rdggs.cells_from_line(0, (179, 10), (-179, 10), plane=False)
        self.assertEqual([str(c) for c in cells], ["R", "Q", "P", "O"])
        # With wrap_antimeridian=True it takes the short way, across the
        # antimeridian -- and traces exactly the same cells as splitting
        # the segment at the antimeridian (RFC 7946 style) and
        # concatenating the halves' traces.
        wrapped = rdggs.cells_from_line(
            0, (179, 10), (-179, 10), plane=False, wrap_antimeridian=True
        )
        self.assertEqual([str(c) for c in wrapped], ["R", "O"])
        halves = rdggs.cells_from_line(
            0, (179, 10), (180, 10), plane=False
        ) + rdggs.cells_from_line(0, (-180, 10), (-179, 10), plane=False)
        deduped = [c for i, c in enumerate(halves) if i == 0 or c != halves[i - 1]]
        self.assertEqual(wrapped, deduped)

        # Consecutive traced cells always touch (edge- or
        # corner-adjacent), across face and region boundaries; endpoints
        # match cell_from_point. Fixed deterministic segments spanning
        # equatorial, polar, and region-crossing geometry.
        segments = [
            (
                2,
                (-11.399091685979357, 64.59420811683768),
                (-57.31373531501049, 87.0761010527302),
            ),
            (
                3,
                (-11.399091685979357, 64.59420811683768),
                (-57.31373531501049, 87.0761010527302),
            ),
            (
                3,
                (165.50972831262268, 26.027432190353828),
                (-12.994806100831056, 60.11220272134576),
            ),
            (2, (-170, -80), (170, -60)),
            (2, (-45, 88), (135, 88)),
        ]
        for res, a, b in segments:
            cells = rdggs.cells_from_line(res, a, b, plane=False)
            self.assertEqual(cells[0], rdggs.cell_from_point(res, a, plane=False))
            self.assertEqual(cells[-1], rdggs.cell_from_point(res, b, plane=False))
            for c0, c1 in itertools.pairwise(cells):
                self.assertTrue(c0.touches(c1), msg=f"{c0} !~ {c1} on {a}->{b}")
                self.assertNotEqual(c0, c1)

        # Exactness regressions: these segments pass through a cell over
        # only a tiny fraction of a cell width (about 1/50 and 1/500 of a
        # cell, confirmed by ultra-dense point location), which
        # sampling-based tracing misses.
        cells = rdggs.cells_from_line(
            3,
            (-11.399091685979357, 64.59420811683768),
            (-57.31373531501049, 87.0761010527302),
            plane=False,
        )
        self.assertIn("N455", [str(c) for c in cells])
        cells = rdggs.cells_from_line(
            3,
            (165.50972831262268, 26.027432190353828),
            (-12.994806100831056, 60.11220272134576),
            plane=False,
        )
        self.assertIn("N241", [str(c) for c in cells])

        # Pinned traces from the scalar bisection implementation this
        # replaced (issue #121): the lockstep ITP solver must find the same
        # crossings to machine precision and hence the same cells.
        import hashlib

        R = rdggs.ellipsoid.R_A
        pinned = [
            (5, (10, 10), (60, 30), False, 195, "d5791b5939d5"),
            (6, (-170, -80), (170, -60), False, 1330, "4925d794e30f"),
            (
                6,
                (165.50972831262268, 26.027432190353828),
                (-12.994806100831056, 60.11220272134576),
                False,
                1572,
                "d2092a18595b",
            ),
            (6, (-45, 88), (135, 88), False, 65, "4de4729f9639"),
            (5, (166.8, -46.4), (178.0, -38.2), False, 59, "7b7180e067d3"),
            (
                4,
                (-11.399091685979357, 64.59420811683768),
                (-57.31373531501049, 87.0761010527302),
                False,
                37,
                "2b6f1e7cf944",
            ),
            (3, (0.9 * R, 0.6 * R), (-0.4 * R, -0.7 * R), True, 45, "b4d9ce28f687"),
        ]
        for res, a, b, plane, count, digest in pinned:
            cells = rdggs.cells_from_line(res, a, b, plane=plane)
            self.assertEqual(len(cells), count, (res, a, b))
            joined = ",".join(str(c) for c in cells).encode()
            self.assertEqual(hashlib.sha1(joined).hexdigest()[:12], digest, (res, a, b))

        # Projection calls are one per inter-crossing interval (the cell
        # lookup) plus a few dozen array calls: one per piece for the scan
        # and one per solver iteration for all crossings together. The
        # scalar bisection this replaced made about 50 per crossing (some
        # 68,000 for this segment).
        class CountingProjection:
            def __init__(self, inner):
                self.inner = inner
                self.calls = 0

            def __call__(self, *args, **kwargs):
                self.calls += 1
                return self.inner(*args, **kwargs)

        counter = CountingProjection(rdggs.rhealpix)
        rdggs.rhealpix = counter
        try:
            cells = rdggs.cells_from_line(6, (-170, -80), (170, -60), plane=False)
        finally:
            del rdggs.__dict__["rhealpix"]
        self.assertEqual(len(cells), 1330)
        self.assertLess(counter.calls, 2 * len(cells))

        # On a grid rotated about the polar axis (lon_0 != 0) the face
        # edges sit on lon_0 + k*90 and the equatorial/polar seam stays
        # at +/-41.9 degrees; segments crossing each must still trace a
        # chain of touching cells that starts and ends where
        # cell_from_point says.
        rotated = RHEALPixDGGS(Ellipsoid(lon_0=129, radians=False))
        for res, a, b in [
            (2, (172, -35), (172, -47)),  # across the southern seam near NZ
            (2, (120, 10), (140, 10)),  # across the face edge at lon_0
            (1, (-30, -30), (-160, 80)),  # face edge at -51, seam, N cap
        ]:
            cells = rotated.cells_from_line(res, a, b, plane=False)
            self.assertEqual(cells[0], rotated.cell_from_point(res, a, plane=False))
            self.assertEqual(cells[-1], rotated.cell_from_point(res, b, plane=False))
            for c0, c1 in itertools.pairwise(cells):
                self.assertTrue(c0.touches(c1), msg=f"{c0} !~ {c1} on {a}->{b}")
                self.assertNotEqual(c0, c1)

    def test_cell_boundaries(self):
        from numpy import allclose

        rdggs = WGS84_003

        # Agreement with each cell's own boundary() -- same count, order,
        # and coordinates -- across cell shapes, across a region
        # boundary, and across mixed resolutions.
        cell_sets = [
            list(rdggs.cell((P, 0)).subcells()),  # quads
            list(rdggs.cell((N, 4)).subcells()),  # cap + darts + skew quads
            [rdggs.cell((Q, i)) for i in (0, 1, 2)]
            + [rdggs.cell((N, i)) for i in (6, 7, 8)],  # region-crossing
            [rdggs.cell((P, 0))] + list(rdggs.cell((P, 0)).subcells()),
        ]
        for cells in cell_sets:
            for n in (2, 3, 7):
                boundaries = rdggs.cell_boundaries(cells, n=n, plane=False)
                for c in cells:
                    expected = c.boundary(n=n, plane=False)
                    self.assertEqual(len(boundaries[c]), len(expected))
                    for got, want in zip(boundaries[c], expected):
                        self.assertTrue(
                            allclose(got, want, rtol=0, atol=1e-9),
                            msg=f"{c} n={n}: {got} != {want}",
                        )

        # The new guarantee: adjacent same-region cells' copies of their
        # shared edge points are identical values, not merely close.
        cells = list(rdggs.cell((P, 0)).subcells())
        n = 5
        boundaries = rdggs.cell_boundaries(cells, n=n, plane=False)
        a = set(map(tuple, boundaries[rdggs.cell((P, 0, 1))]))
        b = set(map(tuple, boundaries[rdggs.cell((P, 0, 2))]))
        self.assertGreaterEqual(len(a & b), n)

        # And the point of it all: strictly fewer projection calls than
        # computing each cell's boundary independently (interior edges
        # projected once, not twice).
        class CountingProjection:
            def __init__(self, inner):
                self.inner = inner
                self.count = 0
                self.calls = 0

            def __call__(self, *args, **kwargs):
                # Count projected points, whether passed singly or as arrays.
                self.count += np.size(args[0])
                self.calls += 1
                return self.inner(*args, **kwargs)

        for face, ratio in ((N, 0.6), (P, 0.15)):
            block = [rdggs.cell((face, i, j)) for i in range(9) for j in range(9)]
            counter = CountingProjection(rdggs.rhealpix)
            rdggs.rhealpix = counter
            try:
                for c in block:
                    c.boundary(n=4, plane=False)
                per_cell_calls = counter.count
                counter.count = 0
                counter.calls = 0
                boundaries = rdggs.cell_boundaries(block, n=4, plane=False)
                batched_calls = counter.count
                batched_invocations = counter.calls
            finally:
                del rdggs.__dict__["rhealpix"]
            self.assertLess(batched_calls, ratio * per_cell_calls, face)
            # One resolution and one region: a single array call. Gathering
            # the planar boundaries projects nothing, since nw_vertex picks a
            # dart's polewards vertex in the plane (issue #122).
            self.assertEqual(batched_invocations, 1, face)

        # In the equatorial region the batch exploits the projection's
        # separability: every point in one lattice column gets one
        # identical longitude and every point in one row one identical
        # latitude, for the whole block, not just where cells touch.
        R = rdggs.ellipsoid.R_A
        pitch = block[0].width(plane=True) / 3
        by_col = {}
        by_row = {}
        for c in block:
            for got, p in zip(boundaries[c], c.boundary(n=4, plane=True)):
                col = round((p[0] + pi * R) / pitch)
                row = round((p[1] + 3 * pi * R / 4) / pitch)
                by_col.setdefault(col, set()).add(got[0])
                by_row.setdefault(row, set()).add(got[1])
        self.assertEqual(len(by_col), 28)
        self.assertEqual(len(by_row), 28)
        self.assertTrue(all(len(v) == 1 for v in by_col.values()))
        self.assertTrue(all(len(v) == 1 for v in by_row.values()))

        # Planar mode is a plain convenience passthrough.
        cells = list(rdggs.cell((P, 0)).subcells())
        boundaries = rdggs.cell_boundaries(cells, n=3, plane=True)
        for c in cells:
            self.assertEqual(boundaries[c], c.boundary(n=3, plane=True))

    def test_index_arrays_match_cell_methods(self):
        # The parsing and geometry layer behind boundary_array and nuclei
        # reproduces Cell.ul_vertex, width, region, the north-west corner
        # choice and nucleus exactly, for every shape, several resolutions
        # and every DGGS configuration that changes the rules.
        from random import Random

        from numpy.testing import assert_array_equal

        from rhealpixdggs.rhp_wrappers import rhp_is_valid

        rng = Random(20260811)
        code = {"equatorial": 0, "north_polar": 1, "south_polar": -1}
        for rdggs in (WGS84_003, WGS84_003_RADIANS, WGS84_123, WGS84_122):
            cells = [c for res in range(3) for c in rdggs.grid(res)]
            for _ in range(300):
                face = rng.choice(["N", "S", "N", "S", "P", "Q"])
                depth = rng.randint(3, 9)
                digits = [rng.randrange(rdggs.N_side**2) for _ in range(depth)]
                cells.append(rdggs.cell([face] + digits))
            indices = [str(c) for c in cells]
            valid, face, digits, resolution = rdggs._parse_indices(indices)
            self.assertTrue(valid.all())
            assert_array_equal(face, [CELLS0.index(i[0]) for i in indices])
            assert_array_equal(resolution, [len(i) - 1 for i in indices])
            x, y, width, region = rdggs._index_geometry(face, digits, resolution)
            ul = np.array([c.ul_vertex(plane=True) for c in cells])
            assert_array_equal(x, ul[:, 0])
            assert_array_equal(y, ul[:, 1])
            assert_array_equal(width, [c.width() for c in cells])
            assert_array_equal(region, [code[c.region()] for c in cells])
            shift = rdggs._nw_corner(face, digits, resolution, x, y, width)
            want = [
                c.vertices(plane=True).index(c.nw_vertex(plane=True)) for c in cells
            ]
            assert_array_equal(shift, want)
            nuclei = rdggs.nuclei(indices, plane=False)
            want_nuclei = np.array([c.nucleus(plane=False) for c in cells])
            self.assertTrue(np.allclose(nuclei, want_nuclei, rtol=0, atol=1e-12))
            assert_array_equal(
                rdggs.nuclei(indices, plane=True),
                np.array([c.nucleus(plane=True) for c in cells]),
            )
        # Validity follows rhp_is_valid; invalid indices give NaN nuclei.
        odd = [
            "N",
            "S8",
            "P44",
            "X1",
            "",
            "N9",
            "Na",
            "P4 ",
            "P-1",
            "n4",
            "N45x",
            "Q" + "0" * 12,
        ]
        valid = WGS84_003._parse_indices(odd)[0]
        self.assertEqual(list(valid), [rhp_is_valid(i, WGS84_003) for i in odd])
        nuclei = WGS84_003.nuclei(odd)
        self.assertTrue(np.isnan(nuclei[~valid]).all())
        self.assertFalse(np.isnan(nuclei[valid]).any())

    def test_rings_match_cell_ring(self):
        # rings() reads the ring number off the index digits for all cells
        # at once, agreeing with Cell.ring. Invalid indices give -1, and
        # mixed resolutions are rejected because ring numbers only compare
        # within one resolution.
        from numpy.testing import assert_array_equal

        for rdggs in (WGS84_003, WGS84_003_RADIANS, WGS84_123, WGS84_122):
            for resolution in range(3):
                cells = list(rdggs.grid(resolution))
                got = rdggs.rings([str(c) for c in cells])
                self.assertEqual(got.dtype.kind, "i")
                assert_array_equal(got, [c.ring() for c in cells])
        assert_array_equal(WGS84_003.rings(["N4", "bad", "P4"]), [0, -1, 3])
        with self.assertRaises(ValueError):
            WGS84_003.rings(["N4", "P"])

    def test_centroids_match_cell_centroid(self):
        # centroids() evaluates Cell.centroid's quadrature rules for all
        # cells of each shape at once; only the summation differs (array
        # sums instead of fsum), so agreement is to rounding.
        from numpy.testing import assert_allclose, assert_array_equal

        for rdggs in (WGS84_003, WGS84_003_RADIANS, WGS84_123, WGS84_122):
            angle = 1.0 if rdggs.ellipsoid.radians else 180 / pi
            cells = [c for res in range(3) for c in rdggs.grid(res)]
            cells += list(rdggs.cell((N, 0, 3)).subcells())
            cells += list(rdggs.cell((S, 3, 1)).subcells())
            indices = [str(c) for c in cells]
            got = rdggs.centroids(indices, plane=False)
            want = np.array([c.centroid(plane=False) for c in cells])
            assert_allclose(got, want, rtol=0, atol=1e-12 * angle)
            assert_array_equal(
                rdggs.centroids(indices, plane=True),
                np.array([c.centroid(plane=True) for c in cells]),
            )
            shapes = {c.ellipsoidal_shape for c in cells}
            self.assertEqual(shapes, {"quad", "cap", "dart", "skew_quad"})
        c = WGS84_003.centroids(["P44", "", "N4"])
        self.assertTrue(np.isnan(c[1]).all())
        self.assertFalse(np.isnan(c[[0, 2]]).any())

    def test_cells_from_points_matches_cell_from_point(self):
        # cells_from_points makes cell_from_point's decisions for every
        # point at once: index strings must agree exactly, including at
        # lattice edges (nudged by one ulp either way), outside the image
        # (empty string for None) and for NaN input.
        from numpy.testing import assert_array_equal

        rng = np.random.default_rng(20260811)

        def scalar(rdggs, res, xs, ys, plane):
            cells = [
                rdggs.cell_from_point(res, (float(a), float(b)), plane=plane)
                for a, b in zip(xs, ys)
            ]
            return np.array(["" if c is None else str(c) for c in cells])

        for rdggs in (WGS84_003, WGS84_003_RADIANS, WGS84_123, WGS84_122):
            R = rdggs.ellipsoid.R_A
            angle = 1.0 if rdggs.ellipsoid.radians else 180 / pi
            lon = rng.uniform(-pi, pi, 400) * angle
            lat = np.arcsin(rng.uniform(-1, 1, 400)) * angle
            w = rdggs.cell_width(2)
            gx = -pi * R + w * np.arange(0, 4 * rdggs.N_side**2 + 1, 3)
            gy = -3 * pi * R / 4 + w * np.arange(0, 6 * rdggs.N_side**2 + 1, 4)
            ex, ey = (a.ravel() for a in np.meshgrid(gx, gy))
            px = np.concatenate(
                [
                    rng.uniform(-pi * R - 1, pi * R + 1, 400),
                    ex,
                    np.nextafter(ex, np.inf),
                    np.nextafter(ex, -np.inf),
                    [np.nan, 0.0],
                ]
            )
            py = np.concatenate(
                [
                    rng.uniform(-3 * pi * R / 4 - 1, 3 * pi * R / 4 + 1, 400),
                    ey,
                    ey,
                    np.nextafter(ey, -np.inf),
                    [0.0, np.nan],
                ]
            )
            for res in (0, 1, 2, 4, 7):
                assert_array_equal(
                    rdggs.cells_from_points(lon, lat, res, plane=False),
                    scalar(rdggs, res, lon, lat, False),
                )
                assert_array_equal(
                    rdggs.cells_from_points(px, py, res, plane=True),
                    scalar(rdggs, res, px, py, True),
                )
        # Shape is preserved and the result is a string array.
        out = WGS84_003.cells_from_points(
            np.zeros((2, 3)), np.full((2, 3), 10.0), 2, plane=False
        )
        self.assertEqual(out.shape, (2, 3))
        self.assertTrue(out.dtype.kind == "U")
        expected = str(WGS84_003.cell_from_point(2, (0.0, 10.0), plane=False))
        self.assertTrue((out == expected).all())

    def test_boundary_array(self):
        import shapely
        from numpy.testing import assert_allclose, assert_array_equal

        rdggs = WGS84_003
        cell_sets = [
            list(rdggs.cell((P, 0)).subcells()),  # quads
            list(rdggs.cell((N, 4)).subcells()),  # cap + darts + skew quads
            [rdggs.cell((Q, i)) for i in (0, 1, 2)]
            + [rdggs.cell((N, i)) for i in (6, 7, 8)],  # region-crossing
            [rdggs.cell((P, 0))] + list(rdggs.cell((P, 0)).subcells()),  # mixed
            [rdggs.cell((S, 4)), rdggs.cell((S, 4)), rdggs.cell((O, 1))],  # repeats
        ]
        for cells in cell_sets:
            indices = [str(c) for c in cells]
            for n in (2, 3, 7):
                # Shape, dtype, order: row k is cells[k].boundary(n) in order.
                b = rdggs.boundary_array(indices, n=n, plane=False)
                self.assertEqual(b.shape, (len(cells), 4 * n - 4, 2))
                self.assertEqual(b.dtype, np.float64)
                for row, c in zip(b, cells):
                    assert_allclose(
                        row, c.boundary(n=n, plane=False), rtol=0, atol=1e-9
                    )
                # Planar mode is the same arithmetic as Cell.boundary(plane=True).
                bp = rdggs.boundary_array(indices, n=n, plane=True)
                for row, c in zip(bp, cells):
                    assert_array_equal(row, np.array(c.boundary(n=n, plane=True)))
                # cell_boundaries is this array as a dictionary of point lists.
                d = rdggs.cell_boundaries(cells, n=n, plane=False)
                for row, c in zip(b, cells):
                    assert_array_equal(row, np.array(d[c]))
        # n below 2 clamps to 2, like boundary(); empty input gives an empty array.
        first = [str(c) for c in cell_sets[0]]
        # A cell's boundary is a function of the cell alone: the same whether
        # it is computed by itself, with its neighbours, or with cells of
        # another resolution, and a shared edge has one set of coordinates.
        block = [f"Q3{a}{b}" for a in range(9) for b in range(9)] + ["Q30", "N", "S8"]
        for n in (2, 4):
            together = rdggs.boundary_array(block, n=n, plane=False)
            for k, index in enumerate(block):
                alone = rdggs.boundary_array([index], n=n, plane=False)[0]
                assert_array_equal(together[k], alone, index)
        self.assertEqual(rdggs.boundary_array(first, n=1).shape, (9, 4, 2))
        self.assertEqual(rdggs.boundary_array([], n=3).shape, (0, 8, 2))
        # Invalid indices give NaN rows in place, valid ones are unaffected.
        b = rdggs.boundary_array(["P0", "", "X1", "P9", "P0", "N44"], n=2)
        self.assertEqual(b.shape, (6, 4, 2))
        self.assertTrue(np.isnan(b[[1, 2, 3]]).all())
        assert_array_equal(b[0], b[4])
        assert_allclose(
            b[5], rdggs.cell((N, 4, 4)).boundary(plane=False), rtol=0, atol=1e-9
        )
        # The rings build valid polygons in one shapely call, except where a
        # ring straddles the antimeridian (or is the cap, which spans every
        # longitude); splitting those is the caller's concern.
        for face in (P, N):
            block = [rdggs.cell((face, i, j)) for i in range(9) for j in range(9)]
            b = rdggs.boundary_array([str(c) for c in block], n=4)
            polygons = shapely.polygons(b)
            self.assertEqual(len(polygons), 81)
            lon_span = b[:, :, 0].max(axis=1) - b[:, :, 0].min(axis=1)
            unwrapped = lon_span < 180
            self.assertTrue(shapely.is_valid(polygons[unwrapped]).all())
            self.assertTrue((shapely.area(polygons[unwrapped]) > 0).all())
            if face == P:
                self.assertTrue(unwrapped.all())
            else:
                self.assertIn(
                    rdggs.cell((N, 4, 4)),
                    [c for c, u in zip(block, unwrapped) if not u],
                )

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
            get = rdggs.cells_from_parallel(1, phi, lam_min, lam_max)
            expect = [rdggs.cell([N, 6])]
            self.assertEqual(get, expect)

            phi = PI / 3
            lam_min = -PI
            lam_max = 0
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
if __name__ == "__main__":
    unittest.main()
