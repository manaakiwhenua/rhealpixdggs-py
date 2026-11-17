"""
Tests for the distance and bearing between cells' nuclei.
"""

import unittest
from random import Random

import numpy as np
import pyproj

from rhealpixdggs.cell import Cell
from rhealpixdggs.dggs import (
    UNIT_003,
    WGS84_003,
    WGS84_003_OGC,
    WGS84_003_RADIANS,
    RHEALPixDGGS,
)


def random_cells(rdggs, count, seed):
    rng = Random(seed)
    cells = []
    for _ in range(count):
        depth = rng.randint(0, 8)
        face = rng.choice("NOPQRS")
        cells.append(rdggs.cell([face, *(rng.randrange(9) for _ in range(depth))]))
    return cells


class NucleusDistanceTestCase(unittest.TestCase):
    def test_sphere_within_bounds_of_geodesic(self):
        # The bounds measured in docs/proximity_error.py (issue #196).
        for rdggs in (WGS84_003, WGS84_003_OGC):
            a = [str(c) for c in random_cells(rdggs, 3000, 1)]
            b = [str(c) for c in random_cells(rdggs, 3000, 2)]
            distance = rdggs.nucleus_distances(a, b)
            exact = rdggs.nucleus_distances(a, b, exact=True)
            apart = exact > 1
            ratio = distance[apart] / exact[apart]
            self.assertLess(np.abs(ratio - 1).max(), 0.00115)
            bearing = rdggs.bearings(a, b)
            azimuth = rdggs.bearings(a, b, exact=True)
            error = np.abs((bearing - azimuth + 180) % 360 - 180)
            near = apart & (exact < 5e6)
            self.assertLess(error[near].max(), 0.1)
            self.assertLess(np.median(error[near]), 0.05)

    def test_exact_is_the_geodesic_between_nuclei(self):
        rdggs = WGS84_003
        cells_a = random_cells(rdggs, 200, 3)
        cells_b = random_cells(rdggs, 200, 4)
        geod = pyproj.Geod(a=rdggs.ellipsoid.a, f=rdggs.ellipsoid.f)
        want = [
            geod.inv(*a.nucleus(plane=False), *b.nucleus(plane=False))
            for a, b in zip(cells_a, cells_b, strict=True)
        ]
        a = [str(c) for c in cells_a]
        b = [str(c) for c in cells_b]
        np.testing.assert_allclose(
            rdggs.nucleus_distances(a, b, exact=True), [w[2] for w in want], rtol=1e-12
        )
        np.testing.assert_allclose(
            rdggs.bearings(a, b, exact=True), [w[0] for w in want], rtol=1e-12
        )

    def test_cell_methods_match_arrays(self):
        rdggs = WGS84_003
        cells_a = random_cells(rdggs, 100, 5)
        cells_b = random_cells(rdggs, 100, 6)
        a = [str(c) for c in cells_a]
        b = [str(c) for c in cells_b]
        for exact in (False, True):
            with self.subTest(exact=exact):
                np.testing.assert_allclose(
                    [
                        x.nucleus_distance(y, exact)
                        for x, y in zip(cells_a, cells_b, strict=True)
                    ],
                    rdggs.nucleus_distances(a, b, exact),
                    rtol=1e-12,
                )
                np.testing.assert_allclose(
                    [
                        x.bearing(y, exact)
                        for x, y in zip(cells_a, cells_b, strict=True)
                    ],
                    rdggs.bearings(a, b, exact),
                    rtol=1e-12,
                )
        p4 = rdggs.cell(("P", 4))
        self.assertEqual(p4.nucleus_distance(p4), 0.0)
        self.assertEqual(p4.nucleus_distance(p4, exact=True), 0.0)
        p1 = rdggs.cell(("P", 1))
        self.assertAlmostEqual(p4.nucleus_distance(p1), p1.nucleus_distance(p4))
        self.assertAlmostEqual(p4.bearing(p1), 0.0)
        self.assertAlmostEqual(p1.bearing(p4), 180.0)
        self.assertAlmostEqual(p4.bearing(rdggs.cell(("P", 5))), 90.0)
        self.assertAlmostEqual(p4.bearing(rdggs.cell(("P", 3))), -90.0)

    def test_every_grid_has_them(self):
        # No index strings needed for the Cell methods.
        rdggs = RHEALPixDGGS(N_side=5)
        a, b = rdggs.cell(("P", 12)), rdggs.cell(("Q", 0, 7))
        geod = pyproj.Geod(a=rdggs.ellipsoid.a, f=rdggs.ellipsoid.f)
        want = geod.inv(*a.nucleus(plane=False), *b.nucleus(plane=False))
        self.assertAlmostEqual(a.nucleus_distance(b, exact=True), want[2])
        self.assertAlmostEqual(a.nucleus_distance(b) / want[2], 1.0, delta=0.00115)

    def test_units(self):
        # The unit sphere gives the central angle; a grid in radians gives
        # the same distance and the bearing in radians.
        self.assertAlmostEqual(
            UNIT_003.nucleus_distances(["P4"], ["P5"])[0], np.pi / 6, places=12
        )
        self.assertAlmostEqual(
            UNIT_003.nucleus_distances(["P4"], ["P5"], exact=True)[0], np.pi / 6
        )
        a, b = ["P1", "N0", "S88"], ["P2", "Q4", "O1"]
        for exact in (False, True):
            np.testing.assert_allclose(
                WGS84_003_RADIANS.nucleus_distances(a, b, exact),
                WGS84_003.nucleus_distances(a, b, exact),
                rtol=1e-12,
            )
            np.testing.assert_allclose(
                WGS84_003_RADIANS.bearings(a, b, exact),
                np.radians(WGS84_003.bearings(a, b, exact)),
                rtol=1e-12,
            )

    def test_pairing_and_invalid(self):
        rdggs = WGS84_003
        one_to_many = rdggs.nucleus_distances(["P4"], ["P1", "P5", "P7"])
        many_to_one = rdggs.nucleus_distances(["P1", "P5", "P7"], ["P4"])
        np.testing.assert_array_equal(one_to_many, many_to_one)
        self.assertEqual(rdggs.nucleus_distances([], []).shape, (0,))
        with self.assertRaises(ValueError):
            rdggs.nucleus_distances(["P1", "P5"], ["P1", "P5", "P7"])
        distance = rdggs.nucleus_distances(["P4", "bad"], ["P1", "P1"])
        self.assertFalse(np.isnan(distance[0]))
        self.assertTrue(np.isnan(distance[1]))
        self.assertTrue(np.isnan(rdggs.bearings(["P4"], ["X"], exact=True)[0]))
        with self.assertRaises(ValueError):
            Cell(WGS84_003, ["P", 4]).nucleus_distance(Cell(WGS84_003_OGC, ["P", 4]))
        with self.assertRaises(ValueError):
            Cell(WGS84_003, ["P", 4]).bearing(Cell(WGS84_003, []))


if __name__ == "__main__":
    unittest.main()
