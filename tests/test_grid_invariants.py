"""Scalar and batched boundaries use the same antimeridian convention."""

import unittest
from math import pi

import numpy as np
import shapely
from numpy.testing import assert_allclose

from rhealpixdggs.dggs import WGS84_002, WGS84_003, RHEALPixDGGS
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID_RADIANS, Ellipsoid


class BoundaryInvariantTestCase(unittest.TestCase):
    def test_antimeridian_edges_keep_the_ring_narrow(self):
        cases = [
            (2, "R11", 1),
            (2, "R13", 1),
            (2, "R31", 1),
            (2, "R33", 1),
            (2, "O00", -1),
            (2, "O02", -1),
            (2, "O20", -1),
            (2, "O22", -1),
            (3, "R2", 1),
            (3, "R22", 1),
            (3, "O0", -1),
            (3, "O00", -1),
        ]
        for side, index, sign in cases:
            for radians in (False, True):
                grid = (
                    RHEALPixDGGS(N_side=side, ellipsoid=WGS84_ELLIPSOID_RADIANS)
                    if radians
                    else (WGS84_002 if side == 2 else WGS84_003)
                )
                cell = grid.cell((index[0], *map(int, index[1:])))
                half = pi if radians else 180.0
                width = half / (2 * side ** (len(index) - 1))
                for n in (2, 3, 6):
                    with self.subTest(index=index, radians=radians, n=n):
                        row = grid.boundary_array([index], n=n)[0]
                        scalar = np.asarray(cell.boundary(n=n, plane=False))
                        assert_allclose(row, scalar, rtol=0, atol=1e-9)
                        for ring in (row, scalar):
                            self.assertTrue((sign * ring[:, 0] > 0).all())
                            self.assertAlmostEqual(np.ptp(ring[:, 0]), width)
                            self.assertIn(sign * half, ring[:, 0])
                            bounds = shapely.polygons(ring).bounds
                            self.assertAlmostEqual(bounds[2] - bounds[0], width)
                        if n == 2:
                            assert_allclose(row, cell.vertices(plane=False), atol=1e-9)

    def test_boundary_array_matches_cell_boundary_on_the_sphere(self):
        for grid in (WGS84_002, WGS84_003, RHEALPixDGGS(ellipsoid=Ellipsoid(lon_0=50))):
            for resolution in (1, 2):
                cells = list(grid.grid(resolution))
                for n in (2, 3):
                    rows = grid.boundary_array([str(cell) for cell in cells], n=n)
                    for cell, row in zip(cells, rows):
                        with self.subTest(side=grid.N_side, cell=str(cell), n=n):
                            # Compare the coordinates themselves, not modulo a full turn.
                            assert_allclose(
                                row, cell.boundary(n=n, plane=False), rtol=0, atol=1e-9
                            )

    def test_antimeridian_roundoff_near_the_poles(self):
        for index in ("N443", "S443"):
            cell = WGS84_003.cell((index[0], *map(int, index[1:])))
            for n in (2, 3, 6):
                with self.subTest(index=index, n=n):
                    row = WGS84_003.boundary_array([index], n=n)[0]
                    assert_allclose(
                        row, cell.boundary(n=n, plane=False), rtol=0, atol=1e-9
                    )
                    self.assertIn(180, row[:, 0])
                    self.assertLess(np.ptp(row[:, 0]), 180)

    def test_caps_and_genuine_crossings_keep_their_longitudes(self):
        for index in ("N", "S"):
            cell = WGS84_003.cell((index,))
            for ring in (
                cell.boundary(plane=False),
                WGS84_003.boundary_array([index])[0],
            ):
                assert_allclose(
                    np.sort(np.asarray(ring)[:, 0]), [-180, -90, 0, 90], atol=1e-9
                )
        rotated = RHEALPixDGGS(ellipsoid=Ellipsoid(lon_0=50))
        cell = rotated.cell(("R", 1))
        for ring in (cell.boundary(plane=False), rotated.boundary_array(["R1"])[0]):
            assert_allclose(
                np.sort(np.asarray(ring)[:, 0]), [-160, -160, 170, 170], atol=1e-9
            )
