"""
Invariants of the cell topology and geometry that must hold on every grid,
checked on N_side 2 and 3 grids with both polar-square placements. The
per-method tests elsewhere mostly use WGS84_003; this is where an even
N_side gets the same scrutiny.
"""

import itertools
import random
import unittest
from math import pi

import numpy as np
import shapely

from rhealpixdggs.cell import RelativePosition as RP
from rhealpixdggs.dggs import WGS84_002, WGS84_003, RHEALPixDGGS
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID, WGS84_ELLIPSOID_RADIANS, Ellipsoid
from rhealpixdggs.rhp_wrappers import k_ring
from rhealpixdggs.zoneset import ZoneSet

GRIDS = {
    "WGS84_002": WGS84_002,
    "WGS84_122": RHEALPixDGGS(
        ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=2, N_side=2
    ),
    "WGS84_003": WGS84_003,
    "WGS84_123": RHEALPixDGGS(
        ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=2, N_side=3
    ),
}
DIAGONALS = ("up_left", "up_right", "down_left", "down_right")
INVERSE = {
    RP.BEFORE: RP.AFTER,
    RP.MEETS: RP.MET_BY,
    RP.OVERLAPS: RP.OVERLAPPED_BY,
    RP.STARTS: RP.STARTED_BY,
    RP.DURING: RP.CONTAINS,
    RP.FINISHES: RP.FINISHED_BY,
    RP.EQUALS: RP.EQUALS,
}
INVERSE.update({v: k for k, v in list(INVERSE.items())})


def lattice_box(cell, depth):
    """
    The cell's planar square on the integer lattice of resolution `depth`
    within its base cell: an exact stand-in for the floating planar square,
    so a shapely oracle sees shared edges as shared rather than a nanometre
    apart. Only for comparing cells of one base cell.
    """
    N = cell.N_side
    rows, cols = cell.suid_rowcol()
    k = len(cell.suid) - 1
    row = sum(d * N ** (k - 1 - i) for i, d in enumerate(rows[1:]))
    col = sum(d * N ** (k - 1 - i) for i, d in enumerate(cols[1:]))
    scale = N ** (depth - k)
    return shapely.box(col * scale, row * scale, (col + 1) * scale, (row + 1) * scale)


class GridInvariantsTestCase(unittest.TestCase):
    def test_neighbours_are_symmetric_and_touch(self):
        for name, rdggs in GRIDS.items():
            for cell in rdggs.grid(2):
                for neighbour in cell.neighbors(plane=True).values():
                    self.assertIn(cell, neighbour.neighbors(plane=True).values(), name)
                    self.assertTrue(cell.touches(neighbour), name)
                for direction in DIAGONALS:
                    diagonal = cell.diagonal_neighbor(direction)
                    if diagonal is not None:
                        self.assertTrue(cell.touches(diagonal), name)

    def test_touching_agrees_with_k_ring_and_distance(self):
        # A cell's 1-ring is exactly the cells it touches, and its size is
        # one more than the count of edge and diagonal neighbours (7 rather
        # than 8 at a cube corner).
        for name, rdggs in GRIDS.items():
            cells = list(rdggs.grid(2))
            for cell in random.Random(1).sample(cells, 24):
                ring = set(k_ring(str(cell), 1, dggs=rdggs))
                expected = 1 + len(cell.neighbors(plane=True))
                expected += sum(
                    cell.diagonal_neighbor(d) is not None for d in DIAGONALS
                )
                self.assertEqual(len(ring), expected, f"{name} {cell}")
                for other in cells:
                    touches = cell.touches(other)
                    self.assertEqual(
                        touches,
                        str(other) in ring and other != cell,
                        f"{name} {cell} {other}",
                    )
                    self.assertEqual(cell.intersects(other), not cell.disjoint(other))
                    if touches:
                        self.assertEqual(
                            cell.distance(other), 0.0, f"{name} {cell} {other}"
                        )

    def test_de9im_and_relative_position_against_exact_lattice(self):
        for name, rdggs in GRIDS.items():
            face = rdggs.cell(("P",))
            cells = [face] + list(face.subcells()) + list(face.subcells(2))
            for a, b in itertools.product(cells, repeat=2):
                want = shapely.relate(lattice_box(a, 2), lattice_box(b, 2))
                self.assertTrue(
                    a.relate(b, want), f"{name} {a} {b}: {want} vs {a._de9im(b)}"
                )
                for direction in ((1, 0), (0, 1), (1, 1)):
                    self.assertEqual(
                        INVERSE[a.relative_position(b, direction)],
                        b.relative_position(a, direction),
                        f"{name} {a} {b} {direction}",
                    )

    def test_tiling_arithmetic(self):
        # region_overlaps and ZoneSet compaction count N_side ** 2 children.
        for name, rdggs in GRIDS.items():
            parent = rdggs.cell(("P", 0))
            children = list(parent.subcells())
            grandchildren = [g for c in children for g in c.subcells()]
            far = rdggs.cell(("S", 3))
            self.assertFalse(parent.region_overlaps(children + [far]), name)
            self.assertTrue(parent.region_overlaps(children[:-1] + [far]), name)
            self.assertFalse(parent.region_overlaps(grandchildren + [far]), name)
            mixed = children[1:] + list(children[0].subcells()) + [far]
            self.assertFalse(parent.region_overlaps(mixed), name)
            union = ZoneSet(rdggs, children[:-1]).union(children[-1], min_res=1)
            self.assertEqual(union.cells, ("P0",), name)

    def test_points_round_trip_and_lines_are_connected(self):
        for name, rdggs in GRIDS.items():
            cells = list(rdggs.grid(2))
            for cell in random.Random(2).sample(cells, 30):
                for point in (cell.nucleus(plane=False), cell.centroid(plane=False)):
                    self.assertEqual(
                        rdggs.cell_from_point(2, point, plane=False),
                        cell,
                        f"{name} {cell}",
                    )
            rng = random.Random(3)
            for _ in range(30):
                start = (rng.uniform(-180, 180), rng.uniform(-80, 80))
                end = (rng.uniform(-180, 180), rng.uniform(-80, 80))
                line = rdggs.cells_from_line(2, start, end, plane=False)
                for a, b in itertools.pairwise(line):
                    self.assertTrue(a.touches(b), f"{name} {a} {b}")

    def test_boundary_array_matches_cell_boundary_on_the_sphere(self):
        for name, rdggs in GRIDS.items():
            cells = list(rdggs.grid(2))
            rows = rdggs.boundary_array([str(c) for c in cells], n=3)
            for row, cell in zip(rows, cells):
                expected = np.array(cell.boundary(n=3, plane=False))
                self.assertTrue(np.allclose(row, expected, atol=1e-9), f"{name} {cell}")

    def test_antimeridian_edge_sign_keeps_ring_span_small(self):
        # A cell whose east or west edge lies on the antimeridian must
        # report that edge with the sign that keeps its ring's longitude
        # span under 180 degrees, identically in both boundary paths and
        # in vertices(), so a planar consumer never sees it as straddling
        # the antimeridian. The span is then exactly the cell's width in
        # longitude, 90 degrees over N_side ** resolution.
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
                rdggs = (
                    RHEALPixDGGS(N_side=side, ellipsoid=WGS84_ELLIPSOID_RADIANS)
                    if radians
                    else (WGS84_002 if side == 2 else WGS84_003)
                )
                cell = rdggs.cell((index[0], *map(int, index[1:])))
                half_turn = pi if radians else 180.0
                width = half_turn / (2 * side ** (len(index) - 1))
                for n in (2, 3, 6):
                    row = rdggs.boundary_array([index], n=n)[0]
                    scalar = np.asarray(cell.boundary(n=n, plane=False))
                    label = f"{index} N_side {side} radians={radians} n={n}"
                    self.assertTrue(np.allclose(row, scalar, atol=1e-9), label)
                    if n == 2:
                        vertices = np.asarray(cell.vertices(plane=False))
                        self.assertTrue(np.allclose(row, vertices, atol=1e-9), label)
                    lons = row[:, 0]
                    self.assertTrue((sign * lons > 0).all(), label)
                    self.assertIn(sign * half_turn, lons, label)
                    self.assertAlmostEqual(np.ptp(lons), width, msg=label)
                    bounds = shapely.polygons(row).bounds
                    self.assertAlmostEqual(bounds[2] - bounds[0], width, msg=label)

    def test_antimeridian_edge_sign_survives_polar_roundoff(self):
        # Near a pole the inverse projection amplifies planar rounding, so
        # an edge meant to lie on the antimeridian can come back a little
        # off it; the rule must still recognise and snap it.
        for index in ("N443", "S443"):
            cell = WGS84_003.cell((index[0], *map(int, index[1:])))
            for n in (2, 3, 6):
                row = WGS84_003.boundary_array([index], n=n)[0]
                self.assertTrue(
                    np.allclose(row, cell.boundary(n=n, plane=False), atol=1e-9),
                    f"{index} n={n}",
                )
                self.assertIn(180.0, row[:, 0], f"{index} n={n}")
                self.assertLess(np.ptp(row[:, 0]), 180, f"{index} n={n}")

    def test_straddling_and_cap_rings_keep_both_signs(self):
        # Cells that genuinely straddle the antimeridian, and cap cells,
        # legitimately mix +180-side and -180-side longitudes; the edge-sign
        # rule must leave them alone.
        for index in ("S0", "N"):
            suid = [ch if not ch.isdigit() else int(ch) for ch in index]
            for lons in (
                WGS84_003.boundary_array([index], n=3)[0, :, 0],
                np.array(WGS84_003.cell(suid).boundary(n=3, plane=False))[:, 0],
            ):
                self.assertTrue((lons > 90).any(), f"{index}: {lons}")
                self.assertTrue((lons < -90).any(), f"{index}: {lons}")
        for lons in (
            WGS84_003.boundary_array(["N"])[0, :, 0],
            np.array(WGS84_003.cell(("S",)).boundary(plane=False))[:, 0],
        ):
            self.assertTrue(np.allclose(sorted(lons), [-180, -90, 0, 90], atol=1e-9))
        # An equatorial cell straddling the antimeridian on a rotated grid.
        rotated = RHEALPixDGGS(ellipsoid=Ellipsoid(lon_0=50))
        for lons in (
            rotated.boundary_array(["R1"])[0, :, 0],
            np.array(rotated.cell(("R", 1)).boundary(plane=False))[:, 0],
        ):
            self.assertTrue(
                np.allclose(sorted(lons), [-160, -160, 170, 170], atol=1e-9)
            )


if __name__ == "__main__":
    unittest.main()
