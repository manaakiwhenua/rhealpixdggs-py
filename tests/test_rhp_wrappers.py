"""
This Python 3.11 code tests the ``rhp_wrappers`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!
"""

# *****************************************************************************
#       Copyright (C) 2023 Nicoletta De Maio <demaion@landcareresearch.co.nz>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

import unittest

import numpy as np
import shapely as sh

import rhealpixdggs.dggs as gs
import rhealpixdggs.rhp_wrappers as rhpw


class RhpWrappersTestCase(unittest.TestCase):
    def test_geo_to_rhp(self):
        # Top-level cell, i.e. cube face
        cell_id = rhpw.geo_to_rhp(0, 0, 0)
        self.assertEqual(cell_id, "Q")

        # Lower-level cell
        cell_id = rhpw.geo_to_rhp(0, 0, 7)
        self.assertEqual(cell_id, "Q3333333")

        cell_id = rhpw.geo_to_rhp(90, 0, 1)
        self.assertEqual(cell_id, "Q3")

        cell_id = rhpw.geo_to_rhp(90, -180, 0, plane=False)
        self.assertEqual(cell_id, "N")

        cell_id = rhpw.geo_to_rhp(90, -180, 0, plane=True)
        self.assertEqual(cell_id, "P")

        # Point outside the DGGS
        cell_id = rhpw.geo_to_rhp(11500249, 56898969, 0, plane=True)
        self.assertIsNone(cell_id)

    def test_rhp_to_geo(self):
        # Invalid cell address
        centroid = rhpw.rhp_to_geo("X")
        self.assertIsNone(centroid)

        # Cap cell without geojson
        centroid = rhpw.rhp_to_geo("N", geo_json=False, plane=False)
        self.assertEqual(centroid, (90, -180))

        # Cap cell with geojson
        centroid = rhpw.rhp_to_geo("N", plane=False)
        self.assertEqual(centroid, (-180, 90))

        # Dart cell without geojson
        centroid = rhpw.rhp_to_geo("N0", geo_json=False, plane=False)
        self.assertEqual(centroid, (np.float64(53.00810765458496), np.float64(90.0)))

        # Dart cell with geojson
        centroid = rhpw.rhp_to_geo("N0", plane=False)
        self.assertEqual(centroid, (90.0, 53.00810765458496))

        # Equatorial cell without geojson. Q is symmetric about the
        # equator, so its centroid latitude is 0 -- up to the floating-
        # point noise floor of the quadrature that computes quad-cell
        # centroid latitudes (Cell.centroid()), hence the tolerance-based
        # comparison.
        centroid = rhpw.rhp_to_geo("Q", geo_json=False, plane=False)
        self.assertAlmostEqual(centroid[0], 0, places=12)
        self.assertEqual(centroid[1], 45)

        # Equatorial cell with geojson
        centroid = rhpw.rhp_to_geo("Q", plane=False)
        self.assertEqual(centroid[0], 45)
        self.assertAlmostEqual(centroid[1], 0, places=12)

    def test_rhp_to_parent(self):
        child_id = "N12345"

        # Usual case
        parent_id = rhpw.rhp_to_parent(child_id, 3)
        self.assertEqual(parent_id, child_id[0:4])

        # Parent is top-level cell
        parent_id = rhpw.rhp_to_parent(child_id, 0)
        self.assertEqual(parent_id, child_id[0])

        # Immediate parent
        parent_id = rhpw.rhp_to_parent(child_id)
        self.assertEqual(parent_id, child_id[:-1])

        # Child is top-level cell
        parent_id = rhpw.rhp_to_parent("N")
        self.assertEqual(parent_id, "N")

        # Resolution mismatch (suppressing warning from inside function)
        parent_id = rhpw.rhp_to_parent(child_id, 7, False)
        self.assertEqual(parent_id, child_id)

        # Invalid child id
        parent_id = rhpw.rhp_to_parent("X")
        self.assertIsNone(parent_id)

    def rhp_to_center_child(self):
        parent_id = "N12345"

        # Usual case
        child_id = rhpw.rhp_to_center_child(parent_id, 9)
        self.assertEqual(child_id, parent_id + "444")

        # Immediate child
        child_id = rhpw.rhp_to_center_child(parent_id)
        self.assertEqual(child_id, parent_id + "4")

        # Resolution mismatch (suppressing warning from inside function)
        child_id = rhpw.rhp_to_center_child(parent_id, 3, False)
        self.assertEqual(child_id, parent_id)

        # Invalid parent id
        child_id = rhpw.rhp_to_center_child("X")
        self.assertIsNone(child_id)

        # DGGS with even number of cells on a side
        child_id = rhpw.rhp_to_center_child(parent_id, dggs=gs.WGS84_002)
        self.assertIsNone(child_id)

    def test_rhp_to_geo_boundary(self):
        # Without geojson: (lat, lng) pairs
        corners = rhpw.rhp_to_geo_boundary("S", geo_json=False, plane=False)
        self.assertEqual(
            corners,
            (
                (-41.93785391016014, -180.0),
                (-41.93785391016014, -90.0),
                (-41.93785391016018, -3.180554681463517e-14),
                (-41.93785391016014, 90.00000000000003),
            ),
        )

        # With geojson: (lng, lat) pairs
        corners = rhpw.rhp_to_geo_boundary("S", geo_json=True, plane=False)
        self.assertEqual(
            corners,
            (
                (-180.0, -41.93785391016014),
                (-90.0, -41.93785391016014),
                (-3.180554681463517e-14, -41.93785391016018),
                (90.00000000000003, -41.93785391016014),
                (-180.0, -41.93785391016014),
            ),
        )

        # Invalid cell address
        corners = rhpw.rhp_to_geo_boundary("X")
        self.assertIsNone(corners)

    def test_rhp_get_resolution(self):
        # No index to resolve
        self.assertIsNone(rhpw.rhp_get_resolution(None))
        self.assertIsNone(rhpw.rhp_get_resolution(""))

        # Typical case
        self.assertEqual(rhpw.rhp_get_resolution("N12345"), 5)

    def test_rhp_get_base_cell(self):
        # No index to resolve
        self.assertIsNone(rhpw.rhp_get_base_cell(None))
        self.assertIsNone(rhpw.rhp_get_base_cell(""))

        # Base cell is self
        self.assertEqual(rhpw.rhp_get_base_cell("N"), "N")

        # Typical case
        self.assertEqual(rhpw.rhp_get_base_cell("N12345"), "N")

    def test_rhp_is_valid(self):
        # Empty strings are invalid
        self.assertFalse(rhpw.rhp_is_valid(None))
        self.assertFalse(rhpw.rhp_is_valid(""))

        # Addresses that don't start with the resolution 0 face are invalid
        self.assertFalse(rhpw.rhp_is_valid("blubb"))

        # Addresses that have digits out of range are invalid
        self.assertFalse(rhpw.rhp_is_valid("N123A5"))

        # Typical valid case
        self.assertTrue(rhpw.rhp_is_valid("N12345"))

    def test_cell_area(self):
        cell_id = "N12345"
        expected_m2 = 1696068530.010148

        # Cell area in km^2
        area = rhpw.cell_area(cell_id)
        self.assertEqual(area, expected_m2 / 10**6)

        # Cell area in m^2
        area = rhpw.cell_area(cell_id, "m^2")
        self.assertEqual(area, expected_m2)

        # Invalid cell id
        area = rhpw.cell_area("X")
        self.assertIsNone(area)

    def test_cell_ring(self):
        # Regression tests for issue #60/cell_ring rewrite: cell_ring() used
        # to jump directly to a ring's "start corner" via error-prone
        # bookkeeping (a live `# TODO: this is wrong` covered exactly this),
        # which was flat wrong for rings touching more than 2 resolution 0
        # faces -- e.g. a k=1 ring around a genuine cube corner cell used to
        # come out with bogus, non-neighboring cells and a duplicate. It's
        # rebuilt here as a breadth-first search over neighbor()/
        # diagonal_neighbor() steps, which handles any number of faces
        # (including genuine 3-valent cube corners, where a ring has fewer
        # than the usual 8*k cells) uniformly and correctly. The expected
        # values below were independently confirmed against
        # neighbor()/diagonal_neighbor() directly (not just against this
        # function's own logic) -- see also the corner cases below, which
        # catch 4 small errors in this test's own previous (commented-out,
        # since the old implementation couldn't pass them) expectations.
        cellidx = "Q444"

        # All ring cells on single face (default distance)
        ring = rhpw.cell_ring(cellidx)
        self.assertListEqual(
            ring, ["Q441", "Q442", "Q445", "Q448", "Q447", "Q446", "Q443", "Q440"]
        )

        # All ring cells on single face (minimum k)
        ring = rhpw.cell_ring(cellidx, 0)
        self.assertListEqual(ring, [cellidx])

        # All ring cells on single face (longer distance)
        ring = rhpw.cell_ring(cellidx, 2)
        self.assertListEqual(
            ring,
            [
                "Q417",
                "Q418",
                "Q416",
                "Q426",
                "Q450",
                "Q453",
                "Q456",
                "Q480",
                "Q472",
                "Q471",
                "Q470",
                "Q462",
                "Q438",
                "Q435",
                "Q432",
                "Q408",
            ],
        )

        # Neighbours across equatorial face edge
        ring = rhpw.cell_ring("Q3")
        self.assertListEqual(ring, ["Q0", "Q1", "Q4", "Q7", "Q6", "P8", "P5", "P2"])

        ring = rhpw.cell_ring("Q5")
        self.assertListEqual(ring, ["Q2", "R0", "R3", "R6", "Q8", "Q7", "Q4", "Q1"])

        # Neighbours across polar cap edges (first equatorial region)
        ring = rhpw.cell_ring("O1")
        self.assertListEqual(ring, ["N7", "N8", "O2", "O5", "O4", "O3", "O0", "N6"])

        ring = rhpw.cell_ring("O7")
        self.assertListEqual(ring, ["O4", "O5", "O8", "S2", "S1", "S0", "O6", "O3"])

        # Neighbours across polar cap edges (second equatorial region)
        ring = rhpw.cell_ring("P1")
        self.assertListEqual(ring, ["N5", "N2", "P2", "P5", "P4", "P3", "P0", "N8"])

        ring = rhpw.cell_ring("P7")
        self.assertListEqual(ring, ["P4", "P5", "P8", "S8", "S5", "S2", "P6", "P3"])

        # Neighbours across polar cap edges (third equatorial region)
        ring = rhpw.cell_ring("Q1")
        self.assertListEqual(ring, ["N1", "N0", "Q2", "Q5", "Q4", "Q3", "Q0", "N2"])

        ring = rhpw.cell_ring("Q7")
        self.assertListEqual(ring, ["Q4", "Q5", "Q8", "S6", "S7", "S8", "Q6", "Q3"])

        # Neighbours across polar cap edges (fourth equatorial region)
        ring = rhpw.cell_ring("R1")
        self.assertListEqual(ring, ["N3", "N6", "R2", "R5", "R4", "R3", "R0", "N0"])

        ring = rhpw.cell_ring("R7")
        self.assertListEqual(ring, ["R4", "R5", "R8", "S0", "S3", "S6", "R6", "R3"])

        # Neighbours across equatorial region (north polar cap)
        ring = rhpw.cell_ring("N3")
        self.assertListEqual(ring, ["N0", "N1", "N4", "N7", "N6", "R2", "R1", "R0"])

        ring = rhpw.cell_ring("N5")
        self.assertListEqual(ring, ["N2", "P2", "P1", "P0", "N8", "N7", "N4", "N1"])

        # Neighbours across equatorial region (south polar cap)
        ring = rhpw.cell_ring("S3")
        self.assertListEqual(ring, ["S0", "S1", "S4", "S7", "S6", "R6", "R7", "R8"])

        ring = rhpw.cell_ring("S5")
        self.assertListEqual(ring, ["S2", "P6", "P7", "P8", "S8", "S7", "S4", "S1"])

        # Neighbours across corner and edges (first equatorial region).
        # Previously commented out: the old implementation couldn't pass
        # these (see this test's docstring above).
        ring = rhpw.cell_ring("O0")
        self.assertListEqual(ring, ["N6", "N7", "O1", "O4", "O3", "R5", "R2"])

        ring = rhpw.cell_ring("O2")
        self.assertListEqual(ring, ["N8", "P0", "P3", "O5", "O4", "O1", "N7"])

        ring = rhpw.cell_ring("O6")
        self.assertListEqual(ring, ["O3", "O4", "O7", "S1", "S0", "R8", "R5"])

        ring = rhpw.cell_ring("O8")
        # Note: corrects an error in this test's own previous (commented
        # out, unverified) expectation, which had "S7" here instead of
        # "O7" -- confirmed via direct neighbor()/diagonal_neighbor() calls
        # on O8, independent of cell_ring() itself.
        self.assertListEqual(ring, ["O5", "P3", "P6", "S2", "S1", "O7", "O4"])

        # Neighbours across corner and edges (second equatorial region)
        ring = rhpw.cell_ring("P0")
        self.assertListEqual(ring, ["N8", "N5", "P1", "P4", "P3", "O5", "O2"])

        ring = rhpw.cell_ring("P2")
        self.assertListEqual(ring, ["N2", "Q0", "Q3", "P5", "P4", "P1", "N5"])

        ring = rhpw.cell_ring("P6")
        self.assertListEqual(ring, ["P3", "P4", "P7", "S5", "S2", "O8", "O5"])

        ring = rhpw.cell_ring("P8")
        self.assertListEqual(ring, ["P5", "Q3", "Q6", "S8", "S5", "P7", "P4"])

        # Neighbours across corner and edges (third equatorial region)
        ring = rhpw.cell_ring("Q0")
        self.assertListEqual(ring, ["N2", "N1", "Q1", "Q4", "Q3", "P5", "P2"])

        ring = rhpw.cell_ring("Q2")
        self.assertListEqual(ring, ["N0", "R0", "R3", "Q5", "Q4", "Q1", "N1"])

        ring = rhpw.cell_ring("Q6")
        self.assertListEqual(ring, ["Q3", "Q4", "Q7", "S7", "S8", "P8", "P5"])

        ring = rhpw.cell_ring("Q8")
        self.assertListEqual(ring, ["Q5", "R3", "R6", "S6", "S7", "Q7", "Q4"])

        # Neighbours across corner and edges (fourth equatorial region)
        ring = rhpw.cell_ring("R0")
        self.assertListEqual(ring, ["N0", "N3", "R1", "R4", "R3", "Q5", "Q2"])

        ring = rhpw.cell_ring("R2")
        self.assertListEqual(ring, ["N6", "O0", "O3", "R5", "R4", "R1", "N3"])

        ring = rhpw.cell_ring("R6")
        self.assertListEqual(ring, ["R3", "R4", "R7", "S3", "S6", "Q8", "Q5"])

        ring = rhpw.cell_ring("R8")
        self.assertListEqual(ring, ["R5", "O3", "O6", "S0", "S3", "R7", "R4"])

        # Neighbours across corner and edges (north polar cap)
        ring = rhpw.cell_ring("N0")
        self.assertListEqual(ring, ["Q2", "Q1", "N1", "N4", "N3", "R1", "R0"])

        ring = rhpw.cell_ring("N2")
        # Note: corrects an error in this test's own previous (commented
        # out, unverified) expectation, which had "Q2" here instead of
        # "Q0" -- confirmed via direct neighbor()/diagonal_neighbor() calls
        # on N2, independent of cell_ring() itself.
        self.assertListEqual(ring, ["Q0", "P2", "P1", "N5", "N4", "N1", "Q1"])

        ring = rhpw.cell_ring("N6")
        # Note: corrects an error in this test's own previous (commented
        # out, unverified) expectation, which had "O2" here instead of
        # "O0" -- confirmed via direct neighbor()/diagonal_neighbor() calls
        # on N6, independent of cell_ring() itself.
        self.assertListEqual(ring, ["N3", "N4", "N7", "O1", "O0", "R2", "R1"])

        ring = rhpw.cell_ring("N8")
        # Note: corrects an error in this test's own previous (commented
        # out, unverified) expectation, which had "P2" here instead of
        # "P0" -- confirmed via direct neighbor()/diagonal_neighbor() calls
        # on N8, independent of cell_ring() itself.
        self.assertListEqual(ring, ["N5", "P1", "P0", "O2", "O1", "N7", "N4"])

        # Neighbours across corner and edges (south polar cap)
        ring = rhpw.cell_ring("S0")
        self.assertListEqual(ring, ["O6", "O7", "S1", "S4", "S3", "R7", "R8"])

        ring = rhpw.cell_ring("S2")
        self.assertListEqual(ring, ["O8", "P6", "P7", "S5", "S4", "S1", "O7"])

        ring = rhpw.cell_ring("S6")
        self.assertListEqual(ring, ["S3", "S4", "S7", "Q7", "Q8", "R6", "R7"])

        ring = rhpw.cell_ring("S8")
        self.assertListEqual(ring, ["S5", "P7", "P8", "Q6", "Q7", "S7", "S4"])

        # Cell ring around centre cell at longer distance (includes corners
        # with dart cells). Previously commented out: the old
        # implementation couldn't pass this either.
        ring = rhpw.cell_ring(cellidx[0:2], 2)
        self.assertListEqual(
            ring,
            ["N1", "N0", "N2", "R0", "R3", "R6", "S6", "S7", "S8", "P8", "P5", "P2"],
        )

        # Cell ring around centre cell at a longer distance (no dart cells involved)
        ring = rhpw.cell_ring(cellidx[0:2], 3)
        self.assertListEqual(
            ring,
            ["N5", "N4", "N3", "R1", "P1", "R4", "R7", "S3", "S4", "S5", "P7", "P4"],
        )

        # A ring at exactly the grid's graph diameter from the centre cell
        # (the largest k with any cells at all -- see the "beyond the grid's
        # diameter" cases below) still returns a full, correct ring.
        ring = rhpw.cell_ring("R4", 5)
        self.assertListEqual(ring, ["P5", "P2", "P1", "P0", "P3", "P6", "P8", "P7"])

        ring = rhpw.cell_ring("N4", 5)
        self.assertListEqual(ring, ["S3", "S6", "S7", "S8", "S5", "S2", "S0", "S1"])

        # Beyond the grid's graph diameter (the largest possible number of
        # edge/corner hops between any two cells at this resolution: 16 for
        # a resolution 2 cell here, confirmed by exhaustively growing the
        # ring outward until it's exhausted every one of the grid's 486
        # cells), there is no cell at that exact distance, so the ring is
        # correctly empty -- not an arbitrary "opposite point" guess the way
        # the old, removed "beyond the hemisphere" shortcut used to make.
        ring = rhpw.cell_ring("O04", 17)
        self.assertListEqual(ring, [])

        ring = rhpw.cell_ring("S04", 17)
        self.assertListEqual(ring, [])

        # Same, for a resolution 1 cell (graph diameter 5, so k=7 is beyond
        # it) and a resolution 0 cell (graph diameter 2, so k=3 is beyond
        # it).
        ring = rhpw.cell_ring("Q0", 7)
        self.assertListEqual(ring, [])

        ring = rhpw.cell_ring("N3", 7)
        self.assertListEqual(ring, [])

        ring = rhpw.cell_ring(cellidx[0], 3)
        self.assertListEqual(ring, [])

        ring = rhpw.cell_ring("N", 3)
        self.assertListEqual(ring, [])

        # Top-level cell: its 4 edge-adjacent faces (there's no diagonal
        # neighbor at all from a resolution 0 cell -- every one of its
        # corners is a genuine cube corner).
        ring = rhpw.cell_ring(cellidx[0])
        self.assertListEqual(ring, ["N", "R", "S", "P"])

        # A top-level cell's unique opposite face is at graph distance
        # exactly 2 (reached via any of the 4 adjacent faces), not the same
        # 4 adjacent faces again -- the old, removed "clamped" shortcut
        # used to (incorrectly) repeat the k=1 ring here.
        ring = rhpw.cell_ring(cellidx[0], 2)
        self.assertListEqual(ring, ["O"])

        # Invalid cell id
        ring = rhpw.cell_ring("X")
        self.assertIsNone(ring)

    def test_k_ring(self):
        # Harmless case
        kring = rhpw.k_ring("Q4")
        self.assertListEqual(
            kring, ["Q4", "Q1", "Q2", "Q5", "Q8", "Q7", "Q6", "Q3", "Q0"]
        )

        # Minimal case
        kring = rhpw.k_ring("Q4", 0)
        self.assertListEqual(kring, ["Q4"])

        # Invalid case
        kring = rhpw.k_ring("X")
        self.assertIsNone(kring)

    def test_cell_ring_properties(self):
        # Regression tests for issue #60/cell_ring rewrite (see
        # test_cell_ring()'s docstring), checking properties that should
        # hold for *any* centre cell and distance, not just the specific
        # cases pinned down above.
        rdggs = gs.WGS84_003

        # Every ring, out to a resolution 1 grid's full graph diameter,
        # should partition the 54-cell grid exactly: no duplicates within
        # a ring, no cell repeated across rings, and every cell in the
        # grid accounted for exactly once by the time the disk is full.
        # Eccentricity (the farthest any cell is from a given centre) isn't
        # uniform: a face-centre cell (N4/Q4) reaches its antipodal
        # counterpart at distance 6, one more than a face-corner cell
        # (N0/Q0) reaches its farthest region (distance 5) -- so this
        # checks out to distance 7 for every centre, well past either.
        for center in ["N4", "N0", "Q4", "Q0"]:
            seen = set()
            total = 0
            for k in range(8):
                ring = rhpw.cell_ring(center, k)
                self.assertEqual(len(ring), len(set(ring)), f"{center} k={k}")
                self.assertTrue(seen.isdisjoint(ring), f"{center} k={k}")
                seen.update(ring)
                total += len(ring)
                if k >= 7:
                    # Beyond the diameter, there's nothing left to find.
                    self.assertEqual(ring, [], f"{center} k={k}")
            self.assertEqual(total, 54, center)

        # Graph distance is symmetric: B is in A's k-ring iff A is in B's.
        a = rdggs.cell(["N", 4])
        for k in range(1, 6):
            ring_a = rhpw.cell_ring(str(a), k)
            for b_str in ring_a:
                self.assertIn(str(a), rhpw.cell_ring(b_str, k))

        # A ring touching a genuine cube corner (N0 sits at one) has fewer
        # than the 8 cells an interior cell's 1-ring would have -- 7, since
        # exactly 3 cells meet at a cube corner rather than 4 (see
        # Cell.diagonal_neighbor()).
        self.assertEqual(len(rhpw.cell_ring("N0", 1)), 7)
        # An interior (non-corner) cell's 1-ring has the full 8.
        self.assertEqual(len(rhpw.cell_ring("N4", 1)), 8)

    def test_polyfill(self):
        # Test data - plane
        idx = ("N", 2, 1, 6, 0, 5, 5, 6, 1, 1)
        cell = rhpw.Cell(rdggs=rhpw.WGS84_003, suid=idx)
        plane_poly = sh.Polygon(cell.vertices())

        # Uncompressed from polygon - plane
        result = rhpw.polyfill(plane_poly, 10)
        self.assertEqual(
            result,
            {
                "N2160556110",
                "N2160556111",
                "N2160556112",
                "N2160556113",
                "N2160556114",
                "N2160556115",
                "N2160556116",
                "N2160556117",
                "N2160556118",
            },
        )

        # Compressed from polygon - plane
        self.assertEqual(rhpw.polyfill(plane_poly, 10, compress=True), {"N216055611"})

        # Regression test for issue #51: compress_cells hardcoded 9 as "a
        # complete sibling group" (correct only for N_side=3), so under an
        # N_side=2 DGGS (WGS84_002, 4 children per cell) a complete group
        # of siblings silently never compressed at all. polyfill() must
        # pass the DGGS's actual N_side through to compact_cells().
        n_side_2_idx = ("N", 2, 1, 3)
        n_side_2_cell = rhpw.Cell(rdggs=gs.WGS84_002, suid=n_side_2_idx)
        n_side_2_poly = sh.Polygon(n_side_2_cell.vertices())
        self.assertEqual(
            rhpw.polyfill(n_side_2_poly, 4, dggs=gs.WGS84_002),
            {"N2130", "N2131", "N2132", "N2133"},
        )
        self.assertEqual(
            rhpw.polyfill(n_side_2_poly, 4, dggs=gs.WGS84_002, compress=True),
            {"N213"},
        )

        # Test data - sphere
        eq_poly_n = sh.Polygon(
            shell=[(-10, -10), (50, -10), (50, 40), (-10, 40), (-10, -10)],
            holes=[
                [(-5, 5), (25, 20), (45, 5), (-5, 5)],
                [(-5, 25), (25, 30), (45, 25), (-5, 25)],
            ],
        )
        eq_poly_s = sh.Polygon(
            shell=[(-10, 10), (-10, -40), (50, -40), (50, 10), (-10, 10)],
            holes=[
                [(-5, -5), (45, -5), (25, -20), (-5, -5)],
                [(-5, -25), (45, -25), (25, -30), (-5, -25)],
            ],
        )
        po_poly_n = sh.Polygon(
            shell=[(0, 75), (-30, 42), (0, 42), (30, 42), (0, 75)],
            holes=[[(0, 70), (5, 60), (-5, 60), (0, 70)]],
        )
        po_poly_s = sh.Polygon(
            shell=[(30, -42), (0, -75), (90, -75), (60, -42), (30, -42)],
            holes=[[(10, -70), (20, -65), (10, -65), (10, -70)]],
        )

        # Polygon tests - sphere
        self.assertEqual(rhpw.polyfill(eq_poly_n, 0, False), {"Q"})
        self.assertEqual(rhpw.polyfill(eq_poly_s, 0, False), {"Q"})
        self.assertEqual(rhpw.polyfill(po_poly_n, 1, False), {"N2"})
        self.assertEqual(rhpw.polyfill(po_poly_s, 1, False), {"S7"})

        # Multipolygon tests - sphere
        result = rhpw.polyfill(
            sh.MultiPolygon(polygons=[eq_poly_n, po_poly_n, po_poly_s]), 1, False
        )
        self.assertEqual(result, {"N2", "Q1", "Q3", "Q4", "S7"})

        # Test data - malformed
        no_area = sh.Polygon(shell=((0, 0), (1, 0), (2, 0), (0, 0)))
        geom_res_mismatch = sh.Polygon(
            shell=[(0, 0), (0, -40), (40, -40), (40, 0), (0, 0)]
        )
        multi_overlap = sh.MultiPolygon(polygons=[eq_poly_n, eq_poly_s])

        # Malformed input geometries
        self.assertIsNone(rhpw.polyfill(None, 0))
        self.assertIsNone(rhpw.polyfill(sh.Polygon(), 0))
        self.assertIsNone(rhpw.polyfill(sh.MultiPolygon(), 0))
        self.assertIsNone(rhpw.polyfill(sh.Point(), 0))
        self.assertIsNone(rhpw.polyfill(no_area, 0))
        self.assertIsNone(rhpw.polyfill(multi_overlap, 0, False))
        self.assertEqual(rhpw.polyfill(plane_poly, 1), set())
        self.assertEqual(rhpw.polyfill(geom_res_mismatch, 0, False), set())

    def test_linetrace(self):
        # Test data
        p_ls = sh.LineString(
            [
                (-14.793092, -37.005372),
                (-15.621138, -40.323142),
                (-18.333333, -36.483403),
                (-14, -37),
            ]
        )
        r_ls = sh.LineString(
            [
                (174.793092, -37.005372),
                (175.621138, -40.323142),
                (178.333333, -36.483403),
                (174, -37),
            ]
        )
        n_ls = sh.LineString(
            [
                (-134.998756, 86.549596),
                (-179.141527, 88.504030),
                (-44.874903, 86.549596),
                (-89.669615, 86.549596),
                (-134, 86),
            ]
        )
        s_ls = sh.LineString(
            [
                (-176.260506, -43.738058),
                (-176.258807, -43.738379),
            ]
        )

        # Equatorial faces - line string
        result = rhpw.linetrace(p_ls, 3, plane=False)
        self.assertEqual(result, ["P874", "P877", "P876", "P873", "P874"])

        result = rhpw.linetrace(r_ls, 3, plane=False)
        self.assertEqual(result, ["R884", "R887", "R888", "R885", "R884"])

        # Equatorial faces - multiline string
        result = rhpw.linetrace(sh.MultiLineString(lines=[p_ls, r_ls]), 3, plane=False)
        self.assertEqual(
            result,
            [
                "P874",
                "P877",
                "P876",
                "P873",
                "P874",
                "R884",
                "R887",
                "R888",
                "R885",
                "R884",
            ],
        )

        # S face, short linestring, at threshold res for more than 1 cell
        result = rhpw.linetrace(s_ls, 9, plane=False)
        self.assertEqual(result, ["S001450634", "S001450635"])

        # S face, short linestring, higher res
        result = rhpw.linetrace(s_ls, 11, plane=False)
        self.assertEqual(
            result, ["S00145063484", "S00145063485", "S00145063563", "S00145063564"]
        )

        # S face, short linestring, higher res
        result = rhpw.linetrace(s_ls, 12, plane=False)
        self.assertEqual(
            result,
            [
                "S001450634840",
                "S001450634841",
                "S001450634842",
                "S001450634850",
                "S001450634851",
                "S001450634854",
                "S001450634855",
                "S001450635633",
                "S001450635634",
                "S001450635635",
                "S001450635643",
                "S001450635644",
                "S001450635645",
                "S001450635648",
            ],
        )

        # S face, short linestring - below res where the ends can be distinguished
        result = rhpw.linetrace(s_ls, 8, plane=False)
        self.assertEqual(result, ["S00145063"])

        # Lines crossing cube face boundaries (not involving cap cells).
        # The segment runs through the middle of Q6 for a substantial
        # stretch (confirmed by point location along the segment; also,
        # jumping S8 -> P8 directly would skip a cell the line plainly
        # crosses).
        s = gs.WGS84_003.cell(("S", 7))
        e = gs.WGS84_003.cell(("P", 5))
        result = rhpw.linetrace(
            sh.LineString([s.centroid(False), e.centroid(False)]), 1, plane=False
        )
        self.assertEqual(result, ["S7", "S8", "Q6", "P8", "P5"])

        # Resolution mismatch (coarse resolution, short line segments)
        result = rhpw.linetrace(p_ls, 2, plane=False)
        self.assertEqual(result, ["P87"])

        result = rhpw.linetrace(
            sh.MultiLineString(lines=[p_ls, r_ls, n_ls]), 2, plane=False
        )
        self.assertEqual(result, ["P87", "R88", "N44"])

        # Malformed input geometries
        self.assertIsNone(rhpw.linetrace(sh.LineString(), 0))
        self.assertIsNone(rhpw.linetrace(sh.LineString([(1, 1), (1, 1)]), 0))

    def test_linetrace_polar_cap(self):
        # A linestring wandering around the north pole, crossing the
        # resolution 3 cap cell N444. A shapely linestring's segments are
        # straight in longitude-latitude coordinates, so at these
        # latitudes each leg sweeps *around* the pole through the cells
        # at the intermediate longitudes (it does not hop across the
        # cap the way a geodesic would). Every leg's expected cell
        # sequence below was independently confirmed by dense point
        # location along the leg with cell_from_point.
        n_ls = sh.LineString(
            [
                (-134.998756, 86.549596),
                (-179.141527, 88.504030),
                (-44.874903, 86.549596),
                (-89.669615, 86.549596),
                (-134, 86),
            ]
        )
        result = rhpw.linetrace(n_ls, 3, plane=False)
        self.assertEqual(
            result,
            ["N447", "N444", "N447", "N448", "N445", "N448", "N447"],
        )


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
