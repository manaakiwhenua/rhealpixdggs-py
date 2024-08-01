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

# Import standard modules.
import unittest

# Import my modules.
import rhealpixdggs.rhp_wrappers as rhpw


class RhpWrappersTestCase(unittest.TestCase):
    def test_geo_to_rhp(self):
        # Top-level cell, i.e. cube face
        cell_id = rhpw.geo_to_rhp(0, 0, 0)
        self.assertEqual(cell_id, "Q")

        # Lower-level cell
        cell_id = rhpw.geo_to_rhp(0, 0, 7)
        self.assertEqual(cell_id, "Q3333333")

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
        self.assertEqual(centroid, (52.948408366310105, 89.99999999999996))

        # Dart cell with geojson
        centroid = rhpw.rhp_to_geo("N0", plane=False)
        self.assertEqual(centroid, (89.99999999999996, 52.948408366310105))

        # Equatorial cell without geojson
        centroid = rhpw.rhp_to_geo("Q", geo_json=False, plane=False)
        self.assertEqual(centroid, (0, 45))

        # Equatorial cell with geojson
        centroid = rhpw.rhp_to_geo("Q", plane=False)
        self.assertEqual(centroid, (45, 0))

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

    def test_rhp_to_geo_boundary(self):
        expected_lat = -41.87385774220941

        # Without geojson: (lat, lng) pairs
        corners = rhpw.rhp_to_geo_boundary("S", geo_json=False, plane=False)
        self.assertEqual(
            corners,
            (
                (expected_lat, -180),
                (expected_lat, -90),
                (expected_lat, 0),
                (expected_lat, 90),
            ),
        )

        # With geojson: (lng, lat) pairs
        corners = rhpw.rhp_to_geo_boundary("S", geo_json=True, plane=False)
        self.assertEqual(
            corners,
            (
                (-180, expected_lat),
                (-90, expected_lat),
                (0, expected_lat),
                (90, expected_lat),
                (-180, expected_lat),
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
        expected_m2 = 1697972133.2114245

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
        cellidx = "Q444"

        # All ring cells on single face (default distance)
        ring = rhpw.cell_ring(cellidx)
        self.assertListEqual(
            ring, ["Q440", "Q441", "Q442", "Q445", "Q448", "Q447", "Q446", "Q443"]
        )

        # All ring cells on single face (minimum k)
        ring = rhpw.cell_ring(cellidx, 0)
        self.assertListEqual(ring, [cellidx])

        # All ring cells on single face (longer distance)
        ring = rhpw.cell_ring(cellidx, 2)
        self.assertListEqual(
            ring,
            [
                "Q408",
                "Q416",
                "Q417",
                "Q418",
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
            ],
        )

        # Neighbours across equatorial face edge
        ring = rhpw.cell_ring("Q3")
        self.assertListEqual(ring, ["P2", "Q0", "Q1", "Q4", "Q7", "Q6", "P8", "P5"])

        ring = rhpw.cell_ring("Q5")
        self.assertListEqual(ring, ["Q1", "Q2", "R0", "R3", "R6", "Q8", "Q7", "Q4"])

        # Neighbours across polar cap edges (first equatorial region)
        ring = rhpw.cell_ring("O1")
        self.assertListEqual(ring, ["N6", "N7", "N8", "O2", "O5", "O4", "O3", "O0"])

        ring = rhpw.cell_ring("O7")
        self.assertListEqual(ring, ["O3", "O4", "O5", "O8", "S2", "S1", "S0", "O6"])

        # Neighbours across polar cap edges (second equatorial region)
        ring = rhpw.cell_ring("P1")
        self.assertListEqual(ring, ["N8", "N5", "N2", "P2", "P5", "P4", "P3", "P0"])

        ring = rhpw.cell_ring("P7")
        self.assertListEqual(ring, ["P3", "P4", "P5", "P8", "S8", "S5", "S2", "P6"])

        # Neighbours across polar cap edges (third equatorial region)
        ring = rhpw.cell_ring("Q1")
        self.assertListEqual(ring, ["N2", "N1", "N0", "Q2", "Q5", "Q4", "Q3", "Q0"])

        ring = rhpw.cell_ring("Q7")
        self.assertListEqual(ring, ["Q3", "Q4", "Q5", "Q8", "S6", "S7", "S8", "Q6"])

        # Neighbours across polar cap edges (fourth equatorial region)
        ring = rhpw.cell_ring("R1")
        self.assertListEqual(ring, ["N0", "N3", "N6", "R2", "R5", "R4", "R3", "R0"])

        ring = rhpw.cell_ring("R7")
        self.assertListEqual(ring, ["R3", "R4", "R5", "R8", "S0", "S3", "S6", "R6"])

        # Neighbours across equatorial region (north polar cap)
        ring = rhpw.cell_ring("N3")
        self.assertListEqual(ring, ["R0", "N0", "N1", "N4", "N7", "N6", "R2", "R1"])

        # ring = rhpw.cell_ring("N5")
        # TODO: assertion

        # TODO: neighbours across equatorial region (south polar cap)
        # ring = rhpw.cell_ring("S3")
        # TODO: assertion

        # ring rhpw.cell_ring("S5")
        # TODO: assertion

        # Neighbours across corner and edges
        # TODO: "O0", "O2", "P0", "P2", "Q0"
        # ring = rhpw.cell_ring("Q2")
        # self.assertListEqual(ring, ["N1", "N0", "R0", "R3", "Q4", "Q3", "Q1"])
        # TODO: "R0", "R2"
        # TODO: "O6", "O8", "P6", "P8", "Q6", "Q8", "R6", "R8"

        # Cell ring at a distance beyond the resolution
        # TODO

        # Top-level cell (regular case)
        ring = rhpw.cell_ring(cellidx[0])
        self.assertListEqual(ring, ["R", "S", "P", "N"])

        # Top-level cell (clamped case)
        ring = rhpw.cell_ring(cellidx[0], 2)
        self.assertListEqual(ring, ["R", "S", "P", "N"])

        # Invalid cell id
        ring = rhpw.cell_ring("X")
        self.assertIsNone(ring)


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
