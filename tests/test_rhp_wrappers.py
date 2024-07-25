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
        self.assertEqual(parent_id, "N123")

        # Parent is top-level cell
        parent_id = rhpw.rhp_to_parent(child_id, 0)
        self.assertEqual(parent_id, "N")

        # Immediate parent
        parent_id = rhpw.rhp_to_parent(child_id)
        self.assertEqual(parent_id, "N1234")

        # Chils is top-level cell
        parent_id = rhpw.rhp_to_parent("N")
        self.assertEqual(parent_id, "N")

        # Resolution mismatch (suppressing warning from inside function)
        parent_id = rhpw.rhp_to_parent(child_id, 7, False)
        self.assertEqual(parent_id, "N12345")

        # Invalid child id
        parent_id = rhpw.rhp_to_parent("X")
        self.assertIsNone(parent_id)

    def rhp_to_center_child(self):
        pass

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
        pass


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
