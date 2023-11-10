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


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
