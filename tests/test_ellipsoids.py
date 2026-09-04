"""
This Python 3.11 code tests the ``ellipsoids`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!

- Alexander Raichev (AR), 2013-01-26: Initial version based on previous test code.
"""

# *****************************************************************************
#       Copyright (C) 2013 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

import unittest

import rhealpixdggs.ellipsoids as ell


class MyTestCase(unittest.TestCase):
    def test_random_point(self):
        # Output should be in the correct range at least.
        for E in [ell.WGS84_ASPHERE, ell.WGS84_ELLIPSOID]:
            lam_min = 0
            lam_max = 35
            phi_min = -47
            phi_max = 58
            for i in range(100):
                p = E.random_point(lam_min, lam_max, phi_min, phi_max)
                self.assertTrue(p[0] >= lam_min and p[0] < lam_max)
                self.assertTrue(p[1] >= phi_min and p[0] < phi_max)

    def test_rejects_nonzero_lat_0(self):
        # Recentring in latitude is not a rotation of the ellipsoid and
        # yields a grid that is neither equal-area nor coherent at the
        # poles (issue #93), so the constructor refuses it in both angle
        # modes; lon_0 recentring and an explicit lat_0=0 are fine.
        for lat_0, radians in [(-37, False), (0.5, False), (0.1, True), (-0.3, True)]:
            with self.assertRaises(ValueError):
                ell.Ellipsoid(lon_0=0, lat_0=lat_0, radians=radians)
        E = ell.Ellipsoid(lon_0=174, lat_0=0, radians=False)
        self.assertEqual((E.lon_0, E.lat_0), (174, 0))


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
