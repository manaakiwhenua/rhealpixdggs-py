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
from random import Random

from numpy import deg2rad, pi

import rhealpixdggs.ellipsoids as ell


class MyTestCase(unittest.TestCase):
    def test_random_point(self):
        # Output should be in the correct range at least.
        for E in [ell.WGS84_ASPHERE, ell.WGS84_ELLIPSOID]:
            lam_min = 0
            lam_max = 35
            phi_min = -47
            phi_max = 58
            rng = Random(5)
            for _ in range(100):
                p = E.random_point(lam_min, lam_max, phi_min, phi_max, rng=rng)
                self.assertTrue(lam_min <= p[0] < lam_max)
                self.assertTrue(phi_min <= p[1] < phi_max)
            # The whole ellipsoid by default, keeping half a degree off the poles.
            for _ in range(100):
                lam, phi = E.random_point(rng=rng)
                self.assertTrue(-180 <= lam < 180)
                self.assertTrue(-89.5 <= phi <= 89.5)
            self.assertEqual(
                E.random_point(rng=Random(9)), E.random_point(rng=Random(9))
            )
        lam, phi = ell.WGS84_ELLIPSOID_RADIANS.random_point(rng=Random(9))
        self.assertTrue(-pi <= lam < pi)
        self.assertTrue(abs(phi) <= pi / 2 - pi / 360)

    def test_comparison_with_other_types(self):
        self.assertNotEqual(ell.WGS84_ELLIPSOID, 3)
        self.assertFalse(ell.WGS84_ELLIPSOID == 3)
        self.assertIs(ell.WGS84_ELLIPSOID.__eq__(3), NotImplemented)

    def test_equality(self):
        E = ell.Ellipsoid
        self.assertEqual(E(lon_0=50), E(lon_0=50))
        self.assertEqual(E(lon_0=50), E(lon_0=410))
        self.assertEqual(E(lon_0=-180), E(lon_0=180))
        self.assertEqual(
            E(lon_0=0.1, radians=True), E(lon_0=0.1 + 2 * pi, radians=True)
        )
        self.assertNotEqual(E(lon_0=50), E(lon_0=0))
        self.assertNotEqual(E(lon_0=50), E(lon_0=50 + 1e-6))
        # The angular unit is part of identity, as for a PROJ CRS.
        self.assertNotEqual(E(), E(radians=True))
        self.assertNotEqual(E(lon_0=50), E(lon_0=deg2rad(50), radians=True))
        self.assertNotEqual(E(), E(f=1 / 298.257222101))

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
