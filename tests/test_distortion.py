"""
Tests for the ``distortion`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!

CHANGELOG:

- Alexander Raichev (AR), 2012-08-28: Initial version.
"""
# *****************************************************************************
#       Copyright (C) 2012 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import third-party modules.
from scipy.spatial.distance import euclidean, norm
from numpy import pi, rad2deg, deg2rad, sqrt, sin, cos

# Import standard modules.
import unittest
from itertools import product

# Import my modules.
import rhealpixdggs.distortion as di
from rhealpixdggs.ellipsoids import (
    WGS84_ASPHERE_RADIANS,
    WGS84_ELLIPSOID,
    WGS84_ELLIPSOID_RADIANS,
)
from rhealpixdggs.projection_wrapper import Proj

# Relative error function.
def rel_err(get, expect):
    a = euclidean(get, expect)
    b = norm(expect)
    if b == 0:
        return a
    else:
        return a / b


sphere = WGS84_ASPHERE_RADIANS
ellipsoid = WGS84_ELLIPSOID_RADIANS
ellipsoid_deg = WGS84_ELLIPSOID
lam = pi / 3
phi = pi / 5
lam_deg, phi_deg = rad2deg([lam, phi])
# Mercator projection:
merc = Proj(ellipsoid=sphere, proj="merc")
# Lambert cylindrical equal area projection:
cea = Proj(ellipsoid=sphere, proj="cea")
cea_e = Proj(ellipsoid=ellipsoid, proj="cea")
cea_ed = Proj(ellipsoid=ellipsoid_deg, proj="cea")


class MyTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_fff_coeffs(self):
        epsilon = 10e-4

        E, F, G = di.fff_coeffs(merc, lam, phi)
        self.assertTrue(rel_err(E, sphere.R ** 2) < epsilon)
        self.assertAlmostEqual(F, 0)
        self.assertTrue(rel_err(G, (sphere.R / cos(phi)) ** 2) < epsilon)

        E, F, G = di.fff_coeffs(cea, lam, phi)
        self.assertTrue(rel_err(E, sphere.R ** 2) < epsilon)
        self.assertAlmostEqual(F, 0)
        self.assertTrue(rel_err(G, (sphere.R * cos(phi)) ** 2) < epsilon)

        # Degrees mode output should agree with radians mode output.
        get = di.fff_coeffs(cea_ed, lam_deg, phi_deg)
        expect = di.fff_coeffs(cea_e, lam, phi)
        for i in range(len(expect)):
            self.assertAlmostEqual(get[i], expect[i])

    def test_scale_factors(self):
        epsilon = 10e-4

        s_M, s_P, s_A, theta = di.scale_factors(merc, lam, phi)
        self.assertTrue(rel_err(s_M, 1 / cos(phi)) < epsilon)
        self.assertTrue(rel_err(s_P, 1 / cos(phi)) < epsilon)
        self.assertTrue(rel_err(s_A, 1 / cos(phi) ** 2) < epsilon)
        self.assertTrue(rel_err(theta, pi / 2) < epsilon)

        s_M, s_P, s_A, theta = di.scale_factors(cea, lam, phi)
        self.assertTrue(rel_err(s_M, cos(phi)) < epsilon)
        self.assertTrue(rel_err(s_P, 1 / cos(phi)) < epsilon)
        self.assertTrue(rel_err(s_A, 1) < epsilon)
        self.assertTrue(rel_err(theta, pi / 2) < epsilon)

        s_M, s_P, s_A, theta = di.scale_factors(cea_e, lam, phi)
        a = ellipsoid.a
        e = ellipsoid.e
        R_A = ellipsoid.R_A
        s_M_correct = a * cos(phi) / (R_A * sqrt(1 - (e * sin(phi)) ** 2))
        s_P_correct = 1 / s_M_correct
        self.assertTrue(rel_err(s_M, s_M_correct) < epsilon)
        self.assertTrue(rel_err(s_P, s_P_correct) < epsilon)
        self.assertTrue(rel_err(s_A, 1) < epsilon)
        self.assertTrue(rel_err(theta, pi / 2) < epsilon)

        # Degrees mode output should agree with radians mode output.
        get = di.scale_factors(cea_ed, lam_deg, phi_deg)
        expect = di.scale_factors(cea_e, lam, phi)
        for i in range(len(expect)):
            self.assertAlmostEqual(get[i], expect[i])

    def test_distortion(self):
        epsilon = 10e-4
        # Conformal projections should have maximum angular distortion
        # equal to 0 and linear distortion equal to 1.0.
        mad, ld, ad = di.distortion(merc, lam, phi)[1:]
        self.assertAlmostEqual(mad, 0.0, places=6)
        self.assertAlmostEqual(ld, 1.0, places=6)

        # Area preserving projections should have area distortion
        # equal to 1.0.
        mad, ld, ad = di.distortion(cea, lam, phi)[1:]
        self.assertAlmostEqual(ad, 1.0, places=6)

        # Degrees mode output should agree with radians mode output.
        get = di.distortion(cea_ed, lam_deg, phi_deg)
        expect = list(di.distortion(cea_e, lam, phi))
        # Entry 1 is an angular measurement.
        expect[1] = rad2deg(expect[1])
        for i in range(len(expect)):
            self.assertAlmostEqual(get[i], expect[i])

    def test_utm_zone(self):
        lam, phi = (-pi + 3 * pi / 180, pi / 3)
        self.assertEqual(di.utm_zone(lam, phi), 1)
        self.assertEqual(di.utm_zone(-lam, -phi), 60)
        lam, phi = deg2rad((-15, 0))
        self.assertEqual(di.utm_zone(lam, phi), 28)
        self.assertEqual(di.utm_zone(0, 85 * pi / 180), None)
        self.assertEqual(di.utm_zone(0, -84 * pi / 180), None)


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
