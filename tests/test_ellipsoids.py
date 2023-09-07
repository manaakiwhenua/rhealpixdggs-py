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

# Import standard modules.
import unittest

# Import my modules.
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


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
