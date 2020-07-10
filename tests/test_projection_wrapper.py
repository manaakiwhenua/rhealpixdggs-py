"""
This Python 3.3 code tests the ``projection_wrapper`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!

CHANGELOG:

- Alexander Raichev (AR), 2013-01-26: Initial version based on a previous tests.
"""
# *****************************************************************************
#       Copyright (C) 2013 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import third-party modules.
from numpy import array, pi, rad2deg, deg2rad
from scipy.spatial.distance import euclidean, norm

# Import standard modules.
import unittest

# Import my modules.
import rhealpixdggs.projection_wrapper as pw
import rhealpixdggs.ellipsoids as ell
from rhealpixdggs.pj_rhealpix import rhealpix_ellipsoid, rhealpix_ellipsoid_inverse

# Relative error function.
def rel_err(get, expect):
    a = euclidean(get, expect)
    b = norm(expect)
    if b == 0:
        return a
    else:
        return a / b


class MyTestCase(unittest.TestCase):
    # ------------------------------------------------------------------------------# Test Proj class.
    def test_call(self):
        # Test scaling.
        ns = 1
        ss = 2
        E = ell.WGS84_ELLIPSOID_RADIANS
        f = pw.Proj(ellipsoid=E, proj="rhealpix", north_square=ns, south_square=ss)
        R_A = E.R_A
        e = E.e
        sample = E.graticule(n=20)
        for p in sample:
            # Forward.
            get = f(*p)
            expect = rhealpix_ellipsoid(*p, north_square=ns, south_square=ss, e=e)
            expect = tuple(R_A * array(expect))
            self.assertEqual(get, expect)
            # Inverse.
            get = f(*get, inverse=True)
            expect = array(expect) / R_A
            expect = rhealpix_ellipsoid_inverse(
                *expect, north_square=ns, south_square=ss, e=e
            )
            expect = tuple(expect)
            self.assertEqual(get, expect)

        # Test degrees mode
        E_DEG = ell.WGS84_ELLIPSOID
        f_deg = pw.Proj(
            ellipsoid=E_DEG, proj="rhealpix", north_square=ns, south_square=ss
        )
        for i in range(len(sample)):
            get = f_deg(*rad2deg(sample[i]))
            expect = f(*sample[i])
            for j in range(len(expect)):
                # Allow for tiny rounding errors.
                self.assertAlmostEqual(get[j], expect[j])

        # Test origin shifting.
        p = (pi / 3, -pi / 6)
        z = (0, 0)
        E = ell.Ellipsoid(R=5, lon_0=p[0], lat_0=p[1])
        f = pw.Proj(ellipsoid=E, proj="rhealpix", north_square=ns, south_square=ss)
        self.assertEqual(f(*p), z)
        self.assertEqual(f(*z, inverse=True), p)

        # Test origin shifting in degrees.
        p = (-20, 50)
        z = (0, 0)
        E_DEG = ell.Ellipsoid(R=5, lon_0=p[0], lat_0=p[1], radians=False)
        f = pw.Proj(ellipsoid=E_DEG, proj="rhealpix", north_square=ns, south_square=ss)
        self.assertEqual(f(*p), z)
        self.assertEqual(f(*z, inverse=True), p)


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
