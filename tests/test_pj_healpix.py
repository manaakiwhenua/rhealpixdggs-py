"""
This Python 3.11 code tests the ``pj_healpix`` module.
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

# Import third-party modules.
from scipy.spatial.distance import euclidean, norm
from numpy import atleast_1d, array, rad2deg, deg2rad, pi, sqrt, sin, arcsin

# Import standard modules.
import unittest

# Import my modules.
import rhealpixdggs.pj_healpix as pjh
from rhealpixdggs.utils import auth_lat, auth_rad


# Relative error function.
def rel_err(get, expect):
    a = euclidean(atleast_1d(get), atleast_1d(expect))
    b = norm(expect)
    if b == 0:
        return a
    else:
        return a / b


# Define lon-lat input points to test.
phi_0 = arcsin(2.0 / 3)
a = (0, pi / 3)
b = (0, -pi / 3)
inputs = [
    (0, 0),
    (0, phi_0),
    (0, -phi_0),
    (pi / 2, 0),
    (-pi / 2, 0),
    (-pi, 0),
    (-pi, pi / 2),
    (-pi, -pi / 2),
    a,
    b,
]


class MyTestCase(unittest.TestCase):
    def test_healpix_sphere(self):
        # Expected outputs of healpix_sphere() applied to inputs.
        sigma_a = sqrt(3 - 3 * sin(a[1]))
        ha = (pi / 4 * (1 - sigma_a), pi / 4 * (2 - sigma_a))
        hb = (ha[0], -ha[1])
        healpix_sphere_outputs = [
            (0, 0),
            (0, pi / 4),
            (0, -pi / 4),
            (pi / 2, 0),
            (-pi / 2, 0),
            (-pi, 0),
            (-3 * pi / 4, pi / 2),
            (-3 * pi / 4, -pi / 2),
            ha,
            hb,
        ]

        # Forward projection should be correct on test points.
        given = inputs
        get = [pjh.healpix_sphere(*p) for p in given]
        expect = healpix_sphere_outputs
        # Fuzz to allow for rounding errors:
        error = 1e-12
        for i in range(len(get)):
            self.assertTrue(rel_err(get[i], expect[i]) < error)

        # Inverse of projection of a point p should yield p.
        given = get
        get = [pjh.healpix_sphere_inverse(*q) for q in given]
        expect = inputs
        for i in range(len(get)):
            self.assertTrue(rel_err(get[i], expect[i]) < error)

        # Inverse projection of p below should return longitude of -pi.
        # Previously, it was returning a number slightly less than pi
        # because of a rounding error, which get magnified by
        # wrap_lon())
        p = (-7 * pi / 8, 3 * pi / 8)
        get = pjh.healpix_sphere_inverse(*p)
        expect = (-pi, arcsin(1 - 1.0 / 12))
        self.assertEqual(get, expect)
        q = (-5 * pi / 6, 5 * pi / 12)
        get = pjh.healpix_sphere_inverse(*q)
        expect = (-pi, arcsin(1 - 1.0 / 27))
        self.assertEqual(get, expect)

    def test_healpix_ellipsoid(self):
        # Expected output of healpix_ellipsoid() applied to inputs.
        e = 0.8
        healpix_ellipsoid_outputs = []
        for p in inputs:
            lam, phi = p
            beta = auth_lat(phi, e=e, radians=True)
            q = pjh.healpix_sphere(lam, beta)
            healpix_ellipsoid_outputs.append(q)

        # Forward projection should be correct on test points.
        given = inputs
        get = [pjh.healpix_ellipsoid(*p, e=e) for p in given]
        expect = healpix_ellipsoid_outputs
        # Fuzz to allow for rounding errors:
        error = 1e-12
        for i in range(len(get)):
            self.assertTrue(rel_err(get[i], expect[i]) < error)

        # Inverse of projection of a point p should yield p.
        given = get
        get = [pjh.healpix_ellipsoid_inverse(*q, e=e) for q in given]
        expect = inputs
        # Fuzz for rounding errors based on the error of the approximation to
        # the inverse authalic latitude function:
        alpha = pi / 4
        alpha_ = auth_lat(
            auth_lat(alpha, e, radians=True), e, radians=True, inverse=True
        )
        error = 10 * rel_err(alpha_, alpha)
        for i in range(len(get)):
            self.assertTrue(rel_err(get[i], expect[i]) < error)

    def test_healpix(self):
        inputs = [(-pi, pi / 3), (0, pi / 4), (pi / 2, -pi / 6)]

        # Should agree with healpix_ellipsoid and healpix_ellipsoid_inverse.
        e = 0.5
        a = 7
        R_A = auth_rad(a, e)
        f = pjh.healpix(a=a, e=e)
        for p in inputs:
            get = f(*p, radians=True)
            expect = tuple(R_A * array(pjh.healpix_ellipsoid(*p, e=e)))
            for i in range(len(expect)):
                self.assertAlmostEqual(get[i], expect[i])
            get = f(*get, radians=True, inverse=True)
            expect = tuple(array(expect) / R_A)
            expect = pjh.healpix_ellipsoid_inverse(*expect, e=e)
            for i in range(len(expect)):
                self.assertAlmostEqual(get[i], expect[i])

        # Should work in degrees mode.
        for p in inputs:
            get = f(*rad2deg(p), radians=False)
            expect = f(*p, radians=True)
            for i in range(len(expect)):
                self.assertAlmostEqual(get[i], expect[i])
            # ------------------------------------------------------------------------------


if __name__ == "__main__":
    unittest.main()
