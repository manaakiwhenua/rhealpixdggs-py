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

import unittest

from numpy import arcsin, array, atleast_1d, pi, rad2deg, sin, sqrt
from scipy.spatial.distance import euclidean, norm

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


_ORACLE_POLYGONS = {}


def shapely_in_healpix_image(x, y):
    # The polygon-based implementation this function replaced (issue #116),
    # kept as the oracle for the arithmetic one.
    from shapely import contains_xy
    from shapely.geometry import Polygon

    if () in _ORACLE_POLYGONS:
        return bool(contains_xy(_ORACLE_POLYGONS[()], x, y))
    eps = 1e-10
    poly = Polygon(
        [
            (-pi - eps, pi / 4 + eps),
            (-3 * pi / 4, pi / 2 + eps),
            (-pi / 2, pi / 4 + eps),
            (-pi / 4, pi / 2 + eps),
            (0, pi / 4 + eps),
            (pi / 4, pi / 2 + eps),
            (pi / 2, pi / 4 + eps),
            (3 * pi / 4, pi / 2 + eps),
            (pi + eps, pi / 4 + eps),
            (pi + eps, -pi / 4 - eps),
            (3 * pi / 4, -pi / 2 - eps),
            (pi / 2, -pi / 4 - eps),
            (pi / 4, -pi / 2 - eps),
            (0, -pi / 4 - eps),
            (-pi / 4, -pi / 2 - eps),
            (-pi / 2, -pi / 4 - eps),
            (-3 * pi / 4, -pi / 2 - eps),
            (-pi - eps, -pi / 4 - eps),
        ]
    )
    _ORACLE_POLYGONS[()] = poly
    return bool(contains_xy(poly, x, y))


def outside_healpix_image(x, y):
    # How far (x, y) lies outside the unfuzzed image, 0 if inside or on it.
    abs_y = abs(y)
    if abs(x) <= pi and abs_y <= pi / 4:
        return 0.0
    cap = min(max(int(2 * x / pi + 2), 0), 3)
    x_c = -3 * pi / 4 + (pi / 2) * cap
    return max(abs(x) - pi, abs(x - x_c) + abs_y - pi / 2, 0.0)


# Offsets from an edge, straddling both fuzz scales (1e-10 and 1e-15).
EDGE_OFFSETS = [
    0.0,
    1e-16,
    -1e-16,
    5e-16,
    -5e-16,
    1e-15,
    -1e-15,
    2e-15,
    -2e-15,
    1e-12,
    -1e-12,
    5e-11,
    -5e-11,
    1e-10,
    -1e-10,
    2e-10,
    -2e-10,
    1e-9,
    -1e-9,
    1e-6,
    -1e-6,
    0.01,
    -0.01,
]
# Offsets along an edge, paired with every EDGE_OFFSETS value across it.
ALONG_OFFSETS = [0.0, 1e-15, -1e-15, 1e-10, -1e-10, 1e-3, -1e-3]


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

    def test_in_healpix_image_matches_polygon_oracle(self):
        # in_healpix_image() is arithmetic (issue #116). It must accept the
        # whole image including its boundary, reject anything more than
        # about 2e-10 outside, and agree with the polygon it replaced
        # everywhere except within that fuzz margin, where the polygon's
        # slightly rotated triangle edges differ from the exact ones.
        from itertools import product
        from math import inf, nan

        from numpy import linspace

        points = [
            (x, y)
            for x in linspace(-pi - 1e-3, pi + 1e-3, 101)
            for y in linspace(-pi / 2 - 1e-3, pi / 2 + 1e-3, 101)
        ]
        for cap in range(4):
            x_c = -3 * pi / 4 + (pi / 2) * cap
            for t in linspace(0, pi / 4, 12):
                for s, dx, dy in product((1, -1), ALONG_OFFSETS, EDGE_OFFSETS):
                    points.append((x_c + t + dx, s * (pi / 2 - t) + dy))
                    points.append((x_c - t + dx, s * (pi / 2 - t) + dy))
        for y in linspace(-pi / 4, pi / 4, 10):
            for dx, dy in product(EDGE_OFFSETS, ALONG_OFFSETS):
                points.append((-pi + dx, y + dy))
                points.append((pi + dx, y + dy))
        for x, y in points:
            got = pjh.in_healpix_image(x, y)
            outside = outside_healpix_image(x, y)
            if outside == 0.0:
                self.assertTrue(got, (x, y))
            elif outside > 2e-10:
                self.assertFalse(got, (x, y))
            else:
                continue
            self.assertEqual(got, shapely_in_healpix_image(x, y), (x, y))
        for x, y in [(nan, 0.0), (0.0, nan), (nan, nan), (inf, 0.0), (0.0, -inf)]:
            self.assertFalse(pjh.in_healpix_image(x, y), (x, y))
            self.assertFalse(shapely_in_healpix_image(x, y), (x, y))

    def test_cap_number_clamped_at_minus_pi(self):
        # Issue #119: the polar cap number is clamped to [0, 3] at both
        # ends. Before, a longitude or planar x a hair below -pi fell in
        # cap -1, and the forward function jumped by a whole
        # (pi/2)(1 - sigma) between lam = -pi and lam = -pi - 1e-15.
        for phi in (0.8, 1.4, -1.0):
            x0, y0 = pjh.healpix_sphere(-pi, phi)
            for d in (1e-16, 1e-12, 1e-10):
                x, y = pjh.healpix_sphere(-pi - d, phi)
                self.assertAlmostEqual(x, x0, delta=1e-9, msg=(phi, d))
                self.assertEqual(y, y0)
                xa, ya = pjh._healpix_sphere_array(array([-pi - d]), array([phi]))
                self.assertAlmostEqual(xa[0], x, places=12, msg=(phi, d))
                self.assertEqual(ya[0], y)
        # Inverse: the fuzzed (-pi, +/-pi/4) corners, admitted by the
        # 1e-10 image margin (offsetting both coordinates by d puts the
        # point 2d outside the triangle edge, so d stays at or below
        # 5e-11), lie in cap 0, land on lam = -pi exactly and take the
        # cap-0 latitude, in the scalar and array paths alike.
        for s in (1, -1):
            for d in (1e-16, 2e-11, 5e-11):
                x, y = -pi - d, s * (pi / 4 + d)
                self.assertTrue(pjh.in_healpix_image(x, y), (x, y))
                lam, phi = pjh.healpix_sphere_inverse(x, y)
                self.assertEqual(lam, -pi, (x, y))
                tau = 2 - 4 * abs(y) / pi
                self.assertEqual(phi, s * arcsin(1 - tau**2 / 3), (x, y))
                la, pa = pjh._healpix_sphere_inverse_array(array([x]), array([y]))
                self.assertEqual(la[0], lam, (x, y))
                self.assertAlmostEqual(pa[0], phi, places=14, msg=(x, y))

    def test_out_of_bounds_raises_value_error(self):
        # Regression test for issue #52: out-of-bounds coordinates used to
        # print an error message and return a sentinel (float("inf"), or
        # bare None for healpix_ellipsoid_inverse) instead of raising --
        # the sentinel case silently propagated inf into any downstream
        # arithmetic, and the None case crashed several frames away with a
        # confusing numpy TypeError the moment healpix()'s closure tried
        # to array()-and-unpack it.
        with self.assertRaises(ValueError):
            pjh.healpix_sphere_inverse(0, 100)
        with self.assertRaises(ValueError):
            pjh.healpix_ellipsoid_inverse(0, 100)
        # Same via the public factory function.
        f = pjh.healpix(a=1, e=0.5)
        with self.assertRaises(ValueError):
            f(0, 100, radians=True, inverse=True)


if __name__ == "__main__":
    unittest.main()
