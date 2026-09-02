"""
This Python 3.11 code tests the ``utils`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!

- Alexander Raichev (AR), 2012-11-21: Initial version.
"""

# *****************************************************************************
#       Copyright (C) 2012 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

import unittest
from math import asin, pi, sin

from scipy.integrate import quad

import rhealpixdggs.utils as ut
from rhealpixdggs.ellipsoids import WGS84_E


def auth_lat_by_defining_integral(phi: float, e: float) -> float:
    """
    Authalic latitude of geodetic latitude `phi` (radians) on an ellipsoid
    of eccentricity `e`, computed by numerically integrating its defining
    integral: beta = arcsin(q(phi)/q(pi/2)), where (after substituting
    u = sin(t) to keep the integrand a smooth rational function)
    q(phi) = (1 - e^2) * integral of du/(1 - e^2 u^2)^2 from 0 to sin(phi).

    This is deliberately independent of auth_lat()'s own closed-form
    formula and power series, so it can serve as ground truth for both.
    """

    def integrand(u):
        return (1 - e**2) / (1 - (e * u) ** 2) ** 2

    # epsabs/epsrel = 1e-13 keeps quad comfortably clear of its internal
    # roundoff limits while remaining orders of magnitude tighter than the
    # assertion tolerances used below.
    q, _ = quad(integrand, 0, sin(phi), epsabs=1e-13, epsrel=1e-13)
    q_pole, _ = quad(integrand, 0, 1, epsabs=1e-13, epsrel=1e-13)
    return asin(q / q_pole)


class UtilsTestCase(unittest.TestCase):
    def test_wrap_longitude(self):
        angles = [0, 90, -90, 180, -180, 360, -360]
        get = list(map(ut.wrap_longitude, angles))
        expect = [0, 90, -90, -180, -180, 0, 0]
        for i in range(len(expect)):
            self.assertAlmostEqual(get[i], expect[i])

    def test_wrap_latitude(self):
        angles = [0, 45, -45, 90, -90, 135, -135]
        get = list(map(ut.wrap_latitude, angles))
        expect = [0, 45, -45, 90, -90, -45, 45]
        for i in range(len(expect)):
            self.assertAlmostEqual(get[i], expect[i])

    def test_auth_lat_radians(self):
        e = WGS84_E
        lat = [0, pi / 6, pi / 4, pi / 3, pi / 2]
        get = []
        for i in lat:
            auth_lat = ut.auth_lat(i, e, inverse=False, radians=True)
            get.append(auth_lat)

        expect = [
            0.0,
            0.5216614083606359,
            0.7831589561070269,
            1.045256493205824,
            1.5707963267948966,
        ]
        for i in range(len(expect)):
            self.assertAlmostEqual(get[i], expect[i])
        # arbitrary-precision results by direct formula for e = 0.08181919104281579
        # 0.5216614083606359112452395455242782683847964175740775168430389631... (lat = pi/6)
        # 0.7831589561070268978794326495267670650498563055857653245023990384... (lat = pi/4)
        # 1.0452564932058239452065937261188103171077779616846897361946711143... (lat = pi/3)
        # 1.5707963267948966192313216916397514420985846996875529104874722961... (lat = pi/2)
        # (1.5707963267948966192313216916397514420985846996875529104874722961... = pi/2)

    def test_auth_lat_degrees(self):
        e = WGS84_E
        lat = [0, 30, 45, 60, 90]
        get = []
        for i in lat:
            auth_lat = ut.auth_lat(i, e, inverse=False, radians=False)
            get.append(auth_lat)

        expect = [
            0.0,
            29.888997033914993,
            44.87170287280393,
            59.888785569338516,
            90.0,
        ]
        for i in range(len(expect)):
            self.assertAlmostEqual(get[i], expect[i])
        # arbitrary-precision results by direct formula for e = 0.08181919104281579
        # 29.88899703391499358370509010777406472843219431935738736016... (lat = 30°)
        # 44.8717028728039283027993984398344502944087940882537429129... (lat = 45°)
        # 59.8887855693385177230524553429705646149296124796160869265... (lat = 60°)
        # 90.0 (lat = 90°)

    def test_auth_lat_radians_inverse(self):
        e = WGS84_E
        lat = [0, pi / 6, pi / 4, pi / 3, pi / 2]
        auth_lat = []
        for i in lat:
            auth_lat.append(ut.auth_lat(i, e, inverse=False, radians=True))

        auth_lat_inverse = []
        for i in auth_lat:
            inverse = ut.auth_lat(i, e, inverse=True, radians=True)
            auth_lat_inverse.append(inverse)

        expect = lat
        for i in range(len(expect)):
            self.assertAlmostEqual(auth_lat_inverse[i], expect[i])

    def test_auth_lat_degrees_inverse(self):
        e = WGS84_E
        lat = [0, 30, 45, 60, 90]
        auth_lat = []
        for i in lat:
            auth_lat.append(ut.auth_lat(i, e, inverse=False, radians=False))

        auth_lat_inverse = []
        for i in auth_lat:
            inverse = ut.auth_lat(i, e, inverse=True, radians=False)
            auth_lat_inverse.append(inverse)

        expect = lat
        for i in range(len(expect)):
            self.assertAlmostEqual(auth_lat_inverse[i], expect[i])

    def test_auth_lat_matches_defining_integral(self):
        # auth_lat() switches implementation on flattening: a closed-form
        # direct formula for f > 1/150, a power series for f <= 1/150. No
        # predefined Ellipsoid has a flattening anywhere near 1/150 (WGS84's
        # f ~ 1/298 uses the series; the sphere presets short-circuit via
        # e == 0), so the direct-formula branch is only reachable through
        # user-constructed high-flattening ellipsoids. Validate BOTH
        # branches against the same independent ground truth -- the
        # numerically-integrated defining integral of authalic latitude --
        # which also confirms the two branches agree with each other across
        # the f = 1/150 switch point.
        cases = [
            # (eccentricity, which branch it exercises)
            (0.8, "direct formula, extreme flattening (f = 0.4)"),
            (0.1155, "direct formula, just above the f = 1/150 threshold"),
            (WGS84_E, "power series, well below the threshold"),
        ]
        lat_degs = [-89, -80, -60, -45, -30, -10, 0, 10, 30, 45, 60, 80, 89]
        for e, label in cases:
            for deg in lat_degs:
                phi = deg * pi / 180
                expected = auth_lat_by_defining_integral(phi, e)
                got_rad = ut.auth_lat(phi, e, radians=True)
                self.assertAlmostEqual(got_rad, expected, places=12, msg=label)
                # Degrees mode must agree with radians mode.
                got_deg = ut.auth_lat(deg, e, radians=False)
                self.assertAlmostEqual(
                    got_deg, expected * 180 / pi, places=10, msg=label
                )

    def test_auth_lat_large_flattening_basic_properties(self):
        # Structural properties of the direct-formula branch (f > 1/150):
        # authalic latitude is an odd function of geodetic latitude, fixes
        # the equator and the poles, pulls every other latitude strictly
        # toward the equator on an oblate ellipsoid, and is strictly
        # increasing.
        e = 0.8
        self.assertEqual(ut.auth_lat(0, e, radians=True), 0)
        self.assertAlmostEqual(ut.auth_lat(pi / 2, e, radians=True), pi / 2, places=12)
        self.assertAlmostEqual(
            ut.auth_lat(-pi / 2, e, radians=True), -pi / 2, places=12
        )
        previous = None
        for deg in range(-89, 90, 2):
            phi = deg * pi / 180
            beta = ut.auth_lat(phi, e, radians=True)
            self.assertAlmostEqual(beta, -ut.auth_lat(-phi, e, radians=True), places=12)
            if deg != 0:
                self.assertLess(abs(beta), abs(phi))
            if previous is not None:
                self.assertGreater(beta, previous)
            previous = beta

    def test_auth_lat_large_flattening_inverse_round_trip(self):
        # The inverse always uses the power series in the third flattening
        # n, truncated at n^6, whatever the flattening -- so for a large
        # flattening the round trip is only approximate: for e = 0.8,
        # n = 0.25 and the truncation error is of order n^7 ~ 6e-5. Assert
        # the round trip lands within 1e-3 radians (~0.06 degrees), which
        # the observed worst case (~5e-4 rad) clears with headroom to spare
        # while still catching any gross regression in either direction of
        # the conversion.
        e = 0.8
        for deg in range(-89, 90, 7):
            phi = deg * pi / 180
            round_trip = ut.auth_lat(
                ut.auth_lat(phi, e, radians=True), e, radians=True, inverse=True
            )
            self.assertAlmostEqual(round_trip, phi, delta=1e-3)
        # For a flattening just above the direct-formula threshold
        # (f ~ 1/149, n ~ 1/298), the same n^7 truncation bound puts the
        # round trip within double-precision noise.
        e = 0.1155
        for deg in range(-89, 90, 7):
            phi = deg * pi / 180
            round_trip = ut.auth_lat(
                ut.auth_lat(phi, e, radians=True), e, radians=True, inverse=True
            )
            self.assertAlmostEqual(round_trip, phi, places=12)


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
