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

import os
import unittest
from math import asin, pi, sin, sqrt

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
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


def auth_lat_series_before_refactor(phi: float, e: float, inverse: bool) -> float:
    """
    The power-series expressions exactly as ``auth_lat`` wrote them before
    the coefficients were extracted into ``_auth_lat_coefficients``; `phi`
    in radians. Kept verbatim so the refactor can be shown bit-identical.
    """
    n = (1 - sqrt(1 - e**2)) / (1 + sqrt(1 - e**2))
    if not inverse:
        authalic_lat = phi + (
            n
            * (
                -4 / 3
                + n
                * (
                    -4 / 45
                    + n
                    * (
                        88 / 315
                        + n
                        * (538 / 4725 + n * (20824 / 467775 + n * (-44732 / 2837835)))
                    )
                )
            )
            * sin(2 * phi)
            + n
            * (
                n
                * (
                    34 / 45
                    + n
                    * (
                        8 / 105
                        + n
                        * (
                            -2482 / 14175
                            + n * (-37192 / 467775 + n * (-12467764 / 212837625))
                        )
                    )
                )
            )
            * sin(4 * phi)
            + n
            * (
                n
                * (
                    n
                    * (
                        -1532 / 2835
                        + n
                        * (
                            -898 / 14175
                            + n * (54968 / 467775 + n * 100320856 / 1915538625)
                        )
                    )
                )
            )
            * sin(6 * phi)
            + n
            * (
                n
                * (
                    n
                    * (
                        n
                        * (
                            6007 / 14175
                            + n * (24496 / 467775 + n * (-5884124 / 70945875))
                        )
                    )
                )
            )
            * sin(8 * phi)
            + n
            * (n * (n * (n * (n * (-23356 / 66825 + n * (-839792 / 19348875))))))
            * sin(10 * phi)
            + n * (n * (n * (n * (n * (n * 570284222 / 1915538625))))) * sin(12 * phi)
        )
        return authalic_lat
    common_lat = phi + (
        n
        * (
            4 / 3
            + n
            * (
                4 / 45
                + n
                * (
                    -16 / 35
                    + n
                    * (-2582 / 14175 + n * (60136 / 467775 + n * 28112932 / 212837625))
                )
            )
        )
        * sin(2 * phi)
        + n
        * (
            n
            * (
                46 / 45
                + n
                * (
                    152 / 945
                    + n
                    * (
                        -11966 / 14175
                        + n * (-21016 / 51975 + n * 251310128 / 638512875)
                    )
                )
            )
        )
        * sin(4 * phi)
        + n
        * (
            n
            * (
                n
                * (
                    3044 / 2835
                    + n
                    * (3802 / 14175 + n * (-94388 / 66825 + n * (-8797648 / 10945935)))
                )
            )
        )
        * sin(6 * phi)
        + n
        * (
            n
            * (
                n
                * (
                    n
                    * (
                        6059 / 4725
                        + n * (41072 / 93555 + n * (-1472637812 / 638512875))
                    )
                )
            )
        )
        * sin(8 * phi)
        + n
        * (n * (n * (n * (n * (768272 / 467775 + n * 455935736 / 638512875)))))
        * sin(10 * phi)
        + n * (n * (n * (n * (n * (n * 4210684958 / 1915538625))))) * sin(12 * phi)
    )
    return common_lat


# RHP_STRICT_BITS=1 asserts the array helpers reproduce the scalar ones to
# the last bit; the default tolerates the ulp-level differences that
# platforms with SIMD transcendental kernels can show.
STRICT_BITS = os.environ.get("RHP_STRICT_BITS") == "1"


def assert_continuous(got, want, atol):
    if STRICT_BITS:
        assert_array_equal(got, want)
    else:
        assert_allclose(got, want, rtol=0, atol=atol)


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

    def test_auth_lat_coefficient_refactor_is_bit_identical(self):
        # The series coefficients come from _auth_lat_coefficients (issue
        # #88); each is the `n * (...)` factor the original expression
        # evaluated before multiplying by the sine, so results must match
        # that expression to the bit.
        phis = [-pi / 2 + k * pi / 4000 for k in range(4001)]
        for e in (WGS84_E, 1e-3, 0.02, 0.1155, 0.8):
            for inverse in (False, True):
                if not inverse and 1 - sqrt(1 - e**2) > 1 / 150:
                    continue  # the direct formula, not the series
                for phi in phis:
                    self.assertEqual(
                        ut.auth_lat(phi, e, inverse=inverse, radians=True),
                        auth_lat_series_before_refactor(phi, e, inverse),
                        (e, inverse, phi),
                    )

    def test_auth_lat_array_matches_scalar(self):
        rng = np.random.default_rng(20260811)
        phis = np.concatenate(
            [
                rng.uniform(-pi / 2, pi / 2, 20000),
                np.linspace(-pi / 2, pi / 2, 2001),
                [0.0, -0.0, pi / 2, -pi / 2, pi / 2 - 1e-9, -pi / 2 + 1e-9],
                pi / 2 - np.logspace(-15, -3, 25),
            ]
        )
        near_pole = np.abs(np.abs(phis) - pi / 2) < 1e-6
        for e in (0.0, WGS84_E, 0.1155, 0.8):
            for inverse in (False, True):
                for radians in (True, False):
                    x = phis if radians else phis * 180 / pi
                    scale = 1 if radians else 180 / pi
                    got = ut._auth_lat_array(x, e, inverse=inverse, radians=radians)
                    want = np.array(
                        [
                            ut.auth_lat(float(v), e, inverse=inverse, radians=radians)
                            for v in x
                        ]
                    )
                    self.assertEqual(got.dtype, np.float64)
                    assert_continuous(got[~near_pole], want[~near_pole], 1e-12 * scale)
                    assert_continuous(got[near_pole], want[near_pole], 1e-9 * scale)

    def test_auth_lat_array_matches_defining_integral(self):
        lat_degs = np.array([-89, -80, -60, -45, -30, -10, 0, 10, 30, 45, 60, 80, 89])
        phis = lat_degs * pi / 180
        for e in (0.8, 0.1155, WGS84_E):
            expected = [auth_lat_by_defining_integral(float(p), e) for p in phis]
            got_rad = ut._auth_lat_array(phis, e, radians=True)
            got_deg = ut._auth_lat_array(lat_degs, e, radians=False)
            for i, want in enumerate(expected):
                self.assertAlmostEqual(got_rad[i], want, places=12)
                self.assertAlmostEqual(got_deg[i], want * 180 / pi, places=10)

    def test_auth_lat_array_shapes(self):
        phi = np.array([[0.1, 0.2], [0.3, 0.4]])
        original = phi.copy()
        out = ut._auth_lat_array(phi, WGS84_E, radians=True)
        self.assertEqual(out.shape, (2, 2))
        assert_array_equal(phi, original)
        # A sphere returns a copy, not the input.
        out = ut._auth_lat_array(phi, 0.0, radians=True)
        assert_array_equal(out, phi)
        self.assertIsNot(out, phi)
        self.assertFalse(np.shares_memory(out, phi))
        # 0-d and integer inputs come back as float64 arrays.
        out = ut._auth_lat_array(np.asarray(45), WGS84_E, radians=False)
        self.assertEqual(out.shape, ())
        self.assertEqual(out.dtype, np.float64)
        self.assertAlmostEqual(float(out), ut.auth_lat(45, WGS84_E, radians=False))

    def test_wrap_arrays_are_bit_identical_to_scalars(self):
        rng = np.random.default_rng(20260811)
        exact = [0.0, -0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 0.5, -0.5, 1.5, -1.5]
        for radians in (True, False):
            half = pi if radians else 180
            values = np.concatenate(
                [
                    rng.uniform(-3 * half, 3 * half, 20000),
                    np.array(exact) * half,
                    np.array(exact) * half + 1e-12,
                    np.array(exact) * half - 1e-12,
                ]
            )
            for array_fn, scalar_fn in (
                (ut._wrap_longitude_array, ut.wrap_longitude),
                (ut._wrap_latitude_array, ut.wrap_latitude),
            ):
                got = array_fn(values, radians=radians)
                want = np.array([scalar_fn(float(v), radians=radians) for v in values])
                assert_array_equal(got, want)
                self.assertEqual(got.dtype, np.float64)
            self.assertEqual(ut._wrap_longitude_array(np.asarray(190.0)).shape, ())


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
