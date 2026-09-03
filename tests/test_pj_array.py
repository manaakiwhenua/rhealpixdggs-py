"""
Tests of the array (vectorised) projection paths against the scalar ones
(issue #88). Discrete outputs and pure-arithmetic functions must match
exactly; outputs that pass through transcendental functions match to
`ATOL` radians by default and exactly under RHP_STRICT_BITS=1.
"""

import os
import unittest
from itertools import product
from math import pi

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from pyproj import Proj

import rhealpixdggs.pj_healpix as pjh
from rhealpixdggs.ellipsoids import WGS84_E
from rhealpixdggs.utils import auth_lat, auth_rad
from tests.test_pj_healpix import ALONG_OFFSETS, EDGE_OFFSETS
from tests.test_pj_healpix_c_reference import geo_gap_deg

STRICT_BITS = os.environ.get("RHP_STRICT_BITS") == "1"
SEED = 20260811
ATOL = 1e-12
POLE_ATOL = 1e-9  # within 1e-6 rad of a pole, sqrt(3(1 - |sin phi|)) amplifies
E_VALUES = (0.0, WGS84_E, 0.1155, 0.8)


def scalar_map(fn, *arrays, **kwargs):
    """
    Apply the scalar function `fn` to every point of the broadcast arrays,
    returning its outputs as arrays of the same shape.
    """
    arrays = np.broadcast_arrays(*(np.asarray(a, dtype=float) for a in arrays))
    shape = arrays[0].shape
    flat = [a.ravel() for a in arrays]
    outs = [fn(*(float(v) for v in point), **kwargs) for point in zip(*flat)]
    if isinstance(outs[0], tuple):
        return tuple(np.array(col, dtype=float).reshape(shape) for col in zip(*outs))
    return np.array(outs).reshape(shape)


def assert_continuous(got, want, atol=ATOL, near_pole=None):
    if STRICT_BITS:
        assert_array_equal(got, want)
    elif near_pole is None:
        assert_allclose(got, want, rtol=0, atol=atol)
    else:
        assert_allclose(got[~near_pole], want[~near_pole], rtol=0, atol=atol)
        assert_allclose(got[near_pole], want[near_pole], rtol=0, atol=POLE_ATOL)


def lonlat_points():
    rng = np.random.default_rng(SEED)
    lam = [rng.uniform(-pi, pi, 10000), rng.uniform(-pi, pi, 1000)]
    phi = [
        np.arcsin(rng.uniform(-1, 1, 10000)),
        rng.choice([-1.0, 1.0], 1000) * rng.uniform(1.4, pi / 2, 1000),
    ]
    exact_phi = [0.0, pjh._PHI0, -pjh._PHI0, pi / 2, -pi / 2]
    seam_phi = [
        s * (pjh._PHI0 + d)
        for s in (1, -1)
        for d in (1e-16, -1e-16, 1e-15, -1e-15, 1e-10, -1e-10)
    ]
    exact_lam = [-pi, -pi / 2, 0.0, pi / 2, -pi + 1e-15, pi - 1e-15]
    grid = list(product(exact_lam, exact_phi + seam_phi))
    lam.append(np.array([g[0] for g in grid]))
    phi.append(np.array([g[1] for g in grid]))
    return np.concatenate(lam), np.concatenate(phi)


LONLAT = lonlat_points()
LONLAT_NEAR_POLE = np.abs(np.abs(LONLAT[1]) - pi / 2) < 1e-6


def healpix_planar_points():
    """
    Planar points inside the HEALPix image (projected random points plus
    every seam with the fuzz-scale offsets) and the offset points that fall
    outside it.
    """
    x, y = scalar_map(pjh.healpix_sphere, *LONLAT)
    seams = []
    for cap in range(4):
        x_c = -3 * pi / 4 + (pi / 2) * cap
        for t in np.linspace(0, pi / 4, 12):
            for s, dx, dy in product((1, -1), ALONG_OFFSETS, EDGE_OFFSETS):
                seams.append((x_c + t + dx, s * (pi / 2 - t) + dy))
                seams.append((x_c - t + dx, s * (pi / 2 - t) + dy))
    for yy in np.linspace(-pi / 4, pi / 4, 10):
        for dx, dy in product(EDGE_OFFSETS, ALONG_OFFSETS):
            seams.append((-pi + dx, yy + dy))
            seams.append((pi + dx, yy + dy))
    for k in range(5):
        for yy in np.linspace(-pi / 2, pi / 2, 9):
            for dx, dy in product(EDGE_OFFSETS, ALONG_OFFSETS):
                seams.append((-pi + k * pi / 2 + dx, yy + dy))
    for s in (1, -1):
        for xx in np.linspace(-pi, pi, 9):
            for dx, dy in product(ALONG_OFFSETS, EDGE_OFFSETS):
                seams.append((xx + dx, s * pi / 2 + dy))
    sx = np.array([p[0] for p in seams])
    sy = np.array([p[1] for p in seams])
    inside = np.array([pjh.in_healpix_image(a, b) for a, b in zip(sx, sy)])
    planar = (np.concatenate([x, sx[inside]]), np.concatenate([y, sy[inside]]))
    return planar, (sx[~inside], sy[~inside])


HEALPIX_PLANAR, HEALPIX_OUTSIDE = healpix_planar_points()
PLANAR_NEAR_POLE = np.abs(np.abs(HEALPIX_PLANAR[1]) - pi / 2) < 1e-6
NOT_FINITE = (
    np.array([np.nan, 0.0, np.nan, np.inf, 0.0]),
    np.array([0.0, np.nan, np.nan, 0.0, -np.inf]),
)


class HealpixArrayTestCase(unittest.TestCase):
    def test_in_healpix_image_array_matches_scalar(self):
        for x, y in (HEALPIX_PLANAR, HEALPIX_OUTSIDE, NOT_FINITE):
            got = pjh._in_healpix_image_array(x, y)
            want = np.array([pjh.in_healpix_image(a, b) for a, b in zip(x, y)])
            assert_array_equal(got, want)
        self.assertTrue(pjh._in_healpix_image_array(*HEALPIX_PLANAR).all())
        self.assertFalse(pjh._in_healpix_image_array(*HEALPIX_OUTSIDE).any())

    def test_cap_number_matches_scalar(self):
        for values in (LONLAT[0], HEALPIX_PLANAR[0], HEALPIX_OUTSIDE[0]):
            want = np.floor(2 * values / pi + 2)
            want = np.where(want >= 4, 3.0, want)
            assert_array_equal(pjh._cap_number(values), want)

    def test_healpix_sphere_array_matches_scalar(self):
        got = pjh._healpix_sphere_array(*LONLAT)
        want = scalar_map(pjh.healpix_sphere, *LONLAT)
        for g, w in zip(got, want):
            assert_continuous(g, w, near_pole=LONLAT_NEAR_POLE)

    def test_healpix_sphere_inverse_array_matches_scalar(self):
        got = pjh._healpix_sphere_inverse_array(*HEALPIX_PLANAR)
        want = scalar_map(pjh.healpix_sphere_inverse, *HEALPIX_PLANAR)
        for g, w in zip(got, want):
            assert_continuous(g, w, near_pole=PLANAR_NEAR_POLE)

    def test_healpix_ellipsoid_arrays_match_scalar(self):
        for e in E_VALUES:
            got = pjh._healpix_ellipsoid_array(*LONLAT, e=e)
            want = scalar_map(pjh.healpix_ellipsoid, *LONLAT, e=e)
            for g, w in zip(got, want):
                assert_continuous(g, w, near_pole=LONLAT_NEAR_POLE)
            got = pjh._healpix_ellipsoid_inverse_array(*HEALPIX_PLANAR, e=e)
            want = scalar_map(pjh.healpix_ellipsoid_inverse, *HEALPIX_PLANAR, e=e)
            for g, w in zip(got, want):
                assert_continuous(g, w, near_pole=PLANAR_NEAR_POLE)

    def test_healpix_factory_arrays_match_scalar(self):
        for a, e, radians in product((1.0, 6378137.0), (0.0, WGS84_E), (True, False)):
            f = pjh.healpix(a=a, e=e)
            R_A = auth_rad(a, e)
            angle = 1.0 if radians else 180 / pi
            lam, phi = LONLAT[0] * angle, LONLAT[1] * angle
            got = f(lam, phi, radians=radians)
            want = scalar_map(f, lam, phi, radians=radians)
            for g, w in zip(got, want):
                self.assertEqual(g.dtype, np.float64)
                assert_continuous(g, w, atol=ATOL * R_A, near_pole=LONLAT_NEAR_POLE)
            # Scaling by R_A and back can push fuzz-edge seam points out of
            # the image; keep the ones both paths accept.
            x, y = HEALPIX_PLANAR[0] * R_A, HEALPIX_PLANAR[1] * R_A
            keep = pjh._in_healpix_image_array(x / R_A, y / R_A)
            x, y, near_pole = x[keep], y[keep], PLANAR_NEAR_POLE[keep]
            got = f(x, y, radians=radians, inverse=True)
            want = scalar_map(f, x, y, radians=radians, inverse=True)
            for g, w in zip(got, want):
                assert_continuous(g, w, atol=ATOL * angle, near_pole=near_pole)
            # Scalars still take the scalar path and come back as numpy scalars.
            out = f(float(lam[0]), float(phi[0]), radians=radians)
            self.assertIsInstance(out, tuple)
            self.assertIsInstance(out[0], np.float64)

    def test_exact_inputs_are_bit_identical(self):
        # Inputs whose transcendental values are exact in any correct kernel
        # must agree exactly on every platform.
        lam = np.array([k * pi / 2 for k in range(-2, 3)] * 3)
        phi = np.repeat([0.0, pi / 2, -pi / 2], 5)
        got = pjh._healpix_sphere_array(lam, phi)
        want = scalar_map(pjh.healpix_sphere, lam, phi)
        for g, w in zip(got, want):
            assert_array_equal(g, w)
        band_x = [k * pi / 2 for k in range(-2, 3)]
        apex_x = [-3 * pi / 4, -pi / 4, pi / 4, 3 * pi / 4]
        x = np.array(band_x * 3 + apex_x * 2 + [-7 * pi / 8, -5 * pi / 6])
        y = np.concatenate(
            [
                np.repeat([0.0, pi / 4, -pi / 4], 5),
                np.repeat([pi / 2, -pi / 2], 4),
                [3 * pi / 8, 5 * pi / 12],
            ]
        )
        got = pjh._healpix_sphere_inverse_array(x, y)
        want = scalar_map(pjh.healpix_sphere_inverse, x, y)
        for g, w in zip(got, want):
            assert_array_equal(g, w)
        # The two seam points invert to exactly -pi, not pi - epsilon.
        self.assertEqual(got[0][-1], -pi)
        self.assertEqual(got[0][-2], -pi)

    def test_round_trips(self):
        lam, phi = LONLAT
        off_pole = np.abs(phi) < pi / 2
        x, y = pjh._healpix_sphere_array(lam, phi)
        lam2, phi2 = pjh._healpix_sphere_inverse_array(x, y)
        assert_allclose(lam2[off_pole], lam[off_pole], rtol=0, atol=1e-12)
        assert_allclose(phi2, phi, rtol=0, atol=1e-12)
        # The inverse authalic series is truncated at n^6, so the round trip
        # is exact to rounding only for small flattenings (the e = 0.8 bound
        # is the one tests/test_utils.py uses).
        # Approaching a pole the polar branch's arcsin(1 - tau^2/3) grows
        # ill-conditioned: the ellipsoidal round trip is ~1e-12 off at 1e-4
        # rad from the pole and ~1e-8 off at the pole itself.
        near = np.abs(np.abs(phi) - pi / 2) < 1e-3
        for e, atol in ((WGS84_E, 1e-12), (0.1155, 1e-12), (0.8, 1e-3)):
            x, y = pjh._healpix_ellipsoid_array(lam, phi, e=e)
            lam2, phi2 = pjh._healpix_ellipsoid_inverse_array(x, y, e=e)
            assert_allclose(lam2[off_pole], lam[off_pole], rtol=0, atol=1e-12)
            assert_allclose(phi2[~near], phi[~near], rtol=0, atol=atol)
            assert_allclose(phi2[near], phi[near], rtol=0, atol=max(atol, 1e-6))

    def test_healpix_array_matches_proj_c_reference(self):
        # PROJ's C HEALPix on the unit sphere, and on the authalic sphere for
        # the ellipsoid, as tests/test_pj_healpix_c_reference.py does.
        lon, lat = LONLAT[0] * 180 / pi, LONLAT[1] * 180 / pi
        x_ref, y_ref = Proj(proj="healpix", R=1)(lon, lat)
        x, y = pjh.healpix(a=1, e=0)(lon, lat)
        assert_allclose(x, x_ref, rtol=1e-12, atol=1e-12)
        assert_allclose(y, y_ref, rtol=1e-12, atol=1e-12)
        a, e = 5.0, WGS84_E
        R_A = auth_rad(a, e)
        beta = np.array([auth_lat(v, e) for v in lat])
        x_ref, y_ref = Proj(proj="healpix", R=R_A)(lon, beta)
        x, y = pjh.healpix(a=a, e=e)(lon, lat)
        assert_allclose(x, x_ref, rtol=1e-12, atol=1e-12 * R_A)
        assert_allclose(y, y_ref, rtol=1e-12, atol=1e-12 * R_A)
        lon2, lat2 = Proj(proj="healpix", a=a, e=e)(x, y, inverse=True)
        gaps = [geo_gap_deg(p, q) for p, q in zip(zip(lon, lat), zip(lon2, lat2))]
        self.assertLess(max(gaps), 1e-5)

    def test_value_error_parity(self):
        rng = np.random.default_rng(SEED)
        good = HEALPIX_PLANAR
        pick = rng.choice(len(HEALPIX_OUTSIDE[0]), 100, replace=False)
        n_out = len(HEALPIX_OUTSIDE[0])
        inverses = [
            pjh._healpix_sphere_inverse_array,
            pjh._healpix_ellipsoid_inverse_array,
            lambda x, y: pjh.healpix(a=1, e=0)(x, y, radians=True, inverse=True),
            lambda x, y: pjh.healpix(a=1, e=0)(x, y, radians=False, inverse=True),
        ]
        for fn in inverses:
            fn(*good)  # a wholly in-image batch does not raise
            for i in pick:
                bx, by = HEALPIX_OUTSIDE[0][i], HEALPIX_OUTSIDE[1][i]
                with self.assertRaises(ValueError):
                    pjh.healpix_sphere_inverse(float(bx), float(by))
                x = np.concatenate([good[0][:5], [bx]])
                y = np.concatenate([good[1][:5], [by]])
                with self.assertRaises(ValueError) as cm:
                    fn(x, y)
                self.assertIn("1 of 6 points", str(cm.exception))
                self.assertIn("index 5", str(cm.exception))
            with self.assertRaises(ValueError) as cm:
                fn(*HEALPIX_OUTSIDE)
            self.assertIn(f"{n_out} of {n_out} points", str(cm.exception))
            with self.assertRaises(ValueError):
                fn(*NOT_FINITE)

    def test_shapes_and_dtypes(self):
        lam = np.array([[0.1, 0.2, 0.3], [1.5, -1.5, 0.0]])
        phi = np.array([[0.4, 1.4, -1.4], [0.0, 1.0, -1.0]])
        lam_copy, phi_copy = lam.copy(), phi.copy()
        x, _ = pjh._healpix_sphere_array(lam, phi)
        self.assertEqual(x.shape, (2, 3))
        self.assertEqual(x.dtype, np.float64)
        assert_array_equal(lam, lam_copy)
        assert_array_equal(phi, phi_copy)
        self.assertFalse(np.shares_memory(x, lam))
        # Broadcasting, integer and empty inputs.
        x, _ = pjh._healpix_sphere_array(lam[0], np.asarray(0.5))
        self.assertEqual(x.shape, (3,))
        x, _ = pjh._healpix_sphere_array(np.array([1, 0]), np.array([1, 0]))
        self.assertEqual(x.dtype, np.float64)
        x, _ = pjh._healpix_sphere_inverse_array(np.array([]), np.array([]))
        self.assertEqual(x.shape, (0,))
        # 0-d array inputs take the array path; as with numpy's own ufuncs,
        # 0-d operands yield numpy scalars.
        f = pjh.healpix(a=2, e=0.1)
        x, _ = f(np.asarray(10.0), np.asarray(80.0))
        self.assertEqual(np.ndim(x), 0)
        self.assertEqual(float(x), f(10.0, 80.0)[0])


if __name__ == "__main__":
    unittest.main()
