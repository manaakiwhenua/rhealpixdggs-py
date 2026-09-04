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
import rhealpixdggs.pj_rhealpix as pjr
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

SQUARE_PAIRS = list(product(range(4), repeat=2))
REGION_CODE = {"equatorial": 0, "north_polar": 1, "south_polar": -1}
# The last 102 LONLAT points are the exact/seam grid; keep it in the smaller
# per-square-pair set.
LONLAT_SMALL = (
    np.concatenate([LONLAT[0][:2000], LONLAT[0][-102:]]),
    np.concatenate([LONLAT[1][:2000], LONLAT[1][-102:]]),
)


def rhealpix_planar_points(north_square, south_square):
    """
    Planar points inside the (north_square, south_square)-rHEALPix image:
    projected random points, the polar squares' diagonals and boundaries
    with fuzz-scale offsets, and the band and square edges; plus the offset
    points that fall outside.
    """
    x, y = scalar_map(
        pjr.rhealpix_sphere,
        *LONLAT_SMALL,
        north_square=north_square,
        south_square=south_square,
    )
    seams = []
    small = (0.0, 1e-16, -1e-16, 1e-15, -1e-15, 2e-15, -2e-15, 1e-10, -1e-10)
    for square, sign in ((north_square, 1), (south_square, -1)):
        u_x, u_y = -3 * pi / 4 + square * pi / 2, sign * pi / 2
        for t in np.linspace(-pi / 4, pi / 4, 25):
            for d in small:
                seams.append((u_x + t, u_y + t + d))
                seams.append((u_x + t, u_y - t + d))
        for cx, cy in product((-pi / 4, 0.0, pi / 4), repeat=2):
            for dx, dy in product(small, repeat=2):
                seams.append((u_x + cx + dx, u_y + cy + dy))
    for x_e in {-pi, pi} | {-pi + k * pi / 2 for k in range(5)}:
        for yy in np.linspace(-3 * pi / 4, 3 * pi / 4, 7):
            for dx, dy in product(EDGE_OFFSETS, ALONG_OFFSETS):
                seams.append((x_e + dx, yy + dy))
    for y_e in (-3 * pi / 4, -pi / 4, pi / 4, 3 * pi / 4):
        for xx in np.linspace(-pi, pi, 7):
            for dx, dy in product(ALONG_OFFSETS, EDGE_OFFSETS):
                seams.append((xx + dx, y_e + dy))
    sx = np.array([q[0] for q in seams])
    sy = np.array([q[1] for q in seams])
    kw = {"north_square": north_square, "south_square": south_square}
    inside = np.array([pjr.in_rhealpix_image(a, b, **kw) for a, b in zip(sx, sy)])

    def accepted(a, b):
        # Inside the rHEALPix image by its 1e-15 fuzz, yet the uncombined
        # point can fall outside the HEALPix image's 1e-10 fuzz.
        try:
            pjr.rhealpix_sphere_inverse(a, b, **kw)
        except ValueError:
            return False
        return True

    ok = inside & np.array([accepted(a, b) for a, b in zip(sx, sy)])
    planar = (np.concatenate([x, sx[ok]]), np.concatenate([y, sy[ok]]))
    rejected = (sx[inside & ~ok], sy[inside & ~ok])
    return planar, (sx[~inside], sy[~inside]), rejected


RHEALPIX_PLANAR = {}
RHEALPIX_OUTSIDE = {}
RHEALPIX_REJECTED = {}
for _pair in SQUARE_PAIRS:
    RHEALPIX_PLANAR[_pair], RHEALPIX_OUTSIDE[_pair], RHEALPIX_REJECTED[_pair] = (
        rhealpix_planar_points(*_pair)
    )
# Random (interior) points only: seam points within the fuzz margins are not
# round-trippable by either path.
LONLAT_INTERIOR = (LONLAT_SMALL[0][:2000], LONLAT_SMALL[1][:2000])
N_INTERIOR_HEALPIX = 11000  # LONLAT's random points; the exact grid follows


def accepts(f, x, y, radians):
    try:
        f(float(x), float(y), radians=radians, inverse=True)
    except ValueError:
        return False
    return True


def triangle_codes(x, y, north_square, south_square, inverse):
    tri, region = [], []
    for a, b in zip(x, y):
        t, r = pjr.triangle(
            float(a),
            float(b),
            north_square=north_square,
            south_square=south_square,
            inverse=inverse,
        )
        tri.append(-1 if t is None else t)
        region.append(REGION_CODE[r])
    return np.array(tri), np.array(region)


class HealpixArrayTestCase(unittest.TestCase):
    def test_in_healpix_image_array_matches_scalar(self):
        for x, y in (HEALPIX_PLANAR, HEALPIX_OUTSIDE, NOT_FINITE):
            got = pjh._in_healpix_image_array(x, y)
            want = np.array([pjh.in_healpix_image(a, b) for a, b in zip(x, y)])
            assert_array_equal(got, want)
        self.assertTrue(pjh._in_healpix_image_array(*HEALPIX_PLANAR).all())
        self.assertFalse(pjh._in_healpix_image_array(*HEALPIX_OUTSIDE).any())

    def test_cap_number_matches_scalar(self):
        # Clamped to [0, 3] at both ends (issue #119): inputs a hair past
        # -pi or pi stay in the outermost cap.
        past_ends = np.array([-pi - 1e-16, -pi - 1e-10, -4.0, pi + 1e-16, 4.0])
        for values in (LONLAT[0], HEALPIX_PLANAR[0], HEALPIX_OUTSIDE[0], past_ends):
            want = np.clip(np.floor(2 * values / pi + 2), 0.0, 3.0)
            assert_array_equal(pjh._cap_number(values), want)
        assert_array_equal(pjh._cap_number(past_ends), [0.0, 0.0, 0.0, 3.0, 3.0])

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


class RhealpixArrayTestCase(unittest.TestCase):
    def test_in_rhealpix_image_array_matches_scalar(self):
        for pair in SQUARE_PAIRS:
            for x, y in (RHEALPIX_PLANAR[pair], RHEALPIX_OUTSIDE[pair], NOT_FINITE):
                got = pjr._in_rhealpix_image_array(x, y, *pair)
                want = np.array(
                    [pjr.in_rhealpix_image(a, b, *pair) for a, b in zip(x, y)]
                )
                assert_array_equal(got, want)
            self.assertTrue(
                pjr._in_rhealpix_image_array(*RHEALPIX_PLANAR[pair], *pair).all()
            )
            self.assertFalse(
                pjr._in_rhealpix_image_array(*RHEALPIX_OUTSIDE[pair], *pair).any()
            )

    def test_triangle_array_matches_scalar(self):
        for pair in SQUARE_PAIRS:
            tri, region = pjr._triangle_array(*HEALPIX_PLANAR, *pair, inverse=False)
            want_tri, want_region = triangle_codes(
                *HEALPIX_PLANAR, *pair, inverse=False
            )
            assert_array_equal(tri, want_tri)
            assert_array_equal(region, want_region)
            tri, region = pjr._triangle_array(
                *RHEALPIX_PLANAR[pair], *pair, inverse=True
            )
            want_tri, want_region = triangle_codes(
                *RHEALPIX_PLANAR[pair], *pair, inverse=True
            )
            assert_array_equal(tri, want_tri)
            assert_array_equal(region, want_region)

    def test_combine_triangles_array_is_bit_identical(self):
        # Pure arithmetic on both sides: exact on every platform.
        for pair in SQUARE_PAIRS:
            got = pjr._combine_triangles_array(*HEALPIX_PLANAR, *pair)
            want = scalar_map(
                pjr.combine_triangles,
                *HEALPIX_PLANAR,
                north_square=pair[0],
                south_square=pair[1],
            )
            for g, w in zip(got, want):
                assert_array_equal(g, w)
            planar = RHEALPIX_PLANAR[pair]
            got = pjr._combine_triangles_array(*planar, *pair, inverse=True)
            want = scalar_map(
                pjr.combine_triangles,
                *planar,
                north_square=pair[0],
                south_square=pair[1],
                inverse=True,
            )
            for g, w in zip(got, want):
                assert_array_equal(g, w)
            # Round trip through the rearrangement on interior points (seam
            # points within the fuzz can come back through another triangle).
            hx, hy = (a[:N_INTERIOR_HEALPIX] for a in HEALPIX_PLANAR)
            fx, fy = pjr._combine_triangles_array(hx, hy, *pair)
            back = pjr._combine_triangles_array(fx, fy, *pair, inverse=True)
            for b, h in zip(back, (hx, hy)):
                assert_allclose(b, h, rtol=0, atol=1e-15)

    def test_rhealpix_sphere_arrays_match_scalar(self):
        for pair in SQUARE_PAIRS:
            kw = {"north_square": pair[0], "south_square": pair[1]}
            got = pjr._rhealpix_sphere_array(*LONLAT_SMALL, **kw)
            want = scalar_map(pjr.rhealpix_sphere, *LONLAT_SMALL, **kw)
            near = np.abs(np.abs(LONLAT_SMALL[1]) - pi / 2) < 1e-6
            for g, w in zip(got, want):
                assert_continuous(g, w, near_pole=near)
            planar = RHEALPIX_PLANAR[pair]
            got = pjr._rhealpix_sphere_inverse_array(*planar, **kw)
            want = scalar_map(pjr.rhealpix_sphere_inverse, *planar, **kw)
            near = np.abs(np.abs(want[1]) - pi / 2) < 1e-6
            for g, w in zip(got, want):
                assert_continuous(g, w, near_pole=near)

    def test_rhealpix_ellipsoid_arrays_match_scalar(self):
        for pair in ((0, 0), (1, 2), (3, 3)):
            kw = {"north_square": pair[0], "south_square": pair[1]}
            for e in E_VALUES:
                got = pjr._rhealpix_ellipsoid_array(*LONLAT_SMALL, e=e, **kw)
                want = scalar_map(pjr.rhealpix_ellipsoid, *LONLAT_SMALL, e=e, **kw)
                near = np.abs(np.abs(LONLAT_SMALL[1]) - pi / 2) < 1e-6
                for g, w in zip(got, want):
                    assert_continuous(g, w, near_pole=near)
                planar = RHEALPIX_PLANAR[pair]
                got = pjr._rhealpix_ellipsoid_inverse_array(*planar, e=e, **kw)
                want = scalar_map(pjr.rhealpix_ellipsoid_inverse, *planar, e=e, **kw)
                near = np.abs(np.abs(want[1]) - pi / 2) < 1e-6
                for g, w in zip(got, want):
                    assert_continuous(g, w, near_pole=near)

    def test_rhealpix_region_hint_arrays(self):
        # region="equatorial" skips the polar rearrangement for every point.
        lam, phi = LONLAT_SMALL
        for pair in SQUARE_PAIRS:
            kw = {
                "north_square": pair[0],
                "south_square": pair[1],
                "region": "equatorial",
            }
            hinted = pjr._rhealpix_sphere_array(lam, phi, **kw)
            plain = pjh._healpix_sphere_array(lam, phi)
            for h, q in zip(hinted, plain):
                assert_array_equal(h, q)
            hinted = pjr._rhealpix_ellipsoid_array(lam, phi, e=WGS84_E, **kw)
            plain = pjh._healpix_ellipsoid_array(lam, phi, e=WGS84_E)
            for h, q in zip(hinted, plain):
                assert_array_equal(h, q)
            # Inverse with the hint on equatorial planar points equals the
            # unhinted inverse there.
            x, y = HEALPIX_PLANAR
            band = (np.abs(y) <= pi / 4) & pjr._in_rhealpix_image_array(x, y, *pair)
            hinted = pjr._rhealpix_sphere_inverse_array(x[band], y[band], **kw)
            plain = pjr._rhealpix_sphere_inverse_array(
                x[band], y[band], north_square=pair[0], south_square=pair[1]
            )
            for h, q in zip(hinted, plain):
                assert_array_equal(h, q)

    def test_rhealpix_factory_arrays_match_scalar(self):
        cases = product(
            (1.0, 6378137.0), (0.0, WGS84_E), ((0, 0), (1, 2)), (True, False)
        )
        for a, e, pair, radians in cases:
            f = pjr.rhealpix(a=a, e=e, north_square=pair[0], south_square=pair[1])
            R_A = auth_rad(a, e)
            angle = 1.0 if radians else 180 / pi
            lam, phi = LONLAT_SMALL[0] * angle, LONLAT_SMALL[1] * angle
            near = np.abs(np.abs(LONLAT_SMALL[1]) - pi / 2) < 1e-6
            got = f(lam, phi, radians=radians)
            want = scalar_map(f, lam, phi, radians=radians)
            for g, w in zip(got, want):
                self.assertEqual(g.dtype, np.float64)
                assert_continuous(g, w, atol=ATOL * R_A, near_pole=near)
            # Scaling by R_A and back moves fuzz-edge seam points by an ulp,
            # which can flip their triangle; keep the points the scalar
            # closure accepts, as the array closure sees the same values.
            x, y = RHEALPIX_PLANAR[pair][0] * R_A, RHEALPIX_PLANAR[pair][1] * R_A
            keep = np.array([accepts(f, a, b, radians) for a, b in zip(x, y)])
            x, y = x[keep], y[keep]
            got = f(x, y, radians=radians, inverse=True)
            want = scalar_map(f, x, y, radians=radians, inverse=True)
            near = np.abs(np.abs(want[1]) - pi / 2 * angle) < 1e-6 * angle
            for g, w in zip(got, want):
                assert_continuous(g, w, atol=ATOL * angle, near_pole=near)
            out = f(float(lam[0]), float(phi[0]), radians=radians)
            self.assertIsInstance(out[0], np.float64)

    def test_rhealpix_round_trips(self):
        lam, phi = LONLAT_INTERIOR
        off_pole = np.abs(phi) < pi / 2
        near = np.abs(np.abs(phi) - pi / 2) < 1e-3
        for pair in ((0, 0), (1, 2), (3, 3)):
            kw = {"north_square": pair[0], "south_square": pair[1]}
            x, y = pjr._rhealpix_sphere_array(lam, phi, **kw)
            lam2, phi2 = pjr._rhealpix_sphere_inverse_array(x, y, **kw)
            assert_allclose(lam2[off_pole], lam[off_pole], rtol=0, atol=1e-12)
            assert_allclose(phi2, phi, rtol=0, atol=1e-12)
            for e, atol in ((WGS84_E, 1e-12), (0.8, 1e-3)):
                x, y = pjr._rhealpix_ellipsoid_array(lam, phi, e=e, **kw)
                lam2, phi2 = pjr._rhealpix_ellipsoid_inverse_array(x, y, e=e, **kw)
                assert_allclose(lam2[off_pole], lam[off_pole], rtol=0, atol=1e-12)
                assert_allclose(phi2[~near], phi[~near], rtol=0, atol=atol)
                assert_allclose(phi2[near], phi[near], rtol=0, atol=max(atol, 1e-6))

    def test_rhealpix_array_matches_proj_c_reference(self):
        lon, lat = LONLAT_SMALL[0] * 180 / pi, LONLAT_SMALL[1] * 180 / pi
        for pair in SQUARE_PAIRS:
            ref = Proj(proj="rhealpix", R=1, north_square=pair[0], south_square=pair[1])
            x_ref, y_ref = ref(lon, lat)
            x, y = pjr.rhealpix(a=1, e=0, north_square=pair[0], south_square=pair[1])(
                lon, lat
            )
            assert_allclose(x, x_ref, rtol=1e-12, atol=1e-12)
            assert_allclose(y, y_ref, rtol=1e-12, atol=1e-12)
        a, e, pair = 5.0, WGS84_E, (1, 2)
        f = pjr.rhealpix(a=a, e=e, north_square=pair[0], south_square=pair[1])
        x, y = f(lon, lat)
        lon2, lat2 = Proj(
            proj="rhealpix", a=a, e=e, north_square=pair[0], south_square=pair[1]
        )(x, y, inverse=True)
        gaps = [geo_gap_deg(p, q) for p, q in zip(zip(lon, lat), zip(lon2, lat2))]
        self.assertLess(max(gaps), 1e-5)

    def test_rhealpix_value_error_parity(self):
        rng = np.random.default_rng(SEED)
        for pair in ((0, 0), (1, 2), (2, 3)):
            good = RHEALPIX_PLANAR[pair]
            outside = RHEALPIX_OUTSIDE[pair]
            pick = rng.choice(len(outside[0]), 40, replace=False)
            kw = {"north_square": pair[0], "south_square": pair[1]}
            f = pjr.rhealpix(a=1, e=0, **kw)
            inverses = [
                lambda x, y, kw=kw: pjr._rhealpix_sphere_inverse_array(x, y, **kw),
                lambda x, y, kw=kw: pjr._rhealpix_ellipsoid_inverse_array(
                    x, y, e=WGS84_E, **kw
                ),
                lambda x, y, f=f: f(x, y, radians=True, inverse=True),
            ]
            for fn in inverses:
                fn(*good)
                for i in pick:
                    bx, by = outside[0][i], outside[1][i]
                    with self.assertRaises(ValueError):
                        pjr.rhealpix_sphere_inverse(float(bx), float(by), **kw)
                    x = np.concatenate([good[0][:5], [bx]])
                    y = np.concatenate([good[1][:5], [by]])
                    with self.assertRaises(ValueError) as cm:
                        fn(x, y)
                    self.assertIn("1 of 6 points", str(cm.exception))
                    self.assertIn(f"({pair[0]}, {pair[1]})-rHEALPix", str(cm.exception))
                with self.assertRaises(ValueError):
                    fn(*NOT_FINITE)
            # Points the rHEALPix check admits but the HEALPix check rejects
            # after uncombining raise in both paths.
            rejected = RHEALPIX_REJECTED[pair]
            self.assertGreater(len(rejected[0]), 0)
            with self.assertRaises(ValueError) as cm:
                pjr._rhealpix_sphere_inverse_array(*rejected, **kw)
            self.assertIn("HEALPix projection", str(cm.exception))
            for a, b in zip(*rejected):
                with self.assertRaises(ValueError):
                    pjr.rhealpix_sphere_inverse(float(a), float(b), **kw)


if __name__ == "__main__":
    unittest.main()
