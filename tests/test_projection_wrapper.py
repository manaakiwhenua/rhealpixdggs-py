"""
This Python 3.11 code tests the ``projection_wrapper`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!
"""

import unittest

import numpy as np
from numpy.testing import assert_allclose

from rhealpixdggs.ellipsoids import (
    WGS84_ELLIPSOID,
    WGS84_ELLIPSOID_RADIANS,
    Ellipsoid,
)
from rhealpixdggs.projection_wrapper import Projection


class MyTestCase(unittest.TestCase):
    def test_homemade_projection_round_trip(self):
        # The homemade healpix/rhealpix projections, called through the
        # wrapper, should invert cleanly: f(p) followed by the inverse
        # should recover p. Exercised in degrees mode (the wrapper reads
        # angle units off the ellipsoid).
        points = [(0, 0), (30, 45), (-120, -60), (170, 80)]
        for proj, kwargs in [
            ("healpix", {}),
            ("rhealpix", {"north_square": 1, "south_square": 2}),
        ]:
            f = Projection(ellipsoid=WGS84_ELLIPSOID, proj=proj, **kwargs)
            for p in points:
                q = f(*p)
                back = f(*q, inverse=True)
                self.assertAlmostEqual(back[0], p[0], places=9, msg=proj)
                self.assertAlmostEqual(back[1], p[1], places=9, msg=proj)

    def test_homemade_projection_round_trip_radians(self):
        # Same as above in radians mode.
        f = Projection(
            ellipsoid=WGS84_ELLIPSOID_RADIANS,
            proj="rhealpix",
            north_square=0,
            south_square=0,
        )
        for p in [(0.5, 0.7), (-2.0, -1.0), (3.0, 1.4)]:
            q = f(*p)
            back = f(*q, inverse=True)
            self.assertAlmostEqual(back[0], p[0], places=12)
            self.assertAlmostEqual(back[1], p[1], places=12)

    def test_pyproj_backed_projection_round_trip(self):
        # A projection not in HOMEMADE_PROJECTIONS dispatches to pyproj.
        # 'cea' (cylindrical equal area) is a stable, universally available
        # PROJ projection. pyproj's own round-trip precision is coarser
        # than the homemade projections' (observed ~5e-9 degrees), hence
        # the looser tolerance.
        f = Projection(ellipsoid=WGS84_ELLIPSOID, proj="cea")
        for p in [(0, 0), (30, 45), (-120, -60)]:
            q = f(*p)
            if p != (0, 0):
                # A real projection to meters happened (not a pass-through).
                self.assertNotEqual(q, p)
            back = f(*q, inverse=True)
            self.assertAlmostEqual(back[0], p[0], delta=1e-6)
            self.assertAlmostEqual(back[1], p[1], delta=1e-6)

    def test_ellipsoid_origin_translation(self):
        # The wrapper recentres on the ellipsoid's (lon_0, lat_0): that
        # point must project to the planar origin, and the planar origin
        # must invert back to it.
        E = Ellipsoid(a=WGS84_ELLIPSOID.a, f=WGS84_ELLIPSOID.f, lon_0=50, lat_0=0)
        f = Projection(ellipsoid=E, proj="rhealpix", north_square=0, south_square=0)
        self.assertEqual(f(50, 0), (0.0, 0.0))
        self.assertEqual(f(0, 0, inverse=True), (50.0, 0.0))

    def test_projection_kwargs_are_passed_through(self):
        # rhealpix's north_square/south_square kwargs reposition the polar
        # squares, so they must change where a polar point lands while
        # leaving equatorial points untouched.
        f0 = Projection(
            ellipsoid=WGS84_ELLIPSOID, proj="rhealpix", north_square=0, south_square=0
        )
        f1 = Projection(
            ellipsoid=WGS84_ELLIPSOID, proj="rhealpix", north_square=1, south_square=0
        )
        self.assertNotEqual(f0(10, 80), f1(10, 80))
        self.assertEqual(f0(10, 20), f1(10, 20))

    def test_unimplemented_homemade_projection_degrades_gracefully(self):
        # 'isea', 'csea', and 'qsc' are listed in HOMEMADE_PROJECTIONS but have
        # no backing rhealpixdggs.pj_* module. Calling them should hit the
        # "Oops! Projection %s is not implemented." branch (returning None),
        # not crash with an uncaught ModuleNotFoundError/AttributeError.
        # Regression test for the `except NameError` bug (see issue #50).
        for proj in ("isea", "csea", "qsc"):
            f = Projection(ellipsoid=WGS84_ELLIPSOID, proj=proj)
            self.assertIsNone(f(0, 0))

    def test_underlying_callable_is_built_once_and_cached(self):
        # Regression test: __call__() used to rebuild its underlying f
        # (re-importing its module and, for homemade projections,
        # recomputing the authalic radius) on every single call, even
        # though it only depends on the ellipsoid's a/e, which never
        # change after construction. See issue #62.
        p = Projection(ellipsoid=WGS84_ELLIPSOID, proj="rhealpix")
        self.assertIsNone(p._f)
        p(0, 0)
        cached = p._f
        self.assertIsNotNone(cached)
        p(0, 0.1)
        p(0.2, -0.1, inverse=True)
        self.assertIs(p._f, cached)

    def test_array_inputs_match_scalar_calls(self):
        # Arrays go through the array path and must agree with a scalar call
        # per point, forward and inverse, in both angle modes.
        for ellipsoid, points, atol_angle in (
            (
                WGS84_ELLIPSOID,
                [(0, 0), (30, 45), (-120, -60), (170, 80), (-180, 0)],
                1e-12,
            ),
            (
                WGS84_ELLIPSOID_RADIANS,
                [(0.5, 0.7), (-2.0, -1.0), (3.0, 1.4), (-3.1, 0.0)],
                1e-14,
            ),
        ):
            for proj, kwargs in [
                ("healpix", {}),
                ("rhealpix", {"north_square": 1, "south_square": 2}),
            ]:
                f = Projection(ellipsoid=ellipsoid, proj=proj, **kwargs)
                u = np.array([p[0] for p in points], dtype=float)
                v = np.array([p[1] for p in points], dtype=float)
                x, y = f(u, v)
                self.assertEqual((x.dtype, y.dtype), (np.float64, np.float64))
                scalar = [f(*p) for p in points]
                assert_allclose(x, [q[0] for q in scalar], rtol=0, atol=1e-8)
                assert_allclose(y, [q[1] for q in scalar], rtol=0, atol=1e-8)
                lam, phi = f(x, y, inverse=True)
                back = [f(*q, inverse=True) for q in scalar]
                assert_allclose(lam, [b[0] for b in back], rtol=0, atol=atol_angle)
                assert_allclose(phi, [b[1] for b in back], rtol=0, atol=atol_angle)
                # A scalar mixed with an array broadcasts.
                x2, _ = f(u, np.float64(points[1][1]) + np.zeros(len(points)))
                self.assertEqual(x2.shape, (len(points),))

    def test_array_round_trip(self):
        rng = np.random.default_rng(20260811)
        lon = rng.uniform(-180, 180, 500)
        lat = np.degrees(np.arcsin(rng.uniform(-1, 1, 500)))
        f = Projection(ellipsoid=WGS84_ELLIPSOID, proj="rhealpix", north_square=2)
        x, y = f(lon, lat)
        lon2, lat2 = f(x, y, inverse=True)
        assert_allclose(lon2, lon, rtol=0, atol=1e-9)
        assert_allclose(lat2, lat, rtol=0, atol=1e-9)

    def test_pyproj_backed_projection_accepts_arrays(self):
        f = Projection(ellipsoid=WGS84_ELLIPSOID, proj="cea")
        lon = np.array([0.0, 30.0, -120.0])
        lat = np.array([0.0, 45.0, -60.0])
        x, y = f(lon, lat)
        self.assertIsInstance(x, np.ndarray)
        lon2, lat2 = f(x, y, inverse=True)
        assert_allclose(lon2, lon, rtol=0, atol=1e-6)
        assert_allclose(lat2, lat, rtol=0, atol=1e-6)

    def test_ellipsoid_origin_translation_arrays(self):
        # The recentring on (lon_0, lat_0) and the wrap back into range apply
        # elementwise; a longitude that crosses the antimeridian after the
        # shift wraps exactly as the scalar call does.
        E = Ellipsoid(a=WGS84_ELLIPSOID.a, f=WGS84_ELLIPSOID.f, lon_0=50, lat_0=0)
        f = Projection(ellipsoid=E, proj="rhealpix", north_square=0, south_square=0)
        x, y = f(np.array([50.0, -170.0]), np.array([0.0, 20.0]))
        self.assertEqual((x[0], y[0]), (0.0, 0.0))
        self.assertEqual((x[1], y[1]), f(-170.0, 20.0))
        lon, lat = f(np.array([0.0, x[1]]), np.array([0.0, y[1]]), inverse=True)
        self.assertEqual((lon[0], lat[0]), (50.0, 0.0))
        self.assertEqual((lon[1], lat[1]), f(x[1], y[1], inverse=True))

    def test_unimplemented_homemade_projection_returns_none_for_arrays(self):
        f = Projection(ellipsoid=WGS84_ELLIPSOID, proj="isea")
        self.assertIsNone(f(np.array([0.0]), np.array([0.0])))

    def test_underlying_callable_is_shared_by_scalar_and_array_calls(self):
        p = Projection(ellipsoid=WGS84_ELLIPSOID, proj="rhealpix")
        p(np.array([0.0]), np.array([0.0]))
        cached = p._f
        p(0, 0.1)
        p(np.array([0.2]), np.array([-0.1]), inverse=True)
        self.assertIs(p._f, cached)


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
