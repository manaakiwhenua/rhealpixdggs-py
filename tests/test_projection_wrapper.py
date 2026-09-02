"""
This Python 3.11 code tests the ``projection_wrapper`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!
"""

import unittest

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


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
