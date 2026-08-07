"""
This Python 3.11 code tests the ``projection_wrapper`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!
"""

# Import standard modules.
import unittest

# Import my modules.
from rhealpixdggs.projection_wrapper import Projection
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID


class MyTestCase(unittest.TestCase):
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
