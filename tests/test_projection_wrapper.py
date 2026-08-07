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


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
