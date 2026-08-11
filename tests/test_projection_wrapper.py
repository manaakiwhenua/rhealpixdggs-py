"""
This Python 3.11 code tests the ``projection_wrapper`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!
"""

# *****************************************************************************
#       Copyright (C) 2026
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import standard modules.
import importlib
import unittest
from unittest.mock import patch

# Import my modules.
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from rhealpixdggs.projection_wrapper import HOMEMADE_PROJECTIONS, Projection


class ProjectionWrapperTestCase(unittest.TestCase):
    def test_implemented_homemade_projection_still_works(self):
        # rhealpix has a real rhealpixdggs.pj_rhealpix module backing it, so
        # this must still project normally and not hit the fallback branch.
        f = Projection(
            ellipsoid=WGS84_ELLIPSOID, proj="rhealpix", north_square=1, south_square=0
        )
        result = f(0, 30)
        self.assertAlmostEqual(float(result[0]), 0.0, places=6)
        self.assertAlmostEqual(float(result[1]), 3740232.8933662786, places=6)

    def test_unimplemented_homemade_projections_degrade_gracefully(self):
        """#50: isea/csea/qsc are declared in HOMEMADE_PROJECTIONS but have no
        backing rhealpixdggs/pj_*.py module. Projection.__call__ used to catch
        only NameError while importlib.import_module() actually raises
        ModuleNotFoundError -- so these crashed uncaught instead of printing
        the intended "not implemented" message and returning None.
        """
        unimplemented = HOMEMADE_PROJECTIONS - {"healpix", "rhealpix"}
        # Confirm the fixture assumption still holds: these really have no
        # backing module, i.e. this test is actually exercising the fallback.
        self.assertTrue(
            unimplemented, "expected at least one unimplemented homemade projection"
        )

        for proj in sorted(unimplemented):
            with self.subTest(proj=proj):
                f = Projection(ellipsoid=WGS84_ELLIPSOID, proj=proj)
                # Must not raise ModuleNotFoundError or any other uncaught
                # exception -- must degrade to None per the documented
                # fallback behaviour.
                result = f(0, 0)
                self.assertIsNone(result)

    def test_module_present_but_missing_attribute_also_degrades_gracefully(self):
        """The AttributeError half of the fix: a homemade-projection module
        that exists but doesn't define the expected same-named callable must
        also degrade gracefully rather than raising.
        """
        real_module = importlib.import_module("rhealpixdggs.pj_rhealpix")
        with patch(
            "rhealpixdggs.projection_wrapper.importlib.import_module",
            return_value=real_module,
        ):
            # pj_rhealpix has no attribute named "rhealpix_missing_fn", so
            # getattr(module, proj) raises AttributeError.
            f = Projection(ellipsoid=WGS84_ELLIPSOID, proj="rhealpix_missing_fn")
            HOMEMADE_PROJECTIONS.add("rhealpix_missing_fn")
            try:
                result = f(0, 0)
            finally:
                HOMEMADE_PROJECTIONS.discard("rhealpix_missing_fn")
            self.assertIsNone(result)


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
