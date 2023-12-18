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

# Import standard modules.
import unittest

# Import my modules.
import rhealpixdggs.utils as ut


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


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
