r"""
This Python 2.7 code tests the ``rhealpix_projection.py`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!

CHANGELOG:

- Alexander Raichev (AR), 2012-11-21: Initial version based on a previous set of tests in test_projections.py.
- AR, 2013-01-15: Added a check that degrees mode yields the same output as radians mode for healpix() and rhealpix().
"""
#*****************************************************************************
#       Copyright (C) 2012 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Import standard modules.
import unittest
# Import my modules.
import projection_tools as pt


class MyTestCase(unittest.TestCase):
    def test_wrap_longitude(self):
        angles = [0, 90, -90, 180, -180, 360, -360]
        get = list(map(pt.wrap_longitude, angles))
        expect = [0, 90, -90, -180, -180, 0, 0]
        for i in range(len(expect)):    
            self.assertAlmostEqual(get[i], expect[i])
        
    def test_wrap_latitude(self):
        angles = [0, 45, -45, 90, -90, 135, -135]
        get = list(map(pt.wrap_latitude, angles))
        expect = [0, 45, -45, 90, -90, -45, 45]
        for i in range(len(expect)):    
            self.assertAlmostEqual(get[i], expect[i])
#------------------------------------------------------------------------------                           
if __name__ == '__main__':
    unittest.main()