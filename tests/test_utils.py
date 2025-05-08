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
from math import pi
import unittest

# Import my modules.
from rhealpixdggs.ellipsoids import WGS84_E
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


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
