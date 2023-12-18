"""
This Python 3.11 code tests the ``PJ_healpix.c`` code of the third party
pyproj module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!

- Alexander Raichev (AR), 2012-11-20: Initial version.
"""
# *****************************************************************************
#       Copyright (C) 2012 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import third-party modules.
from scipy.spatial.distance import euclidean, norm
from scipy import array, sin, arcsin, pi, sqrt, rad2deg, deg2rad
from pyproj import Proj

# Import standard modules.
import unittest

# Import my modules.
from rhealpixdggs.utils import auth_lat, auth_rad


# Relative error function.
def rel_err(get, expect):
    a = euclidean(get, expect)
    b = norm(expect)
    if b == 0:
        return a
    else:
        return a / b


# Define lon-lat input points to test.
RADIANS = False  # Work in radians (True) or degrees (False)?
if RADIANS:
    PI = pi
else:
    PI = 180
phi_0 = arcsin(2.0 / 3)
if not RADIANS:
    phi_0 = rad2deg(phi_0)
a = (0, PI / 3)
b = (0, -PI / 3)
inputs = [
    (0, 0),
    (0, phi_0),
    (0, -phi_0),
    (PI / 2, 0),
    (-PI / 2, 0),
    (-PI, 0),
    (-PI, PI / 2),
    (-PI, -PI / 2),
    a,
    b,
]


class PJHEALPixTestCase(unittest.TestCase):
    def test_healpix_sphere(self):
        # Sphere parameters.
        R = 5
        # Expected outputs of healpix_sphere() applied to inputs.
        if RADIANS:
            sigma_a = sqrt(3 - 3 * sin(a[1]))
        else:
            sigma_a = sqrt(3 - 3 * sin(deg2rad(a[1])))
        ha = (pi / 4 * (1 - sigma_a), pi / 4 * (2 - sigma_a))
        hb = (ha[0], -ha[1])
        healpix_sphere_outputs = [
            (0, 0),
            (0, pi / 4),
            (0, -pi / 4),
            (pi / 2, 0),
            (-pi / 2, 0),
            (-pi, 0),
            (-3 * pi / 4, pi / 2),
            (-3 * pi / 4, -pi / 2),
            ha,
            hb,
        ]
        healpix_sphere_outputs = [tuple(R * array(p)) for p in healpix_sphere_outputs]

        # Forward projection should be correct on test points.
        f = Proj(proj="healpix", R=R)
        given = inputs
        get = [f(*p, radians=RADIANS) for p in given]
        expect = healpix_sphere_outputs
        # Fuzz to allow for rounding errors:
        error = 1e-12
        print()
        print("=" * 80)
        print("HEALPix forward projection, sphere with radius R = %s" % R)
        print("input (radians) / expected output (meters) / received output")
        print("=" * 80)
        for i in range(len(get)):
            print(given[i], expect[i], get[i])
            self.assertTrue(rel_err(get[i], expect[i]) < error)

        # Inverse of projection of a point p should yield p.
        given = get
        get = [f(*q, radians=RADIANS, inverse=True) for q in given]
        expect = inputs
        print("=" * 80)
        print("HEALPix inverse projection, sphere with radius R = %s" % R)
        print("input (meters) / expected output (radians) / received output")
        print("=" * 80)
        for i in range(len(get)):
            print(given[i], expect[i], get[i])
            self.assertTrue(rel_err(get[i], expect[i]) < error)

        # Inverse projection of p below should return longitude of -pi.
        # Previously, it was returning a number slightly less than pi
        # because of a rounding error, which got magnified by
        # wrap_longitude()
        p = R * array((-7 * pi / 8, 3 * pi / 8))
        get = f(*p, radians=RADIANS, inverse=True)
        p1 = arcsin(1 - 1.0 / 12)
        if not RADIANS:
            p1 = rad2deg(p1)
        expect = (-PI, p1)
        self.assertTrue(rel_err(get, expect) < error)

    def test_healpix_ellipsoid(self):
        # Ellipsoid parameters.
        a = 5
        e = 0.8
        R_A = auth_rad(a, e=e)
        # Expected output of healpix_ellipsoid() applied to inputs.
        healpix_ellipsoid_outputs = []
        g = Proj(proj="healpix", R=R_A)
        for p in inputs:
            lam, phi = p
            beta = auth_lat(phi, e=e, radians=RADIANS)
            q = g(lam, beta, radians=RADIANS)
            healpix_ellipsoid_outputs.append(q)

        # Forward projection should be correct on test points.
        f = Proj(proj="healpix", a=a, e=e)
        given = inputs
        get = [f(*p, radians=RADIANS) for p in given]
        expect = healpix_ellipsoid_outputs
        # Fuzz to allow for rounding errors:
        error = 1e-12
        print("=" * 80)
        print(
            "HEALPix forward projection, ellipsoid with major radius a = %s and eccentricity e = %s"
            % (a, e)
        )
        print("input (radians) / expected output (meters) / received output")
        print("=" * 80)
        for i in range(len(get)):
            print(given[i], expect[i], get[i])
            self.assertTrue(rel_err(get[i], expect[i]) < error)

        # Inverse of projection of a point p should yield p.
        given = get
        get = [f(*q, radians=RADIANS, inverse=True) for q in given]
        expect = inputs
        # Fuzz for rounding errors based on the error of the approximation to
        # the inverse authalic latitude function:
        alpha = PI / 4
        alpha_ = auth_lat(
            auth_lat(alpha, e, radians=RADIANS), e, radians=RADIANS, inverse=True
        )
        error = 10 * rel_err(alpha_, alpha)
        print("=" * 80)
        print(
            "HEALPix inverse projection, ellipsoid with major radius a = %s and eccentricity e = %s"
            % (a, e)
        )
        print("input (meters) / expected output (radians) / received output")
        print("=" * 80)
        for i in range(len(get)):
            print(given[i], expect[i], get[i])
            self.assertTrue(rel_err(get[i], expect[i]) < error)

    def test_rhealpix_sphere(self):
        from random import uniform

        # Sphere parameters.
        R = 5
        # Fuzz to allow for rounding errors:
        error = 1e-12
        # Forward projection of random equatorial points should yield the same
        # output as healpix_sphere().
        print("=" * 80)
        print("rHEALPix forward projection, sphere of radius R = %s" % R)
        print("input (radians) / expected output (meters) / received output")
        print("=" * 80)
        eps = 1e-3
        given = [(-PI, PI / 6), (-PI, -PI / 7), (0, 0), (PI / 3, phi_0 - eps)]
        # for i in range(10):
        #     p = (uniform(-pi + eps, pi - eps),
        #          uniform(-phi_0 + eps, phi_0 - eps))
        #     given.append(p)
        h = Proj(proj="healpix", R=R)
        expect = [h(*p, radians=RADIANS) for p in given]
        for ns, ss in product(list(range(4)), repeat=2):
            print("_____ north_square = %s, south_square = %s" % (ns, ss))
            rh = Proj(proj="rhealpix", R=R, north_square=ns, south_square=ss)
            get = [rh(*p, radians=RADIANS) for p in given]
            for i in range(len(get)):
                print(given[i], expect[i], get[i])
                self.assertEqual(get[i], expect[i])

        # Forward projection of polar points should be correct.
        dx = pi / 4
        dy = pi / 12
        north_healpix_output = [
            R * array((-pi + i * pi / 2 + dx, pi / 4 + dy)) for i in range(4)
        ]
        south_healpix_output = [
            R * array((-pi + i * pi / 2 + dx, -pi / 4 - dy)) for i in range(4)
        ]
        h = Proj(proj="healpix", R=R)
        north_given = [
            h(*p, radians=RADIANS, inverse=True) for p in north_healpix_output
        ]
        south_given = [
            h(*p, radians=RADIANS, inverse=True) for p in south_healpix_output
        ]
        for ns, ss in product(list(range(4)), repeat=2):
            print("_____ north_square = %s, south_square = %s" % (ns, ss))
            rh = Proj(proj="rhealpix", R=R, north_square=ns, south_square=ss)
            # Corners of north square.
            ndl = (-pi + ns * pi / 2, pi / 4)
            ndr = (-pi / 2 + ns * pi / 2, pi / 4)
            nur = (-pi / 2 + ns * pi / 2, 3 * pi / 4)
            nul = (-pi + ns * pi / 2, 3 * pi / 4)
            # Corners of south square.
            sdl = (-pi + ss * pi / 2, -3 * pi / 4)
            sdr = (-pi / 2 + ss * pi / 2, -3 * pi / 4)
            sur = (-pi / 2 + ss * pi / 2, -pi / 4)
            sul = (-pi + ss * pi / 2, -pi / 4)
            # north_given should map to the following points,
            # though not necessarily in the order given.
            north_expect = [
                (ndl[0] + dx, ndl[1] + dy),
                (ndr[0] - dy, ndr[1] + dx),
                (nur[0] - dx, nur[1] - dy),
                (nul[0] + dy, nul[1] - dx),
            ]
            north_expect = [tuple(R * array(p)) for p in north_expect]
            # south_given should map to the following points,
            # though not necessarily in the order given.
            south_expect = [
                (sul[0] + dx, sul[1] - dy),
                (sur[0] - dy, sur[1] - dx),
                (sdr[0] - dx, sdr[1] + dy),
                (sdl[0] + dy, sdl[1] + dx),
            ]
            south_expect = [tuple(R * array(p)) for p in south_expect]
            for i, p in enumerate(north_given):
                get = rh(*p, radians=RADIANS)
                expect = north_expect[(i - ns) % 4]
                print(p, expect, get)
                self.assertTrue(rel_err(get, expect) < error)
            for i, p in enumerate(south_given):
                get = rh(*p, radians=RADIANS)
                expect = south_expect[(i - ss) % 4]
                print(p, expect, get)
                self.assertTrue(rel_err(get, expect) < error)

        # The inverse of the projection of a point p should yield p.
        print("=" * 80)
        print("rHEALPix inverse projection, sphere of radius R = %s" % R)
        print("input (meters) / expected output (radians) / received output:")
        print("=" * 80)
        for ns, ss in product(list(range(4)), repeat=2):
            print("_____ north_square = %s, south_square = %s" % (ns, ss))
            f = Proj(proj="rhealpix", R=R, north_square=ns, south_square=ss)
            for p in inputs:
                expect = p
                q = f(*p, radians=RADIANS)
                get = f(*q, radians=RADIANS, inverse=True)
                print(q, expect, get)
                self.assertTrue(rel_err(get, expect) < error)

    def test_rhealpix_ellipsoid(self):
        # Ellipsoid parameters.
        a = 5
        e = 0.8
        R_A = auth_rad(a, e=e)
        # Forward projection should be correct on test points.
        print("=" * 80)
        print(
            "rHEALPix forward projection, ellipsoid with major radius a = %s and eccentricity e = %s"
            % (a, e)
        )
        print("input (radians) / expected output (meters) / received output")
        print("=" * 80)
        given = inputs
        # Fuzz to allow for rounding errors:
        error = 1e-12
        for ns, ss in product(list(range(4)), repeat=2):
            print("_____ north_square = %s, south_square = %s" % (ns, ss))
            expect = []
            g = Proj(proj="rhealpix", R=R_A, north_square=ns, south_square=ss)
            for p in given:
                lam, phi = p
                beta = auth_lat(phi, e=e, radians=RADIANS)
                q = g(lam, beta, radians=RADIANS)
                expect.append(tuple(q))
            f = Proj(proj="rhealpix", a=a, e=e, north_square=ns, south_square=ss)
            get = [f(*p, radians=RADIANS) for p in given]
            for i in range(len(given)):
                print(given[i], expect[i], get[i])
                self.assertTrue(rel_err(get[i], expect[i]) < error)

        # Inverse of projection of point a p should yield p.
        # Fuzz for rounding errors based on the error of the approximation to
        # the inverse authalic latitude function:
        alpha = PI / 4
        alpha_ = auth_lat(
            auth_lat(alpha, e, radians=RADIANS), e, radians=RADIANS, inverse=True
        )
        error = 10 * rel_err(alpha_, alpha)
        print("=" * 80)
        print(
            "HEALPix inverse projection, ellipsoid with major radius a = %s and eccentricity e = %s"
            % (a, e)
        )
        print("input (meters) / expected output (radians) / received output")
        print("=" * 80)
        # The inverse of the projection of a point p should yield p.
        for ns, ss in product(list(range(4)), repeat=2):
            print("_____ north_square = %s, south_square = %s" % (ns, ss))
            f = Proj(proj="rhealpix", a=a, e=e, north_square=ns, south_square=ss)
            for p in inputs:
                expect = p
                q = f(*p, radians=RADIANS)
                get = f(*q, radians=RADIANS, inverse=True)
                print(q, expect, get)
                self.assertTrue(rel_err(get, expect) < error)


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
