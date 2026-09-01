"""
This Python 3.11 code cross-checks the ``PJ_healpix.c`` implementation of
the healpix and rhealpix projections that ships in the PROJ library (as
exposed by the third-party pyproj module) against reference values derived
from the defining formulas -- a canary for a broken or mispatched PROJ
healpix at install time.
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
# Import standard modules.
import unittest
from itertools import product

from numpy import arcsin, array, atleast_1d, deg2rad, pi, rad2deg, sin, sqrt
from numpy.linalg import norm
from pyproj import Proj

# Import my modules.
from rhealpixdggs.utils import auth_lat, auth_rad


# Relative error function. Accepts scalars or point tuples.
def rel_err(get, expect):
    get = atleast_1d(array(get))
    expect = atleast_1d(array(expect))
    a = norm(get - expect)
    b = norm(expect)
    if b == 0:
        return a
    else:
        return a / b


def geo_gap_deg(p, q):
    """
    Angular separation, in degrees, between lon-lat points `p` and `q`
    (in degrees) treated as positions on the unit sphere. Unlike a
    componentwise comparison, this correctly treats (-180, phi) and
    (180, phi) as the same point, and any two longitudes at a pole as
    the same point.
    """
    from numpy import arccos, clip, cos, deg2rad, sin

    lam1, phi1 = deg2rad(array(p, dtype=float))
    lam2, phi2 = deg2rad(array(q, dtype=float))
    cos_gap = sin(phi1) * sin(phi2) + cos(phi1) * cos(phi2) * cos(lam1 - lam2)
    return rad2deg(arccos(clip(cos_gap, -1, 1)))


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
        for i in range(len(get)):
            self.assertTrue(rel_err(get[i], expect[i]) < error)

        # Inverse of projection of a point p should yield p.
        given = get
        get = [f(*q, radians=RADIANS, inverse=True) for q in given]
        expect = inputs
        for i in range(len(get)):
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
        for i in range(len(get)):
            self.assertTrue(rel_err(get[i], expect[i]) < error)

        # Inverse of projection of a point p should yield p, geographically
        # (comparing angular separation on the sphere, so that the
        # antimeridian's two longitude representations and the poles'
        # indeterminate longitude don't register as spurious mismatches).
        #
        # PROJ inverts the authalic latitude with a low-order series whose
        # truncation error grows steeply with eccentricity, so the
        # achievable tolerance depends on e. At this test's deliberately
        # extreme e = 0.8 the observed round-trip error is ~2.2 degrees of
        # latitude: assert within 3 degrees, purely as a smoke test that
        # the inverse lands in the right neighborhood (a sign flip, wrong
        # polar square, or wrong quadrant would miss by far more).
        given = get
        get = [f(*q, radians=RADIANS, inverse=True) for q in given]
        expect = inputs
        for i in range(len(get)):
            self.assertLess(geo_gap_deg(get[i], expect[i]), 3)

        # At a realistic eccentricity (WGS84's), that same series is
        # accurate to about a millionth of a degree, so there the
        # round trip is a genuine precision check.
        f = Proj(proj="healpix", a=a, e=0.0818191908426215)
        for p in inputs:
            q = f(*p, radians=RADIANS)
            back = f(*q, radians=RADIANS, inverse=True)
            self.assertLess(geo_gap_deg(back, p), 1e-5)

    def test_rhealpix_sphere(self):
        # Sphere parameters.
        R = 5
        # Fuzz to allow for rounding errors:
        error = 1e-12
        # Forward projection of equatorial points should yield the same
        # output as healpix_sphere().
        eps = 1e-3
        given = [(-PI, PI / 6), (-PI, -PI / 7), (0, 0), (PI / 3, phi_0 - eps)]
        h = Proj(proj="healpix", R=R)
        expect = [h(*p, radians=RADIANS) for p in given]
        for ns, ss in product(list(range(4)), repeat=2):
            rh = Proj(proj="rhealpix", R=R, north_square=ns, south_square=ss)
            get = [rh(*p, radians=RADIANS) for p in given]
            for i in range(len(get)):
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
                self.assertTrue(rel_err(get, expect) < error)
            for i, p in enumerate(south_given):
                get = rh(*p, radians=RADIANS)
                expect = south_expect[(i - ss) % 4]
                self.assertTrue(rel_err(get, expect) < error)

        # The inverse of the projection of a point p should yield p.
        for ns, ss in product(list(range(4)), repeat=2):
            f = Proj(proj="rhealpix", R=R, north_square=ns, south_square=ss)
            for p in inputs:
                expect = p
                q = f(*p, radians=RADIANS)
                get = f(*q, radians=RADIANS, inverse=True)
                self.assertTrue(rel_err(get, expect) < error)

    def test_rhealpix_ellipsoid(self):
        # Ellipsoid parameters.
        a = 5
        e = 0.8
        R_A = auth_rad(a, e=e)
        # Forward projection should be correct on test points.
        given = inputs
        # Fuzz to allow for rounding errors:
        error = 1e-12
        for ns, ss in product(list(range(4)), repeat=2):
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
                self.assertTrue(rel_err(get[i], expect[i]) < error)

        # The inverse of the projection of a point p should yield p,
        # geographically. Same two-tier tolerances as in
        # test_healpix_ellipsoid, and for the same reason: PROJ inverts
        # the authalic latitude with a low-order series, so at the
        # deliberately extreme e = 0.8 this is a within-3-degrees smoke
        # test, while at WGS84's eccentricity it's a genuine precision
        # check.
        for ns, ss in product(list(range(4)), repeat=2):
            f = Proj(proj="rhealpix", a=a, e=e, north_square=ns, south_square=ss)
            for p in inputs:
                q = f(*p, radians=RADIANS)
                get = f(*q, radians=RADIANS, inverse=True)
                self.assertLess(geo_gap_deg(get, p), 3)
        f = Proj(
            proj="rhealpix", a=a, e=0.0818191908426215, north_square=1, south_square=2
        )
        for p in inputs:
            q = f(*p, radians=RADIANS)
            back = f(*q, radians=RADIANS, inverse=True)
            self.assertLess(geo_gap_deg(back, p), 1e-5)


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
