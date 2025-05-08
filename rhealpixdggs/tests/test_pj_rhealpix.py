"""
This Python 3.11 code tests the ``pj_rhealpix`` module.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!

- Alexander Raichev (AR), 2013-01-26: Initial version based on previous test code.
"""
# *****************************************************************************
#       Copyright (C) 2013 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import third-party modules.
from scipy.spatial.distance import euclidean, norm
from numpy import atleast_1d, array, pi, rad2deg, deg2rad, arcsin

# Import standard modules.
import unittest
from itertools import product

# Import my modules.
import rhealpixdggs.pj_healpix as pjh
import rhealpixdggs.pj_rhealpix as pjr
from rhealpixdggs.utils import auth_lat, auth_rad


# Relative error function.
def rel_err(get, expect):
    a = euclidean(atleast_1d(get), atleast_1d(expect))
    b = norm(expect)
    if b == 0:
        return a
    else:
        return a / b


# Define lon-lat input points.
phi_0 = arcsin(2.0 / 3)
a = (0, pi / 3)
b = (0, -pi / 3)
inputs = [
    (0, 0),
    (0, phi_0),
    (0, -phi_0),
    (pi / 2, 0),
    (-pi / 2, 0),
    (-pi, 0),
    (-pi, pi / 2),
    (-pi, -pi / 2),
    a,
    b,
]


class MyTestCase(unittest.TestCase):
    def test_triangle(self):
        # Create test points in the equatorial, north_square polar, and
        # south_square polar regions.
        dx = pi / 12
        dy = pi / 12
        eq_points = [(-pi + i * pi / 2 + dx, 0) for i in range(4)]
        np_points = [(-pi + i * pi / 2 + dx, pi / 4 + dy) for i in range(4)]
        sp_points = [(-pi + i * pi / 2 + dx, -pi / 4 - dy) for i in range(4)]
        np_points_inv = {}
        for north_square in range(4):
            u = array((-3 * pi / 4 + north_square * pi / 2, pi / 2))
            np_points_inv[north_square] = [
                u + array((0, -dy)),
                u + array((dx, 0)),
                u + array((0, dy)),
                u + array((-dx, 0)),
            ]
        sp_points_inv = {}
        for south_square in range(4):
            u = array((-3 * pi / 4 + south_square * pi / 2, -pi / 2))
            sp_points_inv[south_square] = [
                u + array((0, dy)),
                u + array((dx, 0)),
                u + array((0, -dy)),
                u + array((-dx, 0)),
            ]
        # Test forward.
        for north_square, south_square in product(list(range(4)), repeat=2):
            for p in eq_points:
                cr = pjr.triangle(
                    *p, north_square=north_square, south_square=south_square
                )
                self.assertEqual((None, "equatorial"), cr)
            for i, p in enumerate(np_points):
                cr = pjr.triangle(
                    *p, north_square=north_square, south_square=south_square
                )
                self.assertEqual((i, "north_polar"), cr)
            for i, p in enumerate(sp_points):
                cr = pjr.triangle(
                    *p, north_square=north_square, south_square=south_square
                )
                self.assertEqual((i, "south_polar"), cr)

        # Test inverse.
        for north_square, south_square in product(list(range(4)), repeat=2):
            for p in eq_points:
                cr = pjr.triangle(
                    *p,
                    north_square=north_square,
                    south_square=south_square,
                    inverse=True,
                )
                self.assertEqual((None, "equatorial"), cr)
            for i, p in enumerate(np_points_inv[north_square]):
                cr = pjr.triangle(
                    *p,
                    north_square=north_square,
                    south_square=south_square,
                    inverse=True,
                )
                j = (north_square + i) % 4
                self.assertEqual((j, "north_polar"), cr)
            for i, p in enumerate(sp_points_inv[south_square]):
                cr = pjr.triangle(
                    *p,
                    north_square=north_square,
                    south_square=south_square,
                    inverse=True,
                )
                j = (south_square + i) % 4
                self.assertEqual((j, "south_polar"), cr)

    def test_combine_triangles(self):
        # Create test points in the equatorial, north_square polar, and
        # south_square polar regions.
        dx = pi / 4
        dy = pi / 12
        eq_points = [(-pi + i * pi / 2 + dx, 0) for i in range(4)]
        np_points = [(-pi + i * pi / 2 + dx, pi / 4 + dy) for i in range(4)]
        sp_points = [(-pi + i * pi / 2 + dx, -pi / 4 - dy) for i in range(4)]
        error = 1e-15  # Fuzz to handle small rounding errors.

        # Test forward projection.
        # Equatorial points should stay fixed.
        for north_square, south_square in product(list(range(4)), repeat=2):
            for p in eq_points:
                q = pjr.combine_triangles(
                    *p, north_square=north_square, south_square=south_square
                )
                self.assertEqual(p, q)

        # Test forward projection.
        # Polar points should map to the correct points.
        for north_square, south_square in product(list(range(4)), repeat=2):
            # Corners of north square.
            ndl = (-pi + north_square * pi / 2, pi / 4)
            ndr = (-pi / 2 + north_square * pi / 2, pi / 4)
            nur = (-pi / 2 + north_square * pi / 2, 3 * pi / 4)
            nul = (-pi + north_square * pi / 2, 3 * pi / 4)
            # Corners of south square.
            sdl = (-pi + south_square * pi / 2, -3 * pi / 4)
            sdr = (-pi / 2 + south_square * pi / 2, -3 * pi / 4)
            sur = (-pi / 2 + south_square * pi / 2, -pi / 4)
            sul = (-pi + south_square * pi / 2, -pi / 4)
            # np_points should map to the following points,
            # though not necessarily in the order given.
            np_points_transformed = [
                (ndl[0] + dx, ndl[1] + dy),
                (ndr[0] - dy, ndr[1] + dx),
                (nur[0] - dx, nur[1] - dy),
                (nul[0] + dy, nul[1] - dx),
            ]
            # sp_points should map to the following points,
            # though not necessarily in the order given.
            sp_points_transformed = [
                (sul[0] + dx, sul[1] - dy),
                (sur[0] - dy, sur[1] - dx),
                (sdr[0] - dx, sdr[1] + dy),
                (sdl[0] + dy, sdl[1] + dx),
            ]
            for i, p in enumerate(np_points):
                q = pjr.combine_triangles(
                    *p, north_square=north_square, south_square=south_square
                )
                qq = np_points_transformed[(i - north_square) % 4]
                self.assertTrue(euclidean(q, qq) < error)
            for i, p in enumerate(sp_points):
                q = pjr.combine_triangles(
                    *p, north_square=north_square, south_square=south_square
                )
                qq = sp_points_transformed[(i - south_square) % 4]
                self.assertTrue(euclidean(q, qq) < error)

        # Test inverse projection.
        # The inverse of the projection of a point p should yield p.
        for north_square, south_square in product(list(range(4)), repeat=2):
            for p in eq_points + sp_points + np_points:
                q = pjr.combine_triangles(
                    *p, south_square=south_square, north_square=north_square
                )
                pp = pjr.combine_triangles(
                    *q,
                    south_square=south_square,
                    north_square=north_square,
                    inverse=True,
                )
                self.assertTrue(euclidean(p, pp) < error)
        # Points on the boundary of cells S and N should map to points
        # with y-coordniate -pi/4 and pi/4, respectively.
        for north_square, south_square in product(list(range(4)), repeat=2):
            n_shift = array((north_square * pi / 2, 0))
            s_shift = array((south_square * pi / 2, 0))
            np_boundary = [
                (-pi, pi / 4),
                (-pi, pi / 2),
                (-pi, 3 * pi / 4),
                (-3 * pi / 4, 3 * pi / 4),
                (-pi / 2, 3 * pi / 4),
                (-pi / 2, pi / 2),
                (-pi / 2, pi / 4),
                (-3 * pi / 4, pi / 4),
            ]
            sp_boundary = [(q[0], -q[1]) for q in np_boundary]
            np_boundary = [array(q) + n_shift for q in np_boundary]
            sp_boundary = [array(q) + s_shift for q in sp_boundary]
            for q in np_boundary:
                p = pjr.combine_triangles(
                    *q,
                    south_square=south_square,
                    north_square=north_square,
                    inverse=True,
                )
                self.assertTrue(euclidean(atleast_1d(p[1]), atleast_1d(pi / 4)) < error)
            for q in sp_boundary:
                p = pjr.combine_triangles(
                    *q,
                    south_square=south_square,
                    north_square=north_square,
                    inverse=True,
                )
                self.assertTrue(
                    euclidean(atleast_1d(p[1]), atleast_1d(-pi / 4)) < error
                )

    def test_rhealpix_sphere(self):
        # Test forward projection.
        # Should return the same output as healpix_sphere() of the unit
        # sphere, followed by combine_triangles(), followed by a scaling up.
        for ns, ss in product(list(range(4)), repeat=2):
            for p in inputs:
                q = pjr.rhealpix_sphere(*p, north_square=ns, south_square=ss)
                qq = pjr.combine_triangles(
                    *pjh.healpix_sphere(*p), north_square=ns, south_square=ss
                )
                qq = tuple(qq)
                self.assertEqual(q, qq)

        # Test inverse projection.
        # The inverse of the projection of a point p should yield p.
        error = 1e-15  # Fuzz to handle small rounding errors.
        for ns, ss in product(list(range(4)), repeat=2):
            for p in inputs:
                q = pjr.rhealpix_sphere(*p, north_square=ns, south_square=ss)
                pp = pjr.rhealpix_sphere_inverse(*q, north_square=ns, south_square=ss)
                self.assertTrue(euclidean(p, pp) < error)

    def test_rhealpix_ellipsoid(self):
        # Test map projection.
        # Should return the same output as healpix_ellipsoid(), followed
        # by a scaling down, followed by combine_triangles(),
        # followed by a scaling up.
        e = 0.8
        for ns, ss in product(list(range(4)), repeat=2):
            for p in inputs:
                q = pjr.rhealpix_ellipsoid(*p, north_square=ns, south_square=ss, e=e)
                qq = pjh.healpix_ellipsoid(*p, e=e)
                qq = pjr.combine_triangles(*qq, north_square=ns, south_square=ss)
                self.assertEqual(q, qq)

        # Test inverse projection.
        # The inverse of the projection of a point p should yield p.
        # Fuzz for rounding errors based on the error of the approximation to
        # the inverse authalic latitude function:
        alpha = pi / 4
        alpha_ = auth_lat(auth_lat(alpha, e), e, inverse=True)
        error = 10 * rel_err(alpha_, alpha)
        for ns, ss in product(list(range(4)), repeat=2):
            for p in inputs:
                q = pjr.rhealpix_ellipsoid(*p, north_square=ns, south_square=ss, e=e)
                pp = pjr.rhealpix_ellipsoid_inverse(
                    *q, north_square=ns, south_square=ss, e=e
                )
                self.assertTrue(rel_err(p, pp) < error)

    def test_rhealpix(self):
        inputs = [(-pi, pi / 3), (0, pi / 4), (pi / 2, -pi / 6)]

        # Should agree with rhealpix_ellipsoid and rhealpix_ellipsoid_inverse.
        e = 0.5
        a = 7
        R_A = auth_rad(a, e)
        f = pjr.rhealpix(a=a, e=e)
        for p in inputs:
            get = f(*p, radians=True)
            expect = tuple(R_A * array(pjr.rhealpix_ellipsoid(*p, e=e)))
            for i in range(len(expect)):
                self.assertAlmostEqual(get[i], expect[i])
            get = f(*get, radians=True, inverse=True)
            expect = tuple(array(expect) / R_A)
            expect = pjr.rhealpix_ellipsoid_inverse(*expect, e=e)
            for i in range(len(expect)):
                self.assertAlmostEqual(get[i], expect[i])

        # Should work in degrees mode.
        for p in inputs:
            get = f(*rad2deg(p), radians=False)
            expect = f(*p, radians=True)
            for i in range(len(expect)):
                self.assertAlmostEqual(get[i], expect[i])


if __name__ == "__main__":
    unittest.main()
