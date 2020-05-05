r"""
This Python 2.7 code tests ``projections.py``.
Beware, these tests cover only some functions and only some scenarios.
Keep adding tests!

CHANGELOG:

- Alexander Raichev (AR), 2012-03-22: Initial version.
- AR, 2012-07-10: Added a test in test_healpix_sphere().
"""
#*****************************************************************************
#       Copyright (C) 2011 Alexander Raichev <tortoise.said@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from scipy.spatial.distance import euclidean, norm
import unittest
from itertools import product
from projections import *

# Relative error function.
def rel_err(get, expect):
    a = euclidean(get, expect)
    b = norm(expect)
    if b == 0:
        return a
    else:
        return a/b
        
# Define lon-lat input points to test.
phi_0 = arcsin(2.0/3)
a = (0, pi/3)
b = (0, -pi/3)      
inputs = [
    (0, 0), 
    (0, phi_0), 
    (0, -phi_0), 
    (pi/2, 0), 
    (-pi/2, 0), 
    (-pi, 0),
    (-pi, pi/2), 
    (-pi, -pi/2), 
    a, 
    b,
]


class SCENZGridTestCase(unittest.TestCase):
    def test_healpix_sphere(self):
        # Expected outputs of healpix_sphere() applied to inputs.
        sigma_a = sqrt(3 - 3*sin(a[1]))
        ha = (pi/4*(1 - sigma_a), pi/4*(2 - sigma_a))  
        hb = (ha[0], -ha[1])
        healpix_sphere_outputs = [
            (0, 0), 
            (0, pi/4), 
            (0, -pi/4), 
            (pi/2, 0), 
            (-pi/2, 0), 
            (-pi, 0),
            (-3*pi/4, pi/2), 
            (-3*pi/4, -pi/2), 
            ha, 
            hb
        ]        

        # Forward projection should be correct on test points.
        given = inputs
        get = [healpix_sphere(*p) for p in given]
        expect = healpix_sphere_outputs
        # Fuzz to allow for rounding errors:
        error = 1e-12
        for i in range(len(get)):
            self.assertTrue(rel_err(get[i], expect[i]) < error)
         
        # Inverse of projection of a point p should yield p.
        given = get
        get = [healpix_sphere(*q, inverse=True) for q in given]
        expect = inputs
        for i in range(len(get)):
            self.assertTrue(rel_err(get[i], expect[i]) < error)
            
        # Inverse projection of p below should return longitude of -pi.
        # Previously, it was returning a number slightly less than pi
        # because of a rounding error, which get magnified by
        # wrap_lon())
        p = (-7*pi/8, 3*pi/8)
        get = healpix_sphere(*p, inverse=True)
        expect = (-pi, arcsin(1 - 1.0/12))   
        self.assertEqual(get, expect)
        q = (-5*pi/6, 5*pi/12)
        get = healpix_sphere(*q, inverse=True)
        expect = (-pi, arcsin(1 - 1.0/27))
        self.assertEqual(get, expect)
        
    def test_healpix_ellipsoid(self):
        # Expected output of healpix_ellipsoid() applied to inputs.
        e = 0.8
        healpix_ellipsoid_outputs = []
        for p in inputs:
            lam, phi = p
            beta = auth_lat(phi, e=e)
            q = healpix_sphere(lam, beta)
            healpix_ellipsoid_outputs.append(q)

        # Forward projection should be correct on test points.
        given = inputs
        get = [healpix_ellipsoid(*p, e=e) for p in given]
        expect = healpix_ellipsoid_outputs
        # Fuzz to allow for rounding errors:
        error = 1e-12
        for i in range(len(get)):
            self.assertTrue(rel_err(get[i], expect[i]) < error)
         
        # Inverse of projection of a point p should yield p.
        given = get
        get = [healpix_ellipsoid(*q, e=e, inverse=True) for q in given]
        expect = inputs
        # Fuzz for rounding errors based on the error of the approximation to
        # the inverse authalic latitude function:
        alpha = pi/4
        alpha_ = auth_lat(auth_lat(alpha, e), e, inverse=True)
        error = 10*rel_err(alpha_, alpha)        
        for i in range(len(get)):
            self.assertTrue(rel_err(get[i], expect[i]) < error)
        
    def test_cap(self):
        # Create test points in the equatorial, north polar, and
        # south polar regions.
        dx = pi/12
        dy = pi/12
        eq_points = [(-pi + i*pi/2 + dx, 0 ) for i in range(4)] 
        np_points = [(-pi + i*pi/2 + dx, pi/4 + dy ) for i in range(4)]
        sp_points = [(-pi + i*pi/2 + dx, -pi/4 - dy ) for i in range(4)] 
        np_points_inv = {}
        for north in range(4):
            u = array((-3*pi/4 + north*pi/2, pi/2))
            np_points_inv[north] = [
                u + array((0, -dy)),
                u + array((dx, 0)),
                u + array((0, dy)),
                u + array((-dx, 0))
            ]   
        sp_points_inv = {}
        for south in range(4):
            u = array((-3*pi/4 + south*pi/2, -pi/2))        
            sp_points_inv[south] = [
                u + array((0, dy)),
                u + array((dx, 0)),
                u + array((0, -dy)),
                u + array((-dx, 0))
            ]
        # Test forward.        
        for (north, south) in product(range(4), repeat=2): 
            for p in eq_points:
                cr = cap(*p, north=north, 
                         south=south)
                self.assertEqual((None, 'equatorial'), cr)
            for (i, p) in enumerate(np_points):
                cr = cap(*p, north=north, 
                         south=south)
                self.assertEqual((i, 'north'), cr)
            for (i, p) in enumerate(sp_points):
                cr = cap(*p, north=north, 
                         south=south)
                self.assertEqual((i, 'south'), cr)
    
        # Test inverse.        
        for (north, south) in product(range(4), repeat=2): 
            for p in eq_points:
                cr = cap(*p, north=north, 
                         south=south, inverse=True)
                self.assertEqual((None, 'equatorial'), cr)
            for (i, p) in enumerate(np_points_inv[north]):
                cr = cap(*p, north=north, 
                         south=south, inverse=True)
                j = (north + i) % 4
                self.assertEqual((j, 'north'), cr)
            for (i, p) in enumerate(sp_points_inv[south]):
                cr = cap(*p, north=north, 
                         south=south, inverse=True)
                j = (south + i) % 4
                self.assertEqual((j, 'south'), cr)

    def test_combine_caps(self):
        # Create test points in the equatorial, north polar, and
        # south polar regions.
        dx = pi/4
        dy = pi/12
        eq_points = [(-pi + i*pi/2 + dx, 0 ) for i in range(4)] 
        np_points = [(-pi + i*pi/2 + dx, pi/4 + dy ) for i in range(4)]
        sp_points = [(-pi + i*pi/2 + dx, -pi/4 - dy ) for i in range(4)] 
        error = 1e-15  # Fuzz to handle small rounding errors.
                        
        # Test forward projection.        
        # Equatorial points should stay fixed.
        for (north, south) in product(range(4), repeat=2): 
            for p in eq_points:
                q = combine_caps(*p, north=north, 
                                 south=south)
                self.assertEqual(p, q)
        
        # Test forward projection.
        # Polar points should map to the correct points.
        for (north, south) in product(range(4), repeat=2): 
            # Corners of north square.
            ndl = (-pi + north*pi/2, pi/4)
            ndr = (-pi/2 + north*pi/2, pi/4)
            nur = (-pi/2 + north*pi/2, 3*pi/4)
            nul = (-pi + north*pi/2, 3*pi/4)
            # Corners of south square.
            sdl = (-pi + south*pi/2, -3*pi/4)
            sdr = (-pi/2 + south*pi/2, -3*pi/4)
            sur = (-pi/2 + south*pi/2,  -pi/4)
            sul = (-pi + south*pi/2, -pi/4)
            # np_points should map to the following points,
            # though not necessarily in the order given.
            np_points_transformed = [
                (ndl[0] + dx, ndl[1] + dy), 
                (ndr[0] - dy, ndr[1] + dx),
                (nur[0] - dx, nur[1] - dy),
                (nul[0] + dy, nul[1] - dx)
            ]
            # sp_points should map to the following points,
            # though not necessarily in the order given.
            sp_points_transformed = [
                (sul[0] + dx, sul[1] - dy), 
                (sur[0] - dy, sur[1] - dx),
                (sdr[0] - dx, sdr[1] + dy),
                (sdl[0] + dy, sdl[1] + dx)
            ]
            for (i, p) in enumerate(np_points):
                q = combine_caps(*p, north=north, 
                                 south=south)
                qq = np_points_transformed[(i - north) % 4]
                self.assertTrue(euclidean(q, qq) < error)
            for (i, p) in enumerate(sp_points):
                q = combine_caps(*p, north=north, 
                                 south=south)
                qq = sp_points_transformed[(i - south) % 4]
                self.assertTrue(euclidean(q, qq) < error)
                      
        # Test inverse projection.
        # The inverse of the projection of a point p should yield p.
        for (north, south) in product(range(4), repeat=2): 
            for p in eq_points + sp_points + np_points:
                q = combine_caps(*p, south=south, 
                                 north=north)
                pp = combine_caps(*q, south=south, 
                                  north=north, inverse=True)
                self.assertTrue(euclidean(p, pp) < error)
        # Points on the boundary of cells S and N should map to points
        # with y-coordniate -pi/4 and pi/4, respectively.
        for (north, south) in product(range(4), repeat=2): 
            n_shift = array((north*pi/2, 0))
            s_shift = array((south*pi/2, 0))
            np_boundary = [
                (-pi, pi/4), 
                (-pi, pi/2), 
                (-pi, 3*pi/4), 
                (-3*pi/4, 3*pi/4), 
                (-pi/2, 3*pi/4), 
                (-pi/2, pi/2), 
                (-pi/2, pi/4), 
                (-3*pi/4, pi/4)
            ]
            sp_boundary = [(q[0], -q[1]) for q in np_boundary]
            np_boundary = [array(q) + n_shift for q in np_boundary]
            sp_boundary = [array(q) + s_shift for q in sp_boundary]
            for q in np_boundary:
                p = combine_caps(*q, south=south, 
                                 north=north, inverse=True)
                self.assertTrue(euclidean(p[1], pi/4) < error)
            for q in sp_boundary:
                p = combine_caps(*q, south=south, 
                                 north=north, inverse=True)
                self.assertTrue(euclidean(p[1], -pi/4) < error)
                
    def test_rhealpix_sphere(self):
        # Test forward projection.
        # Should return the same output as healpix_sphere() of the unit
        # sphere, followed by combine_caps(), followed by a scaling up.
        for (ns, ss) in product(range(4), repeat=2): 
            for p in inputs:
                q = rhealpix_sphere(*p, north=ns, south=ss)
                qq = combine_caps(*healpix_sphere(*p), north=ns, 
                                  south=ss)
                qq = tuple(qq)
                self.assertEqual(q, qq)
    
        # Test inverse projection.
        # The inverse of the projection of a point p should yield p.
        error = 1e-15   # Fuzz to handle small rounding errors.
        for (ns, ss) in product(range(4), repeat=2): 
            for p in inputs:
                q = rhealpix_sphere(*p, north=ns, south=ss)
                pp = rhealpix_sphere(*q, north=ns, south=ss, 
                                     inverse=True)                
                self.assertTrue(euclidean(p, pp) < error)
                
    def test_rhealpix_ellipsoid(self):
        # Test map projection.
        # Should return the same output as healpix_ellipsoid(), followed
        # by a scaling down, followed by combine_caps(), 
        # followed by a scaling up.
        e = 0.8
        for (ns, ss) in product(range(4), repeat=2): 
            for p in inputs:
                q = rhealpix_ellipsoid(*p, north=ns, south=ss, 
                                       e=e)
                qq = healpix_ellipsoid(*p, e=e)
                qq = combine_caps(*qq, north=ns, south=ss)
                self.assertEqual(q, qq)
                
        # Test inverse projection.
        # The inverse of the projection of a point p should yield p.
        # Fuzz for rounding errors based on the error of the approximation to
        # the inverse authalic latitude function:
        alpha = pi/4
        alpha_ = auth_lat(auth_lat(alpha, e), e, inverse=True)
        error = 10*rel_err(alpha_, alpha)        
        for (ns, ss) in product(range(4), repeat=2): 
            for p in inputs:
                q = rhealpix_ellipsoid(*p, north=ns, south=ss, 
                                       e=e)
                pp = rhealpix_ellipsoid(*q, north=ns, south=ss, 
                                        e=e, inverse=True)                
                self.assertTrue(rel_err(p, pp) < error)
#------------------------------------------------------------------------------                           
if __name__ == '__main__':
    unittest.main()