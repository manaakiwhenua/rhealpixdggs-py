r"""
This Python 2.7 code tests ``grids.py`.
Beware, while these tests cover most functions and scenarios,
they don't cover them all.
Keep adding tests!

CHANGELOG:

- Alexander Raichev (AR), 2012-03-01: Initial version
- AR, 2012-03-30: Added more tests in the interim.
- AR, 2012-07-04: Refactored cell family tests.
- AR, 2012-07-06: Added another test for cell nuclei.

"""
#*****************************************************************************
#       Copyright (C) 2011 Alexander Raichev <tortoise.said@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from scipy import array, sin, cos, tan, arctan, arcsin
import unittest
from random import randint, uniform 
import projections
from grids import *

class SCENZGridTestCase(unittest.TestCase):

    def setUp(self):
        self.E = Earth('WGS84', degrees=False)
        
        # Euclidean distance between two points.
        from scipy.spatial.distance import euclidean
        self.distance = euclidean 
        #------------------------------------------------------------------------------            
    # Test that the projections scale properly.
    def test_healpix_and_rhealpix(self):
        unit_sphere = Earth('unit_sphere', degrees=False)
        WGS84_asphere = Earth('WGS84_asphere', degrees=False)
        WGS84 = Earth('WGS84', degrees=False)
        R = projections.WGS84_R
        e = projections.WGS84_E
        
        # Test forward projections.
        # Arbritrary (lon, lat) point and its authalic latitude version.
        p = (pi/5, pi/3)
        pa = (p[0], projections.auth_lat(p[1], e))
        # Projection of p from WGS84_asphere should be R times bigger
        # than the projection of p from unit_sphere.
        self.assertEqual(WGS84_asphere.healpix(*p), 
                         scale(unit_sphere.healpix(*p), R))
        self.assertEqual(WGS84_asphere.rhealpix(*p), 
                         scale(unit_sphere.rhealpix(*p), R))
        # Projection of p from WGS84 should be R times bigger
        # than the projection of pa from unit_sphere.
        self.assertEqual(WGS84.healpix(*p), 
                         scale(unit_sphere.healpix(*pa), R))
        self.assertEqual(WGS84.rhealpix(*p), 
                         scale(unit_sphere.rhealpix(*pa), R))
        
        # Test inverse projections.
        # Arbitrary (x, y) point and the authalic latitude version of 
        # its inverse projection onto the unit sphere.
        q = (-pi/7, -pi/5)
        pp = unit_sphere.healpix(*q, inverse=True)
        pp = (pp[0], projections.auth_lat(pp[1], e, inverse=True))                  
        # The inverse projection of q onto unit_sphere should equal the 
        # inverse projection of R*q onto WGS84_asphere.
        self.assertEqual(WGS84_asphere.healpix(*scale(q, R), inverse=True), 
                         unit_sphere.healpix(*q, inverse=True))
        self.assertEqual(WGS84_asphere.rhealpix(*scale(q, R), inverse=True), 
                         unit_sphere.rhealpix(*q, inverse=True))
        # The inverse projection of q*R onto WGS84 should equal ppp.
        self.assertEqual(WGS84.healpix(*scale(q, R), inverse=True), pp)
        self.assertEqual(WGS84.rhealpix(*scale(q, R), inverse=True), pp)
                         #------------------------------------------------------------------------------            
    # Test Cell methods
    def test_suid_from_index(self):
        # Test order='level'.
        # The suid with index 0 should be N.
        suid = ('N',)
        got = Cell.suid_from_index(self.E, 0, order='level')
        self.assertEqual(got, suid)
        
        # Should be correct on cell P1.
        suid = ('P', 1)
        i = 25
        got = Cell.suid_from_index(self.E, i, order='level')
        self.assertEqual(got, suid)
        
        # Test order='post'.
        # The suid with index 0 should be N00...0.
        suid = ['N'] + [0 for i in range(self.E.max_level)]
        suid = tuple(suid)
        got = Cell.suid_from_index(self.E, 0, order='post')
        self.assertEqual(got, suid)
        
        # Should be correct on cell P1.
        suid = ('P', 1)
        def num(k):
            return self.E.num_cells(level_a=k, subcells=True)
        i = 2 * num(0) + 1 * num(1) + num(1) - 1
        got = Cell.suid_from_index(self.E, i, order='post')
        self.assertEqual(got, suid)
                
    def test_Cell_init(self):
        # Should set the earth, suid, and level attributes correctly.
        suid = ('S', 0, 1, 2, 6, 3, 1)
        level = 6
        C = Cell(self.E, suid)
        self.assertEqual(C.earth, self.E)
        self.assertEqual(C.suid, suid)
        self.assertEqual(C.level, level)
        
        # Should not create invalid cells.
        suid = 'G0'
        self.assertRaises(AssertionError, Cell, self.E, suid)
                    
        # Should be correct on cell P1.
        suid = 'P1'
        i = 25
        got = Cell(self.E, level_order_index=i)
        self.assertEqual(str(got), suid)
        
        suid = 'P1'
        def num(k):
            return self.E.num_cells(level_a=k, subcells=True)
        i = 2 * num(0) + 1 * num(1) + num(1) - 1
        got = Cell(self.E, post_order_index=i)
        self.assertEqual(str(got), suid)
        
    def test_Cell_str9(self):
        C = Cell(self.E, 'O01452')
        self.assertEqual(C.str9(), '101452') 
            
    def test_Cell_suid_colrow(self):
        # Should work for level 0 cells.
        for suid in Earth.CELLS0:
            self.assertEqual(suid, Cell(self.E, suid).suid_colrow()[0][0])
            self.assertEqual(suid, Cell(self.E, suid).suid_colrow()[1][0])

        # Should work on an arbitrary cell.
        c = Cell(self.E, 'N138686')
        col, row = c.suid_colrow()
        for i in range(1, c.level + 1):
            self.assertEqual(c.suid[i], col[i] + 3*row[i])    
    
    def test_Cell_eq(self):
        A = Cell(self.E, 'O123')
        AA = Cell(self.E, 'O123')
        B = Cell(self.E, 'N2875')
        self.assertTrue(A == AA)
        self.assertFalse(A == B)
    
    def test_Cell_ne(self):
        A = Cell(self.E, 'O123')
        AA = Cell(self.E, 'O123')
        B = Cell(self.E, 'N2875')
        self.assertFalse(A != AA)
        self.assertTrue(A != B)
        
    def test_Cell_lt(self):
        a = Cell(self.E, 'N7681')
        b = Cell(self.E, 'O123')
        c = Cell(self.E, 'O12')
        self.assertTrue(a < b)
        self.assertTrue(a < c)
        self.assertFalse(a < a)
        self.assertFalse(b < a)
        self.assertFalse(c < a)
    
    def test_Cell_le(self):
        a = Cell(self.E, 'N7681')
        b = Cell(self.E, 'O123')
        bb = Cell(self.E, 'O123')
        c = Cell(self.E, 'O12')
        self.assertTrue(b <= bb)
        self.assertTrue(bb <= b)
        self.assertTrue(a <= b)
        self.assertTrue(a <= c)
        self.assertFalse(b <= a)
        self.assertFalse(c <= a)
        
    def test_Cell_gt(self):
        a = Cell(self.E, 'N7681')
        b = Cell(self.E, 'O123')
        c = Cell(self.E, 'O12')
        self.assertFalse(a > b)
        self.assertFalse(a > c)
        self.assertFalse(a > a)
        self.assertTrue(b > a)
        self.assertTrue(c > a)
        
    def test_Cell_ge(self):
        a = Cell(self.E, 'N7681')
        b = Cell(self.E, 'O123')
        bb = Cell(self.E, 'O123')
        c = Cell(self.E, 'O12')
        self.assertTrue(b >= bb)
        self.assertTrue(bb >= a)
        self.assertFalse(a >= b)
        self.assertFalse(a >= c)
        self.assertTrue(b >= a)
        self.assertTrue(c >= a)
        
    def test_Cell_index(self):
        # Empty cell should have index None
        empty = self.E.cell()
        self.assertEqual(empty.index(order='level'), None)
        self.assertEqual(empty.index(order='post'), None)

        # Test order='level'.
        # A cell should have index one greater than its predecessor
        # at the same level.
        c = self.E.cell('N2748')
        b = c.predecessor()
        self.assertEqual(c.index(order='level'), b.index(order='level') + 1)
        
        # It should invert suid_from_index().
        N = self.E.num_cells(0, self.E.max_level)
        k = 3616048 #randint(0, N - 1)
        got = Cell(self.E, Cell.suid_from_index(self.E, k, order='level')).\
              index(order='level')
        self.assertEqual(got, k)

        # Test order='post'.
        # The cell N00...0 should have index 0.
        suid = 'N' + '0' * self.E.max_level
        a = Cell(self.E, suid)
        self.assertEqual(a.index(order='post'), 0)
        
        # A cell should have index one greater than its last child.
        cc = list(c.subcells())[8]
        self.assertEqual(c.index(order='post'),       
                         cc.index(order='post') + 1)
        
        # It should invert suid_from_index().
        k = 3616048 #randint(0, N - 1)
        got = Cell(self.E, Cell.suid_from_index(self.E, k, order='post')).\
              index(order='post')
        self.assertEqual(got, k)

    def test_Cell_successor(self):
        A = Cell(self.E, 'N81')
        As0 = Cell(self.E, 'O')
        As1 = Cell(self.E, 'O0')
        As2 = Cell(self.E, 'N82')
        As3 = Cell(self.E, 'N820')
        self.assertTrue(A.successor(0) == As0)
        self.assertTrue(A.successor(1) == As1)
        self.assertTrue(A.successor(2) == As2)
        self.assertTrue(A.successor(3) == As3)
        A = Cell(self.E, 'S88')
        self.assertTrue(A.successor(1) == None)
    
    def test_Cell_predecessor(self):
        A = Cell(self.E, 'O01')
        Ap0 = Cell(self.E, 'N')
        Ap1 = Cell(self.E, 'N8')
        Ap2 = Cell(self.E, 'O00')
        Ap3 = Cell(self.E, 'O018')
        self.assertTrue(A.predecessor(0) == Ap0)
        self.assertTrue(A.predecessor(1) == Ap1)
        self.assertTrue(A.predecessor(2) == Ap2)
        self.assertTrue(A.predecessor(3) == Ap3)
        A = Cell(self.E, 'N00')
        self.assertTrue(A.predecessor(1) == None)
    
    def test_Cell_subcell(self):
        C = Cell(self.E, 'S120')
        D = Cell(self.E, 'S1205')
        self.assertTrue(D.subcell(C))
        self.assertFalse(C.subcell(D))
    
    def test_Cell_subcells(self):
        l = 6        
        C = Cell(self.E, 'S105773')
        s = list(C.subcells(l + 1))
        t = [Cell(self.E, str(C) + str(i)) for i in range(9)]
        for i in range(9):
            self.assertTrue(s[i] == t[i])    
            
    def test_Cell_ul_vertex(self):
        # Should work for level 0 cells.
        for X in self.E.grid(0):
            got = X.ul_vertex()
            should_get = self.E.ul_vertex[str(X)]
            self.assertEqual(got, should_get)
        
        for X in self.E.grid(0):
            got = X.ul_vertex(surface='ellipsoid')
            p = self.E.ul_vertex[str(X)]
            should_get = self.E.project(p, source='rhealpix', target='lonlat')
            self.assertEqual(got, should_get)
            
        # Should work on children cells.
        a = Cell(self.E, 'S23')
        l = a.level
        x, y = a.ul_vertex()        
        w = self.E.cell_width(l + 1)
        error = 1e-15 # Error tolerance.
        for (i,j) in product(range(3), repeat=2):
            s = str(i + j*3)
            b = Cell(self.E, str(a) + s)
            xx, yy = b.ul_vertex()
            xp, yp = (x + i*w, y - j*w) 
            self.assertTrue(self.distance([xx, yy], [xp, yp]) < error)
        
        a = Cell(self.E, 'S23')
        l = a.level
        x, y = a.ul_vertex()        
        w = self.E.cell_width(l + 1)
        error = 1e-15 # Error tolerance.
        for (i,j) in product(range(3), repeat=2):
            s = str(i + j*3)
            b = Cell(self.E, str(a) + s)
            xx, yy = b.ul_vertex(surface='ellipsoid')
            xp, yp = self.E.project((x + i*w, y - j*w), \
                                    source='rhealpix', target='lonlat') 
            self.assertTrue(self.distance([xx, yy], [xp, yp]) < error)
          
    def test_Cell_nucleus_and_vertices(self):
        # The nucleus and vertices of a cell should equal the nucleus 
        # and vertices of its appropriate children.
        a = self.E.cell('P367')
        a_nv = a.nucleus_and_vertices()
        appropriate_child = {0: 4, 1: 0, 2: 6, 3: 8, 4: 2}
        error = 1e-5  # Error tolerance in meters.
        for (k, v) in appropriate_child.items():
            aa = list(a.subcells())[v]
            aa_nv = aa.nucleus_and_vertices()
            self.assertTrue(self.distance(a_nv[k], aa_nv[k]) < error)

        # Nuclei of children should be in correct position 
        # relative to parent cell in rHEALPix projection.
        a = Cell(self.E, 'S74121')   # Arbitrary cell.
        w = a.width()
        (x, y) = a.ul_vertex()        
        error = 1e-5  # Error tolerance in meters.
        for (row, col) in product(range(3), repeat=2):
            s = str(row*3 + col)
            # Child cell in (row, column) position relative to a:
            b = Cell(self.E, str(a) + s)   
            nucleus = b.nucleus_and_vertices()[0]
            correct_nucleus = (x + w/6 + (w/3)*col, y - w/6 - (w/3)*row)
            self.assertTrue(self.distance(nucleus, correct_nucleus) < error)
        
    def test_Cell_atomic_rotate(self):
        # Should return correct values.
        s = []  # Function values.
        t = []  # Correct values.
        for q in range(4):
            s.append ([Cell.atomic_rotate(q, x) for x in range(9)])
        t.append([0, 1, 2, 3, 4, 5, 6, 7, 8])
        t.append([2, 5, 8, 1, 4, 7, 0, 3, 6])
        t.append([8, 7, 6, 5, 4, 3, 2, 1, 0])
        t.append([6, 3, 0, 7, 4, 1, 8, 5, 2])
        for q in range(4):
            self.assertEqual(s[q], t[q])
            
    def test_Cell_neighbor(self):
        X = Cell(self.E, 'S0')
        Xl = X.neighbor('left')
        Y = Cell(self.E, 'R8')
        self.assertTrue(Xl == Y)
        
    def test_Cell_neighbors(self):
        # Should work on level 1 cells with suids of the form s0
        neighbors = []
        for s in Earth.CELLS0:
           nb = Cell(self.E, s + '0').neighbors()
           for d in ['left', 'right', 'up', 'down']:
               neighbors.append(nb[d])
        correct_neighbors = [
        Cell(self.E, 'R0'), Cell(self.E, 'N1'), \
        Cell(self.E, 'Q2'), Cell(self.E, 'N3'), \
        Cell(self.E, 'R2'), Cell(self.E, 'O1'), \
        Cell(self.E, 'N6'), Cell(self.E, 'O3'), \
        Cell(self.E, 'O2'), Cell(self.E, 'P1'), \
        Cell(self.E, 'N8'), Cell(self.E, 'P3'), \
        Cell(self.E, 'P2'), Cell(self.E, 'Q1'), \
        Cell(self.E, 'N2'), Cell(self.E, 'Q3'), \
        Cell(self.E, 'Q2'), Cell(self.E, 'R1'), \
        Cell(self.E, 'N0'), Cell(self.E, 'R3'), \
        Cell(self.E, 'R8'), Cell(self.E, 'S1'), \
        Cell(self.E, 'O6'), Cell(self.E, 'S3')]
        for i in range(len(neighbors)):
            self.assertTrue(neighbors[i] == correct_neighbors[i])
            
    def test_region(self):
        X = self.E.cell('P0')
        self.assertEqual(X.region(), 'equatorial')
        Y = self.E.cell('S78')
        self.assertEqual(Y.region(), 'polar')
        
    def test_ellipsoidal_shape(self):
        for s in self.E.CELLS0[1:5]:
            X = self.E.cell(s)
            self.assertEqual(X.ellipsoidal_shape(), 'rectangle')
        for s in [self.E.CELLS0[0], self.E.CELLS0[5]]:
            for t in ['', '4', '44']:
                X = self.E.cell(s + t)
                self.assertEqual(X.ellipsoidal_shape(), 'cap')
            for t in ['408', '0', '44446', '24']:
                X = self.E.cell(s + t)
                self.assertEqual(X.ellipsoidal_shape(), 'dart')
            for t in ['418', '1', '465', '23']:
                X = self.E.cell(s + t)
                self.assertEqual(X.ellipsoidal_shape(), 'trapezoid')
                
    def test_centroid(self):
        print 'Testing centroid method, which will take a minute...'
        # The centroid of a planar cell is its nucleus. 
        for suid in ['Q7', 'S22']:
            X = self.E.cell(suid)
            centroid = X.centroid()
            nucleus = X.nucleus_and_vertices()[0]
            self.assertEqual(centroid, nucleus)
            
        # The centroid of a ellipsoidal cap cell is also its nucleus.
        for suid in ['N444', 'S']:
            X = self.E.cell(suid)
            centroid = X.centroid(surface='ellipsoid')
            nucleus = X.nucleus_and_vertices(surface='ellipsoid')[0]
            self.assertEqual(centroid, nucleus)
        
        # For other ellipsoidal cells, test centroid() against a Monte Carlo 
        # approximation of the centroid.
        def monte_carlo_centroid(cell):
            nv = cell.nucleus_and_vertices()
            lam_nucleus = cell.earth.project(nv[0], 'rhealpix', 'lonlat')[0]
            vertices = nv[1:]
            x1, x2 = vertices[0][0], vertices[3][0]
            y1, y2 = vertices[1][1], vertices[0][1]
            N = 10000
            sample_points = []
            for i in range(N):
                x, y = uniform(x1, x2), uniform(y1, y2)
                lam, phi = self.E.project((x, y), source='rhealpix', 
                                          target='lonlat')
                sample_points.append(array((lam, phi)))
            lam_bar, phi_bar = sum(sample_points)/N
            sample_var = sum([array((
                                     (p[0] - lam_bar)**2,                                   
                                     (p[1] - phi_bar)**2))
                              for p in sample_points])/(N - 1)
            lam_bar_err = sqrt(sample_var[0]/N) # Approximately 
            phi_bar_err = sqrt(sample_var[1]/N) # Approximately
            if lam_nucleus in [-pi, -pi/2, 0, pi/2]:
                lam_bar = lam_nucleus
            return lam_bar, phi_bar, abs(lam_bar_err), abs(phi_bar_err)
        
        print
        for suid in ['Q7', 'N7', 'N77', 'N6', 'N66', 'S1',
                     'S11', 'S111', 'S2', 'S24']:
            X = self.E.cell(suid)
            print "Testing centroid(surface='ellipsoid') for %s cell %s..."\
              % (X.ellipsoidal_shape(), X)
            lam_bar, phi_bar = X.centroid(surface='ellipsoid')
            lam_bar_approx, phi_bar_approx, lam_bar_err, phi_bar_err = \
            monte_carlo_centroid(X)
            # print 'lam_bar and approx', lam_bar, lam_bar_approx, lam_bar_err
            # print 'phi_bar and approx', phi_bar, phi_bar_approx, phi_bar_err
            self.assertTrue(self.distance(lam_bar, lam_bar_approx) <\
                                                      10*lam_bar_err)
            self.assertTrue(self.distance(phi_bar, phi_bar_approx) <\
                                                      10*phi_bar_err)
                                    
#------------------------------------------------------------------------------        
    # Test CellFamily methods    
    def test_CellFamily_init(self):
        X = self.E.cell('S0')
        Y = self.E.cell('P33')
        Z = self.E.cell('S')
        F = CellFamily([X, Y, X, Z, Y])
        cells = [Y, Z]
        min_level = 0
        max_level = 2 
        self.assertEqual(len(F.cells), len(cells))
        for i in range(len(F.cells)):
            self.assertEqual(F.cells[i].suid, cells[i].suid)
        self.assertEqual(F.min_level, min_level)
        self.assertEqual(F.max_level, max_level)
     
    def test_CellFamily_contains(self):
        f = CellFamily(list(self.E.grid(1)))
        c = self.E.cell('O5')
        d = self.E.cell('O55')
        self.assertTrue(c in f)
        self.assertTrue(d not in f)
        
    def test_CellFamily_delitem(self):
        f = CellFamily(list(self.E.grid(1)))
        c = self.E.cell('N0')
        del(f[0])
        self.assertTrue(c not in f)
        
    def test_CellFamily_add(self):
        f = CellFamily(list(self.E.grid(0)))
        del f[2]    # Remove cell C.
        
        # Should add cell C7 to f in correct spot.
        c = self.E.cell('P7')
        f.add(c)
        self.assertTrue(str(f[2]), str(c))

        # Should not add cell B7 to f since it's a subcell of a cell in f.
        d = self.E.cell('O7')
        f.add(d)
        self.assertTrue(d not in f)
   
    def test_CellFamily_minimize(self):
        # The minimizeion of the empty cell family should be itself.
        f = CellFamily()
        fc = f.minimize()
        self.assertTrue(f is fc)
        self.assertEqual(fc.min_level, None)
        self.assertEqual(fc.max_level, None)
        
        # Minimize.
        A00 = self.E.cell('S00')
        B = self.E.cell('O')
        C00 = self.E.cell('P00')
        C02 = self.E.cell('P02')
        f = CellFamily([A00, C00, C02] + list(B.subcells(3)))
        fc = f.minimize()
        g = CellFamily([A00, B, C00, C02])
        self.assertTrue(fc == g)
        self.assertEqual(fc.min_level, 0)
        self.assertEqual(fc.max_level, 2)
            
        # This family can't be minimizeed.
        cells = list(B.subcells(1))
        cells.pop(7)    
        f = CellFamily(cells)
        fc = f.minimize()
        self.assertTrue(f == fc)
        self.assertEqual(fc.min_level, 1)
        self.assertEqual(fc.max_level, 1)
        
    def test_CellFamily_union(self):
        # Union with the empty cell family.
        f = CellFamily()
        g = CellFamily(list(self.E.grid(0)))
        self.assertTrue(f.union(g) == g)
        
        # Union.
        p17 = self.E.cell('P17')
        p23 = self.E.cell('P23')
        p345 = self.E.cell('P345')
        p346 = self.E.cell('P346')
        p381 = self.E.cell('P381')
        p2 = self.E.cell('P2')
        p3 = self.E.cell('P3')
        f1 = CellFamily([p17, p23, p345, p346, p381])
        f2 = CellFamily([p2, p3])
        u = f1.union(f2)
        uu = CellFamily([p17, p2, p3])
        self.assertTrue(u == uu)
        
        # Union with filtering.
        u = f1.union(f2, filter_level=1)
        uu = CellFamily([p2, p3])
        self.assertTrue(u == uu)
        
        u = f1.union(f2, filter_level=0)
        uu = CellFamily()
        self.assertTrue(u == uu)

    def test_CellFamily_union_all(self):
        # Union.
        p17 = self.E.cell('P17')
        p23 = self.E.cell('P23')
        p345 = self.E.cell('P345')
        p346 = self.E.cell('P346')
        p381 = self.E.cell('P381')
        p2 = self.E.cell('P2')
        p3 = self.E.cell('P3')
        p34 = self.E.cell('P34')
        p38 = self.E.cell('P38')
        f1 = CellFamily([p17, p23, p345, p346, p381])
        f2 = CellFamily([p2, p3])
        f3 = CellFamily([p2, p34, p38])
        u = f1.union_all([f2, f3])
        uu = CellFamily([p17, p2, p3])
        self.assertTrue(u == uu)
        
        # Union with filtering.
        u = f1.union_all([f2, f3], filter_level=1)
        uu = CellFamily([p2, p3])
        self.assertTrue(u == uu)
        
        u = f1.union_all([f2, f3], filter_level=0)
        uu = CellFamily()
        self.assertTrue(u == uu)
        
    def test_CellFamily_intersect(self):
        # Intersect with an empty cell family.
        f = CellFamily()
        g = CellFamily(list(self.E.grid(0)))
        self.assertTrue(f.intersect(g) == f)

        # Intersect.
        p17 = self.E.cell('P17')
        p23 = self.E.cell('P23')
        p345 = self.E.cell('P345')
        p346 = self.E.cell('P346')
        p381 = self.E.cell('P381')
        p2 = self.E.cell('P2')
        p3 = self.E.cell('P3')
        f1 = CellFamily([p17, p23, p345, p346, p381])
        f2 = CellFamily([p2, p3])
        v = f1.intersect(f2)
        vv = CellFamily([p23, p345, p346, p381])
        self.assertTrue(v == vv)

        # Intersect with filtering.
        v = f1.intersect(f2, filter_level=2)
        vv = CellFamily([p23])
        self.assertTrue(v == vv)

        v = f1.intersect(f2, filter_level=1)
        vv = CellFamily()
        self.assertTrue(v == vv)

    def test_CellFamily_intersect_all(self):
        # Intersect.
         p17 = self.E.cell('P17')
         p23 = self.E.cell('P23')
         p345 = self.E.cell('P345')
         p346 = self.E.cell('P346')
         p381 = self.E.cell('P381')
         p2 = self.E.cell('P2')
         p3 = self.E.cell('P3')
         p34 = self.E.cell('P34')
         p38 = self.E.cell('P38')
         f1 = CellFamily([p17, p23, p345, p346, p381])
         f2 = CellFamily([p2, p3])
         f3 = CellFamily([p2, p34, p38])
         v = f1.intersect_all([f2, f3])
         vv = CellFamily([p23, p345, p346, p381])
         self.assertTrue(v == vv)
         
         # Intersect with filtering.
         v = f1.intersect_all([f2, f3], filter_level=2)
         vv = CellFamily([p23])
         self.assertTrue(v == vv)

         v = f1.intersect_all([f2, f3], filter_level=1)
         vv = CellFamily()
         self.assertTrue(v == vv)

#------------------------------------------------------------------------------        
    # Test Earth methods.
    def test_Earth_init(self):
        # Atomic neighbors of integers 0,...,8 should be consistent. 
        # For example, the down neighbor of the up neighbor of k should be k.        
        an = self.E.atomic_neighbors
        for k in range(9):
            self.assertEqual(an[an[k]['up']]['down'], k)
            self.assertEqual(an[an[k]['right']]['left'], k)
        
    def test_Earth_grid(self):
        # Should output the correct number of cells.
        l = 3
        cells = list(self.E.grid(l))
        self.assertEqual(len(cells), 6*9**l)
        
        # Pick a random cell of index k in cells and check that cells[k+1]
        # is the successor cell.
        k = randint(0, 6*9**l - 2)
        a = cells[k]
        b = cells[k + 1]
        c = a.successor()
        self.assertEqual(str(b), str(c))
            
    def test_Earth_interval(self):
        # Should produce the correct number of cells
        # and each cell should have level order index one greater 
        # than the previous.
        count = 0
        A = self.E.cell('N1')
        B = self.E.cell('S28')
        start_index = self.E.cell('N20').index(order='level')
        end_index = self.E.cell('S28').index(order='level')
        correct_count = end_index - start_index + 1
        count = 0
        old_index = start_index - 1
        for X in self.E.interval(A, B):
            new_index = X.index(order='level')
            self.assertEqual(new_index, old_index + 1)
            old_index = new_index
            count += 1
        self.assertEqual(count, correct_count)
        
        # Should produce the correct number of cells
        # and each cell should have level order index one greater 
        # than the previous.
        count = 0
        A = self.E.cell('N20')
        B = self.E.cell('S2')
        start_index = self.E.cell('N20').index(order='level')
        end_index = self.E.cell('S28').index(order='level')
        correct_count = end_index - start_index + 1
        count = 0
        old_index = start_index - 1
        for X in self.E.interval(A, B):
            new_index = X.index(order='level')
            self.assertEqual(new_index, old_index + 1)
            old_index = new_index
            count += 1
        self.assertEqual(count, correct_count)

    def test_Earth_cell_from_point(self):
         # The nucleus of a cell should yield the cell.
         for surface in ['plane', 'ellipsoid']:
             for x in self.E.grid(1):
                 c = x.nucleus_and_vertices(surface=surface)[0]
                 y = self.E.cell_from_point(x.level, c, surface=surface)
                 self.assertEqual(x.suid, y.suid)
                
    def test_Earth_cell_from_region(self):
        # For any cell X in the image of the rhealpix projecton
        # with nucleus c and width w, a square centered at c 
        # with width w/2 should yield X.
        X = self.E.cell('S47')
        c = X.nucleus_and_vertices()[0]
        w = X.width()
        ul = (c[0] - w/4, c[1] + w/4)
        dr = (c[0] + w/4, c[1] - w/4)
        Y = self.E.cell_from_region(ul, dr)
        self.assertEqual(X.suid, Y.suid)
    
        # The nuclei of two cells at the same level that have 
        # the same parent cell should yield the parent cell.
        A = self.E.cell('O301')
        B = self.E.cell('O305')
        should_get = 'O30'
        Ac = A.nucleus_and_vertices()[0]
        Bc = B.nucleus_and_vertices()[0]
        got = str(self.E.cell_from_region(Ac, Bc))
        self.assertEqual(got, should_get)
    
        # Should work on defunct cases.
        # Rectangle with only two vertices.
        A = self.E.cell('O301')
        B = self.E.cell('O302')
        should_get = 'O30'
        Ac = A.nucleus_and_vertices()[0]
        Bc = B.nucleus_and_vertices()[0]
        got = str(self.E.cell_from_region(Ac, Bc))
        self.assertEqual(got, should_get)
        # Rectangle with only one vertex.
        A = self.E.cell('O301')
        B = self.E.cell('O301')
        should_get = 'O301' + '4'*(self.E.max_level - 3)
        Ac = A.nucleus_and_vertices()[0]
        Bc = B.nucleus_and_vertices()[0]
        got = str(self.E.cell_from_region(Ac, Bc))
        self.assertEqual(got, should_get)  
        
        # Test ellipsoidal regions.
        # An ellipsoidal cap contained in S4 but not S44.
        X = self.E.cell('S4')
        vertices = X.nucleus_and_vertices(surface='ellipsoid')[1:]
        phi = vertices[0][1] - pi/100
        ul = (-pi, -pi/2)
        dr = (-pi, phi)
        got = str(self.E.cell_from_region(ul, dr, surface='ellipsoid'))         
        should_get = 'S4'
        self.assertEqual(got, should_get)
        
        # An ellipsoidal rectangle sharing two vertices with an
        # ellipsoidal cell.
        X = self.E.cell('N7')   # A four-edged polar ellipsoidal cell
        vertices = X.nucleus_and_vertices(surface='ellipsoid')[1:]
        fuzz = 1e-3 # To avoid rounding errors
        ul = max(vertices[0][0], vertices[1][0]) + fuzz, \
             max(vertices[0][1], vertices[3][1]) - fuzz
        dr = min(vertices[2][0], vertices[3][0]) - fuzz, \
             min(vertices[1][1], vertices[2][1]) + fuzz
        got = str(self.E.cell_from_region(ul, dr, surface='ellipsoid'))     
        should_get = str(X)
        self.assertEqual(got, should_get)
#------------------------------------------------------------------------------                           
if __name__ == '__main__':
    unittest.main()