r"""
This Python 2.7 module implements the SCENZ-Grid discrete global grid system (DGGS), which is an instance of the rHEALPix DGGS.

CHANGELOG:

- Alexander Raichev (AR), 2011-07-11: Initial version
- AR, 2011-07-26: Added new functions. 
  Improved suids.
- AR, 2011-08-05: Extended function capababilities 
  to the WGS84 ellipsoid. Improved suids.
- AR, 2011-08-31: Added optional central meridian shift in 
  HEALPix and rHEALPix projections. 
  Added optional polar squares shift in rHEALPix projection.
  Added Earth class and restructured and tidied everything.
  Fixed rounding bug in HEALPix projection. 
- AR, 2011-09-01: Added option for entering angles in degrees
  in projections and set degrees as default mode in Earth class.
  Fixed a rounding bug in auth_lat().
- AR, 2011-09-27: Split projection functions and grid 
  functions into separate files. 
- AR, 2011-10-07: Improved variable and method suids. 
  Tested more. Documented more.
- AR, 2012-03-06: Refactored to update attribute and method suids and
  redo cell ordering.
- AR, 2012-03-12: Refactored some more and added CellFamily class.
- AR, 2012-03-28: Improved CellFamily methods.
- AR, 2012-04-05: Reformatted code to comform to 
  `standard Python style <http://www.python.org/dev/peps/pep-0008/>`_. 
  Fixed minimize(), union(), and intersect() to handle empty cell families properly. 
- AR, 2012-04-11: Improved index() and it's inverse.
- AR, 2012-04-13: Rewrote CellFamily to inherit from collections.MutableSequence. Added and optimized CellFamily methods.
- AR, 2012-04-16: Removed wildcard imports to avoid suidspace pollution.
- AR, 2012-04-23: Added a `filter_level` option to cell family intersectionand union operations.
- AR, 2012-04-25: Changed level 0 cell suids to 'N', 'O', ..., 'S' labeled left to right and top to bottom and changed subcell suids to 0, 1, ..., 8 labeled left to right and top to bottom. Consequently had to change the code that converts between cell suid and cell location.        
- AR, 2012-05-09: Changed 'name' to 'suid' and improved intersect_all() slightly.
- AR, 2012-06-08: Changed 'suid_xy' to 'suid_colrow' and cleaned up location(), cell_from_point(), and cell_from_region(). 
- AR, 2012-07-04: Simplified minimize() and eliminated preprocess().
- AR, 2012-07-05: Simplified total ordering definition in Cell class.
- AR, 2012-07-06: Added a width() method to Cell class.
- AR, 2012-09-17: Fixed rounding errors in nucleus_and_vertices() which occurred when interpolation > 0.
- AR, 2012-10-15: Changed names. Changed 'location' to 'ul_vertex', 'center' to 'nucleus', 'corners' to 'vertices'. Introduced the 'surface' keyword for cell operations that depend on the distinction between planar and ellipsoidal cells. - AR, 2012-10-24: Added ellipsoidal_shape() and centroid().
    
EXAMPLES::

Create a WGS84 ellipsoid model of the Earth::

    >>> E = Earth(ellps='WGS84', south=1, north=2)
    >>> print E
    ____Earth model____
    lengths measured in meters and angles measured in degrees
    ellipsoid: WGS84
        major radius: 6378137
        flattening factor: 0.00335281066475
        minor radius: 6356752.31425
        eccentricity: 0.0818191908426
        authalic sphere radius: 6371007.18092
    central meridian: 0
    north pole square position: 2
    south pole square position: 1
    max areal resolution: 1.0
    max cell level: 15
        
Use the HEALPix projection for this Earth to project the
(longitude, latitude) point (45, 60) onto a plane and back again::

    >>> p = E.healpix(45, 60); p
    (5003777.338885325, 6823798.3004998406)
    >>> E.healpix(*p, inverse=True)
    (45.0, 59.999999998490921)
 
Do the same but use the rHEALPix projection::

    >>> p = E.rhealpix(45, 60); p
    (5003777.338885325, 6823798.3004998406)
    >>> E.rhealpix(*p, inverse=True)
    (45.0, 59.999999998490921)
    
Create the 'Q3' cell of this Earth and find its four neighbors::

    >>> c = E.cell('Q3')
    >>> for k, v in c.neighbors().items():
    ...     print k, v
    ... 
    down Q6
    right Q4
    up Q0
    left P5
    
Find the nucleus and vertices of this cell and the level 2 cell 
that contains the nucleus::

    >>> nucleus = c.nucleus_and_vertices()[0]; nucleus
    (1667925.7796284417, 0.0)
    >>> print E.cell_from_point(3, nucleus)
    Q344
        
NOTES:

All lengths are measured in meters and all angles are measured in radians 
unless indicated otherwise. 
Points lying on the plane are given in rectangular (horizontal, vertical) coordinates, and points lying on the ellipsoid are given in geodetic (longitude, latitude) coordinates unless indicated otherwise.
Below, GGS abbreviates global grid system.

Except when manipulating positive integers, I avoid the modulo function '%' 
and insted write everything in terms of 'floor()'.
This is because Python interprets the sign of '%' differently than
Java or C, and I don't want to confuse people who are translating this code
to those languages.  
"""
#*****************************************************************************
#       Copyright (C) 2012 Alexander Raichev <tortoise.said@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http: //www.gnu.org/licenses/
#*****************************************************************************

# Third-party module imports
from scipy import array, base_repr, ceil, log, pi, sqrt, integrate, floor
# Python module imports
from string import join
from itertools import product
# My module imports
import projections as projs 

def scale(x, factor, inverse=False):
    r"""
    Multiply (all the elements of) `x` by the number `factor` or by 1/`factor` 
    if `inverse=True`.
    Here `x` can be a number or a tuple of numbers.
    
    EXAMPLES::
    
        >>> scale(5, 2)
        10
        >>> scale(5, 2,inverse=True)
        2.5
        >>> scale((5, 3.1),2)
        (10, 6.2)
        >>> scale((5, 3.1),2, inverse=True)
        (2.5, 1.55)
        
    """
    #print "Scale %f" % (factor);
    #print x;
    if inverse:
        c = 1/float(factor)
    else:
        c = factor        
    if isinstance(x, tuple):
        y = tuple([c*xx for xx in x])
    else:
        y = c*x
    return y


class Earth(object):
    r"""
    Represents a model of the Earth its planar and ellipsoidal GGSs.
    """
    # Earth ellipsoids to choose from.
    unit_sphere = {'figure': 'sphere', 'R': 1}
    WGS84 = {'figure': 'ellipsoid', 'a': projs.WGS84_A, 'f': projs.WGS84_F, 
             'b': projs.WGS84_B, 'e': projs.WGS84_E, 'R': projs.WGS84_R}
    WGS84_asphere = {'figure': 'sphere', 'R': WGS84['R']}
    ellipsoids = {'unit_sphere': unit_sphere, 'WGS84': WGS84, 
                  'WGS84_asphere': WGS84_asphere}
    
    # Level 0 cell suids, which are anamolous.
    CELLS0 = ['N','O','P','Q','R','S'] 
    
    # Dictionary of the ordering (Morton order) of subcells within a cell 
    # in terms of the horizontal-vertical (column-row) coordinates of 
    # the subcells with the origin placed at the downleft corner. 
    # Here's a diagram of the ordering and coordinates:
    #
    # 0 | 0 1 2
    # 1 | 3 4 5
    # 2 | 6 7 8
    #   --------
    #     0 1 2  
    #
    # The subcells are being labeled from left to right and top to
    # bottom to match up with how files are read.
    child_order = {}   
    for (col, row) in product(range(3), repeat=2):    
        order = col + row*3   
        # Handy to have both coord and order as dictionary keys.
        child_order[(col,row)] = order
        child_order[order] = (col, row)

    def __init__(self, ellps='WGS84', degrees=True, lon0=0, north=0, south=0,
                 max_areal_resolution=1.0):
        # A little error checking.
        assert ellps in Earth.ellipsoids, \
        "ellps must be one of " + str(Earth.ellipsoids.keys())
        self.ellps = ellps
        
        # Set the shape attributes of this Earth model based on ellps 
        # and Earth.shapes.
        for attr in Earth.ellipsoids[ellps]:
            setattr(self, attr, Earth.ellipsoids[ellps][attr])
        
        # Set the degrees and prime meridian attributes.    
        self.degrees = degrees
        self.lon0 = lon0

        # Set the grid attributes.
        # Where rHEALPix assembles the north polar caps:
        self.north = north % 4  
        # Where rHEALPix assembles the south polar caps:
        self.south = south % 4  
        self.max_areal_resolution = max_areal_resolution
        # Find the maximum grid level needed to have cells of area
        # at most max_areal_resolution.
        self.max_level = int(ceil(
          log(self.R**2*(2*pi/3)/max_areal_resolution)/\
          (2*log(3)) ))
        
        # Store the upper left vertices of the level 0 cells 
        # in the rHEALPix grid hierarchy for this Earth.
        # The default layout by CELLS0 index is
        #
        # 0
        # 1 2 3 4
        # 5.
        #
        CELLS0 = Earth.CELLS0
        ul_vertex = {  # Location for radius = 1
            CELLS0[0]: (-pi + self.north * pi / 2, 3 * pi / 4),
            CELLS0[1]: (-pi, pi / 4),
            CELLS0[2]: (-pi / 2, pi / 4),
            CELLS0[3]: (0, pi / 4),
            CELLS0[4]: (pi / 2, pi / 4),
            CELLS0[5]: (-pi + self.south * pi / 2, -pi / 4)
        }
        # Scale up ul_vertex for radius = self.R
        self.ul_vertex = dict([(x, scale(ul_vertex[x], self.R)) 
                              for x in ul_vertex.keys()])
        
        # Initialize atomic neighbor relationships among cells.    
        # Dictionary of up, right, down, and left neighbors of 
        # level 0 cells and subcells 0--8, aka the atomic neighbors. 
        # Based on the layouts
        #                   
        #   0                             
        #   1 2 3 4   (but folded into a cube) and
        #   5                              
        #                   
        #   0 1 2 
        #   3 4 5 
        #   6 7 8.
        # 
        an = {} 
        # Neighbors of CELLS0[1], ..., CELLS0[4]    
        an[CELLS0[1]] = {'left': CELLS0[4], 'right': CELLS0[2], 
                         'down': CELLS0[5], 'up': CELLS0[0]}
        an[CELLS0[2]] = {'left': CELLS0[1], 'right': CELLS0[3], 
                         'down': CELLS0[5], 'up': CELLS0[0]}
        an[CELLS0[3]] = {'left': CELLS0[2], 'right': CELLS0[4], 
                         'down': CELLS0[5], 'up': CELLS0[0]}
        an[CELLS0[4]] = {'left': CELLS0[3], 'right': CELLS0[1], 
                         'down': CELLS0[5], 'up': CELLS0[0]}
        # Neighbors of CELLS0[0] and CELLS0[5] depend on 
        # volues of north and south, respectively.
        nn = self.north
        an[CELLS0[0]] = {'down': CELLS0[(nn + 0) % 4 + 1],   
                         'right': CELLS0[(nn + 1) % 4 + 1],
                         'up': CELLS0[(nn + 2) % 4 + 1],
                         'left': CELLS0[(nn + 3) % 4 + 1]}
        ss = self.south
        an[CELLS0[5]] = {'up': CELLS0[(ss + 0) % 4 + 1],   
                         'right': CELLS0[(ss + 1) % 4 + 1],
                         'down': CELLS0[(ss + 2) % 4 + 1],
                         'left': CELLS0[(ss + 3) % 4 + 1]}

        # Neighbors of 0, 1, ..., 8.
        for i in range(9):
            an[i] = {'left': i - 1, 'right': i + 1, 
                     'up': (i - 3) % 9, 'down': (i + 3) % 9}
        # Adjust left and right edge cases.
        for i in [0, 3, 6]:
            an[i]['left'] = an[i]['left'] + 3
        for i in [2, 5, 8]:
            an[i]['right'] = an[i]['right'] - 3
        self.atomic_neighbors = an          

    def __str__(self):
        if self.degrees:
            angles = 'degrees'
        else:
            angles = 'radians'
        result =["____Earth model____"]
        result.append("lengths measured in meters and angles measured in " 
                      + angles)
        result.append("ellipsoid: " + self.ellps)
        if self.figure == 'sphere':
            result.append("    radius: %s" % self.R)
        else:
            result.append("    major radius: %s" % self.a)
            result.append("    flattening factor: %s" % self.f)
            result.append("    minor radius: %s" % self.b)
            result.append("    eccentricity: %s" % self.e)
            result.append("    authalic sphere radius: %s" % self.R)
        result.append("central meridian: %s" % self.lon0)
        result.append("north pole square position: %s" % self.north)
        result.append("south pole square position: %s" % self.south)
        result.append("max areal resolution: %s" % self.max_areal_resolution)
        result.append("max cell level: %s" % self.max_level)
        return "\n".join(result)
        
    def healpix(self, u, v, inverse=False):
        r"""
        Return the HEALPix projection of point `(u, v)` (or its inverse if 
        `inverse=True`) for this Earth model.
        
        EXAMPLES::
        
            >>> E = Earth('unit_sphere')
            >>> p = (0, 0)
            >>> q = (-180, 90)
            >>> print E.healpix(*p) == (0, 0)
            True
            >>> print E.healpix(*q) == (-3 * pi / 4, pi / 2)
            True
            
        """
        if self.figure == 'sphere':
            if not inverse:
                return scale(projs.healpix_sphere(u, v, lon0=self.lon0, 
                                                 degrees=self.degrees),
                             self.R)
            else:
                return projs.healpix_sphere(*scale((u, v), self.R,
                                                  inverse=True),
                                           lon0=self.lon0, 
                                           degrees=self.degrees,
                                           inverse=True)
        else:
            if not inverse:
                return scale(projs.healpix_ellipsoid(u, v, e=self.e, 
                                                    lon0=self.lon0,
                                                    degrees=self.degrees),
                             self.R)
            else:
                return projs.healpix_ellipsoid(*scale((u, v), self.R, 
                                              inverse=True),
                                              e=self.e, lon0=self.lon0, 
                                              degrees=self.degrees,
                                              inverse=True)
            
            
    def rhealpix(self, u, v, inverse=False):
        r"""
        Return the rHEALPix projection of the point `(u, v)` (or its inverse if 
        `inverse=True`) appropriate to this Earth.
        
        EXAMPLES::
        
            >>> E = Earth('unit_sphere')
            >>> p = (0, 0)
            >>> q = (-180, -90)
            >>> print E.rhealpix(*p) == (0, 0)
            True
            >>> print E.rhealpix(*q) == (-3 * pi / 4, -pi / 2)
            True
        
        """
        if self.figure == 'sphere':
            if not inverse:
                return scale(projs.rhealpix_sphere(u, v, lon0=self.lon0, 
                                                  degrees=self.degrees,
                                                  north=self.north, 
                                                  south=self.south),
                             self.R)
            else:
                return projs.rhealpix_sphere(*scale((u, v), self.R,
                                                   inverse=True),
                                            lon0=self.lon0, 
                                            degrees=self.degrees,
                                            north=self.north, 
                                            south=self.south,
                                            inverse=True)
        else:
            if not inverse:
                return scale(projs.rhealpix_ellipsoid(u, v, e=self.e, 
                                                     lon0=self.lon0,
                                                     degrees=self.degrees,
                                                     north=self.north, 
                                                     south=self.south),
                             self.R)
            else:
                return projs.rhealpix_ellipsoid(*scale((u, v), self.R,
                                               inverse=True),
                                               e=self.e, lon0=self.lon0, 
                                               degrees=self.degrees,
                                               north=self.north, 
                                               south=self.south,
                                               inverse=True)
            
    def combine_caps(self, u, v, inverse=False):
        r"""
        Return the combine_caps() transformation of the point `(u, v)` 
        (or its inverse if `inverse=True`) appropriate to this Earth.
        It maps the HEALPix projection to the rHEALPix projection.
        
        EXAMPLES::
        
            >>> E = Earth('unit_sphere')
            >>> p = (0, 0)
            >>> q = (-pi / 4, pi / 2)
            >>> print E.combine_caps(*p) == (0, 0)
            True
            >>> print E.combine_caps(*q) == (-3 * pi / 4, pi / 2)
            True
            
        """
        return scale(projs.combine_caps(*scale((u, v), self.R, inverse=True), 
                                       north=self.north, south=self.south,              
                                       inverse=inverse),
                     self.R)
                     
    def project(self, p, source='lonlat', target='rhealpix'):
        r"""
        Transform a point `p` from the image of the `source` 
        projection to the image of the `target` projection.   
    
        INPUT:

         - `p` - A point in geographic space, in the HEALPix projection,
           or in the rHEALPix projection. 
         - `source, target` - One of 'lonlat', 'healpix', or 'rhealpix'.
         
        EXAMPLES::
        
            >>> from numpy import rad2deg, arcsin
            >>> E = Earth('unit_sphere')
            >>> p = (-pi / 2, 3 * pi / 4)
            >>> phi0 = rad2deg(arcsin(2.0 / 3))
            >>> print phi0
            41.8103148958
            >>> print E.project(p, source='rhealpix', target='lonlat') == (0, phi0)
            True
            
        """
        if source == 'lonlat':
            if target == 'healpix':
                return self.healpix(*p)
            elif target == 'rhealpix':
                return self.rhealpix(*p)
            else:
                return p
        elif source == 'healpix':
            if target == 'lonlat':
                return self.healpix(*p, inverse=True)
            elif target == 'rhealpix':
                return self.combine_caps(*p)
            else:
                return p
        elif source == 'rhealpix':
            if target == 'lonlat':
                return self.rhealpix(*p, inverse=True)
            elif target == 'healpix':
                return self.combine_caps(*p, inverse=True)
            else:
                return p
        else:
            return p

    def cell(self, suid=None, level_order_index=None, post_order_index=None):
        r"""
        Return a cell of this Earth model either from its SUID or 
        from its level and index. The cell can be interpreted as lying on the 
        plane or on the ellipsoid, depending on the `surface` option given 
        in the instance methods below.
        
        EXAMPLES::
        
            >>> E = Earth()
            >>> c = E.cell('N45')
            >>> print isinstance(c, Cell)
            True
            >>> print c
            N45
            
        """
        return Cell(self, suid, level_order_index, post_order_index)
                
    def grid(self, level):
        r"""
        Generator function for all the cells at level `level`.
        
        EXAMPLES::
        
            >>> E = Earth()
            >>> g = E.grid(0)
            >>> print [str(x) for x in g]
            ['N', 'O', 'P', 'Q', 'R', 'S']
            
        """
        suid = [Earth.CELLS0[0]] + [0 for i in range(level)]
        c = self.cell(suid)
        yield c
        cs = c.successor(level)
        while cs:
            yield cs
            cs = cs.successor(level)
    
    def num_cells(self, level_a, level_b=None, subcells=False):
        r"""
        Return the number of cells at levels `level_a`
        to `level_b` (inclusive).
        Assume `level_a <= level_b`.
        If `subcells=True`, then return the number of subcells at levels
        `level_a` to `level_b` (inclusive) of a cell at level `level_a`.
        If `level_b=None` and `subcells=False, then return the number of
        cells at level `level_a`.
        If `level_b=None` and `subcells=True`, then return the number of
        subcells from level `level_a` to level `self.max_level`.
        
        EXAMPLES::
            
            >>> E = Earth()
            >>> E.num_cells(0)
            6
            >>> E.num_cells(0, 1)
            60
            >>> E.num_cells(0, subcells=True)
            231627523606480
            >>> E.num_cells(0, 1, subcells=True)
            10
            >>> E.num_cells(5, 6, subcells=True)
            10
                        
        """
        if subcells:
            if (level_b is None) or (level_b < level_a):
                level_b = self.max_level
            num = (9**(level_b - level_a + 1) - 1) / 8
        else:
            if (level_b is None) or (level_b < level_a):
                level_b = level_a
            num = 6 * (9**(level_b + 1) - 9**(level_a)) / 8
        return num
          
    def sample(self, cells, k):
        r"""
        Return a `k` length list of unique elements chosen at random without
        replacement from the cell population sequence `cells`. 
        """
        import random
        return random.sample(cells, k)
            
    def random_cell(self, level=None):
        r"""
        Return a cell chosen of level `level` chosen uniformly at random among
        all level `level` cells.
        If `level=None`, then `level` is first chosen uniformly randomly in
        [0,..,self.max_level]
        """
        from random import randint
        
        if level == None:
            level = randint(0, self.max_level)
        suid = []
        suid.append(Earth.CELLS0[randint(0, 6)])
        for i in range(1, level + 1):
            suid.append(randint(0, 8))
        return Cell(self, suid)    
        
    def cell_width(self, level, surface='plane'):
        r"""
        For this Earth model, return the width of a planar cell at level 
        `level`.
        If `surface != 'plane'`, then return None, because ellipsoidal cells
        don't have constant width.
        
        EXAMPLES::
        
            >>> E = Earth('unit_sphere')
            >>> E.cell_width(0)
            1.5707963267948966
            >>> E.cell_width(1)
            0.5235987755982988
        """
        if surface == 'plane':
            return self.R *(pi/2)*3**(-level)
    
    def cell_area(self, level, surface='plane'):
        r"""
        For this Earth model, return the area of a planar or ellipsoiadl cell
        at level `level`.
        """
        w = self.cell_width(level)
        if surface == 'plane':
            return w**2
        else:
            return 8/(3*pi)*w**2
            
    def interval(self, a, b):
        r"""
        Generator function for all the level `max(a.level, a.level)` 
        cells between cell `a` and cell `b` (inclusive and with respect to the     
        postorder ordering on cells).
        Note that `a` and `b` don't have to lie at the same level.
        
        EXAMPLES::
        
            >>> E = Earth()
            >>> a = E.cell('N10')
            >>> b = E.cell('N3')
            >>> print [str(z) for z in E.interval(a, b)]
            ['N10', 'N11', 'N12', 'N13', 'N14', 'N15', 'N16', 'N17', 'N18', 'N20', 'N21', 'N22', 'N23', 'N24', 'N25', 'N26', 'N27', 'N28', 'N30', 'N31', 'N32', 'N33', 'N34', 'N35', 'N36', 'N37', 'N38']
        
        """
        # Choose the starting cell, which might not be A.
        level = max(a.level, b.level)
        if a.level < level:
            cell = a.successor(level)
        else: 
            cell = Cell(self, a.suid[:level + 1])
        while cell <= b:
            yield cell
            cell = cell.successor(level)
            
    def cell_from_point(self, level, p, surface='plane'):
        r"""          
        For this Earth model, return the planar or ellipsoidal cell at level 
        `level` that contains the point `p`.
        
        INPUT:
        
        - `level` - Cell level.
        - `p` - A point on the plane or ellipsoid given in rectangular or 
          geodetic coordinates, respectively.
        - `surface` - (Optional; default='plane') One of 'plane' or 
          'ellipsoid', indicating whether the given point and returned cell 
          both lie on the plane or on the ellipsoid.
          
        EXAMPLES::
        
            >>> E = Earth('WGS84')
            >>> p = (0, 0)
            >>> c = E.cell_from_point(2, p)
            >>> print c
            Q33
            
        """
        # Get the rectangular coordinates of p.
        if surface == 'plane':
            x, y = p
        else:
            x, y = self.project(p, 'lonlat', 'rhealpix')

        # Determine the level 0 cell c0 that (x, y) lies in, 
        # since level 0 cells are anamolous.
        R = self.R
        if y > R * pi / 4:
            suid = [Earth.CELLS0[0]]
        elif y < -R * pi / 4:
            suid = [Earth.CELLS0[5]]            
        elif x < -R * pi / 2:
            suid = [Earth.CELLS0[1]]
        elif x >= -R * pi / 2 and x < 0:
            suid = [Earth.CELLS0[2]]
        elif x >= 0 and x < R * pi / 2:
            suid = [Earth.CELLS0[3]]
        else:
            suid = [Earth.CELLS0[4]]
        
        if level == 0:
            # Done.
            return Cell(self, suid)
            
        # Compute the horizontal and vertical distances between (x, y) and
        # the ul_vertex of c0 as fractions of the width of c0.
        w = self.cell_width(0)
        dx = abs(x - self.ul_vertex[suid[0]][0])/w
        dy = abs(y - self.ul_vertex[suid[0]][1])/w
        if dx == 1:
            # Border case. Take a smidgen off dx.
            dx -= self.max_areal_resolution/(2*R*w)
        if dy == 1:
            # Border case. Take a smidgen off yd.
            dy -= self.max_areal_resolution/(2*R*w)
        
        # Compute the ternary expansions of dx and dy and truncate them
        # at index level to get the column and row SUIDs of 
        # the level level cell c containing (x,y).   
        suid_col = base_repr(int(float(str(dx * 3**level))), 3)
        suid_row = base_repr(int(float(str(dy * 3**level))), 3) 
        # Using int(float(str(.))) instead of the straightforward int(.),
        # because the latter was giving me rounding errors.
        # Prefix with the appropriate amount of zeros.
        suid_col = '0'*(level - len(suid_col)) + suid_col
        suid_row = '0'*(level - len(suid_row)) + suid_row
        
        # Use the column and row SUIDs of c to get the SUID of c.
        for i in range(level):
            suid.append(Earth.child_order[(int(suid_col[i]), 
                                           int(suid_row[i]))])
        return Cell(self, suid) 
        
    def cell_from_region(self, ul, dr, surface='plane'):
        r"""
        Return the smallest planar or ellipsoidal cell wholly containing 
        the region bounded by the axis-aligned rectangle with upper left 
        and lower right vertices given by the the points `ul` and `dr`, 
        respectively.  
        If such as cell does not exist, then return None.
        If `surface='plane'`, then `ul` and `dr` and the returned cell
        are interpreted as lying on the plane.
        If `surface='ellipsoid'`, then they are interpreted as lying on
        the ellipsoid and the bounding rectangle is interpreted as being
        longitude-latitude aligned.    
        To specify an ellipsoidal cap region (whose centroid is the north
        or south pole and whose boundary is a parallel of latitude),
        set ul = (-pi, pi/2) or (-pi, -pi/2), indicating a north or south pole,
        respectively, and set dr = (-pi, phi), where phi is the latitude of
        the boundary.
        Works on defunct (point and line) planar rectangles too.
        
        EXAMPLES::
        
            >>> E = Earth('unit_sphere')
            >>> p = (0, 0)
            >>> q = (pi/4, pi/4)
            >>> c = E.cell_from_region(p, q)
            >>> print c
            Q
            
        """    
        if surface != 'plane':
            # Compute planar ul and dr as follows.
            # Get all four vertices of the ellipsoidal rectangle.
            if abs(ul[1]) == pi/2:
                # Ellipsoidal cap case.
                vertices = [(-3*pi/4, dr[1]), (-pi/4, dr[1]), (pi/4, dr[1]), 
                           (3*pi/4, dr[1])]
            else:
                vertices = [ul, (ul[0], dr[1]), dr, (dr[0], ul[1])]
            # Project the vertices onto the plane.
            vertices = [self.project(p, 'lonlat', 'rhealpix') for p in
                        vertices]
            # Find the upper left and lower right vertices of the 
            # planar bounding rectangle.
            ul = (min([p[0] for p in vertices]), max([p[1] for p in vertices]))
            dr = (max([p[0] for p in vertices]), min([p[1] for p in vertices]))
        
        # Find the level max_level cells containing ul and dr.
        level = self.max_level
        ul_cell = self.cell_from_point(level, ul)
        dr_cell = self.cell_from_point(level, dr) 
        ul_suid = ul_cell.suid
        dr_suid = dr_cell.suid
            
        # Find the longest common prefix of ul_suid and dr_suid.
        least = level + 1   # Default if the suids agree everywhere 
        for i in range(level + 1):
            if (ul_suid[i] != dr_suid[i]): 
                least = i
                break
        if least == 0:
            # No one cell contains R.
            return None
        else:
            return self.cell(ul_suid[:least])

from functools import total_ordering
@total_ordering     
# Defines a total ordering on Cell instances from the definitions of __eq__(), __ne__(), and __le__() below.
class Cell(object):
    r"""
    Represents a cell of the planar or ellipsoidal global grid system (GGS)
    relative to a given Earth model.
    Cell identifiers are of the form (p_0, p_1,...,p_l), where p_0 is one of 
    the characters 'A', 'B', 'C', 'D', 'E', 'F' and p_i for i > 0 is one of
    the characters '0', '1',..., '8'.
    """
    @staticmethod
    def suid_from_index(earth, index, order='level'):
        r"""
        Return the suid of a cell from its index.
        The index is according to the cell ordering `order`,
        which can be 'level' (default) or 'post'.
        See the `index()` docstring for more details on orderings.
        For internal use.
        """
        from math import log
        if order == 'post':
            # Compute suid from post order index one character at a time.
            suid = []
            p = index
            def num(k):
                return earth.num_cells(level_a=k, subcells=True)
            # Consider the tree T of all cells.
            # The indices of the cells in the six subtrees rooted at the 
            # level 0 cells lie in the intervals 
            # [0, num(0)), [num(0), 2 * num(0)),..., [5 * num(0), 6 * num(0)),
            # respectively.
            # So computing p // num(0) gives us the first character of
            # the suid of our cell c. 
            # Setting p = p % num(0) gives us the post order index
            # of c relative to the subtree rooted at cell suid[0].
            # The relative indices of the cells in the nine subtree rooted at
            # the nine children of suid[0] lie in the intervals
            # [0, num(1)), [num(1), 2 * num(1)),..., [8 * num(1), 9 * num(1)),
            # respectively.
            # So computing p // num(1) gives us suid[1].
            # Repeating this procedure until p == num(i) - 1 for some i 
            # (which will happen when i = max_level at the latest) 
            # gives us all the characters of suid.
            for i in range(earth.max_level + 1):
                n = num(i)
                q, r = divmod(p, n)
                suid.append(q)
                p = r
                if p == n - 1:
                    break
            # Replace first digit with appropriate letter.
            suid[0] = Earth.CELLS0[suid[0]]
            suid = tuple(suid)
        else:
            # Compute suid from level order index.
            def ind(k):
                r"""
                Return the level order index of the first cell at
                level k.
                """
                return 6 * ((9**k - 1)/8) 
            # The cells at level L have indices in the interval
            # [ind(L), ind(L + 1)).      
            k = int(log(8*(index/6.0) + 1, 9))  
            # k = L  or L + 1. Find out which one.
            remainder = index - ind(k) 
            if  remainder >= 0:
                # Then k = L
                L = k
            else:
                # Then k = L + 1
                L = k - 1
                remainder = index - ind(L)
            # Now compute cell suid from remainder.
            suid = base_repr(remainder, 9)
            # If necessary, prepend with zeros to get a length L + 1 string.
            suid = '0'*(L + 1 - len(suid)) + suid
            suid = [int(s) for s in suid]
            # Replace first digit with appropriate letter.
            suid[0] = Earth.CELLS0[suid[0]]
            suid = tuple(suid)
        return suid
         
    def __init__(self, earth, suid=None, level_order_index=None, 
                 post_order_index=None):
        r"""
        Create a cell either from its suid or from its level order or
        post order index.  The cell can be interpreted as lying on the plane
        or on the ellipsoid, depending on the `surface` option given in the
        instance methods below.   
        
        EXAMPLES::
        
            >>> E = Earth(ellps='unit_sphere')
            >>> c = Cell(E, 'N3')
            >>> print c
            N3
            >>> c = Cell(E, level_order_index=9)
            >>> print c
            N3
            >>> c = Cell(E, post_order_index=9)
            >>> print c
            N
            
        """
        self.earth = earth  # Earth model that self is a part of.
        max_level = earth.max_level
        self.suid = ()      # Spatially unique identifier of self.
        self.level = None   # Level of self in grid hierarchy.
        if suid:
            # A little error checking.
            assert (isinstance(suid, list) or isinstance(suid, tuple) or 
            isinstance(suid, str)), \
            'Cell suid must be a list, tuple, or string. Got %s.' % suid
            assert suid[0] in Earth.CELLS0, \
            'suid[0] must lie in %s. Got %s.' % (Earth.CELLS0, suid[0])
            assert (len(suid) > 0) and (len(suid) - 1 <= earth.max_level), \
            'Need 0 < len(suid) - 1 <= %s. Got %s.' % (earth.max_level, suid)
            
            self.suid = [suid[0]] + [int(n) for n in suid[1:]] 
            self.suid = tuple(self.suid)                
        elif level_order_index:          
            self.suid = Cell.suid_from_index(self.earth, level_order_index,       
                                             order='level') 
        elif post_order_index:
            self.suid = Cell.suid_from_index(self.earth, post_order_index,       
                                             order='post') 
        self.level = len(self.suid) - 1 
        
    def __nonzero__(self):
        r"""
        Return my truth value (True or False).
        """
        return bool(self.suid)
        
    def __str__(self):
        # Convert self.suid to a string without commas or parentheses.
        return join(map(str, self.suid), '')

    def __eq__(self, other):
        r"""
        The equality relation on cells.
        Cells are equal if they have the same earth and suid.
        """   
        if (other is not None) and (self.earth == other.earth) and\
          (self.suid == other.suid):
            return True
        else:
            return False

    def __ne__(self, other):
        r"""
        The inequality relation on cells. 
        Since Python 2.6 doesn't automatically create reverse relations
        from given ones, I must define this seemingly redundant relation.
        """   
        return not self == other

    def __le__(self, other):
        r"""
        The (strictly) less-than relation on cells.
        Derived from the post order traversal of the tree T of all cells
        defined in the `index()` docstring. 
        Return True if (`self.suid < other.suid` and 
        `self.suid` is not a prefix of `other.suid`) or 
        `self` is a subcell of `other`.
        Here < is the lexicographic order.
        Returns False otherwise.
        """
        s = str(self)
        t = str(other)
        if (s <= t and not t.startswith(s)) or s.startswith(t):
            return True
        else:
            return False
   
    def str9(self):
        r"""
        Return the suid of `self` as a base 9 numeral.
        Just need to convert the leading character to a digit between 0 and 8.
        
        EXAMPLES::
        
            >>> E = Earth()
            >>> c = Cell(E, 'N130534')
            >>> print c.str9()
            0130534
            
        """
        suid = list(self.suid)
        suid[0] = Earth.CELLS0.index(suid[0])
        return join(map(str, suid), '')
            
    def index(self, order='level'):
        r"""
        Return the index of `self` when it's ordered according to `order`.
        Here `order` can be 'level' (default) or 'post'.        
        Indices start at 0.
        The empty cell has index None.

        The ordering comes from the way of traversing the tree T 
        of all cells defined as follows.
        The root of T is a non-cell place holder. 
        The children of the root are the cells A < B < ... < F.
        The children of a cell in T with suid s are s0 < s1 < ... < s8.
        
        The level order index of a nonempty cell is its position 
        (starting from 0) in the level order traversal of T starting at cell A.
        
        The post order index of a nonempty cell is its position 
        (starting from 0) in the post order traversal of T.
         
        EXAMPLES::
        
            >>> E = Earth(ellps='unit_sphere')
            >>> c = Cell(E, 'N2')
            >>> print c.index(order='level')
            8
            >>> print c.index(order='post')
            2
            
        """   
        L = self.level     
        mL = self.earth.max_level
        if not self.suid:
            return None
        if order == 'post':
            s = self.str9()
            def num(k):
                return self.earth.num_cells(level_a=k, subcells=True)
            result = sum(int(s[i]) * num(i) for i in range(L + 1)) +\
                     num(L) - 1
        else:
            # Calculate level order index.
            result = self.earth.num_cells(level_a=0, level_b = L - 1) +\
                     int(self.str9(), 9)
        return result

    def suid_colrow(self):
        r"""
        Return the pair of row- and column-suids of `self`, each as tuples.
        
        EXAMPLES::
        
            >>> E = Earth()
            >>> c = Cell(E, 'N73')
            >>> xn, yn = c.suid_colrow()
            >>> print xn == ('N', 1, 0)
            True
            >>> print yn == ('N', 2, 1)
            True
            
        """
        suid_col = [self.suid[0]]
        suid_row = [self.suid[0]]
        for n in self.suid[1:]:
            col, row = Earth.child_order[n]
            suid_col.append(col)
            suid_row.append(row)
        return tuple(suid_col), tuple(suid_row) 
    
    def width(self, surface='plane'):
        r"""
        Return the width of this cell. If `surface='ellipsoid'`, then
        return None, because ellipsoidal cells don't have a fixed
        width.
        
        EXAMPLES::
        
            >>> c = Cell(Earth('unit_sphere'), 'N8')
            >>> print c
            N8
            >>> c.width() == pi/2*3**(-1)
            True
            
        """
        return self.earth.cell_width(self.level, surface=surface)      
              
    def area(self, surface='plane'):
        r"""
        Return the area of this cell.
        """
        return self.earth.cell_area(self.level, surface=surface)      
        
    def successor(self, level=None):
        r"""
        Return the least level `level` cell greater than `self`.
        Note: `self` need not be a level `level` cell.
        
        EXAMPLES::
        
            >>> c = Cell(Earth(), 'N82')
            >>> print str(c.successor())
            N83
            >>> print str(c.successor(0))
            O
            >>> print str(c.successor(1))
            O0
            >>> print str(c.successor(3))
            N830
            
        """
        suid = list(self.suid)
        if level is None:
            level = self.level
        if level < self.level:
            # Truncate suid at level level and return its successor. 
            return Cell(self.earth, suid[:level + 1]).successor()
        elif level > self.level:
            # Find the level self.level successor of suid  
            # and pad it with zeros.
            suid = list(self.successor().suid) + \
                   [0 for i in range(level - self.level)]
            return Cell(self.earth, suid)
        
        # Can now assume level = self.level.
        # First, find the greatest index i such that suid[i] != 8.
        greatest = 0
        for i in reversed(range(1, level + 1)):
            if suid[i] != 8:
                greatest = i
                break
                
        # Second, increment suid[greatest] and append all zeros
        # if possible.
        if greatest == 0:
            # suid[greatest] is a letter.
            i = Earth.CELLS0.index(suid[greatest])
            if i == 5:
                # End of the line. No successor.
                return None
            else:
                suid = [Earth.CELLS0[i + 1]] + [0 for j in range(level)]
        else:
            # suid[greatest] is a number < 8.
            suid = suid[0:greatest] + [suid[greatest] + 1] + \
                   [0 for j in range(level - greatest)]
        return Cell(self.earth, suid)
              
    def predecessor(self, level=None):
        r"""
        Return the greatest level `level` cell less than `self`.
        Note: `self` need not be a level `level` cell.
        
        EXAMPLES::
        
            >>> c = Cell(Earth(), 'N08')
            >>> print str(c.predecessor())
            N07
            >>> print str(c.predecessor(0))
            None
            >>> print str(c.predecessor(1))
            None
            >>> print str(c.predecessor(3))
            N088
            
        """
        suid = list(self.suid)
        if level is None:
            level = self.level
        if level < self.level:
            # Return predecessor of suid[:level + 1] 
            return Cell(self.earth, suid[:level + 1]).predecessor()
        elif level > self.level:
            # Return suid padded with eights.
            suid = suid + [8 for i in range(level - self.level)]
            return Cell(self.earth, suid)
        
        # Can now assume level = self.level.
        # Find the predecessor of suid.
        # First, find the greatest index i such that suid[i] != 0.
        greatest = 0
        for i in reversed(range(1, level + 1)):
            if suid[i] != 0:
                greatest = i
                break

        # Second, decrement suid[greatest] and append all eights
        # if possible.
        if greatest == 0:
            # suid[greatest] is a letter.
            i = Earth.CELLS0.index(suid[greatest])
            if i == 0:
                # End of the line. No predecessor.
                return None
            else:
                suid = [Earth.CELLS0[i - 1]] + [8 for i in range(level)]
        else:
            # nome[greatest] is a number > 0.
            suid = suid[0:greatest] + [suid[greatest] - 1] + \
                   [8 for i in range(level - greatest)]
        return Cell(self.earth, suid)        

    def subcell(self, other):
        r"""
        Subcell (subset) relation on cells.
        
        EXAMPLES::
        
            >>> a = Cell(Earth(), 'N1')
            >>> b = Cell(Earth(), 'N')
            >>> print a.subcell(b)
            True
            >>> print b.subcell(a)
            False
            
        """
        return str(self).startswith(str(other))
        
    def subcells(self, level=None):
        r"""
        Generator function for the set of level `level` subcells of `self`.  
        If `level=None`, then return a generator function for the children
        of `self`.
        
        EXAMPLES::
        
            >>> c = Cell(Earth(), 'N')
            >>> print [str(z) for z in c.subcells(2)]
            ['N00', 'N01', 'N02', 'N03', 'N04', 'N05', 'N06', 'N07', 'N08', 'N10', 'N11', 'N12', 'N13', 'N14', 'N15', 'N16', 'N17', 'N18', 'N20', 'N21', 'N22', 'N23', 'N24', 'N25', 'N26', 'N27', 'N28', 'N30', 'N31', 'N32', 'N33', 'N34', 'N35', 'N36', 'N37', 'N38', 'N40', 'N41', 'N42', 'N43', 'N44', 'N45', 'N46', 'N47', 'N48', 'N50', 'N51', 'N52', 'N53', 'N54', 'N55', 'N56', 'N57', 'N58', 'N60', 'N61', 'N62', 'N63', 'N64', 'N65', 'N66', 'N67', 'N68', 'N70', 'N71', 'N72', 'N73', 'N74', 'N75', 'N76', 'N77', 'N78', 'N80', 'N81', 'N82', 'N83', 'N84', 'N85', 'N86', 'N87', 'N88']
              
        """
        L = self.level
        if level is None:
            level = L + 1
        if level < L:
            return  # Stop iteration
        if level == L:
            yield self
            return 
        for t in product(range(9), repeat=level - L):
            yield Cell(self.earth, list(self.suid)+list(t))
                    
    def ul_vertex(self, surface='plane'):
        r"""
        Return the upper left (northwest) vertex of this planar or 
        ellipsoidal cell.
        
        WARNING: The upper left vertex of a cell might not lie in the cell, 
        because not all cells contain their boundary.
                  
        EXAMPLES::
        
            >>> c = Cell(Earth('unit_sphere'), 'N0')
            >>> print c.ul_vertex() == (-pi, 3*pi/4)
            True
            
        """
        # Call this cell c.
        # Find the location of the level 0 cell c0 containing c.
        x0, y0 = self.earth.ul_vertex[self.suid[0]]
        level = self.level
            
        # The column and row SUIDs of c give the the horizontal and vertical
        # distances, respectively, between the ul_vertex of c0 and 
        # the ul_vertex of c as fractions of the width of c0.
        suid_col, suid_row = self.suid_colrow()
        dx = sum(3**(level - i) * suid_col[i] for i in range(1, level + 1)) *\
             3**(-level)
        dy = sum(3**(level - i) * suid_row[i] for i in range(1, level + 1)) *\
             3**(-level)
        # Multiplied by 3**(-level) at the end instead of at each term to
        # avoid rounding errors.
        
        # Use dx, dy, and the width of c0 to calculate the ul_vertex of c.
        E = self.earth
        x = x0 + E.cell_width(0) * dx
        y = y0 - E.cell_width(0) * dy 

        # Project onto ellipsoid if necessary.
        if surface != 'plane':
            x, y = self.earth.project((x, y), source='rhealpix', 
                                      target='lonlat')        
        return x, y
        
    def nucleus_and_vertices(self, surface='plane', interpolation=0):
        r"""
        Return the nucleus of this planar or ellipsoidal cell along 
        with its corner points interpolated by `3**(interpolation)` points 
        along each edge.
        The nucleus is the first point in the output list and the remaining 
        points are the boundary points enumerated from the upleft corner
        to the downleft corner to the downright corner to the upright corner
        to the upleft corner (but not including the downright corner in
        the list again).
        The output for ellipsoidal cells is the projection onto the ellipsoid
        of the output for planar cells.  In particular, while the 
        nucleus of a planar cell is its centroid, the nucleus 
        of an ellipsoidal cell is not its centroid.
        To compute the centroid of a cell, use centroid() below.
                
        NOTE:
                
        The interpolation option is intended mostly for drawing the
        boundary of `self` with varying degrees of smoothness.
        
        EXAMPLES::
        
            >>> E = Earth('unit_sphere')
            >>> c = E.cell('N')
            >>> c.nucleus_and_vertices()
            [(-2.3561944901923448, 1.5707963267948966), (-3.141592653589793, 2.356194490192345), (-3.1415926535897931, 0.78539816339744828), (-1.5707963267948966, 0.78539816339744828), (-1.5707963267948966, 2.3561944901923448)]
        """
        # Compute the nucleus and boundary points in the plane.
        ul = self.ul_vertex()
        w = self.earth.cell_width(self.level)
        result = [tuple(ul + array((w/2, -w/2)))]   # nucleus
        point = ul
        result.append(point)
        # Number of points to interpolate along each edge of the square cell.    
        if interpolation < 0:
            interpolation = 0
        n = 3**interpolation
        for direction in [(0, -1), (1, 0), (0, 1), (-1, 0)]:
            for j in range(n):
                temp = array(point) + (float(j + 1)/n)*w*array(direction)
                result.append(tuple(temp))
            point = result[-1]
        # Remove the last point because it's the first point.
        result.pop()
        
        # Return the appropriate projection of result.
        if surface == 'plane':
            return result
        else:
            return [self.earth.project(p, 'rhealpix', 'lonlat') 
                    for p in result]
    
    def region(self):
        r"""
        Return the region, 'equatorial' or 'polar', of this cell.
        """
        if self.suid[0] in Earth.CELLS0[1:5]:
            return 'equatorial'
        else:
            return 'polar'
            
    def ellipsoidal_shape(self):
        r"""
        Return the shape of this cell ('rectangle', 'cap', 'dart', or 
        'trapezoid') when viewed on the ellipsoid.
        """
        import re   # For regular expressions.
        
        suid = str(self)
        if suid[0] in Earth.CELLS0[1:5]:
            return 'rectangle'
        if re.match(r"4*$", suid[1:]):
            return 'cap'
        if re.match(r"[048]+$", suid[1:]) or re.match(r"[246]+$", suid[1:]):
            return 'dart'
        return 'trapezoid'
        
    def centroid(self, surface='plane'):
        r"""
        Return the centroid of this planar or ellipsoidal cell.
        """
        if surface == 'plane':
            # Then this cell's centroid is its nucleus.
            return self.nucleus_and_vertices(surface='plane')[0]
        
        # This cell is ellipsoidal.
        # So we have to do some work.            
        nv = self.nucleus_and_vertices(surface='ellipsoid')
        nucleus = nv[0]
        vertices = nv[1:]
        shape = self.ellipsoidal_shape()
        if shape == 'cap':
            return nucleus
        if shape == 'rectangle':
            lam_bar = nucleus[0]
            phi_bar = sum([v[1] for v in vertices])/4    
            return lam_bar, phi_bar   
        planar_vertices = self.nucleus_and_vertices()[1:]
        x1 = min([v[0] for v in planar_vertices])
        x2 = max([v[0] for v in planar_vertices])
        y1 = min([v[1] for v in planar_vertices])
        y2 = max([v[1] for v in planar_vertices])
        area = (x2 - x1)**2
        lam = lambda x, y: self.earth.project((x, y), 'rhealpix', 'lonlat')[0]
        phi = lambda x, y: self.earth.project((x, y), 'rhealpix', 'lonlat')[1]
        if shape == 'dart':
            lam_bar = nucleus[0]
            phi_bar = (1/area)*\
                      integrate.dblquad(phi, y1, y2, \
                                        lambda x: x1, lambda x: x2)[0]
            return lam_bar, phi_bar 
        # Now shape == 'trapezoid'.
        # Option 1: Compute lam_bar by numerical integration:
        # lam_bar = (1/area)*\
        #       integrate.dblquad(lam, y1, y2, lambda x: x1, lambda x: x2)[0]
        # Option 2: Evaluate the integral symbolically 
        # and then plug in values, which might speed up the calculation:
        w = x2 - x1
        R_A = self.earth.R
        hx0, hy0 = self.earth.project(nucleus, 'lonlat', 'healpix') 
        # x and y extremes of the HEALPix projection of this cell's interior:
        hx1 = hx0 - w/2
        hx2 = hx0 + w/2
        # Without loss of generality, force HEALPix y coordinates into 
        # the northern hemisphere:
        hy1 = abs(hy0) - w/2 
        hy2 = abs(hy0) + w/2   
        # Compute xc.
        cap_number = floor(2*hx0/(pi*R_A) + 2)
        if cap_number >= 4: 
            # Rounding error. 
            cap_number = 3
        xc = -3*pi/4 + (pi/2)*cap_number    
        integral = lambda x, y: (pi/8)*x*(2*R_A*xc - x)*\
                   log(1 - 2*y/(pi*R_A)) + xc*x*y
        lam_bar = (1/area)*\
                  (integral(hx2, hy2) - integral(hx1, hy2) -\
                   integral(hx2, hy1) + integral(hx1, hy1))
        phi_bar = (1/area)*\
                  integrate.dblquad(phi, y1, y2, lambda x: x1, lambda x: x2)[0]
        return lam_bar, phi_bar
        
    @staticmethod
    def atomic_rotate(quarter_turns, x):
        r"""
        Return the function g that represents the table ::

            0 1 2
            3 4 5
            6 7 8

        rotated anticlockwise by `quarter_turns` quarter turns.
        The table is then read from left to right and top to bottom
        to give the values g(0), g(1), ..., g(8).

        INPUT:

        - `quarter_turns` - 0, 1, 2, or 3.
        - `x` - 0, 1, 2, ..., 8, Earth.CELLS0[5], Earth.CELLS0[1], ..., Earth.CELLS0[0].
        
        EXAMPLES::

            >>> Cell.atomic_rotate(0, 0)
            0
            >>> Cell.atomic_rotate(1, 0)
            2
            >>> Cell.atomic_rotate(2, 0)
            8
            >>> Cell.atomic_rotate(3, 0)
            6
            >>> Cell.atomic_rotate(4, 0)
            0
            
        """
        # Subcell numbers move like this:
        f = {
        0: 2,
        1: 5,
        2: 8,
        3: 1,
        4: 4,
        5: 7,
        6: 0,
        7: 3,
        8: 6
        }

        # Level 0 cells stay fixed.
        for n in Earth.CELLS0:
            f[n] = n

        quarter_turns = quarter_turns % 4
        if quarter_turns == 1:
            return f[x]
        elif quarter_turns == 2:
            return f[f[x]]
        elif quarter_turns == 3:
            return f[f[f[x]]]
        else:
            return x
        
    def rotate(self, quarter_turns):
        r"""
        Return the cell that is the result of rotating this cell 
        `quarter_turns` quarter_turns.
        Used in neighbor().
    
        EXAMPLES::

            >>> c = Cell(Earth(), 'N0')
            >>> print c.rotate(0)
            N0
            >>> print c.rotate(1)
            N2
            >>> print c.rotate(2)
            N8
            >>> print c.rotate(3)
            N6
            >>> print c.rotate(4)
            N0
            
        """    
        quarter_turns = quarter_turns % 4
        suid = [Cell.atomic_rotate(quarter_turns, n) for n in self.suid]  
        return Cell(self.earth, suid)
        
    def neighbor(self, direction):
        r"""
        Return the neighboring cell of this cell in the 
        direction `direction`.
        The direction is relative to this planar neighbor diagram,
        where `self` is the middle cell ::

                            up
                          *-----*
                          |     |
                          |     |
                          |     |
                    *-----*-----*-----*
                    |     | 012 |     |
               left |     | 345 |     | right 
                    |     | 678 |     |
                    *-----*-----*-----*
                          |     |
                          |     |
                          |     |
                          *-----*
                           down
        
        The tricky part is that the neighbor relationships of the 
        six level 0 cells is determined by the positions of those cells
        on the surface of a cube, one on each face (and not on a plane).
        So sometimes rotating cells in needed to compute neighbors.
        
        INPUT:
        
        - `direction` - One of the strings 'up', 'right', 'down', or 'left'.

        EXAMPLES::

            >>> c = Cell(Earth(), 'N0')
            >>> print c.neighbor('down')
            N3
            
        """        
        dirs = ['left','right','down','up']
        assert direction in dirs, \
        "direction must be one of " + str(dirs)
        an = self.earth.atomic_neighbors

        # First, compute the neighbor of self  naively, that is, 
        # without considering rotations.
        self_suid = self.suid
        neighbor_suid = []
        border = {'left': [0, 3, 6], 'right': [2, 5, 8], \
        'up': [0, 1, 2], 'down': [6, 7, 8]}
        crossed_all_borders = False
        # Scan from the back to the front of Xsuid.
        for i in reversed(range(len(self_suid))):
            n = self_suid[i]
            if crossed_all_borders:
                neighbor_suid.append(n)
            else:
                neighbor_suid.append(an[n][direction])
                if n not in border[direction]:
                    crossed_all_borders = True
        neighbor_suid.reverse()
        neighbor = Cell(self.earth, neighbor_suid)
        
        # Second, rotate the neighbor if necessary.
        # If self is a polar cell and neighbor is not, or vice versa, 
        # then rotate neighbor accordingly.  
        self0 = self_suid[0]
        CELLS0 = Earth.CELLS0
        neighbor0 = neighbor_suid[0]
        if (self0 == CELLS0[5] and neighbor0 == an[self0]['left']) \
        or (self0 == an[CELLS0[5]]['right'] and neighbor0 == CELLS0[5]) \
        or (self0 == CELLS0[0] and neighbor0 == an[self0]['right']) \
        or (self0 == an[CELLS0[0]]['left'] and neighbor0 == CELLS0[0]):
            neighbor = neighbor.rotate(1)
        elif (self0 == CELLS0[5] and neighbor0 == an[self0]['down']) \
        or (self0 == an[CELLS0[5]]['down'] and neighbor0 == CELLS0[5]) \
        or (self0 == CELLS0[0] and neighbor0 == an[self0]['up']) \
        or (self0 == an[CELLS0[0]]['up'] and neighbor0 == CELLS0[0]):
            neighbor = neighbor.rotate(2)
        elif (self0 == CELLS0[5] and neighbor0 == an[self0]['right']) \
        or (self0 == an[CELLS0[5]]['left'] and neighbor0 == CELLS0[5]) \
        or (self0 == CELLS0[0] and neighbor0 == an[self0]['left']) \
        or (self0 == an[CELLS0[0]]['right'] and neighbor0 == CELLS0[0]):
            neighbor = neighbor.rotate(3)
        return neighbor
        
    def neighbors(self):
        r"""
        Return a dictionary of the left, right, down, and up neighbors of this 
        cell. 
        Direction is relative to this planar neighbor diagram, where `self` is 
        the middle cell ::

                             up
                           *-----*
                           |     |
                           |     |
                           |     |
                     *-----*-----*-----*
                     |     | 012 |     |
                left |     | 345 |     | right 
                     |     | 678 |     |
                     *-----*-----*-----*
                           |     |
                           |     |
                           |     |
                           *-----*
                            down

         The tricky part is that the neighbor relationships of the 
         six level 0 cells is determined by the positions of those cells
         on the surface of a cube, one on each face (and not a plane).
         So sometimes rotating cells in needed to compute neighbors.
        
        EXAMPLES::

            >>> c = Cell(Earth(), 'N0')
            >>> for k, v in c.neighbors().items():
            ...     print k, v
            ... 
            down N3
            right N1
            up Q2
            left R0
            
        """ 
        result = {}       
        for d in ['left', 'right', 'down','up']:
            result[d] = self.neighbor(d)
        return result


import collections
from bisect import bisect_left, insort_left

class CellFamily(collections.MutableSequence):
    r"""
    A subcell-free ordered list of cells along with the minimum level of 
    the cells in family the maximum level of the cells in family.
    """
     
    def __init__(self, cells=None, min_level=None, max_level=None,
                 sort=True, eliminate_subcells=True):
        r"""
        Returns a cell family given a list of cells.
        Do not set `sort=False` or `eliminate_subcells=False` unless you know
        that `cells` is sorted or subcell-free, respectively.
        
        If `sort=True`, then uses O(n * log(n)) cell operations, where n = 
        len(cells).
        If `sort=False` and (`eliminate_subcells=True` or `min_level=None` or 
        `max_level=None`), then uses O(n) cell operations.
        Otherwise uses O(1) cell operations.
              
        EXAMPLES::
            
            >>> E = Earth()
            >>> cells = [E.cell('N8'), E.cell('Q77')] + list(E.grid(0))
            >>> f = CellFamily(cells); print f
            ['N', 'O', 'P', 'Q', 'R', 'S']
            
        """
        if not cells:
            # The cells of the family:
            self.cells = []
            # The Earth model of the family:
            self.earth = None
            # The level of the biggest cells in the family:        
            self.min_level = None
            # The level of the smallest cells in the family:
            self.max_level = None
        else:
            self.cells = cells
            self.earth = cells[0].earth  
            self.min_level = min_level  
            self.max_level = max_level  
            if sort:
                # O(n*log(n)), where n = len(cells). 
                self.cells = sorted(cells)
            if eliminate_subcells:
                # O(n).
                # Assume self.cells is sorted.
                # Elimated subcells by scanning in reverse order.
                temp = [self.cells.pop()]    # Last cell.
                for c in reversed(self.cells):
                    if c.subcell(temp[0]):
                        next
                    else:
                        temp.insert(0, c)
                self.cells = temp                
            if (min_level is None) or (max_level is None):
                # O(n)
                # Compute the levels.
                self.min_level = self.cells[0].level
                self.max_level = self.min_level
                for c in self.cells:
                    self.min_level = min(self.min_level, c.level)
                    self.max_level = max(self.max_level, c.level)
                
    # First define the methods necessary to inherit from
    # collections.MutableSequence.      
    def __len__(self):
        r"""
        Return the length of self. 
        """
        return len(self.cells)
        
    def __iter__(self):
        r"""
        Return an iterator object for self.
        """
        return iter(self.cells)
        
    def __contains__(self, key):
        r"""
        Return True if `key` is in `self` and False otherwise. 
        Uses bisection to search for `key` in O(log(n)) comparisons,
        where n is the number of cells in `self`.
        """
        cells = self.cells
        i = bisect_left(self.cells, key, 0)
        return i < len(cells) and cells[i] == key
        
    def __getitem__(self, index):
        r"""
        Return `self[index]`.
        """
        return self.cells[index]
        
    def __setitem__(self, index, value):
        r"""
        Set `self[index] = value`.
        """
        self.cells[index] = value
        return
        
    def __delitem__(self, index):
        r"""
        Remove `self[index]`.    
        """
        del self.cells[index]
        
    def insert(self, index, value):
        r"""
        Insert the `value` into `self` at index `index`.
        Warning: for internal use only; does not maintain sorted order.
        Use `add()` to add cells to `self`.
        """
        self.cells.insert(index, value)
        
    # Now define the rest of the methods i need.
    def __eq__(self, other):
        r"""
        Return True if `self` and `other` are equal and False otherwise.
        """
        return (self.earth == other.earth) and (self.cells == other.cells)

    def __ne__(self, other):
        r"""
        Return True if `self` and `other` are unequal and False otherwise.
        """
        return not self == other

    def __str__(self):
        r"""
        Return a string representation of `self`.
        """
        return str([str(cell) for cell in self.cells])
        
    def add(self, key):
        r"""
        Add `key` to `self` if `key` is not a subcell of any cell of `self`.
        Uses bisection to find the insertion point in O(log(len(self)))
        comparisons.
        
        EXAMPLES::
        
            >>> E = Earth()
            >>> f = CellFamily(list(E.grid(0)))
            >>> del f[2]; print f
            ['N', 'O', 'Q', 'R', 'S']
            >>> c = E.cell('P3')
            >>> f.add(c)
            >>> print f
            ['N', 'O', 'P3', 'Q', 'R', 'S']
            >>> c = E.cell('O3')
            >>> f.add(c)    # Should do nothing.
            >>> print f
            ['N', 'O', 'P3', 'Q', 'R', 'S']
            
        """
        cells = self.cells
        # Find insertion point for key.
        i = bisect_left(cells, key, 0)
        if i < len(cells) and not key.subcell(cells[i]):
            self.insert(i, key)

    def minimize(self):
        r"""
        Return the minimal cell family that has the same region 
        as `self`.

        Uses O(len(self)) cell operations.

        EXAMPLES::

            >>> E = Earth()
            >>> c = Cell(E, 'N')
            >>> f = CellFamily(list(c.subcells(1)))
            >>> print f
            ['N0', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8']
            >>> print f.minimize()
            ['N']

        """
        import re   # Use regular expressions for this function.

        # Deal with empty family case.
        if not self.cells:
            return self

        # Initialize variables.
        # Processed cells.
        black = []
        # Encountered but unprocesed cells:
        gray = []
        # Unencountered cells (or rather SUIDs) that might get built up from 
        # smaller gray cells:
        white = []

        # Loop through all cells in self (in ascending order) 
        # and group cells into bigger cells when possible.
        for c in self:             
            gray.append(c)
            suid = str(c)

            # Is c the first corner piece of bigger cells?
            # That is, does c's suid end in zeros?
            s = re.split(r'(0+)$', suid)[:2]
            if s[0] != suid:
                # Then c can start filling up bigger cells.
                # Color these bigger cells white.
                s, zeros = s[0], s[1]
                white.extend([s + '0'*i for i in range(len(zeros))])
            else:
                # c is not the first corner piece of bigger cells.
                # Is c continuing an interval of gray cells? 
                # That is, is c's predecessor gray?
                if c.predecessor() in gray:
                    # Is c the last corner piece of a white cell?
                    # That is, does c's suid end in eights?
                    t = re.split(r'(8+)$', suid)[:2]
                    if t[0] != suid:
                        # Then color gray the biggest white cell that c is 
                        # the last corner piece of.
                        t, eights = t[0], t[1]
                        for i in range(len(eights)):
                            prefix = t + '8'*i
                            if prefix in white:
                                big = Cell(self.earth, prefix)
                                gray.append(big)
                                break
                else:
                    # c is not continuing an interval of gray cells.
                    # So no more grouping of gray cells can be done.
                    # So color all the gray cells black and 
                    # reset gray and white cell lists.
                    black.extend(gray)
                    gray = []
                    white = []                                
        # Color final lot of gray cells black.
        black.extend(gray)

        # Create cell family from black cell list.        
        return CellFamily(black, sort=False)    
        
    def union(self, other, filter_level=None):
        r"""
        Return the cell family u that is the  minimizeed union of `self` 
        and `other`, that is, the smallest cell family whose region is the 
        union of the regions of `self` and `other`.
        If a nonnegative integer `filter_level` is given, then return the cell
        family comprised of the cell's of u at levels <= `filter_level`.
        Do not set `minimize_first=False` unless you know that the input
        cell families are already minimizeed.
        
        Uses O(max(m, n)) cell operations, where m = `len(self)` and
        n = `len(other)`.
        
        EXAMPLES::
        
            >>> a = Cell(Earth(), 'N')
            >>> b1 = Cell(Earth(), 'P1')
            >>> f = CellFamily([b1] + list(a.subcells(1))[:4])
            >>> print f
            ['N0', 'N1', 'N2', 'N3', 'P1']
            >>> g = CellFamily(list(a.subcells(1))[3:])
            >>> print g
            ['N3', 'N4', 'N5', 'N6', 'N7', 'N8']
            >>> print f.union(g)
            ['N', 'P1']
            >>> print f.union(g, filter_level=0)
            ['N']

        """
        f = self
        g = other
        
        # Deal with empty family case.
        if not f.cells:
            return g
        if not g.cells:
            return f
                    
        # Union the regions of f and g.
        # We can't ignore cells of level > filter_level,
        # because they may aggregate to bigger (smaller level) cells.    
        u = []  # The union  
        i = 0   # Index of current cell in f under investigation.
        j = 0   # Index of current cell in g under investigation.
        while i < len(f) and j < len(g):
            a = f[i]
            b = g[j]
            if a.subcell(b):
                # Discard/skip a.
                i += 1
            elif b.subcell(a):
                # Discard/skip b.
                j += 1
            elif a < b:
                # Add a to u.
                u.append(a)
                i += 1
            else:
                # b < a
                # Add b to u    
                u.append(b)
                j += 1
                    
        # At this point, we've processed all of f's cells or 
        # all of g's cells.
        # Add the remaining unprocessed cells (possibly none) to u.
        if i >= len(f):
            u.extend(g[j:])
        else:
            # j >= len(g)
            u.extend(f[i:])
                
        # Minimize u, which is already sorted and subcell-free.
        u = CellFamily(u, sort=False, eliminate_subcells=False).minimize()
 
        # Filter by level if desired.
        if isinstance(filter_level, int) and filter_level >= 0 \
        and filter_level < u.max_level:
            u = CellFamily([x for x in u if x.level <= filter_level], 
                           sort=False, eliminate_subcells=False)
        return u
    
    def union_all(self, others, filter_level=None, minimize_first=True):
        r"""
        Return the cell family u that is the union of the regions of this cell 
        family and cell families in the list `others`. 
        If a nonnegative integer `filter_level` is given, then return the cell
        family comprised of the cell's of u at levels <= `filter_level`.

        Uses O(n) cell operations, where n is the maximum of `len(self)` and
        the lengths of the cell families in `others`.
        
        EXAMPLES::
        
            >>> E = Earth()
            >>> cell_families = []
            >>> for c in E.grid(1):
            ...     cell_families.append(CellFamily([c]))
            >>> print len(cell_families)
            54
            >>> f = cell_families.pop(0)
            >>> print f
            ['N0']
            >>> print f.union_all(cell_families)
            ['N', 'O', 'P', 'Q', 'R', 'S']
        """
        families = [self] + others
            
        # Build the union iteratively.
        # Do not filter yet, because small cells may aggregate to big cells
        # in the course of unioning.
        u = families.pop(0)
        for f in families:
            u = u.union(f)
            
        # Filter by level if desired.
        if isinstance(filter_level, int) and filter_level >= 0 \
        and filter_level < u.max_level:
            u = CellFamily([x for x in u if x.level <= filter_level],
                           sort=False, eliminate_subcells=False)
        return u
        
    def intersect(self, other, filter_level=None):
        r"""
        Return the minimal cell family v for the intersection of the regions 
        of `self` and `other`. 
        If a nonnegative integer `filter_level` is given, then return the cell
        family comprised of the cell's of v at levels <= `filter_level`. 
 
        Uses O(max(m, n)) cell operations, where m = `len(self)` and
        n = `len(other)`.
        
        EXAMPLES::
        
            >>> a = Cell(Earth(), 'N')
            >>> b = Cell(Earth(), 'P')
            >>> f = CellFamily([b] + list(a.subcells(1))[:4])
            >>> print f
            ['N0', 'N1', 'N2', 'N3', 'P']
            >>> g = CellFamily([b] + list(a.subcells(1))[3:])
            >>> print g
            ['N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'P']
            >>> print f.intersect(g)
            ['N3', 'P']
            >>> print f.intersect(g, filter_level=0)
            ['P']
            
        """
        f = self
        g = other
        
        # If f or g is empty, return the empty cell family.
        if not (f and g):
            return CellFamily()       
            
        # Set default filter level if none given.
        if not (isinstance(filter_level, int) and filter_level >= 0):
            filter_level = f.earth.max_level
             
        # Intersect the regions of f and g,
        # ignoring cells at level > filter_level.
        v = []  # The intersection.
        i = 0   # Index of current cell in f under investigation.
        j = 0   # Index of current cell in g under investigation.
        while i < len(f) and j < len(g):
            a = f[i]
            b = g[j]
            if a.level > filter_level:
                # Discard/skip a.
                i += 1
                continue
            if b.level > filter_level:
                # Discard/skip b.
                j += 1
                continue
            if a.subcell(b):
                # Add a to v.
                v.append(a)
                i += 1
            elif b.subcell(a):
                # Add b to v.
                v.append(b)
                j += 1
            elif a < b:
                # Discard/skip a.
                i += 1
            else:
                # b < a
                # Discard/skip  b.
                j += 1
                        
        # Minimize v, which is already sorted and subcell-free. 
        return CellFamily(v, sort=False, eliminate_subcells=False).minimize()
        
    def intersect_all(self, others, filter_level=None):
        r"""
        Return the minimal cell family v for the intersection of the regions
        of `self` and the cell families in the list `others`.
        If a nonnegative integer `filter_level` is given, then return the cell
        family comprised of the cell's of v at levels <= `filter_level`.

        Uses O(n) cell operations, where n is the maximum of `len(self)` and
        the lengths of the cell families in `others`.

        EXAMPLES::
        
            >>> E = Earth()
            >>> cell_families = []
            >>> for c in E.grid(0):
            ...     cell_families.append(CellFamily([c]))
            ... 
            >>> f = cell_families.pop(0)
            >>> print f.intersect_all(cell_families)
            []
            
        """
        families = [self] + others

        # Build the intersection iteratively.
        # Filter along the way.
        v = families.pop(0)
        for f in families:
            v = v.intersect(f, filter_level=filter_level)
            if not v:
                # v is empty and so we can stop now.
                break
        return v