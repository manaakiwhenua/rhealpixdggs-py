r"""
The Python 2.7 module implements the cartographic projections used in SCENZ-Grid.

CHANGELOG:

- Alexander Raichev (AR), 2011-09-27: Split projection functions and pixelation 
  functions into separate files. 
- AR, 2011-10-04: Tidied a little. Improved variable names.
- AR, 2011-10-06: Added more examples and more default arguments.
- AR, 2011-10-23: Added more error-checking, such as the function in_image().
- AR, 2012-03-07: Changed 'rot' to 'rotate'.
- AR, 2012-03-23: Fixed a conceptual bug in combine_caps() and get_cap(). Now polar caps actually assemble on top of the caps numbered north and south.
- AR, 2012-03-26: Changed in_image() to accept a radius argument.
- AR, 2012-04-16: Imported NumPy instead of SciPy.
- AR, 2012-06-13: Simplified get_cap() to cap() and simplified combine_caps().
- AR, 2012-07-10: Corrected a subtle longitude rounding error in healpix_sphere(inverse=True).
- AR, 2012-08-31: Redefined the projections for only the unit sphere and 
  the ellipsoid with eccentricity e and authalic sphere equal to the unit
  sphere. Will handle scaling of spheres/ellipsoids in grids.py  
  instead. This simplifies the code and puts it more in line with proj4
  design, which also deals with scaling as a pre-/post-processing step outside
  of the projection definitions.
- AR, 2012-09-03: Corrected WGS84_f, the WGS84 ellipsoid flattening factor.

NOTE:

All lengths are measured in meters and all angles are measured in radians 
unless indicated otherwise. 

These projections work only for ellipsoids and not for general (triaxial) ellipsoids.

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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# scipy is a third-party module
from numpy import pi, floor, ceil, log, sqrt, sin, cos, tan, \
arcsin, arccos, arctan, deg2rad, rad2deg, sign, array, identity, dot 


# Matrix for anticlockwise rotation by pi/2: 
rotate1 = array([[0, -1], [1, 0]])
 # Matrix for anticlockwise rotation by pi:  
rotate2 = dot(rotate1, rotate1) 
# Matrix for anticlockwise rotation by 3*pi/2.
rotate3 = dot(rotate2, rotate1)  
# Dictionary of all powers of rotate1 and its inverse (rotate3).
rotate = {0: identity(2, int), 1: rotate1, 2: rotate2, 3: rotate3, -1: rotate3,
          -2: rotate2, -3: rotate1}

def healpix_sphere(u, v, lon0=0, degrees=False, inverse=False):
    r"""
    Compute the forward and inverse signature functions of the HEALPix 
    projection of the unit sphere with central meridian `lon0`. 
        
    INPUT:
    
    - `u, v` - If `inverse=False`, then these are geographic coordinates.
      Any input angles are accepted and get standardized to lie in the 
      intervals -pi <= `u` < pi and -pi/2 <= `v` <= pi/2 
      via wrap_longitude() and wrap_latitude().
      If `inverse=True`, then these are planar coordinates in the image
      of the HEALPix projection of the unit sphere.
    - `lon0` - (Optional; default = 0) Central meridian. Any angle is accepted
      and gets standardized to lie in the interval -pi <= `lon0` < pi.
    - `degrees` - (Optional; default = False) Boolean. If True, then use 
      degrees for input/output instead of radians.
    - `inverse` - (Optional; default = False) Boolean. If False, then compute
      forward function. If True, then compute inverse function.
    
    EXAMPLES::
    
        >>> healpix_sphere(0, 0, lon0=-45, degrees=True); pi/4
        (0.78539816339744828, 0.0)
        0.7853981633974483
        >>> healpix_sphere(-45, 0, lon0=-45, degrees=True)
        (0.0, 0.0)
        >>> healpix_sphere(180, 0, lon0=-45, degrees=True)
        (-2.3561944901923448, 0.0)
        >>> healpix_sphere(-180, 0, lon0=-45, degrees=True)
        (-2.3561944901923448, 0.0)
        >>> u, v = -50, 82
        >>> x, y = healpix_sphere(u, v, degrees=True)
        >>> healpix_sphere(x, y, degrees=True, inverse=True); u, v
        (-49.999999999999993, 82.000000000000028)
        (-50, 82)
        
    """
    if not inverse:
        # Compute forward function.
        lam, phi = u, v
        if degrees:
            # Convert to radians.
            lam, phi, lon0 = deg2rad((lam, phi, lon0))
        lam = wrap_longitude(lam - lon0)   # Shift central meridian to `lon0`.
        phi = wrap_latitude(phi)  
        phi0 = arcsin(2.0/3)
        
        # Equatorial region.
        if abs(phi) <= phi0:
            x = lam 
            y = 3*pi/8*sin(phi)
        
        # Polar region.
        else:
            sigma = sqrt(3*(1 - abs(sin(phi))))
            cap_number = floor(2*lam/pi + 2)
            if cap_number >= 4: 
                # Rounding error
                cap_number = 3
            lamc = -3*pi/4 + (pi/2)*cap_number  
            x = lamc + (lam - lamc)*sigma
            y = sign(phi)*pi/4*(2 - sigma)
        return (x, y)
    else:
        # Compute inverse function.
        if degrees:
            # Convert to radians.
            lon0 = deg2rad(lon0)
        x, y = u, v
        # Throw error if input coordinates are out of bounds.
        if not in_image(x, y):
            print "Error: input coordinates (%f, %f) are out of bounds" % \
                  (x, y)
            return
        y0 = pi/4
        # Equatorial region.
        if abs(y) <= y0:
            lam = x
            phi = arcsin(8*y/(3*pi))
        # Polar region but not the poles.
        elif abs(y) < pi/2:
            cap_number = floor(2*x/pi + 2)
            if cap_number >= 4: 
                # Rounding error.
                cap_number = 3
            xc = -3*pi/4 + (pi/2)*cap_number    
            tau = 2 - 4*abs(y)/pi
            lam = xc + (x - xc)/tau
            phi = sign(y)*arcsin(1 - tau**2/3)
            # Handle rounding errors in longitude.
            if lam < -pi:
                lam = -pi
            elif lam > pi:
                lam = pi
        # Poles.
        else:
            lam = -pi - lon0
            phi = sign(y)*pi/2
        # Shift central meridian to 0.
        lam = wrap_longitude(lam + lon0)   
        if degrees:
            # Convert to degrees.
            lam, phi = rad2deg((lam, phi))
        return lam, phi
        
def healpix_ellipsoid(u, v, e, lon0=0, degrees=False, inverse=False):
    r"""
    Compute the forward and inverse signature functions of the HEALPix 
    projection of an oblate ellipsoid with eccentricity `e` whose authalic
    sphere is the unit sphere with central meridian `lon0`. 
        
    INPUT:
    
    - `u, v` - If `inverse=False`, then these are geographic coordinates.
      Any input angles are accepted and get standardized to lie in the intervals 
      -pi <= `u` < pi and -pi/2 <= `v` <= pi/2 via wrap_longitude() and 
      wrap_latitude().
      If `inverse=True`, then these are planar coordinates in the image
      of the HEALPix projection of a sphere of radius 1.
    - `e` - Eccentricity.
    - `lon0` - (Optional; default = 0) Central meridian. Any angle is accepted
      and gets standardized to lie in the interval -pi <= `lon0` < pi.
    - `degrees` - (Optional; default = False) Boolean. If True, then use 
      degrees for input/output instead of radians.
    - `inverse` - (Optional; default = False) Boolean. If False, then compute
      forward function. If True, then compute inverse function.
    
    EXAMPLES::
    
        >>> u, v = 82,-45
        >>> x, y = healpix_ellipsoid(u, v, e=WGS84_E, degrees=True)
        >>> healpix_ellipsoid(x, y,e=WGS84_E, degrees=True, inverse=True); u, v
        (81.999999999999986, -44.999999995319058)
        (82, -45)
    """
    if not inverse:
        # Compute forward function.
        lam, phi = u, v
        if degrees:
            # Convert to radians.
            lam, phi, lon0 = deg2rad((lam, phi, lon0))
        bet = auth_lat(phi, e)
        return healpix_sphere(lam, bet, lon0=lon0, degrees=False)  
    else:
        # Compute inverse function.
        # Throw error if input coordinates are out of bounds.
        if not in_image(u, v):
            print "Error: input coordinates (%f, %f) are out of bounds" % \
                  (u, v)
            return
        if degrees:
            # Convert to radians.
            lon0 = deg2rad(lon0)
        lam, bet = healpix_sphere(u, v, lon0=lon0, degrees=False, inverse=True)
        phi = auth_lat(bet, e, inverse=True)    
        if degrees:
            # Convert to degrees
            lam, phi = rad2deg((lam, phi))
        return lam, phi
        
def rhealpix_sphere(u, v, lon0=0, north=0, south=0, degrees=False,
                    inverse=False):
    r"""
    Compute the forward and inverse signature functions of the rHEALPix 
    projection of the unit sphere with central meridian `lon0`. 
    The north pole square is put in position `north`, and the south pole
    square is put in position `south`.
        
    INPUT:
    
    - `u, v` - If `inverse` = False, then these are geographic coordinates.
      Any input angles are accepted and get standardized to lie in the intervals 
      -pi <= `u` < pi and -pi/2 <= `v` <= pi/2 via wrap_longitude() and 
      wrap_latitude().
      If `inverse` = True, then these are planar coordinates in the image
      of the rHEALPix projection of the unit sphere.
    - `lon0` - (Optional; default = 0) Central meridian. Any angle is accepted
      and gets standardized to lie in the interval -pi <= `lon0` < pi.
    - `north, south` - (Optional; defaults = 0, 0) Integers between 0 and 3 
      indicating positions of north polar and south polar squares, respectively.
    - `degrees` - (Optional; default = False) Boolean. If True, then use 
      degrees for input/output instead of radians.
    - `inverse` - (Optional; default = False) Boolean. If False, then compute
      forward function. If True, then compute inverse function.
    
    EXAMPLES::
    
        >>> u, v = -50, 82
        >>> x, y = rhealpix_sphere(u, v, north=0, south=0, degrees=True)
        >>> rhealpix_sphere(x, y, degrees=True, inverse=True); u, v
        (-50.000000000000036, 82.000000000000028)
        (-50, 82)
                
    NOTE:
    
    The polar squares are labeled 0, 1, 2, 3 from east to west like this::
        
        east    lon0    west        
        *---*---*---*---*
        | 0 | 1 | 2 | 3 |
        *---*---*---*---*
        |   |   |   |   |
        *---*---*---*---*
        | 0 | 1 | 2 | 3 |
        *---*---*---*---*    
    """
    # Ensure north and south lie in {0, 1, 2, 3}.
    north = north % 4
    south = south % 4
    if not inverse:
        # Compute forward function.
        x, y = healpix_sphere(u, v, lon0=lon0, degrees=degrees)
        return combine_caps(x, y, north=north, south=south)
    else:
        # Compute inverse function.
        # Throw error if input coordinates are out of bounds.
        if not in_image(u, v, proj='rhealpix_sphere', south=south,
                        north=north):
            print "Error: input coordinates (%f,%f) are out of bounds" % (u, v)
            return
        x, y = combine_caps(u, v, north=north, south=south, inverse=True)
        return healpix_sphere(x, y, lon0=lon0, degrees=degrees, inverse=True)

def rhealpix_ellipsoid(u, v, e, lon0=0, north=0, south=0, degrees=False,
                      inverse=False):
    r"""
    Compute the forward and inverse signature functions of the rHEALPix 
    projection of an oblate ellipsoid with eccentricity `e` whose 
    authalic sphere is the unit sphere with central meridian `lon0`. 
    The north pole square is put in position `north`, and the south pole
    square is put in position `south`.
        
    INPUT:
    
    - `u, v` - If `inverse=False`, then these are geographic coordinates.
      Any input angles are accepted and get standardized to lie in the 
      intervals 
      -pi <= `u` < pi and -pi/2 <= `v` <= pi/2 via wrap_longitude() and 
      wrap_latitude().
      If `inverse= True`, then these are planar coordinates in the
      image of the rHEALPix projection of the unit sphere.
    - `e` - Eccentricity.
    - `lon0` - (Optional; default = 0) Central meridian. Any angle is accepted
      and gets standardized to lie in the interval -pi <= `lon0` < pi.
    - `north, south` - (Optional; defaults = 0, 0) Integers between 0 and 3 
      indicating positions of north polar and south polar squares, 
      respectively.
      See rhealpix_sphere() docstring for a diagram.
    - `degrees` - (Optional; default = False) Boolean. If True, then use 
      degrees for input/output instead of radians.
    - `inverse` - (Optional; default = False) Boolean. If False, then compute
      forward function. If True, then compute inverse function.

    EXAMPLES::
    
        >>> u, v = 90, 72
        >>> x, y = healpix_ellipsoid(u, v, e=WGS84_E, degrees=True)
        >>> healpix_ellipsoid(x, y,e=WGS84_E, degrees=True, inverse=True); u, v
        (90.000000000000028, 71.99999999914445)
        (90, 72)
    """
    # Ensure north and south lie in {0, 1,2, 3}.
    north = north % 4
    south = south % 4
    if not inverse:
        # Compute forward function.
        x, y = healpix_ellipsoid(u, v, e=e, lon0=lon0, degrees=degrees) 
        return combine_caps(x, y, north=north, south=south)  
    else:
        # Compute inverse function.
        # Throw error if input coordinates are out of bounds.
        if not in_image(u, v, proj='rhealpix_sphere', south=south, 
                        north=north):
            print "Error: input coordinates (%f,%f) are out of bounds" % (u, v)
            return
        x, y = combine_caps(u, v, north=north, south=south, inverse=True)
        return healpix_ellipsoid(x, y, e=e, lon0=lon0, degrees=degrees, 
                                 inverse=True)   
        
def combine_caps(x, y, north=0, south=0, inverse=False):
    r"""
    Rearrange point `(x, y)` in the HEALPix projection by 
    combining the polar caps into two polar squares.
    Put the north polar square in position `north` and 
    the south polar square in position `south`.
    If `inverse=True`, uncombine the polar caps.

    INPUT:

    - `x, y` - Coordinates in the HEALPix projection of the unit sphere.
    - `north, south` - Integers between 0 and 3 indicating the positions
      of the north polar square and south polar square respectively.
      See rhealpix_sphere() docstring for a diagram.
    - `inverse` - (Optional; default = False) Boolean. If False, then compute
      forward function. If True, then compute inverse function.
        
    EXAMPLES::
    
        >>> u, v = -pi/4, pi/3
        >>> x, y = combine_caps(u, v); x, y
        (-1.8325957145940459, 1.5707963267948966)
        >>> combine_caps(x, y, inverse=True); u, v
        (-0.78539816339744828, 1.0471975511965976)
        (-0.7853981633974483, 1.0471975511965976)
        
    """
    c, region = cap(x, y, north=north, south=south, inverse=inverse)
    if region == 'equatorial':
        # (x,y) remains fixed
        return x, y
        
    xy = array((x, y))
    tc = array((-3*pi/4 + c*pi/2, sign(y)*pi/2))
    if not inverse:
        # Forward function.
        # Rotate (x, y) about tc and then translate it to 
        # the tip u of the polar cap it will be assembled upon.
        if region == 'north':
            u = array((-3*pi/4 + north*pi/2, pi/2))
            x, y = dot(rotate[c - north], xy - tc) + u     
        elif region == 'south':
            u = array((-3*pi/4 + south*pi/2, -pi/2))
            x, y = dot(rotate[-(c - south)], xy - tc) + u
    else:
        # Inverse function.
        # Unrotate (x, y) about u and then translate it to tc.
        if region == 'north':
            u = array((-3*pi/4 + north*pi/2, pi/2))
            x, y = dot(rotate[-(c - north)], xy - u) + tc
        elif region == 'south':
            u = array((-3*pi/4 + south*pi/2, -pi/2))
            x, y = dot(rotate[c - south], xy - u) + tc
        
    return x, y
        
def cap(x, y, north=0, south=0, inverse=False):
    r"""
    Return the number of the polar cap and region that `(x, y)` lies in.
    If `inverse=False`, then assume `(x,y)` lies in the image of the HEALPix 
    projection of the unit sphere.
    If `inverse=True`, then assume `(x,y)` lies in the image of the 
    `(north, south)`-rHEALPix projection of the unit sphere.

    INPUT:
    
    - `x, y` - Coordinates in the HEALPix or rHEALPix (if `inverse=True`) 
      projection of the unit sphere.
    - `north, south` - Integers between 0 and 3 indicating the positions
      of the north pole square and south pole square respectively.
      See rhealpix_sphere() docstring for a diagram.
    - `inverse` - (Optional; default = False) Boolean. If False, then compute
      forward function. If True, then compute inverse function.
      
    OUTPUT:
    
    The pair (cap_number, region).
    Here region equals 'north' (polar), 'south' (polar), or 'equatorial', 
    indicating where `(x, y)` lies.
    If region = 'equatorial', then cap_number = None.
    Suppose now that region != 'equatorial'.
    If `inverse=False`, then cap_number is the number (= 0, 1, 2, or 3) of 
    the HEALPix polar cap Z that `(x, y)` lies in.
    If `inverse=True`, then cap_number is the number (= 0, 1, 2, or 3) of 
    the HEALPix polar cap that `(x, y)` will get moved into.
    
    EXAMPLES::
    
        >>> cap(-pi/4, pi/4 + 0.1)
        (1, 'north')
        >>> cap(-3*pi/4 + 0.1, pi/2, inverse=True)
        (1, 'north')
        
    
    NOTE:
    
    In the HEALPix projection, the polar caps are labeled 0--3 from 
    east to west like this::

            *       *       *       *       
          * 0 *   * 1 *   * 2 *   * 3 *        
        *-------*-------*-------*-------*
        |       |       |       |       |
        |       |       |       |       |
        |       |       |       |       |
        *-------*-------*-------*-------*
          * 0 *   * 1 *   * 2 *   * 3 *
            *       *       *       *
    
    In the rHEALPix projection these polar caps get rearranged
    into a square with the caps numbered `north` and `south` remaining fixed.
    For example, if `north` = 1 and `south` = 3, then the caps get
    rearranged this way:: 

        North polar square:     *-------*       
                                | * 3 * |    
                                | 0 * 2 |    
                                | * 1 * |    
                            ----*-------*----
                            
        South polar square: ----*-------*----
                                | * 3 * |
                                | 2 * 0 |
                                | * 1 * |
                                *-------*    
        
    """
    if not inverse:
        # Forward function.
        # Find the region (x, y) lies in.
        if y > pi/4:
            region = 'north'
        elif y < -pi/4:
            region = 'south'
        else:
            region = 'equatorial'
        
        # Find the cap number of (x, y) in the image of the HEALPix projection.
        if region == 'equatorial':
            cap_number = None
        else:
            if x < -pi/2:
                cap_number = 0
            elif x >= -pi/2 and x < 0:
                cap_number = 1
            elif x >= 0 and x < pi/2:
                cap_number = 2
            else:
                cap_number = 3
    else:
        # Inverse function.
        # Find the region (x, y) lies in.
        if y > pi/4:
            region = 'north'
        elif y < -pi/4:
            region = 'south'
        else:
            region = 'equatorial'
        
        # Find HEALPix polar cap number that (x, y) moves to when 
        # the rHEALPix polar square is disassembled.
        eps = 1e-15     # Fuzz to avoid some rounding errors.
        if region == 'equatorial':
            cap_number = None
        elif region == 'north':
            L1 = x - (-3*pi/4 + (north - 1)*pi/2)
            L2 = -x + (-3*pi/4 + (north + 1)*pi/2)
            if y < L1 - eps and y >= L2 - eps:
                cap_number = (north + 1) % 4
            elif y >= L1 - eps and y > L2 + eps:       
                cap_number = (north + 2) % 4
            elif y > L1 + eps and y <= L2 + eps:       
                cap_number = (north + 3) % 4
            else:
                cap_number = north
        else:
            # region == 'south':
            L1 = x - (-3*pi/4 + (south + 1)*pi/2)
            L2 = -x + (-3*pi/4 + (south - 1)*pi/2)
            if y <= L1 + eps and y > L2 + eps:       
                cap_number = (south + 1) % 4
            elif y < L1 - eps and y <= L2 + eps:       
                cap_number = (south + 2) % 4
            elif y >= L1 - eps and y < L2 - eps:       
                cap_number = (south + 3) % 4
            else:
                cap_number = south
                
    return cap_number, region   
        
def in_image(x, y, proj='healpix_sphere', north=0, south=0):
    r"""
    Return True if `(x, y)` lies in the image of the projection `proj`
    of the unit sphere, where `proj` is either
    'healpix_sphere' or 'rhealpix_sphere' and `north` and `south` indicate 
    the positions of the polar squares in case `proj` = 'rhealpix_sphere'.
    Return False otherwise.
        
    EXAMPLES::
    
        >>> eps = 0     # Test boundary points.
        >>> hp = [
        ... (-pi - eps, pi/4),
        ... (-3*pi/4, pi/2 + eps),
        ... (-pi/2, pi/4 + eps),
        ... (-pi/4, pi/2 + eps),
        ... (0, pi/4 + eps),
        ... (pi/4, pi/2 + eps),
        ... (pi/2, pi/4 + eps),
        ... (3*pi/4, pi/2 + eps),
        ... (pi + eps, pi/4),
        ... (pi + eps,-pi/4),
        ... (3*pi/4,-pi/2 - eps),
        ... (pi/2,-pi/4 - eps),
        ... (pi/4,-pi/2 - eps),
        ... (0,-pi/4 - eps),
        ... (-pi/4,-pi/2 - eps),
        ... (-pi/2,-pi/4 - eps),
        ... (-3*pi/4,-pi/2 - eps),
        ... (-pi - eps,-pi/4)
        ... ]
        >>> for p in hp:
        ...     if not in_image(*p):
        ...             print 'Fail'
        ... 
        >>> in_image(0, 0)
        True
        >>> in_image(0, pi/4 + 0.1)
        False
        >>> eps = 0     # Test boundary points.
        >>> north, south = 0, 0
        >>> rhp = [
        ... (-pi - eps, pi/4 + eps),
        ... (-pi + north*pi/2 - eps, pi/4 + eps),
        ... (-pi + north*pi/2 - eps, 3*pi/4 + eps),
        ... (-pi + (north + 1)*pi/2 + eps, 3*pi/4 + eps),
        ... (-pi + (north + 1)*pi/2 + eps, pi/4 + eps),
        ... (pi + eps, pi/4 + eps),
        ... (pi + eps,-pi/4 - eps),
        ... (-pi + (south + 1)*pi/2 + eps,-pi/4 - eps),
        ... (-pi + (south + 1)*pi/2 + eps,-3*pi/4 - eps),
        ... (-pi + south*pi/2 - eps,-3*pi/4 - eps),
        ... (-pi + south*pi/2 -eps,-pi/4 - eps),
        ... (-pi - eps,-pi/4 - eps)
        ... ]
        >>> for p in rhp:
        ...     if not in_image(*p, proj='rhealpix_sphere'):
        ...             print 'Fail'
        ... 
        >>> in_image(0, 0,proj='rhealpix_sphere')
        True
        >>> in_image(0, pi/4 + 0.1, proj='rhealpix_sphere')
        False
        
    """    
    # matplotlib is a third-party module.
    from matplotlib.nxutils import pnpoly
    
    # Fuzz to slightly expand (r)HEALPix image boundary so that pnpoly() 
    # counts points on the original boundary as lying in the image.
    eps = 1e-12
    if proj == 'healpix_sphere':
        vertices = [
        (-pi - eps, pi/4),
        (-3*pi/4, pi/2 + eps),
        (-pi/2, pi/4 + eps),
        (-pi/4, pi/2 + eps),
        (0, pi/4 + eps),
        (pi/4, pi/2 + eps),
        (pi/2, pi/4 + eps),
        (3*pi/4, pi/2 + eps),
        (pi + eps, pi/4),
        (pi + eps,-pi/4),
        (3*pi/4,-pi/2 - eps),
        (pi/2,-pi/4 - eps),
        (pi/4,-pi/2 - eps),
        (0,-pi/4 - eps),
        (-pi/4,-pi/2 - eps),
        (-pi/2,-pi/4 - eps),
        (-3*pi/4,-pi/2 - eps),
        (-pi - eps,-pi/4)
        ]
    else:
        # Assume proj == 'rhealpix_sphere'
        vertices = [
        (-pi - eps, pi/4 + eps),
        (-pi + north*pi/2 - eps, pi/4 + eps),
        (-pi + north*pi/2 - eps, 3*pi/4 + eps),
        (-pi + (north + 1)*pi/2 + eps, 3*pi/4 + eps),
        (-pi + (north + 1)*pi/2 + eps, pi/4 + eps),
        (pi + eps, pi/4 + eps),
        (pi + eps, -pi/4 - eps),
        (-pi + (south + 1)*pi/2 + eps, -pi/4 - eps),
        (-pi + (south + 1)*pi/2 + eps, -3*pi/4 - eps),
        (-pi + south*pi/2 - eps, -3*pi/4 - eps),
        (-pi + south*pi/2 -eps, -pi/4 - eps),
        (-pi - eps, -pi/4 - eps)
        ]
    return bool(pnpoly(x, y, vertices))    
        
def wrap_longitude(lam, degrees=False):
    r"""
    Given a point p on the unit circle at angle `lam` from the positive x-axis, 
    return its angle theta in the range -pi <= theta < pi.
    If `radians` = True, then `lam` and the output are given in radians.
    Otherwise, they are given in degrees.
    
    EXAMPLES::
    
        >>> wrap_longitude(2*pi + pi)
        -3.1415926535897931
        >>> wrap_longitude(-185, degrees=True)
        175.0
        >>> wrap_longitude(-180, degrees=True)
        -180.0
        >>> wrap_longitude(185, degrees=True)
        -175.0
    """
    if degrees:
        # Convert to radians.
        lam = deg2rad(lam)
    if lam < -pi or lam >= pi:
        result = lam - 2*pi*floor(lam/(2*pi))    # x mod 2*pi
        if result >= pi:
            result = result - 2*pi
    else:
        result = lam
    if degrees:
        # Convert to degrees.
        result = rad2deg(result)
    return result
        
def wrap_latitude(phi, degrees=False):
    r"""
    Given a point p on the unit circle at angle `phi` from the positive x-axis,
    if p lies in the right half of the circle, then return its angle that lies 
    in the interval [-pi/2, pi/2].
    If p lies in the left half of the circle, then reflect it through the 
    origin, and return the angle of the reflected point that lies in the 
    interval [-pi/2, pi/2].
    If `radians` = True, then `phi` and the output are given in radians.
    Otherwise, they are given in degrees.
    
    EXAMPLES::

        >>> wrap_latitude(45, degrees=True)
        45.0
        >>> wrap_latitude(-45, degrees=True)
        -45.0
        >>> wrap_latitude(90, degrees=True)
        90.0
        >>> wrap_latitude(-90, degrees=True)
        -90.0
        >>> wrap_latitude(135, degrees=True)
        -45.0
        >>> wrap_latitude(-135, degrees=True)
        45.0
    """
    if degrees:
        # Convert to radians.
        phi = deg2rad(phi)
    # Put phi in range -pi <= phi < pi.
    phi = wrap_longitude(phi)
    if abs(phi) <= pi/2:
        result = phi
    else:
        result = phi - sign(phi)*pi
    if degrees:
        # Convert to degrees.
        result = rad2deg(result)
    return result
    
def auth_lat(phi, e, inverse=False):
    r"""
    Given a point of geographic latitude `phi` on an ellipse of 
    eccentricity `e`, return the authalic latitude of the point.
    If `inverse=True`, then compute its inverse approximately.
    
    EXAMPLES::
    
        >>> beta = auth_lat(pi/4, WGS84_E); beta
        0.78315895611802222
        >>> auth_lat(beta, WGS84_E, inverse=True); pi/4
        0.78539816331575041
        0.7853981633974483
    """
    if not inverse:
        # Return authalic latitude.
        q = ((1 - e**2)*sin(phi))/(1 - (e*sin(phi))**2) - \
            (1 - e**2)/(2*e)*log((1 - e*sin(phi))/(1 + e*sin(phi)))
        qp = 1 - (1 - e**2)/(2*e)*log((1 - e)/(1 + e))
        ratio = q / qp
        # Avoid rounding errors.
        if abs(ratio) > 1:
            # Make abs(ratio) = 1
            ratio = sign(ratio) 
        return arcsin(ratio)
    else:
        # Otherwise return good approximation of inverse of authalic latitude.
        return phi + (e**2 / 3 + 31*e**4/180 + 517*e**6/5040) *\
               sin(2*phi) + (23*e**4/360 + 251*e**6/3780) *\
               sin(4*phi) + (761*e**6/45360)*sin(6*phi)
        
def para_lat(phi, e, inverse=False):
    r"""
    Given a point at geographic latitude `phi` on an ellipse of 
    eccentricity `e`, return the parametric latitude of the point.
    If `inverse=True`, then compute its inverse.
    
    EXAMPLES::
    
        >>> eta = para_lat(pi/3, WGS84_E); eta
        1.0457420826841251
        >>> para_lat(eta, WGS84_E, inverse=True); pi/3
        1.0471975511965976
        1.0471975511965976
        
    """
    if not inverse:
        # Return parametric latitude.
        if phi == pi/2 or phi == -pi/2:
            return phi
        else:
            return arctan(sqrt(1-e**2)*tan(phi))
    else:
        # Otherwise return inverse of parametric latitude.
        if phi == pi/2 or phi == -pi/2:
            return phi
        else:
            return arctan(tan(phi)/sqrt(1-e**2))
            
def ellipsoid_parameters(a, f):
    r"""
    Given the major radius `a` and flattening factor `f` of a
    ellipsoid, return its minor radius, eccentricity, and authalic 
    sphere radius. 
    
    EXAMPLES::
    
        >>> ellipsoid_parameters(WGS84_A, WGS84_F)
        (6356752.314245179, 0.081819190842621486, 6371007.1809184756)
    """
    e = sqrt(f*(2 - f))
    qp = 1 - (1 - e**2)/(2*e)*log((1 - e)/(1 + e)) 
    R_A = a*sqrt(qp/2)
    b = a*(1 - f)
    return b, e, R_A

# WGS84 ellipsoid parameters: major radius, flattening factor, minor radius,
# eccentricity, authalic radius.
WGS84_A = 6378137
WGS84_F = 1/298.257223563   
WGS84_B, WGS84_E, WGS84_R = ellipsoid_parameters(WGS84_A, WGS84_F)