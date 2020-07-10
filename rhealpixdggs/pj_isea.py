"""
This Python 3.3 module implements the icosohedral Snyder equal area map projection (ISEA) for ellipsoids of revolution by using the spherical version of the projection (from the PROJ.4 library) on the ellipsoid's authalic sphere.

CHANGELOG:

- Alexander Raichev (AR), 2012-01-25: Initial version.
- AR, 2013-07-23: Ported to Python 3.3.

NOTE:

All lengths are measured in meters and all angles are measured in radians 
unless indicated otherwise. 
By 'ellipsoid' below, i mean an oblate ellipsoid of revolution.
"""
# *****************************************************************************
#       Copyright (C) 2012 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import third-party modules.
import pyproj
from numpy import rad2deg, deg2rad, array, pi

# Import my modules.
from rhealpixdggs.utils import my_round, auth_lat, auth_rad


def isea_sphere(lam, phi):
    """
    Compute the signature function of the icosohedral Snyder equal area
    (ISEA) map projection of the unit sphere. 
        
    INPUT:
    
    - `lam, phi` - Geodetic longitude-latitude coordinates in radians.
      Assume -pi <= `lam` < pi and -pi/2 <= `phi` <= pi/2. 
    
    EXAMPLES::
    
    """
    f = pyproj.Proj(proj="isea", R=1)
    return f(lam, phi, radians=True)


def isea_sphere_inverse(x, y):
    """
    Compute the inverse of isea_sphere().

    INPUT:
    
    - `x, y` - Planar coordinates in meters in the image of the
      ISEA projection of the unit sphere.
    
    NOTES:
    
    Not implemented here nor in PROJ.4.
    """
    pass


def isea_ellipsoid(lam, phi, e=0):
    """
    Compute the signature functions of the ISEA may projection of an oblate
    ellipsoid with eccentricity `e` whose authalic sphere is the unit sphere.
        
    INPUT:
    
    - `lam, phi` - Geodetic longitude-latitude coordinates in radians.
      Assume -pi <= `lam` < pi and -pi/2 <= `phi` <= pi/2. 
    - `e` - Eccentricity of the oblate ellipsoid.
    """
    beta = auth_lat(phi, e, radians=True)
    return isea_sphere(lam, beta)


def isea_ellipsoid_inverse(x, y, e=0):
    """
    Compute the inverse of isea_ellipsoid().
    
    INPUT:
    
    - `x, y` - Planar coordinates in meters in the image of the
      ISEA projection of the ellipsoid with eccentricity `e` whose
      authalic sphere is the unit sphere.
    
    NOTES:
    
    Not implemented here nor in PROJ.4.        
    """
    pass


def isea(a=1, e=0):
    """
    Return a function object that wraps the ISEA projection and its inverse
    of an ellipsoid with major radius `a` and eccentricity `e`.
    
    EXAMPLES::
    
        >>> f = isea(a=2, e=0)
        >>> print(my_round(f(0, pi/3, radians=True), 15))
        (-4.7097959155097699, 2.9195761776404252)
        >>> g = isea(a=2, e=0.1)
        >>> print(my_round(g(0, 60, radians=False), 15))
        (-4.6978889550868486, 2.9179977608222689)
        
    OUTPUT:
    
    - A function object of the form f(u, v, radians=False, inverse=False).
    """
    R_A = auth_rad(a, e)

    def f(u, v, radians=False, inverse=False):
        if not inverse:
            lam, phi = u, v
            if not radians:
                # Convert to radians.
                lam, phi = deg2rad([lam, phi])
            return tuple(R_A * array(isea_ellipsoid(lam, phi, e=e)))
        else:
            # Scale down to R_A = 1.
            x, y = array((u, v)) / R_A
            lam, phi = array(isea_ellipsoid_inverse(x, y, e=e))
            if not radians:
                # Convert to degrees.
                lam, phi = rad2deg([lam, phi])
            return lam, phi

    return f
