"""
This Python 3.3 module implements the cubic Snyder equal area map projection (CSEA) as described in [Snyd1992]_.

.. [Snyd1992] J. P. Snyder, An equal-area map projection for polyhedral globes, Cartographica 29 (1992), 10--21.

CHANGELOG:

- Alexander Raichev (AR), 2013-01-31: Initial version.
- AR, 2013-02-08: Fixed a bug in the xy-positions of the square faces.
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
from numpy import (
    pi,
    rad2deg,
    deg2rad,
    floor,
    sin,
    cos,
    arccos,
    arctan,
    arctan2,
    sqrt,
    array,
)

# Import my modules.
from rhealpixdggs.utils import auth_lat, auth_rad

g = arctan(sqrt(2))  # Snyder approximates this as 54.73561032 degrees.
tang = sqrt(2)
G = pi / 3
theta = pi / 4  # = 2*pi divided by twice the number of sides of a square
cottheta = 1
alpha = pi / 2  # = 2*(pi/2 - theta)
R_prime = sqrt(pi / 6)  # = (1/tang)*sqrt(2*(G - theta)/(sin(theta)*cos(theta)))
A = pi / 2 - g
face_centers_lonlat = [
    (-2 * pi / 3, A),
    (-pi / 3, -A),
    (0, A),
    (pi / 3, -A),
    (2 * pi / 3, A),
    (pi, -A),
]
B = sqrt(pi / 3)  # = R_prime*tang
face_centers_xy = [
    (-2 * B, B / 2),
    (-B, -B / 2),
    (0, B / 2),
    (B, -B / 2),
    (2 * B, B / 2),
    (3 * B, -B / 2),
]
EPS = 1e-12  # Fuzz to avoid rounding errors.


def csea_sphere(lam, phi):
    """
    Compute the signature function of the icosohedral Snyder equal area
    (ISEA) map projection of the unit sphere. 
        
    INPUT:
    
    - `lam, phi` - Geodetic longitude-latitude coordinates in radians.
      Assume -pi <= `lam` < pi and -pi/2 <= `phi` <= pi/2.     
    """
    # (lam, phi) lies in one of these two squares:
    face_1 = int(floor((lam + pi) / (pi / 3)))
    face_2 = (face_1) - 1 % 6
    for i in [face_1, face_2]:
        lam_0, phi_0 = face_centers_lonlat[i]
        z = arccos(sin(phi_0) * sin(phi) + cos(phi_0) * cos(phi) * (cos(lam - lam_0)))
        if z > g:
            # (lam, phi) does not lie in this square
            continue
        # Might have found the correct square face.
        Az = arctan2(
            cos(phi) * sin(lam - lam_0),
            cos(phi_0) * sin(phi) - sin(phi_0) * cos(phi) * cos(lam - lam_0),
        )
        j = 0
        while Az <= 0 or Az > alpha:
            if Az <= 0:
                Az += alpha
                j += 1
            else:
                Az -= alpha
                j -= 1
        q = arctan(tang / (cos(Az) + sin(Az) * cottheta))
        if z > q + EPS:
            # (lam, phi) does not lie in this square
            continue
        # Found the correct square face.
        H = arccos(sin(Az) * sin(G) * cos(g) - cos(Az) * cos(G))
        A_G = Az + G + H - pi
        Az_prime = arctan2(2 * A_G, (R_prime * tang) ** 2 - 2 * A_G * cottheta)
        d_prime = R_prime * tang / (cos(Az_prime) + sin(Az_prime) * cottheta)
        f = d_prime / (2 * R_prime * sin(q / 2))
        rho = 2 * R_prime * f * sin(z / 2)
        Az_prime -= j * alpha
        x = rho * sin(Az_prime) + face_centers_xy[i][0]
        y = rho * cos(Az_prime) + face_centers_xy[i][1]
        return x, y


def csea_sphere_inverse(x, y):
    """
    Compute the inverse of csea_sphere().

    INPUT:
    
    - `x, y` - Planar coordinates in meters in the image of the
      ISEA projection of the unit sphere.
    
    NOTES:
    
    Not implemented here nor in PROJ.4.
    """
    pass


def csea_ellipsoid(lam, phi, e=0):
    """
    Compute the signature functions of the ISEA may projection of an oblate
    ellipsoid with eccentricity `e` whose authalic sphere is the unit sphere.
        
    INPUT:
    
    - `lam, phi` - Geodetic longitude-latitude coordinates in radians.
      Assume -pi <= `lam` < pi and -pi/2 <= `phi` <= pi/2. 
    - `e` - Eccentricity of the oblate ellipsoid.
    """
    beta = auth_lat(phi, e, radians=True)
    return csea_sphere(lam, beta)


def csea_ellipsoid_inverse(x, y, e=0):
    """
    Compute the inverse of csea_ellipsoid().
    
    INPUT:
    
    - `x, y` - Planar coordinates in meters in the image of the
      ISEA projection of the ellipsoid with eccentricity `e` whose
      authalic sphere is the unit sphere.
    
    NOTES:
    
    Not implemented here nor in PROJ.4.        
    """
    pass


def csea(a=1, e=0):
    """
    Return a function object that wraps the ISEA projection and its inverse
    of an ellipsoid with major radius `a` and eccentricity `e`.
                
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
            return tuple(R_A * array(csea_ellipsoid(lam, phi, e=e)))
        else:
            # Scale down to R_A = 1.
            x, y = array((u, v)) / R_A
            lam, phi = array(csea_ellipsoid_inverse(x, y, e=e))
            if not radians:
                # Convert to degrees.
                lam, phi = rad2deg([lam, phi])
            return lam, phi

    return f
