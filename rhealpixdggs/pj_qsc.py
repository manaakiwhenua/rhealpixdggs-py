"""
This Python 3.3 module implements the quadrilateralized spherical cube (QSC) map projection as described in [CaGr2002]_.

.. [CaGr2002] M. R. Calabretta and E. W. Greisen, Representations of celestial coordinates in FITS, Astronomy & Astrophysics 395 (2002), 1077--1122.

CHANGELOG:

- Alexander Raichev (AR), 2012-02-04: Initial version.
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
from numpy import pi, rad2deg, deg2rad, sin, cos, argmax, sqrt, arctan, arcsin, array

# Import my modules.
from rhealpixdggs.utils import auth_lat, auth_rad

face_centers_lonlat = [
    (0, pi / 2),
    (0, 0),
    (pi / 2, 0),
    (pi, 0),
    (3 * pi / 2, 0),
    (0, -pi / 2),
]


def qsc_sphere(lam, phi):
    """
    Compute the signature function of the icosohedral Snyder equal area
    (ISEA) map projection of the unit sphere. 
        
    INPUT:
    
    - `lam, phi` - Geodetic longitude-latitude coordinates in radians.
      Assume -pi <= `lam` < pi and -pi/2 <= `phi` <= pi/2.     
    """
    # Find out what face (lam, phi) lies on.
    L = cos(lam) * cos(phi)
    M = sin(lam) * cos(phi)
    N = sin(phi)
    face = argmax([N, L, M, -L, -M, -N])
    if face == 0:
        xi, eta, zeta = M, -L, N
    elif face == 1:
        xi, eta, zeta = M, N, L
    elif face == 2:
        xi, eta, zeta = -L, N, M
    elif face == 3:
        xi, eta, zeta = -M, N, -L
    elif face == 4:
        xi, eta, zeta = L, N, -M
    else:
        # face = 5
        xi, eta, zeta = M, L, -N
    lam_0, phi_0 = face_centers_lonlat[face]
    if xi > abs(eta) or eta > abs(xi):
        S = 1
    else:
        S = -1
    if abs(xi) > abs(eta):
        omega = eta / xi
    else:
        if eta == 0:
            omega = 0
        else:
            omega = xi / eta
    u = (pi / 4) * S * sqrt((1 - zeta) / (1 - 1 / sqrt(2 + omega ** 2)))
    v = u / (pi / 12) * (arctan(omega) - arcsin(omega / (sqrt(2 * (1 + omega ** 2)))))
    if abs(xi) > abs(eta):
        x, y = array([lam_0, phi_0]) + array([u, v])
    else:
        x, y = array([lam_0, phi_0]) + array([v, u])
    return x, y


def qsc_sphere_inverse(x, y):
    """
    Compute the inverse of qsc_sphere().

    INPUT:
    
    - `x, y` - Planar coordinates in meters in the image of the
      ISEA projection of the unit sphere.
    
    NOTES:
    
    Not implemented here nor in PROJ.4.
    """
    pass


def qsc_ellipsoid(lam, phi, e=0):
    """
    Compute the signature functions of the ISEA may projection of an oblate
    ellipsoid with eccentricity `e` whose authalic sphere is the unit sphere.
        
    INPUT:
    
    - `lam, phi` - Geodetic longitude-latitude coordinates in radians.
      Assume -pi <= `lam` < pi and -pi/2 <= `phi` <= pi/2. 
    - `e` - Eccentricity of the oblate ellipsoid.
    """
    beta = auth_lat(phi, e, radians=True)
    return qsc_sphere(lam, beta)


def qsc_ellipsoid_inverse(x, y, e=0):
    """
    Compute the inverse of qsc_ellipsoid().
    
    INPUT:
    
    - `x, y` - Planar coordinates in meters in the image of the
      ISEA projection of the ellipsoid with eccentricity `e` whose
      authalic sphere is the unit sphere.
    
    NOTES:
    
    Not implemented here nor in PROJ.4.        
    """
    pass


def qsc(a=1, e=0):
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
            return tuple(R_A * array(qsc_ellipsoid(lam, phi, e=e)))
        else:
            # Scale down to R_A = 1.
            x, y = array((u, v)) / R_A
            lam, phi = array(qsc_ellipsoid_inverse(x, y, e=e))
            if not radians:
                # Convert to degrees.
                lam, phi = rad2deg([lam, phi])
            return lam, phi

    return f
