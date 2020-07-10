"""
This Python 3.3 code implements the zenithial orthotriangular (ZOT) map projection as described in [Dutt1991]_.

.. [Dutt1991] Geoffrey Dutton, Zenithial orthotriangular projection, 
Auto-Carto 10 (1991), ACSM-ASPRS Annual Convention, Baltimore, Technical Papers, v. 6, 77--95. 

CHANGELOG:

- Alexander Raichev (AR), 2013-01-30: Initial version.
- AR, 2013-07-23: Ported to Python 3.3.

NOTE:

All lengths are measured in meters and all angles are measured in radians 
unless indicated otherwise. 
By 'ellipsoid' below, i mean an oblate ellipsoid of revolution.
"""
# *****************************************************************************
#       Copyright (C) 2013 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import standard modules.
from math import pi

FLOPS = [1, 1, -1, -1, -1, 1, 1, -1]
P2 = pi / 2


def zot_sphere(lam, phi, diam=1, flops=FLOPS):
    """
    Compute the signature function of the ZOT projection of the unit sphere. 
        
    INPUT:
    
    - `lam, phi` - Geodetic longitude-latitude coordinates in radians.
      Assume -pi <= `lam` < pi and -pi/2 <= `phi` <= pi/2. 
    
    EXAMPLES::
    """
    s = diam / pi
    org = int((P2 - phi) // P2)
    octant = int((org + 1) * (lam // P2))
    x1 = 2 - ((octant + org - 1) % 2)
    x2 = 3 - x1
    hs = 1 - (2 * org)
    r = (s * FLOPS[octant], -s * hs * FLOPS[9 - octant])
    c = (-org * r[x1], -org * r[x2])
    clp = P2 - (hs * phi)
    olp = clp * (lam // P2) / P2
    px = c[x1] + r[x1] * abs(clp - olp)
    py = c[x2] + r[x2] * olp
    return px, py


def zot_sphere_inverse(x, y):
    """
    Compute the inverse of the zot_sphere().
    
    INPUT:
    
    - `x, y` - Planar coordinates in meters in the image of the
      ZOT projection of the unit sphere.
    
    EXAMPLES::
    """
    pass


def zot_ellipsoid(lam, phi, e=0.5):
    pass


def zot_ellipsoid_inverse(x, y, e=0.5):
    pass
