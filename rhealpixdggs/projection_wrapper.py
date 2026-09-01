"""
This Python 3.11 module implements a wrapper for map projections.

- Alexander Raichev (AR), 2013-01-25: Refactored code from release 0.3.

NOTE:

All lengths are measured in meters and all angles are measured in radians
unless indicated otherwise.
By 'ellipsoid' below, I mean an oblate ellipsoid of revolution.
"""

# *****************************************************************************
#       Copyright (C) 2013 Alexander Raichev <alex.raichev@gmail.com>
#
#  Distributed under the terms of the GNU Lesser General Public License (LGPL)
#                  http://www.gnu.org/licenses/
# *****************************************************************************

# Import third-party modules.
# Import standard modules.
import importlib

import pyproj

from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID

# Import my modules.
# my_round is doctest-only: the doctests use it from the module globals.
from rhealpixdggs.utils import my_round, wrap_latitude, wrap_longitude  # noqa: F401

# Homemade map projections, as opposed to those in the PROJ.4 library.
# Remove 'healpix' and 'rhealpix' to use the PROJ.4 versions instead,
# assuming you have the *correct/patched* PROJ.4 versions.
HOMEMADE_PROJECTIONS = {"healpix", "rhealpix", "isea", "csea", "qsc"}
# HOMEMADE_PROJECTIONS = {"isea", "csea", "qsc"}


class Projection:
    """
    Represents a map projection of a given ellipsoid.

    INSTANCE ATTRIBUTES:

    - `ellipsoid` - An ellipsoid (Ellipsoid instance) to project.
    - `proj` - The name (string) of the map projection, either a valid PROJ.4
      projection name or a valid homemade projection name.
    - `kwargs` - Keyword arguments (dictionary) needed for the projection's
      definition, but not for the definition of the ellipsoid.  For example,
      these could be {'north_square':1, 'south_square': 2} for the rhealpix
      projection.

    EXAMPLES::

        >>> from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
        >>> f = Projection(ellipsoid=WGS84_ELLIPSOID, proj='rhealpix', north_square=1, south_square=0)
        >>> print(tuple(x.tolist() for x in my_round(f(0, 30), 15)))
        (0.0, 3740232.8933662786)
        >>> f = Projection(ellipsoid=WGS84_ELLIPSOID, proj='cea')
        >>> print(my_round(f(0, 30), 15))
        (0.0, 3171259.315518537)

    NOTES:

    When accessing a homemade map projection assume that it can be called via
    a function g(a, e), where a is the major radius of the ellipsoid to be
    projected and e is its eccentricity.
    The output of g should be a function object of the form
    f(u, v, radians=False, inverse=False).
    For example, see the healpix() function in ``pj_healpix.py``.
    """

    def __init__(self, ellipsoid=WGS84_ELLIPSOID, proj=None, **kwargs):
        self.proj = proj
        # Keyword arguments related to the projection but not to its
        # underlying ellipsoid, e.g. for rhealpix these could be
        # {'north_square':1, 'south_square': 2}:
        self.kwargs = kwargs
        self.ellipsoid = ellipsoid
        # Lazily built and cached by _get_f(). `a`/`e` (the only ellipsoid
        # attributes this depends on) never change after construction, so
        # the underlying projection callable only ever needs to be built
        # once, no matter how many times __call__() is invoked.
        self._f = None

    def __str__(self):
        result = ["map projection:"]
        result.append("    proj = %s" % self.proj)
        result.append("    kwargs = %s" % self.kwargs)
        result.append("    ellipsoid:")
        for k, v in sorted(self.ellipsoid.__dict__.items()):
            result.append(" " * 8 + k + " = " + str(v))
        return "\n".join(result)

    def _get_f(self):
        """
        Return the underlying f(u, v, radians=False, inverse=False)
        callable for this projection, building it on first use and
        caching it thereafter. Building it (re-importing its module and,
        for homemade projections, recomputing the authalic radius) is not
        free, and __call__() may invoke it dozens of times per cell for
        callers like Cell.boundary().
        """
        if self._f is None:
            a = self.ellipsoid.a
            e = self.ellipsoid.e
            if self.proj in HOMEMADE_PROJECTIONS:
                try:
                    # Import projection module for proj.
                    module = importlib.import_module("rhealpixdggs.pj_" + self.proj)
                    self._f = getattr(module, self.proj)(a=a, e=e, **self.kwargs)
                except (AttributeError, ModuleNotFoundError):
                    print("Oops! Projection %s is not implemented." % self.proj)
                    return None
            else:
                # Use a projection from the PROJ library.
                self._f = pyproj.Proj(proj=self.proj, a=a, e=e, **self.kwargs)
        return self._f

    def __call__(self, u, v, inverse=False):
        ellipsoid = self.ellipsoid
        lon_0 = ellipsoid.lon_0
        lat_0 = ellipsoid.lat_0
        radians = ellipsoid.radians
        f = self._get_f()
        if f is None:
            return
        if not inverse:
            # Translate longitudes and latitudes so that
            # (lon_0, lat_0) maps to (0, 0) in the plane.
            lam = wrap_longitude(u - lon_0, radians=radians)
            phi = wrap_latitude(v - lat_0, radians=radians)
            return f(lam, phi, radians=radians)
        else:
            lam, phi = f(u, v, radians=radians, inverse=True)
            # Translate longitudes and latitudes so that
            # (0, 0) in the plane maps to (lon_0, lat_0) on the ellipsoid.
            lam = wrap_longitude(lam + lon_0, radians=radians)
            phi = wrap_latitude(phi + lat_0, radians=radians)
            return lam, phi
