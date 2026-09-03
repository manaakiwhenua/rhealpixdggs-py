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

import importlib
from typing import Any, overload

import numpy as np
import pyproj

from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID, Ellipsoid

# my_round is doctest-only: the doctests use it from the module globals.
from rhealpixdggs.utils import (  # noqa: F401
    FloatArray,
    ProjectionFunction,
    _wrap_latitude_array,
    _wrap_longitude_array,
    my_round,
    wrap_latitude,
    wrap_longitude,
)

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

    def __init__(
        self,
        ellipsoid: Ellipsoid = WGS84_ELLIPSOID,
        proj: str | None = None,
        **kwargs: Any,
    ) -> None:
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
        self._f: ProjectionFunction | pyproj.Proj | None = None

    def __str__(self) -> str:
        result = ["map projection:"]
        result.append(f"    proj = {self.proj}")
        result.append(f"    kwargs = {self.kwargs}")
        result.append("    ellipsoid:")
        for k, v in sorted(self.ellipsoid.__dict__.items()):
            result.append(" " * 8 + k + " = " + str(v))
        return "\n".join(result)

    def _get_f(self) -> ProjectionFunction | pyproj.Proj | None:
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
                    print(f"Oops! Projection {self.proj} is not implemented.")
                    return None
            else:
                # Use a projection from the PROJ library.
                self._f = pyproj.Proj(proj=self.proj, a=a, e=e, **self.kwargs)
        return self._f

    @overload
    def __call__(
        self, u: float, v: float, inverse: bool = ...
    ) -> tuple[float, float] | None: ...

    @overload
    def __call__(
        self, u: FloatArray, v: FloatArray, inverse: bool = ...
    ) -> tuple[FloatArray, FloatArray] | None: ...

    @overload
    def __call__(
        self, u: float | FloatArray, v: float | FloatArray, inverse: bool = ...
    ) -> tuple[float, float] | tuple[FloatArray, FloatArray] | None: ...

    def __call__(
        self, u: float | FloatArray, v: float | FloatArray, inverse: bool = False
    ) -> tuple[float, float] | tuple[FloatArray, FloatArray] | None:
        """
        Project the point `(u, v)`, or invert it if `inverse` = True. `u` and
        `v` may be floats or numpy arrays of a common shape; arrays are
        projected in one pass and come back as a pair of float64 arrays.
        """
        ellipsoid = self.ellipsoid
        lon_0 = ellipsoid.lon_0
        lat_0 = ellipsoid.lat_0
        radians = ellipsoid.radians
        f = self._get_f()
        if f is None:
            return None
        if isinstance(u, np.ndarray) or isinstance(v, np.ndarray):
            u_arr, v_arr = np.broadcast_arrays(
                np.asarray(u, dtype=np.float64), np.asarray(v, dtype=np.float64)
            )
            if not inverse:
                lam_arr = _wrap_longitude_array(u_arr - lon_0, radians=radians)
                phi_arr = _wrap_latitude_array(v_arr - lat_0, radians=radians)
                x, y = f(lam_arr, phi_arr, radians=radians)
                return np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
            lam_arr, phi_arr = f(u_arr, v_arr, radians=radians, inverse=True)
            lam_arr = _wrap_longitude_array(
                np.asarray(lam_arr, dtype=np.float64) + lon_0, radians=radians
            )
            phi_arr = _wrap_latitude_array(
                np.asarray(phi_arr, dtype=np.float64) + lat_0, radians=radians
            )
            return lam_arr, phi_arr
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
