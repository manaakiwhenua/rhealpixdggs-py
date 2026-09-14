"""
Spherical harmonic transforms on the nuclei of an rHEALPix grid.

The nuclei of a resolution lie on isolatitude rings with equally spaced
longitudes (see :doc:`isolatitude`), which is the layout ducc0's ring-based
spherical harmonic routines accept. This module builds ducc0's description
of the rings from ``RHEALPixDGGS.ring_table``, lays cell values out in
``cells_on_ring`` order, and wraps ``ducc0.sht.synthesis`` (coefficients to
nucleus values) and ``ducc0.sht.pseudo_analysis`` (nucleus values to
coefficients, a least-squares fit, since the nuclei are not a quadrature).

Coefficients follow the healpy convention: a real field has complex
coefficients :math:`a_{\\ell m}` for :math:`0 \\le m \\le \\ell \\le
\\ell_{max}`, stored m-major with :math:`(\\ell, m)` at ``alm_index(lmax, l,
m)``, and its value is :math:`\\sum_\\ell a_{\\ell 0} Y_{\\ell 0} + 2
\\operatorname{Re} \\sum_{m > 0} a_{\\ell m} Y_{\\ell m}`.

The transforms work on the authalic sphere, where the cells are equal-area
and the rings sit at the closed-form latitudes. With ``n`` cells per
base-cell side, a field band-limited up to about degree ``1.8 n`` is
recovered by ``analysis`` in a dozen iterations; towards degree ``2 n`` the
fit is ill-conditioned and slow, and at ``2 n`` the short polar rings alias
and it fails.

Requires ducc0: ``pip install rhealpixdggs[harmonics]``.
"""

from collections.abc import Mapping
from typing import Any, NamedTuple

import numpy as np

try:
    import ducc0
except ImportError as error:  # pragma: no cover
    raise ImportError(
        "rhealpixdggs.harmonics needs ducc0; install it with "
        "`pip install rhealpixdggs[harmonics]`"
    ) from error

from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs.utils import FloatArray


def alm_size(lmax: int) -> int:
    """
    Return the number of coefficients of a real field band-limited at
    degree `lmax`: one for each ``0 <= m <= l <= lmax``.

    EXAMPLES::

        >>> alm_size(2)
        6

    """
    return (lmax + 1) * (lmax + 2) // 2


def alm_index(lmax: int, l: int, m: int) -> int:
    """
    Return the position of the coefficient of degree `l` and order `m` in
    the m-major layout healpy and ducc0 use for a field band-limited at
    degree `lmax`.

    EXAMPLES::

        >>> [alm_index(2, l, m) for m in range(3) for l in range(m, 3)]
        [0, 1, 2, 3, 4, 5]

    """
    return m * (2 * lmax + 1 - m) // 2 + l


def ring_geometry(rdggs: RHEALPixDGGS, resolution: int) -> dict[str, np.ndarray]:
    """
    Return the description of the resolution `resolution` nuclei that
    ducc0's ring-based routines take as keyword arguments: `theta`, the
    colatitude of each ring on the authalic sphere in radians; `nphi`, the
    cells on it; `phi0`, the longitude in radians of its first nucleus in
    ``cells_on_ring`` order (0 on a single-cell pole ring); and
    `ringstart`, its offset in the flat array ``to_ring_order`` builds.

    EXAMPLES::

        >>> from rhealpixdggs.dggs import WGS84_003
        >>> geometry = ring_geometry(WGS84_003, 1)
        >>> geometry['nphi'].tolist()
        [1, 8, 12, 12, 12, 8, 1]
        >>> geometry['ringstart'].tolist()
        [0, 1, 9, 21, 33, 45, 53]

    """
    table = rdggs.ring_table(resolution)
    scale = 1.0 if rdggs.ellipsoid.radians else np.pi / 180
    quarter = np.pi / 2
    return {
        "theta": quarter - table.authalic_latitude * scale,
        "nphi": table.population.astype(np.uint64),
        "phi0": np.nan_to_num(table.first_longitude) * scale,
        "ringstart": (np.cumsum(table.population) - table.population).astype(np.uint64),
    }


def ring_order(rdggs: RHEALPixDGGS, resolution: int) -> np.ndarray:
    """
    Return the index strings of every resolution `resolution` cell in the
    order of the flat array the transforms use: ring by ring from the north
    pole, each ring in ``cells_on_ring`` order.

    EXAMPLES::

        >>> from rhealpixdggs.dggs import WGS84_003
        >>> ring_order(WGS84_003, 1)[:10].tolist()
        ['N4', 'N6', 'N7', 'N8', 'N5', 'N2', 'N1', 'N0', 'N3', 'O0']

    """
    rings = len(rdggs.ring_table(resolution).population)
    return np.concatenate([rdggs.cells_on_ring(resolution, i) for i in range(rings)])


def to_ring_order(
    rdggs: RHEALPixDGGS, resolution: int, values: Mapping[str, float]
) -> FloatArray:
    """
    Return the values of every resolution `resolution` cell, keyed by index
    string in `values`, as the flat float64 array the transforms use (see
    ``ring_order``). Every cell of the resolution must be present: raise a
    ValueError naming the first missing one otherwise.

    EXAMPLES::

        >>> from rhealpixdggs.dggs import WGS84_003
        >>> values = {str(c): c.nucleus(plane=False)[1] for c in WGS84_003.grid(1)}
        >>> to_ring_order(WGS84_003, 1, values).round(3)[:9].tolist()
        [90.0, 58.528, 58.528, 58.528, 58.528, 58.528, 58.528, 58.528, 58.528]

    """
    order = ring_order(rdggs, resolution)
    missing = [index for index in order if index not in values]
    if missing:
        raise ValueError(
            f"values lack {len(missing)} of the {len(order)} resolution "
            f"{resolution} cells, starting with {missing[0]!r}"
        )
    return np.array([values[index] for index in order], dtype=np.float64)


def from_ring_order(
    rdggs: RHEALPixDGGS, resolution: int, values: FloatArray
) -> dict[str, float]:
    """
    Return the flat array `values`, laid out as ``ring_order`` describes,
    as a dictionary keyed by index string: the inverse of ``to_ring_order``.

    EXAMPLES::

        >>> from rhealpixdggs.dggs import WGS84_003
        >>> import numpy as np
        >>> back = from_ring_order(WGS84_003, 1, np.arange(54.0))
        >>> back['N4'], back['O0']
        (0.0, 9.0)

    """
    order = ring_order(rdggs, resolution)
    flat = np.asarray(values, dtype=np.float64).reshape(-1)
    if flat.shape != order.shape:
        raise ValueError(
            f"expected {len(order)} values for resolution {resolution}, got {len(flat)}"
        )
    return {index: float(value) for index, value in zip(order, flat)}


def _check_alm(alm: Any, lmax: int) -> np.ndarray:
    coefficients = np.asarray(alm, dtype=np.complex128).reshape(-1)
    if len(coefficients) != alm_size(lmax):
        raise ValueError(
            f"lmax {lmax} needs {alm_size(lmax)} coefficients, got {len(coefficients)}"
        )
    return coefficients


def synthesis(
    rdggs: RHEALPixDGGS,
    resolution: int,
    alm: Any,
    lmax: int,
    nthreads: int = 1,
) -> FloatArray:
    """
    Return the value at every resolution `resolution` nucleus, in
    ``ring_order``, of the real field with coefficients `alm` (healpy
    layout, ``alm_size(lmax)`` complex numbers), evaluated on the authalic
    sphere by ``ducc0.sht.synthesis``.

    EXAMPLES::

        >>> from rhealpixdggs.dggs import WGS84_003
        >>> import numpy as np
        >>> alm = np.zeros(alm_size(1), dtype=complex)
        >>> alm[alm_index(1, 1, 0)] = 1.0  # Y_10, proportional to cos(colatitude)
        >>> values = synthesis(WGS84_003, 1, alm, lmax=1)
        >>> values.round(4)[[0, 9, 21, 33]].tolist()
        [0.4886, 0.2172, 0.0, -0.2172]

    """
    coefficients = _check_alm(alm, lmax)
    values = ducc0.sht.synthesis(
        alm=coefficients[None, :],
        lmax=lmax,
        spin=0,
        nthreads=nthreads,
        **ring_geometry(rdggs, resolution),
    )
    return np.asarray(values[0], dtype=np.float64)


class AnalysisInfo(NamedTuple):
    """
    How ``analysis`` stopped: `stop_reason` as ducc0 reports it (1: the
    solution converged to the requested tolerance, 2: a least-squares
    solution was found, 3: the problem is too ill-conditioned, 7: `maxiter`
    was reached), the `iterations` taken, the `residual` norm relative to
    the norm of the values, and ducc0's `quality` measure of the
    least-squares solution.
    """

    stop_reason: int
    iterations: int
    residual: float
    quality: float


def analysis(
    rdggs: RHEALPixDGGS,
    resolution: int,
    values: Mapping[str, float] | FloatArray,
    lmax: int,
    maxiter: int = 100,
    epsilon: float = 1e-10,
    nthreads: int = 1,
) -> tuple[np.ndarray, AnalysisInfo]:
    """
    Return the coefficients (healpy layout, ``alm_size(lmax)`` complex
    numbers) of the real field band-limited at degree `lmax` that best fits
    `values` at the resolution `resolution` nuclei, and how the fit stopped.
    `values` is either a flat array in ``ring_order`` or a mapping from
    index string to value, which goes through ``to_ring_order``.

    The nuclei are not a quadrature, so this is ``ducc0.sht.pseudo_analysis``,
    an iterative least-squares solve with at most `maxiter` iterations to
    relative tolerance `epsilon`. With ``n = N_side ** resolution``, a
    field band-limited up to about degree ``1.8 n`` is recovered to
    rounding in a dozen iterations; between that and ``2 n`` the fit is
    ill-conditioned and may need `maxiter` raised; at ``2 n`` and above the
    short polar rings alias and the result is not the coefficients you
    started from. Check ``info.stop_reason`` and ``info.residual`` rather
    than trusting the result.

    EXAMPLES::

        >>> from rhealpixdggs.dggs import WGS84_003
        >>> import numpy as np
        >>> alm = np.zeros(alm_size(2), dtype=complex)
        >>> alm[alm_index(2, 2, 1)] = 0.5 - 0.25j
        >>> values = synthesis(WGS84_003, 1, alm, lmax=2)
        >>> recovered, info = analysis(WGS84_003, 1, values, lmax=2)
        >>> bool(np.abs(recovered - alm).max() < 1e-9), info.stop_reason
        (True, 1)

    """
    if isinstance(values, Mapping):
        flat = to_ring_order(rdggs, resolution, values)
    else:
        flat = np.asarray(values, dtype=np.float64).reshape(-1)
    geometry = ring_geometry(rdggs, resolution)
    cells = int(geometry["nphi"].sum())
    if len(flat) != cells:
        raise ValueError(
            f"expected {cells} values for resolution {resolution}, got {len(flat)}"
        )
    alm, stop_reason, iterations, residual, quality = ducc0.sht.pseudo_analysis(
        map=flat[None, :],
        lmax=lmax,
        spin=0,
        maxiter=maxiter,
        epsilon=epsilon,
        nthreads=nthreads,
        **geometry,
    )
    info = AnalysisInfo(
        int(stop_reason), int(iterations), float(residual), float(quality)
    )
    return np.asarray(alm[0], dtype=np.complex128), info
