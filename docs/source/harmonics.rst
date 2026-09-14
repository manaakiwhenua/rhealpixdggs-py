Spherical Harmonics
===================

The nuclei of a resolution lie on isolatitude rings with equally spaced
longitudes (see :doc:`isolatitude`), which is the layout that ducc0_'s
ring-based spherical harmonic routines accept. ``rhealpixdggs.harmonics``
feeds them the ring table and handles the bookkeeping: it needs the
optional ``harmonics`` extra, which brings ducc0 (wheels are published for
Linux, macOS and Windows), and the rest of the package works without it.

.. _ducc0: https://pypi.org/project/ducc0/

Conventions
-----------

Coefficients follow healpy: a real field band-limited at degree
:math:`\ell_{max}` has one complex coefficient :math:`a_{\ell m}` for each
:math:`0 \le m \le \ell \le \ell_{max}`, stored m-major with
:math:`(\ell, m)` at ``alm_index(lmax, l, m)`` in an array of
``alm_size(lmax)`` entries, and its value is

.. math::

   f(\theta, \lambda) = \sum_{\ell} a_{\ell 0} Y_{\ell 0}(\theta, \lambda)
   + 2 \operatorname{Re} \sum_{\ell} \sum_{m > 0} a_{\ell m}
   Y_{\ell m}(\theta, \lambda),

so the :math:`m = 0` coefficients of a real field are real. Values are
attached to nuclei on the authalic sphere, where the cells are equal-area
and the rings sit at their closed-form latitudes, and they travel as a
flat array in ``ring_order``: ring by ring from the north pole, each ring
in ``cells_on_ring`` order. ``to_ring_order`` and ``from_ring_order``
convert between that array and a dictionary keyed by cell index.

Synthesis and analysis
----------------------

``synthesis`` evaluates a field given by its coefficients at every
nucleus; ``analysis`` finds the coefficients that best fit values at the
nuclei. The nuclei are not a quadrature (their rings are not at
Gauss-Legendre latitudes), so analysis is ducc0's iterative least-squares
solve, ``pseudo_analysis``, and it returns how the fit stopped alongside
the coefficients. A round trip on resolution 2::

    >>> import numpy as np
    >>> from rhealpixdggs.dggs import WGS84_003
    >>> from rhealpixdggs import harmonics
    >>> rdggs, resolution, lmax = WGS84_003, 2, 12
    >>> geometry = harmonics.ring_geometry(rdggs, resolution)
    >>> sorted(geometry), int(geometry['nphi'].sum())
    (['nphi', 'phi0', 'ringstart', 'theta'], 486)
    >>> rng = np.random.default_rng(0)
    >>> size = harmonics.alm_size(lmax)
    >>> alm = rng.normal(size=size) + 1j * rng.normal(size=size)
    >>> alm[:lmax + 1] = alm[:lmax + 1].real
    >>> values = harmonics.synthesis(rdggs, resolution, alm, lmax)
    >>> by_cell = harmonics.from_ring_order(rdggs, resolution, values)
    >>> len(by_cell), harmonics.ring_order(rdggs, resolution)[0]
    (486, np.str_('N44'))
    >>> recovered, info = harmonics.analysis(rdggs, resolution, by_cell, lmax)
    >>> info.stop_reason, info.iterations < 20
    (1, True)
    >>> bool(np.abs(recovered - alm).max() < 1e-8)
    True

``ring_geometry`` is the description ducc0 takes as keyword arguments,
straight from ``ring_table``: colatitude per ring on the authalic sphere,
cells per ring, longitude of the first nucleus, and each ring's offset in
the flat array. It is exposed so that other ducc0 routines, such as
``adjoint_synthesis`` or the spin-weighted transforms, can be called on
the same layout.

Analysing a field sampled at the nuclei works the same way. The field
:math:`\sin\theta \cos\lambda` of the :doc:`isolatitude` page is a single
spherical harmonic, and the fit finds exactly that::

    >>> ids = harmonics.ring_order(rdggs, resolution)
    >>> theta = np.repeat(geometry['theta'], geometry['nphi'].astype(int))
    >>> lon = np.radians(rdggs.nuclei(ids)[:, 0])
    >>> alm, info = harmonics.analysis(rdggs, resolution, np.sin(theta) * np.cos(lon), lmax=4)
    >>> nonzero = np.abs(alm) > 1e-8
    >>> [(l, m) for m in range(5) for l in range(m, 5) if nonzero[harmonics.alm_index(4, l, m)]]
    [(1, 1)]

How high a degree a resolution supports
---------------------------------------

An equatorial ring of :math:`4n` cells (:math:`n = N_{side}^k`) resolves
longitudinal wavenumbers up to :math:`2n`, and the polar rings far fewer,
so a band-limited field is only recoverable below degree :math:`2n`, and
the least-squares problem grows ill-conditioned as the degree approaches
that limit. Measured on WGS84 with :math:`N_{side} = 3`, synthesising
random coefficients and analysing them back with the default 100
iterations:

.. list-table::
   :header-rows: 1
   :widths: 12 12 12 14 14 36

   * - Resolution
     - Cells
     - :math:`n`
     - Degree
     - Degree / :math:`n`
     - Iterations and relative error
   * - 2
     - 486
     - 9
     - 9 / 13 / 17 / 18
     - 1.0 / 1.4 / 1.9 / 2.0
     - 6, 7, 20 iterations at :math:`10^{-10}` or better; at degree 18 the solver reports convergence but the coefficients are wrong (error 0.8)
   * - 3
     - 4374
     - 27
     - 27 / 40 / 48 / 53 / 54
     - 1.0 / 1.5 / 1.8 / 2.0 / 2.0
     - 5, 6, 12 iterations at :math:`10^{-10}`; 77 iterations at :math:`10^{-9}`; at degree 54 the iteration limit is hit and the coefficients are wrong (error 0.3)
   * - 4
     - 39366
     - 81
     - 81 / 121 / 145 / 153 / 157
     - 1.0 / 1.5 / 1.8 / 1.9 / 1.9
     - 5, 5, 12 iterations at :math:`10^{-9}` or better; 93 iterations at :math:`10^{-8}`; degree 157 needs 666 iterations for :math:`10^{-7}`

At degree :math:`2n` several different sets of coefficients produce the
same values at the nuclei, so the solver returns one of them, and it is
not the one the field was made from: the values are reproduced, the
residual is tiny, and the stop reason can still say converged. The
output does not warn you, so the degree has to be chosen in advance. A
safe rule is degrees up to about :math:`1.8 n`, where a dozen iterations
suffice; between there and :math:`2n`, raise ``maxiter`` and watch
``info.stop_reason`` (1 means converged, 7 means the iteration limit was
hit) and ``info.residual``, knowing that a clean stop is necessary but
not sufficient. These numbers are for exact
band-limited input; a real field sampled at the nuclei has power above
any chosen degree, and that power aliases into the fit, so choose the
degree from the field's smoothness, not from this table alone.

healpy's ``map2alm`` and ``alm2map`` cannot be used instead of this
module: they assume HEALPix's ring populations, which differ from
rHEALPix's.

Module reference
----------------

.. automodule:: rhealpixdggs.harmonics
    :members:
    :undoc-members:
    :show-inheritance:
