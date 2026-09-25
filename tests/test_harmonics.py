"""
Tests for rhealpixdggs.harmonics, the ducc0-backed spherical harmonic
transforms on the nuclei of a grid. Skipped when ducc0 is not installed.
"""

import unittest
from importlib.util import find_spec

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.special import sph_harm_y

from rhealpixdggs.dggs import WGS84_003, RHEALPixDGGS
from rhealpixdggs.ellipsoids import WGS84_ASPHERE, WGS84_ELLIPSOID

if find_spec("ducc0"):
    from rhealpixdggs import harmonics

WGS84_122 = RHEALPixDGGS(
    ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=2, N_side=2
)
SPHERE_003 = RHEALPixDGGS(ellipsoid=WGS84_ASPHERE, N_side=3)


def nucleus_authalic_angles(rdggs, resolution):
    """Colatitude and longitude in radians, on the authalic sphere, of every
    nucleus in ring order."""
    table = rdggs.ring_table(resolution)
    ids = np.concatenate(
        [rdggs.cells_on_ring(resolution, i) for i in range(len(table.population))]
    )
    theta = np.radians(90 - np.repeat(table.authalic_latitude, table.population))
    phi = np.radians(rdggs.nuclei(ids)[:, 0])
    return ids, theta, phi


def random_alm(rng, lmax):
    """Coefficients of a real band-limited field: complex, with real m = 0."""
    alm = rng.normal(size=harmonics.alm_size(lmax)) + 1j * rng.normal(
        size=harmonics.alm_size(lmax)
    )
    alm[: lmax + 1] = alm[: lmax + 1].real
    return alm


@unittest.skipUnless(find_spec("ducc0"), "ducc0 is not installed")
class HarmonicsTestCase(unittest.TestCase):
    def test_alm_layout_is_healpy_order(self):
        # m-major triangular storage: (l, m) at m (2 lmax + 1 - m) / 2 + l,
        # a bijection onto range(alm_size(lmax)).
        for lmax in (0, 1, 5, 12):
            size = harmonics.alm_size(lmax)
            self.assertEqual(size, (lmax + 1) * (lmax + 2) // 2)
            indices = sorted(
                harmonics.alm_index(lmax, l, m)
                for m in range(lmax + 1)
                for l in range(m, lmax + 1)
            )
            self.assertEqual(indices, list(range(size)))
            self.assertEqual(harmonics.alm_index(lmax, lmax, 0), lmax)
            if lmax:
                self.assertEqual(harmonics.alm_index(lmax, 1, 1), lmax + 1)

    def test_ring_geometry_matches_ring_table(self):
        # ducc0 wants colatitude in radians, pixels per ring, the azimuth of
        # the first pixel and each ring's offset in the flat map, as unsigned
        # integers where it says so. Pole rings have azimuth 0.
        for rdggs, resolution in ((WGS84_003, 2), (WGS84_122, 3), (SPHERE_003, 1)):
            table = rdggs.ring_table(resolution)
            geometry = harmonics.ring_geometry(rdggs, resolution)
            self.assertEqual(set(geometry), {"theta", "nphi", "phi0", "ringstart"})
            assert_allclose(geometry["theta"], np.radians(90 - table.authalic_latitude))
            assert_array_equal(geometry["nphi"], table.population)
            self.assertEqual(geometry["nphi"].dtype, np.uint64)
            self.assertEqual(geometry["ringstart"].dtype, np.uint64)
            assert_array_equal(
                geometry["ringstart"], np.cumsum(table.population) - table.population
            )
            single = table.population == 1
            assert_array_equal(geometry["phi0"][single], 0.0)
            assert_allclose(
                geometry["phi0"][~single], np.radians(table.first_longitude[~single])
            )

    def test_ring_order_packing_round_trips(self):
        # to_ring_order lays cell values out as ducc0's flat map, ring by
        # ring in cells_on_ring order; from_ring_order reads them back.
        rdggs, resolution = WGS84_003, 1
        cells = [str(c) for c in rdggs.grid(resolution)]
        values = {c: float(i) for i, c in enumerate(sorted(cells))}
        flat = harmonics.to_ring_order(rdggs, resolution, values)
        self.assertEqual(flat.shape, (6 * 9,))
        geometry = harmonics.ring_geometry(rdggs, resolution)
        for i, (start, count) in enumerate(
            zip(geometry["ringstart"], geometry["nphi"], strict=True)
        ):
            ring_ids = rdggs.cells_on_ring(resolution, i)
            assert_array_equal(
                flat[int(start) : int(start + count)], [values[c] for c in ring_ids]
            )
        self.assertEqual(harmonics.from_ring_order(rdggs, resolution, flat), values)
        del values["P4"]
        with self.assertRaises(ValueError):
            harmonics.to_ring_order(rdggs, resolution, values)

    def test_synthesis_matches_scipy_at_the_nuclei(self):
        # A single unit coefficient synthesises to the real field it stands
        # for under the healpy convention: Re Y_l0 for m = 0, 2 Re Y_lm
        # otherwise, evaluated at the nuclei's authalic coordinates.
        lmax = 6
        for rdggs, resolution in ((WGS84_003, 2), (WGS84_122, 3), (SPHERE_003, 2)):
            _, theta, phi = nucleus_authalic_angles(rdggs, resolution)
            for l, m in ((0, 0), (1, 0), (1, 1), (3, 2), (5, 0), (6, 6)):
                alm = np.zeros(harmonics.alm_size(lmax), dtype=complex)
                alm[harmonics.alm_index(lmax, l, m)] = 1.0
                got = harmonics.synthesis(rdggs, resolution, alm, lmax)
                want = (1 if m == 0 else 2) * np.real(sph_harm_y(l, m, theta, phi))
                self.assertEqual(got.shape, theta.shape)
                assert_allclose(got, want, rtol=0, atol=1e-12)

    def test_synthesis_rejects_the_wrong_number_of_coefficients(self):
        with self.assertRaises(ValueError):
            harmonics.synthesis(WGS84_003, 1, np.zeros(5, dtype=complex), lmax=3)

    def test_analysis_recovers_coefficients_below_degree_2n(self):
        # The nuclei are not a quadrature, so analysis is a least-squares
        # fit. Up to about degree 1.8n it recovers a band-limited field's
        # coefficients to rounding in a dozen iterations; towards 2n the
        # fit is ill-conditioned and slow, and at 2n aliasing on the short
        # polar rings defeats it. This pins the usable degree.
        rng = np.random.default_rng(20260914)
        for rdggs, resolution in (
            (WGS84_003, 2),
            (WGS84_003, 3),
            (WGS84_122, 3),
            (WGS84_122, 4),
        ):
            n = rdggs.N_side**resolution
            for lmax in (n, int(1.8 * n)):
                alm = random_alm(rng, lmax)
                values = harmonics.synthesis(rdggs, resolution, alm, lmax)
                recovered, info = harmonics.analysis(rdggs, resolution, values, lmax)
                self.assertEqual(recovered.shape, alm.shape)
                self.assertEqual(info.stop_reason, 1)
                self.assertLessEqual(info.iterations, 20)
                self.assertLess(info.residual, 1e-8)
                assert_allclose(recovered, alm, rtol=0, atol=1e-8 * np.abs(alm).max())
            lmax = 2 * n
            alm = random_alm(rng, lmax)
            values = harmonics.synthesis(rdggs, resolution, alm, lmax)
            recovered, info = harmonics.analysis(rdggs, resolution, values, lmax)
            self.assertGreater(np.abs(recovered - alm).max() / np.abs(alm).max(), 1e-3)

    def test_analysis_of_cell_values_by_index(self):
        # The dict form: values keyed by index string go through
        # to_ring_order, and the result is the same as the flat form.
        rdggs, resolution, lmax = WGS84_003, 2, 8
        rng = np.random.default_rng(3)
        alm = random_alm(rng, lmax)
        flat = harmonics.synthesis(rdggs, resolution, alm, lmax)
        by_index = harmonics.from_ring_order(rdggs, resolution, flat)
        recovered, _ = harmonics.analysis(rdggs, resolution, by_index, lmax)
        assert_allclose(recovered, alm, rtol=0, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
