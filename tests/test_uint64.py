"""
Tests for the uint64 zone identifiers and level-row-col positions of cells.
"""

import unittest

import numpy as np

from rhealpixdggs.dggs import WGS84_002, WGS84_003, WGS84_003_OGC, RHEALPixDGGS
from rhealpixdggs.rhp_wrappers import rhp_from_uint64, rhp_to_uint64

# Grids with the polar squares in other columns and N_side 2.
SQUARES_12 = RHEALPixDGGS(N_side=3, north_square=1, south_square=2)
SQUARES_31 = RHEALPixDGGS(N_side=2, north_square=3, south_square=1)


class Uint64TestCase(unittest.TestCase):
    # uint64 identifiers from Ecere's DGGAL 0.0.7, which implements the
    # registered DGGRS, at levels 1 to 14 on every kind of face.
    DGGAL = (
        ("R5", 576460756598390795),
        ("N4", 576460753377165313),
        ("S4", 576460759819616257),
        ("N080", 1729382263352721414),
        ("N456", 1729382271942656015),
        ("S453", 1729382328850972687),
        ("R553110", 3458765673461713695),
        ("N438261", 3458764926137401660),
        ("S444513", 3458766466956919156),
        ("P210222", 3458765296578331946),
        ("O5753828785", 5764709888211076940),
        ("N4417831111", 5764637876641887222),
        ("S4706361228", 5764771556425230009),
        ("N45633856162200", 8073216189467775304),
        ("S45654453234507", 8083479663435502261),
        ("Q58805442070635", 8078850503671873639),
    )

    def cell(self, rdggs, index):
        return rdggs.cell(rdggs.parse_index(index))

    def test_matches_dggal(self):
        g = WGS84_003_OGC
        for index, value in self.DGGAL:
            with self.subTest(index=index):
                self.assertEqual(self.cell(g, index).uint64, value)
                self.assertEqual(str(g.cell_from_uint64(value)), index)
        indices = [index for index, _ in self.DGGAL]
        values = [value for _, value in self.DGGAL]
        encoded = g.uint64_from_indices(indices)
        self.assertEqual(encoded.dtype, np.uint64)
        self.assertEqual(encoded.tolist(), values)
        self.assertEqual(g.indices_from_uint64(encoded), indices)

    def test_level_row_col(self):
        g = WGS84_003
        # Level 0: N above O in column 0, O to R across the middle band.
        expected = {
            "N": (0, 0, 0),
            "O": (0, 1, 0),
            "P": (0, 1, 1),
            "Q": (0, 1, 2),
            "R": (0, 1, 3),
            "S": (0, 2, 0),
        }
        for index, lrc in expected.items():
            with self.subTest(index=index):
                self.assertEqual(self.cell(g, index).level_row_col, lrc)
        self.assertEqual(self.cell(g, "P4").level_row_col, (1, 4, 4))
        self.assertEqual(self.cell(g, "P44").level_row_col, (2, 13, 13))
        self.assertEqual(self.cell(SQUARES_12, "N").level_row_col, (0, 0, 1))
        self.assertEqual(self.cell(SQUARES_12, "S").level_row_col, (0, 2, 2))

    def test_round_trips(self):
        for rdggs in (WGS84_003_OGC, WGS84_002, SQUARES_12, SQUARES_31):
            for resolution in range(4):
                with self.subTest(N_side=rdggs.N_side, resolution=resolution):
                    cells = list(rdggs.grid(resolution))
                    indices = [str(c) for c in cells]
                    values = [c.uint64 for c in cells]
                    self.assertEqual(len(set(values)), len(values))
                    for cell, value in zip(cells, values, strict=True):
                        self.assertEqual(rdggs.cell_from_uint64(value), cell)
                        self.assertEqual(
                            rdggs.cell_from_level_row_col(*cell.level_row_col), cell
                        )
                        self.assertEqual(
                            rdggs.cell_from_level_row_col_text(cell.level_row_col_text),
                            cell,
                        )
                    self.assertEqual(
                        rdggs.uint64_from_indices(indices).tolist(), values
                    )
                    self.assertEqual(rdggs.indices_from_uint64(values), indices)

    def test_grid_without_index_strings(self):
        # N_side 5 has no index strings; the uint64 is its compact identifier.
        rdggs = RHEALPixDGGS(N_side=5)
        for cell in rdggs.grid(2):
            self.assertEqual(rdggs.cell_from_uint64(cell.uint64), cell)
        with self.assertRaises(ValueError):
            rdggs.uint64_from_indices(["N"])
        with self.assertRaises(ValueError):
            rdggs.indices_from_uint64([0])

    def test_sorted_values_are_scanline_over_the_frame(self):
        cells = list(WGS84_003.grid(2))
        by_value = sorted(cells, key=lambda c: c.uint64)
        by_position = sorted(cells, key=lambda c: c.level_row_col)
        self.assertEqual(by_value, by_position)

    def test_text(self):
        cell = self.cell(WGS84_003_OGC, "R553110")
        self.assertEqual(cell.level_row_col_text, "6-438-B1F")
        self.assertEqual(format(cell.uint64, "016X"), "3000010E00000B1F")
        self.assertEqual(WGS84_003_OGC.cell_from_level_row_col_text("6-438-b1f"), cell)
        for text in ("6-438", "6-438-B1G", "6_438_B1F", " 6-438-B1F", "-6-438-B1F"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                WGS84_003_OGC.cell_from_level_row_col_text(text)

    def test_rejects_values_naming_no_cell(self):
        g = WGS84_003
        p = 3  # cells per face side at level 1
        cases = [
            (0, 0, 1),  # beside N, where no face is
            (0, 2, 3),  # beside S
            (1, 3 * p, 0),  # below the frame
            (1, 0, 4 * p),  # right of the frame
            (g.max_resolution + 1, 0, 0),  # finer than the grid goes
            (-1, 0, 0),
        ]
        for level, row, col in cases:
            with self.subTest(level=level, row=row, col=col):
                with self.assertRaises(ValueError):
                    g.cell_from_level_row_col(level, row, col)
                if level >= 0:
                    value = level << 59 | row << 30 | col
                    with self.assertRaises(ValueError):
                        g.cell_from_uint64(value)
                    with self.assertRaises(ValueError):
                        g.indices_from_uint64([value])
        for value in (-1, 1 << 64):
            with self.subTest(value=value), self.assertRaises(ValueError):
                g.cell_from_uint64(value)
        with self.assertRaises(ValueError):
            g.indices_from_uint64([-1])
        with self.assertRaises(ValueError):
            g.uint64_from_indices(["N", "X1"])

    def test_too_fine_for_the_fields(self):
        # N_side 7 fills the 29 row bits after level 9 (3 * 7**10 > 2**29).
        rdggs = RHEALPixDGGS(N_side=7, max_areal_resolution=1e-3)
        self.assertGreater(rdggs.max_resolution, 9)
        fits = rdggs.cell(["S", *[48] * 9])
        self.assertEqual(rdggs.cell_from_uint64(fits.uint64), fits)
        with self.assertRaises(ValueError):
            _ = rdggs.cell(["S", *[48] * 10]).uint64

    def test_wrappers(self):
        self.assertEqual(rhp_to_uint64("R553110"), 3458765673461713695)
        self.assertEqual(rhp_from_uint64(3458765673461713695), "R553110")
        self.assertIsNone(rhp_to_uint64("INVALID"))
        self.assertIsNone(rhp_from_uint64(1))
        # The identifier depends on N_side and the polar squares, not the
        # ellipsoid's central meridian.
        self.assertEqual(
            rhp_to_uint64("R553110", WGS84_003),
            rhp_to_uint64("R553110", WGS84_003_OGC),
        )


if __name__ == "__main__":
    unittest.main()
