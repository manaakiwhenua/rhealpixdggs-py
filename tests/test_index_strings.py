"""
Index strings, a face letter followed by one digit per resolution, are
defined for N_side 2 and 3 and read and written by one parser and one
formatter on the grid (issue #146). These tests check the round trips on
those grids and that every string entry point raises for any other grid.
"""

import json
import unittest

import numpy as np
from shapely.geometry import LineString, Polygon

from rhealpixdggs import export, rhp_wrappers
from rhealpixdggs.cell import CELLS0, Cell
from rhealpixdggs.conversion import compact_cells, compress_order_cells
from rhealpixdggs.dggs import WGS84_002, WGS84_003, RHEALPixDGGS
from rhealpixdggs.zoneset import ZoneSet

GRIDS = [WGS84_002, WGS84_003]
LARGE = [RHEALPixDGGS(N_side=4), RHEALPixDGGS(N_side=5)]
SQUARE = Polygon([(0.0, 0.0), (0.0, 20.0), (20.0, 20.0), (20.0, 0.0)])
LINE = LineString([(0.0, 0.0), (20.0, 20.0)])


class IndexStringsTestCase(unittest.TestCase):
    def test_round_trip_on_every_cell_to_resolution_2(self):
        for rdggs in GRIDS:
            with self.subTest(N_side=rdggs.N_side):
                self.assertTrue(rdggs.has_index_strings)
                for resolution in range(3):
                    for cell in rdggs.grid(resolution):
                        index = str(cell)
                        self.assertEqual(rdggs.format_index(cell.suid), index)
                        self.assertEqual(rdggs.parse_index(index), cell.suid)
                        self.assertEqual(rdggs.index_resolution(index), resolution)
                        self.assertEqual(rdggs.index_ancestor(index, resolution), index)
                        self.assertEqual(len(index), resolution + 1)
                        self.assertTrue(rhp_wrappers.rhp_is_valid(index, rdggs))

    def test_parse_rejects_malformed_strings(self):
        for rdggs in GRIDS:
            too_big = str(rdggs.N_side**2)
            bad = ["", "X", "n1", "N-1", "N 1", "N1.", "Na", "N²", "N" + too_big]
            for index in bad + ["1N"]:
                with self.subTest(N_side=rdggs.N_side, index=index):
                    self.assertIsNone(rdggs.parse_index(index))
                    self.assertFalse(rhp_wrappers.rhp_is_valid(index, rdggs))

    def test_parse_rejects_strings_beyond_max_resolution(self):
        for rdggs in GRIDS:
            deepest = "N" + "0" * rdggs.max_resolution
            too_deep = deepest + "0"
            with self.subTest(N_side=rdggs.N_side):
                suid = rdggs.parse_index(deepest)
                self.assertEqual(rdggs.cell(suid).suid, suid)
                self.assertIsNone(rdggs.parse_index(too_deep))
                self.assertFalse(rhp_wrappers.rhp_is_valid(too_deep, rdggs))
                with self.assertRaises(ValueError):
                    rdggs.cell(("N",) + (0,) * (rdggs.max_resolution + 1))
                valid, _, _, _ = rdggs._parse_indices([deepest, too_deep])
                self.assertEqual(valid.tolist(), [True, False])
                with self.assertRaises(ValueError):
                    rdggs.level_row_col_from_indices([too_deep])

    def test_codespace_and_uri(self):
        from rhealpixdggs.dggs import WGS84_003_OGC, WGS84_003_RADIANS

        self.assertEqual(WGS84_003_OGC.codespace, WGS84_003_OGC.uri)
        self.assertIsNone(WGS84_003.codespace)
        own = "https://example.org/dggrs/n2"
        mine = RHEALPixDGGS(N_side=2, uri=own)
        self.assertEqual(mine.uri, own)
        self.assertEqual(mine.codespace, own)
        self.assertEqual(mine.sub_zone_order, "mortonCurve")
        self.assertEqual(mine, WGS84_002)
        # The registered URI names the registered parameters, whatever a
        # constructor is told.
        ogc = WGS84_003_OGC.uri
        same = RHEALPixDGGS(WGS84_003_OGC.ellipsoid, uri=ogc)
        self.assertEqual(same.uri, ogc)
        named = RHEALPixDGGS(WGS84_003_OGC.ellipsoid, uri=own)
        self.assertEqual(named.uri, ogc)
        with self.assertRaises(ValueError):
            RHEALPixDGGS(N_side=2, uri=ogc)
        with self.assertRaises(ValueError):
            RHEALPixDGGS(N_side=2, uri=ogc.replace("https", "http"))
        self.assertIsNone(WGS84_003_RADIANS.codespace)

    def test_dggrs_definition_agrees_with_ogc(self):
        # The fields OGC's published rHEALPix definition fixes, at
        # https://www.opengis.net/def/dggrs/OGC/1.0/rHEALPix; the prose
        # and links are not compared.
        from rhealpixdggs.dggs import WGS84_003_OGC

        d = WGS84_003_OGC.to_dggrs_definition()
        self.assertEqual(d["uri"], WGS84_003_OGC.uri)
        definition = d["dggh"]["definition"]
        self.assertEqual(definition["spatialDimensions"], 2)
        self.assertEqual(definition["temporalDimensions"], 0)
        self.assertEqual(
            definition["projString"], "+proj=rhealpix +lon_0=50 +ellps=WGS84"
        )
        self.assertEqual(definition["basePolyhedron"], "cube")
        self.assertEqual(definition["refinementRatio"], 9)
        self.assertEqual(
            definition["refinementStrategy"],
            ["nestedChildCell", "nodeCentredChildCell"],
        )
        self.assertEqual(
            definition["constraints"], {"cellAxisAligned": True, "cellEqualSized": True}
        )
        self.assertEqual(definition["zoneTypes"], ["square"])
        parameters = d["dggh"]["parameters"]
        self.assertEqual(parameters["ellipsoid"], "[EPSG:7030]")
        self.assertAlmostEqual(
            parameters["orientation"]["latitude"], -41.937853910160136, places=12
        )
        self.assertEqual(parameters["orientation"]["longitude"], 50)
        self.assertEqual(parameters["orientation"]["azimuth"], 0)
        self.assertEqual(d["zirs"]["textZIRS"]["type"], "hierarchicalConcatenation")
        self.assertEqual(d["zirs"]["uint64ZIRS"]["type"], "ogc2DTMSHexLevelRowCol")
        self.assertEqual(d["subZoneOrder"]["type"], "scanline")
        json.dumps(d)

    def test_dggrs_definition_of_other_grids(self):
        from rhealpixdggs.ellipsoids import Ellipsoid

        d = WGS84_002.to_dggrs_definition()
        self.assertNotIn("uri", d)
        self.assertEqual(d["dggh"]["definition"]["refinementRatio"], 4)
        self.assertEqual(
            d["dggh"]["definition"]["refinementStrategy"], ["nestedChildCell"]
        )
        self.assertEqual(d["subZoneOrder"]["type"], "mortonCurve")
        own = "https://example.org/dggrs/n4"
        n4 = RHEALPixDGGS(N_side=4, north_square=1, south_square=2, uri=own)
        d = n4.to_dggrs_definition(title="N4", description="mine")
        self.assertEqual((d["title"], d["description"], d["uri"]), ("N4", "mine", own))
        self.assertEqual(
            d["dggh"]["definition"]["projString"],
            "+proj=rhealpix +lon_0=0 +ellps=WGS84 +north_square=1 +south_square=2",
        )
        self.assertIn("(N, 10, 15)", d["zirs"]["textZIRS"]["description"])
        radians = RHEALPixDGGS(
            ellipsoid=Ellipsoid(R=6371000, lon_0=np.pi / 2, radians=True)
        )
        d = radians.to_dggrs_definition()
        self.assertEqual(
            d["dggh"]["definition"]["projString"],
            "+proj=rhealpix +lon_0=90 +R=6371000.0",
        )
        self.assertEqual(d["dggh"]["parameters"]["orientation"]["longitude"], 90)
        self.assertAlmostEqual(
            d["dggh"]["parameters"]["orientation"]["latitude"], -41.8103149, places=6
        )

    def test_digit_range_follows_n_side(self):
        self.assertEqual(WGS84_003.parse_index("P48"), ("P", 4, 8))
        self.assertIsNone(WGS84_002.parse_index("P48"))
        self.assertEqual(WGS84_002.parse_index("P31"), ("P", 3, 1))

    def test_ancestors_and_children(self):
        for rdggs in GRIDS:
            with self.subTest(N_side=rdggs.N_side):
                cell = rdggs.cell(("P", 1, 0, 2))
                index = str(cell)
                ancestors = [rdggs.index_ancestor(index, r) for r in range(4)]
                self.assertEqual(ancestors, ["P", "P1", "P10", "P102"])
                self.assertEqual(rdggs.index_ancestor(index, 9), index)
                children = rdggs.index_children("P1")
                expected = [str(c) for c in rdggs.cell(("P", 1)).subcells()]
                self.assertEqual(children, expected)
                self.assertEqual(len(children), rdggs.N_side**2)
                self.assertEqual(rhp_wrappers.rhp_to_parent(index, dggs=rdggs), "P10")
                self.assertEqual(
                    rhp_wrappers.rhp_to_parent(index, res=1, dggs=rdggs), "P1"
                )
                self.assertEqual(rhp_wrappers.rhp_get_resolution(index, rdggs), 3)
                self.assertEqual(rhp_wrappers.rhp_get_base_cell(index, rdggs), "P")

    def test_center_child_uses_the_formatter(self):
        self.assertEqual(rhp_wrappers.rhp_to_center_child("P1", dggs=WGS84_003), "P14")
        self.assertEqual(
            rhp_wrappers.rhp_to_center_child("P1", res=3, dggs=WGS84_003), "P144"
        )

    def test_compaction_groups_by_parent(self):
        children = WGS84_002.index_children("P1")
        self.assertEqual(compact_cells(children, N_side=2), {"P1"})
        self.assertEqual(compress_order_cells(children, N_side=2), ["P1"])
        self.assertEqual(compact_cells(children[:-1], N_side=2), set(children[:-1]))

    def test_invalid_strings_raise_where_they_used_to_be_grouped(self):
        for call in (
            lambda: compact_cells(["N0", "ba"], N_side=3),
            lambda: compress_order_cells(["N0", "ba"], N_side=3),
            lambda: ZoneSet(WGS84_003)._cell("bad"),
        ):
            with self.assertRaisesRegex(ValueError, "invalid cell index"):
                call()

    def test_large_n_side_has_no_index_strings(self):
        for rdggs in LARGE:
            with self.subTest(N_side=rdggs.N_side):
                self.assertFalse(rdggs.has_index_strings)
                cell = rdggs.cell(("N", 10, 15))
                self.assertEqual(str(cell), "(N, 10, 15)")
                self.assertEqual(rdggs.cell(cell.suid[:-1]), rdggs.cell(("N", 10)))
                self.assertEqual(len(list(cell.subcells())), rdggs.N_side**2)

    def test_every_string_entry_point_raises_for_large_n_side(self):
        for rdggs in LARGE:
            self._check_every_string_entry_point_raises(rdggs)

    def _check_every_string_entry_point_raises(self, rdggs):
        cell = rdggs.cell(("N", 10, 15))
        n = rdggs.N_side
        calls = {
            "parse_index": lambda: rdggs.parse_index("N1015"),
            "format_index": lambda: rdggs.format_index(cell.suid),
            "index_resolution": lambda: rdggs.index_resolution("N1015"),
            "index_ancestor": lambda: rdggs.index_ancestor("N1015", 1),
            "index_children": lambda: rdggs.index_children("N"),
            "rhp_is_valid": lambda: rhp_wrappers.rhp_is_valid("N1015", rdggs),
            "geo_to_rhp": lambda: rhp_wrappers.geo_to_rhp(10.0, 10.0, 1, dggs=rdggs),
            "rhp_to_geo": lambda: rhp_wrappers.rhp_to_geo("N1", dggs=rdggs),
            "rhp_to_parent": lambda: rhp_wrappers.rhp_to_parent("N1", dggs=rdggs),
            "rhp_to_center_child": lambda: rhp_wrappers.rhp_to_center_child(
                "N1", dggs=rdggs
            ),
            "rhp_to_geo_boundary": lambda: rhp_wrappers.rhp_to_geo_boundary(
                "N1", dggs=rdggs
            ),
            "rhp_get_resolution": lambda: rhp_wrappers.rhp_get_resolution("N1", rdggs),
            "rhp_get_base_cell": lambda: rhp_wrappers.rhp_get_base_cell("N1", rdggs),
            "cell_area": lambda: rhp_wrappers.cell_area("N1", dggs=rdggs),
            "cell_ring": lambda: rhp_wrappers.cell_ring("N1", dggs=rdggs),
            "k_ring": lambda: rhp_wrappers.k_ring("N1", dggs=rdggs),
            "polyfill": lambda: rhp_wrappers.polyfill(SQUARE, 1, dggs=rdggs),
            "polyfill_array": lambda: rhp_wrappers.polyfill_array(
                SQUARE, 1, dggs=rdggs
            ),
            "linetrace": lambda: rhp_wrappers.linetrace(LINE, 1, dggs=rdggs),
            "shapes": lambda: rdggs.shapes(["N1"]),
            "nuclei": lambda: rdggs.nuclei(["N1"]),
            "centroids": lambda: rdggs.centroids(["N1"]),
            "boundary_array": lambda: rdggs.boundary_array(["N1"]),
            "cells_from_points": lambda: rdggs.cells_from_points(
                np.array([10.0]), np.array([10.0]), 1
            ),
            "cells_in_box": lambda: rdggs.cells_in_box(1, (0.0, 10.0), (10.0, 0.0)),
            "ZoneSet": lambda: ZoneSet(rdggs, [cell]),
            "cell_table": lambda: export.cell_table([cell], rdggs),
            "geometries": lambda: export.geometries([cell], rdggs),
            "to_geojson": lambda: export.to_geojson([cell], rdggs),
            "to_csv": lambda: export.to_csv([cell], rdggs),
            "compact_cells": lambda: compact_cells(["N1"], N_side=n),
            "compress_order_cells": lambda: compress_order_cells(["N1"], N_side=n),
        }
        for name, call in calls.items():
            with (
                self.subTest(N_side=n, entry_point=name),
                self.assertRaisesRegex(ValueError, "N_side 2 and 3"),
            ):
                call()

    def test_face_letters_are_the_resolution_0_index_strings(self):
        for rdggs in GRIDS:
            with self.subTest(N_side=rdggs.N_side):
                self.assertEqual([str(c) for c in rdggs.grid(0)], list(CELLS0))
                for face in CELLS0:
                    self.assertEqual(rdggs.parse_index(face), (face,))
                self.assertIsInstance(Cell(rdggs, ["N"]), Cell)


if __name__ == "__main__":
    unittest.main()
