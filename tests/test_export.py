"""
Tests for rhealpixdggs.export: GeoJSON and delimited-text export of cells.
"""

import csv
import io
import json
import unittest

import numpy as np
import shapely
from numpy.testing import assert_allclose, assert_array_equal
from shapely.geometry import shape

from rhealpixdggs import export
from rhealpixdggs.dggs import WGS84_003, WGS84_003_RADIANS, RHEALPixDGGS
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID, Ellipsoid

WGS84_122 = RHEALPixDGGS(
    ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=2, N_side=2
)
ROTATED = RHEALPixDGGS(ellipsoid=Ellipsoid(lon_0=50), N_side=3)
PROPERTIES = [
    "index",
    "resolution",
    "shape",
    "area",
    "centroid_lon",
    "centroid_lat",
    "nucleus_lon",
    "nucleus_lat",
]


def grid_indices(rdggs, resolution):
    return [str(c) for c in rdggs.grid(resolution)]


class ExportTestCase(unittest.TestCase):
    def test_feature_collection_structure(self):
        indices = ["P0", "N0", "N1", "N44", "N46", "S00", "N"]
        fc = export.to_geojson(indices)
        self.assertEqual(fc["type"], "FeatureCollection")
        self.assertEqual(set(fc), {"type", "features"})
        self.assertEqual([f["id"] for f in fc["features"]], indices)
        for feature in fc["features"]:
            self.assertEqual(feature["type"], "Feature")
            self.assertEqual(set(feature), {"type", "id", "geometry", "properties"})
            self.assertEqual(list(feature["properties"]), PROPERTIES)
            props = feature["properties"]
            self.assertIsInstance(props["index"], str)
            self.assertIsInstance(props["resolution"], int)
            self.assertIsInstance(props["shape"], str)
            for key in PROPERTIES[3:]:
                self.assertIs(type(props[key]), float, key)
            self.assertIn(feature["geometry"]["type"], ("Polygon", "MultiPolygon"))
            coords = np.asarray(shapely.get_coordinates(shape(feature["geometry"])))
            self.assertTrue((np.abs(coords[:, 0]) <= 180).all())
            self.assertTrue((np.abs(coords[:, 1]) <= 90).all())
        # Plain JSON in and out: builtin types only.
        self.assertEqual(json.loads(json.dumps(fc)), fc)

    def test_geometries_follow_rfc7946_over_whole_grids(self):
        grids = [
            (WGS84_003, 0),
            (WGS84_003, 1),
            (WGS84_003, 2),
            (WGS84_122, 1),
            (WGS84_122, 2),
            (ROTATED, 1),
            (ROTATED, 2),
        ]
        for rdggs, resolution in grids:
            indices = grid_indices(rdggs, resolution)
            geoms = export.geometries(indices, dggs=rdggs)
            self.assertEqual(geoms.shape, (len(indices),))
            self.assertTrue(shapely.is_valid(geoms).all(), (rdggs, resolution))
            polygons = shapely.get_parts(geoms)
            self.assertTrue(
                all(shapely.is_ccw(p.exterior) for p in polygons), (rdggs, resolution)
            )
            coords = shapely.get_coordinates(geoms)
            self.assertTrue(
                (coords[:, 0] >= -180).all() and (coords[:, 0] <= 180).all()
            )
            self.assertTrue((coords[:, 1] >= -90).all() and (coords[:, 1] <= 90).all())
            # The features tile the longitude-latitude rectangle exactly
            # (planar areas), so nothing is lost or doubled at the cuts.
            self.assertAlmostEqual(shapely.area(geoms).sum(), 360 * 180, places=6)
            self.assertAlmostEqual(shapely.union_all(geoms).area, 360 * 180, places=6)
            # Cells crossing the antimeridian are MultiPolygons with an
            # eastern and a western part; everything else is a Polygon.
            # A ring whose longitudes jump by more than half a turn wraps;
            # unwrapped eastwards it genuinely crosses if it reaches past
            # 180 (a cell whose east edge is the antimeridian does not).
            # A vertex at a pole has an arbitrary longitude and is ignored.
            rings = rdggs.boundary_array(indices, n=2)
            lon = np.where(np.abs(rings[:, :, 1]) > 90 - 1e-9, np.nan, rings[:, :, 0])
            wraps = np.nanmax(lon, axis=1) - np.nanmin(lon, axis=1) > 180
            unwrapped = np.where(lon < 0, lon + 360, lon)
            caps = rdggs.shapes(indices) == "cap"
            # Genuine crossing needs extent on both sides of 180: a cell
            # whose west edge is reported as +180 lies wholly to the west.
            crossing = (
                ~caps
                & wraps
                & (np.nanmax(unwrapped, axis=1) > 180 + 1e-9)
                & (np.nanmin(unwrapped, axis=1) < 180 - 1e-9)
            )
            for index, geom, crosses in zip(indices, geoms, crossing):
                if crosses:
                    self.assertEqual(geom.geom_type, "MultiPolygon", index)
                    self.assertEqual(len(geom.geoms), 2, index)
                    east, west = geom.geoms
                    self.assertAlmostEqual(east.bounds[2], 180.0, places=9)
                    self.assertAlmostEqual(west.bounds[0], -180.0, places=9)
                else:
                    self.assertEqual(geom.geom_type, "Polygon", index)
            expected_crossing = {
                (id(WGS84_003), 1): ["N6", "S0"],
                (id(ROTATED), 1): ["N3", "R1", "R4", "R7", "S3"],
            }
            if (id(rdggs), resolution) in expected_crossing:
                self.assertEqual(
                    [i for i, c in zip(indices, crossing) if c],
                    expected_crossing[(id(rdggs), resolution)],
                )
            # Caps are polygons reaching the pole along both sides of the
            # antimeridian.
            for index, geom, cap in zip(indices, geoms, caps):
                if cap:
                    pole = 90.0 if index.startswith("N") else -90.0
                    points = {tuple(c) for c in shapely.get_coordinates(geom)}
                    self.assertIn((180.0, pole), points, index)
                    self.assertIn((-180.0, pole), points, index)

    def test_sampling(self):
        rdggs = WGS84_003
        for n in (2, 3, 6):
            quad, dart, skew = export.geometries(["P4", "N0", "N1"], n=n)
            self.assertEqual(len(quad.exterior.coords), 5)
            self.assertEqual(len(dart.exterior.coords), 4 * n - 3)
            self.assertEqual(len(skew.exterior.coords), 4 * n - 3)
            # A polar ring is the cell's boundary, counter-clockwise, still
            # starting at the north-west vertex.
            boundary = np.array(rdggs.cell(("N", 1)).boundary(n=n, plane=False))
            expected = np.roll(boundary[::-1], 1, axis=0)
            assert_allclose(
                np.array(skew.exterior.coords)[:-1], expected, rtol=0, atol=1e-9
            )
        # n below 2 clamps to 2.
        (dart,) = export.geometries(["N0"], n=1)
        self.assertEqual(len(dart.exterior.coords), 5)

    def test_properties_match_the_cell_api(self):
        rdggs = WGS84_003
        indices = grid_indices(rdggs, 1) + ["P44", "N44", "S08"]
        table = export.cell_table(indices, dggs=rdggs)
        assert_array_equal(table.indices, indices)
        assert_array_equal(table.resolution, [len(i) - 1 for i in indices])
        assert_array_equal(table.shape, rdggs.shapes(indices))
        assert_array_equal(
            table.area, [rdggs.cell_area(len(i) - 1, plane=False) for i in indices]
        )
        assert_array_equal(table.centroid, rdggs.centroids(indices))
        assert_array_equal(table.nucleus, rdggs.nuclei(indices))
        fc = export.to_geojson(indices, dggs=rdggs)
        for k, feature in enumerate(fc["features"]):
            props = feature["properties"]
            cell = rdggs.cell([indices[k][0]] + [int(d) for d in indices[k][1:]])
            self.assertEqual(props["shape"], cell.ellipsoidal_shape)
            self.assertEqual(props["area"], cell.area(plane=False))
            # The array forms differ from the per-cell methods in the last
            # bits (array sums rather than fsum), as their docstrings say.
            assert_allclose(
                (props["centroid_lon"], props["centroid_lat"]),
                cell.centroid(plane=False),
                rtol=0,
                atol=1e-9,
            )
            assert_allclose(
                (props["nucleus_lon"], props["nucleus_lat"]),
                cell.nucleus(plane=False),
                rtol=0,
                atol=1e-9,
            )

    def test_grids_in_radians_and_other_ellipsoids(self):
        indices = grid_indices(WGS84_003, 1)
        degrees = export.to_geojson(indices, dggs=WGS84_003)
        radians = export.to_geojson(indices, dggs=WGS84_003_RADIANS)
        for a, b in zip(degrees["features"], radians["features"]):
            assert_allclose(
                shapely.get_coordinates(shape(a["geometry"])),
                shapely.get_coordinates(shape(b["geometry"])),
                rtol=0,
                atol=1e-9,
            )
            for key in PROPERTIES[3:]:
                self.assertAlmostEqual(a["properties"][key], b["properties"][key])
        sphere = RHEALPixDGGS(ellipsoid=Ellipsoid(R=1), N_side=3)
        fc = export.to_geojson(["P0"], dggs=sphere)
        self.assertEqual(
            fc["features"][0]["properties"]["area"], sphere.cell_area(1, plane=False)
        )

    def test_inputs(self):
        rdggs = WGS84_003
        cells = list(rdggs.cell(("P", 0)).subcells())
        as_cells = export.to_geojson(cells, dggs=rdggs)
        as_strings = export.to_geojson([str(c) for c in cells], dggs=rdggs)
        as_array = export.to_geojson(np.array([str(c) for c in cells]), dggs=rdggs)
        self.assertEqual(as_cells, as_strings)
        self.assertEqual(as_array, as_strings)
        # Duplicates and order are preserved.
        fc = export.to_geojson(["P1", "P0", "P1"], dggs=rdggs)
        self.assertEqual([f["id"] for f in fc["features"]], ["P1", "P0", "P1"])
        self.assertEqual(
            export.to_geojson([]), {"type": "FeatureCollection", "features": []}
        )
        self.assertEqual(export.to_csv([]), ",".join(PROPERTIES) + "\n")
        for call in (
            lambda: export.to_geojson(["P0", "X9", ""]),
            lambda: export.to_csv(["P0", "X9"]),
            lambda: export.geometries(["X9"]),
            lambda: export.cell_table(["P0", "X9"]),
        ):
            with self.assertRaises(ValueError) as raised:
                call()
            self.assertIn("'X9'", str(raised.exception))

    def test_csv(self):
        indices = ["P0", "N0", "N", "S44"]
        text = export.to_csv(indices)
        self.assertTrue(text.endswith("\n") and not text.endswith("\n\n"))
        rows = list(csv.reader(io.StringIO(text)))
        self.assertEqual(rows[0], PROPERTIES)
        self.assertEqual(len(rows), len(indices) + 1)
        fc = export.to_geojson(indices)
        for row, feature in zip(rows[1:], fc["features"]):
            props = feature["properties"]
            self.assertEqual(row[0], props["index"])
            self.assertEqual(int(row[1]), props["resolution"])
            self.assertEqual(row[2], props["shape"])
            for k, key in enumerate(PROPERTIES[3:], start=3):
                self.assertEqual(float(row[k]), props[key])
        tabbed = export.to_csv(indices, delimiter="\t")
        self.assertEqual(tabbed.splitlines()[0], "\t".join(PROPERTIES))
        self.assertEqual(len(tabbed.splitlines()), len(indices) + 1)

    def test_geometries_agree_with_geojson(self):
        indices = grid_indices(WGS84_003, 1)
        geoms = export.geometries(indices)
        fc = export.to_geojson(indices)
        for geom, feature in zip(geoms, fc["features"]):
            self.assertTrue(shapely.equals_exact(geom, shape(feature["geometry"]), 0))


if __name__ == "__main__":
    unittest.main()
