"""
Tests for rhealpixdggs.export: GeoJSON and delimited-text export of cells.
"""

import csv
import io
import json
import unittest
import xml.etree.ElementTree as ET
from importlib import resources

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
            for index, geom, crosses in zip(indices, geoms, crossing, strict=True):
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
                    [i for i, c in zip(indices, crossing, strict=True) if c],
                    expected_crossing[(id(rdggs), resolution)],
                )
            # Caps are polygons reaching the pole along both sides of the
            # antimeridian.
            for index, geom, cap in zip(indices, geoms, caps, strict=True):
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
        for a, b in zip(degrees["features"], radians["features"], strict=True):
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
        for row, feature in zip(rows[1:], fc["features"], strict=True):
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
        for geom, feature in zip(geoms, fc["features"], strict=True):
            self.assertTrue(shapely.equals_exact(geom, shape(feature["geometry"]), 0))


if __name__ == "__main__":
    unittest.main()


GML = "{http://www.opengis.net/gml/3.2}"
RHP = f"{{{export.GML_NAMESPACE}}}"


def zones(text):
    root = ET.fromstring(text)
    return root, root.findall(f"{RHP}member/{RHP}Zone")


class GMLTestCase(unittest.TestCase):
    def test_document_structure(self):
        cells = ["N", "P0", "N46", "S8"]
        root, members = zones(
            export.to_gml(cells, ROTATED, codespace="https://example.org/g")
        )
        self.assertEqual(root.tag, f"{RHP}ZoneCollection")
        self.assertEqual(root.get(f"{GML}id"), "zones")
        self.assertEqual(
            root.get("{http://www.w3.org/2001/XMLSchema-instance}schemaLocation"),
            f"{export.GML_NAMESPACE} {export.GML_SCHEMA_LOCATION}",
        )
        self.assertEqual(root.find(f"{RHP}dggrs").text, "https://example.org/g")
        self.assertEqual([z.get(f"{GML}id") for z in members], cells)
        table = export.cell_table(cells, ROTATED)
        for k, zone in enumerate(members):
            identifier = zone.find(f"{GML}identifier")
            self.assertEqual(
                (identifier.text, identifier.get("codeSpace")),
                (cells[k], "https://example.org/g"),
            )
            self.assertEqual(
                zone.find(f"{RHP}resolution").text, str(table.resolution[k])
            )
            self.assertEqual(zone.find(f"{RHP}shape").text, table.shape[k])
            area = zone.find(f"{RHP}area")
            self.assertEqual((float(area.text), area.get("uom")), (table.area[k], "m2"))
            self.assertEqual(
                [b.text for b in zone.findall(f"{RHP}boundaryType")],
                list(export.boundary_types(cells[k], ROTATED)),
            )
            # GML writes EPSG:4326 latitude first.
            pos = zone.find(f"{RHP}position/{GML}Point/{GML}pos")
            self.assertEqual(
                [float(v) for v in pos.text.split()],
                [table.centroid[k, 1], table.centroid[k, 0]],
            )
            point = zone.find(f"{RHP}position/{GML}Point")
            self.assertEqual(point.get("srsName"), "urn:ogc:def:crs:EPSG::4326")
            nucleus = zone.find(f"{RHP}nucleus/{GML}Point/{GML}pos")
            self.assertEqual(
                [float(v) for v in nucleus.text.split()],
                [table.nucleus[k, 1], table.nucleus[k, 0]],
            )

    def test_geometry_matches_geojson(self):
        cells = ["N", "P0", "N46", "S8", "N4", "O"]
        _, members = zones(export.to_gml(cells, ROTATED, n=4, codespace="x"))
        features = export.to_geojson(cells, ROTATED, n=4)["features"]
        for zone, feature in zip(members, features, strict=True):
            geometry = zone.find(f"{RHP}geometry")[0]
            if feature["geometry"]["type"] == "Polygon":
                self.assertEqual(geometry.tag, f"{GML}Polygon")
                polygons = [geometry]
            else:
                self.assertEqual(geometry.tag, f"{GML}MultiSurface")
                polygons = geometry.findall(f"{GML}surfaceMember/{GML}Polygon")
            self.assertEqual(geometry.get("srsName"), "urn:ogc:def:crs:EPSG::4326")
            self.assertEqual(geometry.get(f"{GML}id"), f"{feature['id']}.geometry")
            rings = feature["geometry"]["coordinates"]
            if feature["geometry"]["type"] == "Polygon":
                rings = [rings]
            for polygon, ring in zip(polygons, rings, strict=True):
                values = [
                    float(v)
                    for v in polygon.find(
                        f"{GML}exterior/{GML}LinearRing/{GML}posList"
                    ).text.split()
                ]
                expected = [v for lon, lat in ring[0] for v in (lat, lon)]
                self.assertEqual(values, expected)

    def test_boundary_types_match_the_edges(self):
        # Sample each planar edge and project it: a smallCircle keeps its
        # latitude, a greatEllipse its longitude, a projectedLine neither.
        for rdggs in (WGS84_003, WGS84_122, RHEALPixDGGS(N_side=4)):
            for resolution in range(3):
                for cell in rdggs.grid(resolution):
                    types = export.boundary_types(cell, rdggs)
                    (x0, x1), (y0, y1) = cell.xy_range()
                    t = np.linspace(0.1, 0.9, 5)
                    edges = [
                        [(x0 + s * (x1 - x0), y1) for s in t],
                        [(x1, y1 - s * (y1 - y0)) for s in t],
                        [(x1 - s * (x1 - x0), y0) for s in t],
                        [(x0, y0 + s * (y1 - y0)) for s in t],
                    ]
                    if len(types) == 1:
                        self.assertEqual(cell.ellipsoidal_shape, "cap")
                        edges = [[point for edge in edges for point in edge]]
                    for kind, points in zip(types, edges, strict=True):
                        ll = np.array(
                            [rdggs.rhealpix(x, y, inverse=True) for x, y in points]
                        )
                        same_lat = np.ptp(ll[:, 1]) < 1e-9
                        turns = (np.diff(ll[:, 0]) + 180) % 360 - 180
                        same_lon = bool(np.all(np.abs(turns) < 1e-9))
                        expected = {
                            "smallCircle": same_lat,
                            "greatEllipse": same_lon,
                            "projectedLine": not same_lat and not same_lon,
                        }[kind]
                        self.assertTrue(expected, msg=f"{rdggs.N_side} {cell} {kind}")

    def test_codespace_and_other_ellipsoids(self):
        from rhealpixdggs.dggs import WGS84_003_OGC

        _, members = zones(export.to_gml(["P0"], WGS84_003_OGC))
        self.assertEqual(
            members[0].find(f"{GML}identifier").get("codeSpace"), WGS84_003_OGC.uri
        )
        with self.assertRaises(ValueError):
            export.to_gml(["P0"], WGS84_003)
        with self.assertRaises(ValueError):
            export.to_gml(["P0", "nope"], WGS84_003, codespace="x")
        _, members = zones(
            export.to_gml(["P0"], RHEALPixDGGS(N_side=2, uri="https://example.org/n2"))
        )
        self.assertEqual(
            members[0].find(f"{GML}identifier").get("codeSpace"),
            "https://example.org/n2",
        )
        # Not WGS84: geodetic coordinates, no srsName claimed.
        sphere = RHEALPixDGGS(
            ellipsoid=Ellipsoid(R=6371000), N_side=3, uri="https://example.org/s"
        )
        _, members = zones(export.to_gml(["P0"], sphere))
        self.assertIsNone(members[0].find(f"{RHP}geometry/{GML}Polygon").get("srsName"))
        # Radians grids are written in degrees, like the GeoJSON.
        _, members = zones(export.to_gml(["P0"], WGS84_003_RADIANS, codespace="x"))
        self.assertEqual(
            members[0].find(f"{RHP}nucleus/{GML}Point/{GML}pos").text.split()[1],
            repr(float(export.cell_table(["P0"], WGS84_003).nucleus[0, 0])),
        )

    def test_planar_document(self):
        from rhealpixdggs.dggs import WGS84_003_OGC

        cells = ["P0", "N46", "N", "S8"]
        root, members = zones(export.to_gml(cells, WGS84_003_OGC, planar=True))
        crs = root.find(f"{RHP}crs/{GML}ProjectedCRS")
        self.assertEqual(crs.get(f"{GML}id"), "rhealpix")
        values = [v.text for v in crs.iter(f"{GML}value")]
        self.assertEqual(values, ["50", "0", "0"])
        self.assertEqual(
            crs.find(f"{GML}baseGeodeticCRS").get("{http://www.w3.org/1999/xlink}href"),
            "urn:ogc:def:crs:EPSG::4326",
        )
        for index, zone in zip(cells, members, strict=True):
            cell = WGS84_003_OGC.cell(WGS84_003_OGC.parse_index(index))
            polygon = zone.find(f"{RHP}geometry/{GML}Polygon")
            self.assertEqual(polygon.get("srsName"), "#rhealpix")
            values = [
                float(v)
                for v in polygon.find(
                    f"{GML}exterior/{GML}LinearRing/{GML}posList"
                ).text.split()
            ]
            ring = list(zip(values[0::2], values[1::2], strict=True))
            self.assertEqual(ring[0], ring[-1])
            self.assertEqual(sorted(ring[:-1]), sorted(cell.vertices(plane=True)))
            nucleus = zone.find(f"{RHP}nucleus/{GML}Point/{GML}pos").text.split()
            self.assertEqual(
                [float(v) for v in nucleus], list(cell.nucleus(plane=True))
            )
        with self.assertRaises(ValueError):
            export.to_gml(
                ["P0"], RHEALPixDGGS(ellipsoid=Ellipsoid(R=1), uri="x"), planar=True
            )

    def test_validates_against_the_schema(self):
        # The schema imports GML 3.2.1 from schemas.opengis.net; skip when
        # xmlschema or the network is missing.
        try:
            import xmlschema
        except ImportError:
            self.skipTest("xmlschema is not installed")
        path = resources.files("rhealpixdggs") / "schemas" / "zone.xsd"
        try:
            schema = xmlschema.XMLSchema(str(path))
        except Exception as error:  # noqa: BLE001 - any fetch failure
            self.skipTest(f"GML schemas unavailable: {error}")
        cells = ["N", "P0", "N46", "S8", "N4", "O", "N0"]
        schema.validate(
            export.to_gml(cells, ROTATED, n=3, codespace="https://example.org/g")
        )
        schema.validate(
            export.to_gml(["N0", "S3"], WGS84_122, codespace="https://example.org/h")
        )
        schema.validate(
            export.to_gml(
                cells, ROTATED, codespace="https://example.org/g", planar=True
            )
        )
