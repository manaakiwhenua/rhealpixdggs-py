"""
The CRS a grid publishes (OGC Topic 21 A.7): the rHEALPix projected CRS
on WGS 84 or the grid's ellipsoid, which PROJ's own rhealpix must apply
to the same planar coordinates the package computes.
"""

import json
import random
import unittest
import warnings
from math import radians

import pyproj

from rhealpixdggs.dggs import (
    UNIT_003,
    WGS84_002,
    WGS84_003,
    WGS84_003_OGC,
    WGS84_003_RADIANS,
    RHEALPixDGGS,
)
from rhealpixdggs.ellipsoids import Ellipsoid

GRIDS = {
    "ogc": WGS84_003_OGC,
    "lon_0 0": WGS84_003,
    "radians": WGS84_003_RADIANS,
    "N_side 2": WGS84_002,
    "unit sphere": UNIT_003,
    "GRS80, squares 1 and 2": RHEALPixDGGS(
        ellipsoid=Ellipsoid(a=6378137, f=1 / 298.257222101, lon_0=10),
        north_square=1,
        south_square=2,
    ),
}


class CRSTestCase(unittest.TestCase):
    def test_proj_applies_the_definition(self):
        # The package's plane is the CRS: PROJ's rhealpix, given the
        # published definition, lands within micrometres of the package's
        # projection everywhere.
        rng = random.Random(0)
        for name, rdggs in GRIDS.items():
            with self.subTest(grid=name):
                crs = rdggs.crs
                transformer = pyproj.Transformer.from_crs(
                    crs.geodetic_crs, crs, always_xy=True
                )
                scale = rdggs.ellipsoid.R_A
                for _ in range(300):
                    lon, lat = rng.uniform(-180, 180), rng.uniform(-89.9, 89.9)
                    if rdggs.ellipsoid.radians:
                        x1, y1 = rdggs.rhealpix(radians(lon), radians(lat))
                    else:
                        x1, y1 = rdggs.rhealpix(lon, lat)
                    x2, y2 = transformer.transform(lon, lat)
                    self.assertAlmostEqual(x1 / scale, x2 / scale, places=12)
                    self.assertAlmostEqual(y1 / scale, y2 / scale, places=12)

    def test_wgs84_grids_are_on_epsg_4326(self):
        for rdggs in (WGS84_003_OGC, WGS84_003, WGS84_003_RADIANS, WGS84_002):
            crs = rdggs.crs
            self.assertEqual(crs.geodetic_crs.to_authority(), ("EPSG", "4326"))
            self.assertTrue(crs.is_projected)
            self.assertEqual(crs.axis_info[0].unit_name, "metre")
            self.assertIn("ENSEMBLEACCURACY[2.0]", crs.to_wkt())
            self.assertEqual(crs.coordinate_operation.method_name, "PROJ rhealpix")

    def test_parameters_name_the_grid(self):
        crs = GRIDS["GRS80, squares 1 and 2"].crs
        params = {p.name: p.value for p in crs.coordinate_operation.params}
        self.assertEqual(params, {"lon_0": 10, "north_square": 1, "south_square": 2})
        self.assertEqual(crs.ellipsoid.inverse_flattening, 298.257222101)
        self.assertNotEqual(crs.geodetic_crs.to_authority(), ("EPSG", "4326"))
        # The radians grid still has lon_0 in degrees in the CRS.
        params = {
            p.name: p.value for p in WGS84_003_RADIANS.crs.coordinate_operation.params
        }
        self.assertEqual(params["lon_0"], 0)
        self.assertEqual(UNIT_003.crs.ellipsoid.semi_major_metre, 1)

    def test_definitions_round_trip(self):
        for name, rdggs in GRIDS.items():
            with self.subTest(grid=name):
                crs = rdggs.crs
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    self.assertEqual(pyproj.CRS.from_wkt(crs.to_wkt()), crs)
                    self.assertEqual(pyproj.CRS.from_json(crs.to_json()), crs)
                definition = rdggs.to_dggrs_definition()
                projjson = definition["dggh"]["definition"]["crs"]
                self.assertEqual(projjson["type"], "ProjectedCRS")
                self.assertEqual(pyproj.CRS.from_json_dict(projjson), crs)
                json.dumps(definition)
