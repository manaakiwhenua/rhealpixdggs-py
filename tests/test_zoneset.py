"""
Tests for rhealpixdggs.zoneset and the hierarchy predicates on Cell: the
remaining operations of OGC Topic 21 v2.0 Table 53.
"""

import random
import unittest

import numpy as np
import shapely

from rhealpixdggs.cell import CELLS0, Cell
from rhealpixdggs.dggs import WGS84_003, RHEALPixDGGS
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from rhealpixdggs.rhp_wrappers import k_ring
from rhealpixdggs.zoneset import BoundaryType, ZoneSet

N, O, P, Q, R, S = CELLS0
WGS84_122 = RHEALPixDGGS(
    ellipsoid=WGS84_ELLIPSOID, north_square=1, south_square=2, N_side=2
)


def cell(index, rdggs=WGS84_003):
    return rdggs.cell([index[0]] + [int(d) for d in index[1:]])


def planar_geometry(zones):
    boxes = []
    for c in zones:
        (x0, x1), (y0, y1) = c.xy_range()
        boxes.append(shapely.box(x0, y0, x1, y1))
    return shapely.union_all(boxes) if boxes else shapely.Polygon()


class HierarchyPredicateTestCase(unittest.TestCase):
    def test_parent_of_and_child_of(self):
        p4, p40, p400 = cell("P4"), cell("P40"), cell("P400")
        self.assertTrue(p4.parent_of(p40))
        self.assertTrue(p40.child_of(p4))
        self.assertFalse(p4.parent_of(p400))  # a grandchild, not a child
        self.assertFalse(p400.child_of(p4))
        self.assertFalse(p40.parent_of(p4))
        self.assertFalse(p4.parent_of(p4))
        self.assertFalse(p4.parent_of(cell("P")))
        # inherit_id makes no difference in a nested hierarchy.
        self.assertTrue(p4.parent_of(p40, inherit_id=True))
        self.assertTrue(p40.child_of(p4, inherit_id=True))
        for method in (p4.parent_of, p4.child_of, p4.sibling_of):
            with self.assertRaises(ValueError):
                method(WGS84_003.cell())
            with self.assertRaises(ValueError):
                method(Cell(RHEALPixDGGS(N_side=4), (P, 0)))

    def test_sibling_of(self):
        p40, p41, p44, p32, p20 = (cell(i) for i in ("P40", "P41", "P44", "P32", "P20"))
        # Same parent: siblings either way. P32 sits left of P40 in the
        # neighbouring parent P3: a sibling only when inheritance is
        # ignored. P20 is not adjacent to P40 at all.
        self.assertTrue(p40.sibling_of(p41))
        self.assertTrue(p40.sibling_of(p41, inherit_id=True))
        self.assertTrue(p40.sibling_of(p44, inherit_id=True))  # corner-adjacent
        self.assertTrue(p40.sibling_of(p32))
        self.assertFalse(p40.sibling_of(p32, inherit_id=True))
        self.assertFalse(p40.sibling_of(p20))
        self.assertFalse(p40.sibling_of(p20, inherit_id=True))
        # A cell is its own sibling, as the spec's sibling() sets include
        # the cell; different resolutions are never siblings.
        self.assertTrue(p40.sibling_of(p40))
        self.assertFalse(p40.sibling_of(cell("P4")))
        self.assertFalse(p40.sibling_of(cell("P411")))


class ZoneSetTestCase(unittest.TestCase):
    def test_construction(self):
        rdggs = WGS84_003
        zs = ZoneSet(rdggs, ["P40", cell("P41"), np.str_("P4")])
        self.assertEqual(zs.cells, ("P4", "P40", "P41"))  # sorted, nested kept
        self.assertEqual(len(zs), 3)
        self.assertEqual([str(c) for c in zs], ["P4", "P40", "P41"])
        self.assertIn("P40", zs)
        self.assertIn(cell("P40"), zs)
        self.assertNotIn("P42", zs)
        self.assertEqual(zs.resolutions, (1, 2))
        self.assertEqual(zs, ZoneSet(rdggs, ["P41", "P40", "P4", "P4"]))  # a set
        self.assertNotEqual(zs, ZoneSet(rdggs, ["P4"]))
        self.assertEqual(len(ZoneSet(rdggs)), 0)
        with self.assertRaises(ValueError):
            ZoneSet(rdggs, ["P4", "X9"])
        with self.assertRaises(ValueError):
            ZoneSet(rdggs, [Cell(RHEALPixDGGS(N_side=4), (P, 0))])

    def test_parent_child_sibling_follow_the_spec_examples(self):
        # Clause 8.3.3's examples for zone 40, transposed to P40.
        p4, p40, p400 = cell("P4"), cell("P40"), cell("P400")
        children = [str(c) for c in p4.subcells()]
        grandchildren = [str(g) for c in p4.subcells() for g in c.subcells()]
        self.assertEqual(p40.parent().cells, ("P4",))
        self.assertEqual(p400.parent(levels=2).cells, ("P4", "P40"))
        self.assertEqual(p400.parent(levels=5).cells, ("P", "P4", "P40"))
        self.assertEqual(cell("P").parent().cells, ())
        self.assertEqual(set(p40.sibling(inherit_id=True).cells), set(children))
        self.assertEqual(set(p40.sibling().cells), set(k_ring("P40", 1)))
        self.assertEqual(set(p40.sibling(levels=2).cells), set(k_ring("P40", 2)))
        self.assertEqual(
            set(p400.sibling(inherit_id=True, levels=2).cells), set(grandchildren)
        )
        self.assertEqual(set(p4.child(inherit_id=True).cells), set(children))
        self.assertEqual(
            set(p4.child(inherit_id=True, levels=2).cells),
            set(children) | set(grandchildren),
        )
        # The same operations on a set apply to every member.
        zs = ZoneSet(WGS84_003, ["P40", "P41"])
        self.assertEqual(zs.parent().cells, ("P4",))
        self.assertEqual(len(zs.child()), 18)
        self.assertEqual(
            set(zs.sibling().cells), set(k_ring("P40", 1)) | set(k_ring("P41", 1))
        )
        with self.assertRaises(ValueError):
            p40.parent(levels=0)

    def test_set_algebra_matches_shapely_on_the_plane(self):
        rdggs = WGS84_003
        rng = random.Random(20260917)
        pool = [str(c) for f in (P, Q) for c in rdggs.cell((f,)).subcells()]
        pool += [str(g) for i in pool[:6] for g in cell(i).subcells()]
        for _ in range(12):
            a = ZoneSet(rdggs, rng.sample(pool, rng.randint(1, 6)))
            b = ZoneSet(rdggs, rng.sample(pool, rng.randint(1, 6)))
            ga, gb = planar_geometry(a), planar_geometry(b)
            for method, expected in (
                (a.union, ga.union(gb)),
                (a.intersection, ga.intersection(gb)),
                (a.difference, ga.difference(gb)),
                (a.sym_difference, ga.symmetric_difference(gb)),
            ):
                result = method(b)
                self.assertIsInstance(result, ZoneSet)
                # The result's cells are disjoint and, together, the
                # shapely region: equal area, and covered by it.
                cells = list(result)
                for x in cells:
                    for y in cells:
                        self.assertTrue(x is y or x.disjoint(y) or x.touches(y))
                got = planar_geometry(cells)
                self.assertAlmostEqual(
                    got.area / ga.area, expected.area / ga.area, places=9
                )
                self.assertLess(got.difference(expected).area, 1e-9 * ga.area)
                # Resolutions within the default range, the operands' own.
                low = min(a.resolutions + b.resolutions)
                high = max(a.resolutions + b.resolutions)
                self.assertTrue(all(low <= r <= high for r in result.resolutions))

    def test_set_algebra_details(self):
        rdggs = WGS84_003
        p4 = ZoneSet(rdggs, ["P4"])
        p40 = ZoneSet(rdggs, ["P40"])
        # Intersection with a descendant is the descendant.
        self.assertEqual(p4.intersection(p40).cells, ("P40",))
        self.assertEqual(p40.intersection(p4).cells, ("P40",))
        # Difference refines down to the removed cell's resolution.
        eight = p4.difference(p40)
        self.assertEqual(len(eight), 8)
        self.assertTrue(all(len(i) == 3 and i != "P40" for i in eight.cells))
        # The union of the eight children and the ninth compacts back to
        # P4 when the range allows resolution 1; by default the range is
        # the operands' own, resolution 2, so it stays as nine cells.
        self.assertEqual(eight.union(p40, min_res=1).cells, ("P4",))
        self.assertEqual(len(eight.union(p40)), 9)
        # The resolution range is a filter: cells finer than max_res drop.
        self.assertEqual(p4.difference(ZoneSet(rdggs, ["P400"]), max_res=1).cells, ())
        self.assertEqual(len(p4.difference(ZoneSet(rdggs, ["P400"]), max_res=2)), 8)
        # Disjoint operands.
        p5 = ZoneSet(rdggs, ["P5"])
        self.assertEqual(p4.intersection(p5).cells, ())
        self.assertEqual(p4.union(p5).cells, ("P4", "P5"))
        self.assertEqual(p4.difference(p5).cells, ("P4",))
        self.assertEqual(p4.sym_difference(p4).cells, ())
        # A whole face compacts to itself.
        self.assertEqual(
            ZoneSet(rdggs, [str(c) for c in rdggs.cell((P,)).subcells()])
            .union(ZoneSet(rdggs), min_res=0)
            .cells,
            ("P",),
        )
        # Cells and sets mix as operands; grids must agree.
        self.assertEqual(p4.union(cell("P5")).cells, ("P4", "P5"))
        with self.assertRaises(ValueError):
            p4.union(ZoneSet(WGS84_122, ["P0"]))

    def test_buffer(self):
        rdggs = WGS84_003
        p44 = cell("P44")
        w = p44.width()
        # Touching cells are at distance 0: buffer(0) is the 1-ring.
        self.assertEqual(set(p44.buffer(0).cells), set(k_ring("P44", 1)))
        # Ring 2 cells lie at planar distance w (edge) to sqrt(2) w (corner).
        self.assertEqual(set(p44.buffer(1.5 * w).cells), set(k_ring("P44", 2)))
        ring_two = set(k_ring("P44", 2)) - set(k_ring("P44", 1))
        self.assertEqual(
            set(p44.buffer(1.2 * w).cells) - set(k_ring("P44", 1)),
            {i for i in ring_two if p44.distance(cell(i)) <= 1.2 * w},
        )
        # On the ellipsoid every kept cell is within the distance and every
        # dropped neighbour of the result is not.
        d = 400_000.0
        kept = p44.buffer(d, plane=False)
        for c in kept:
            self.assertLessEqual(p44.distance(c, plane=False), d)
        fringe = {i for c in kept for i in k_ring(str(c), 1)} - set(kept.cells)
        for i in fringe:
            self.assertGreater(p44.distance(cell(i), plane=False), d)
        # A set buffers every member.
        zs = ZoneSet(rdggs, ["P44", "Q44"])
        self.assertEqual(
            set(zs.buffer(0).cells), set(k_ring("P44", 1)) | set(k_ring("Q44", 1))
        )

    def test_geometry_attributes(self):
        rdggs = WGS84_003
        zs = ZoneSet(rdggs, ["P40", "P41", "P43"])
        geometry = zs.geometry
        self.assertTrue(geometry.is_valid)
        self.assertAlmostEqual(
            geometry.area, sum(shapely.area(g) for g in zs.geometries), places=9
        )
        self.assertEqual(zs.boundary.geom_type, "LineString")
        self.assertTrue(zs.convex_hull.covers(geometry))
        self.assertIs(zs.boundary_type, BoundaryType.PROJECTED_LINE)
        self.assertEqual(BoundaryType.PROJECTED_LINE.value, "projectedLine")
        self.assertEqual(len(BoundaryType), 6)
        self.assertTrue(ZoneSet(rdggs).geometry.is_empty)


if __name__ == "__main__":
    unittest.main()
