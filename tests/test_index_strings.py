"""
Index strings, a face letter followed by one digit per resolution, are
defined for N_side 2 and 3 and read and written by one parser and one
formatter on the grid (issue #146). These tests check the round trips on
those grids and that every string entry point raises for any other grid.
"""

import numpy as np
import pytest
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


@pytest.mark.parametrize("rdggs", GRIDS)
def test_round_trip_on_every_cell_to_resolution_2(rdggs):
    assert rdggs.has_index_strings
    for resolution in range(3):
        for cell in rdggs.grid(resolution):
            index = str(cell)
            assert rdggs.format_index(cell.suid) == index
            assert rdggs.parse_index(index) == cell.suid
            assert rdggs.index_resolution(index) == resolution
            assert rdggs.index_ancestor(index, resolution) == index
            assert len(index) == resolution + 1
            assert rhp_wrappers.rhp_is_valid(index, rdggs)


@pytest.mark.parametrize("rdggs", GRIDS)
def test_parse_rejects_malformed_strings(rdggs):
    too_big = str(rdggs.N_side**2)
    bad = ["", "X", "n1", "N-1", "N 1", "N1.", "Na", "N²", "N" + too_big, "1N"]
    for index in bad:
        assert rdggs.parse_index(index) is None, index
        assert not rhp_wrappers.rhp_is_valid(index, rdggs), index


def test_digit_range_follows_n_side():
    assert WGS84_003.parse_index("P48") == ("P", 4, 8)
    assert WGS84_002.parse_index("P48") is None
    assert WGS84_002.parse_index("P31") == ("P", 3, 1)


@pytest.mark.parametrize("rdggs", GRIDS)
def test_ancestors_and_children(rdggs):
    cell = rdggs.cell(("P", 1, 0, 2))
    index = str(cell)
    ancestors = [rdggs.index_ancestor(index, r) for r in range(4)]
    assert ancestors == ["P", "P1", "P10", "P102"]
    assert rdggs.index_ancestor(index, 9) == index
    children = rdggs.index_children("P1")
    assert children == [str(c) for c in rdggs.cell(("P", 1)).subcells()]
    assert len(children) == rdggs.N_side**2
    assert rhp_wrappers.rhp_to_parent(index, dggs=rdggs) == "P10"
    assert rhp_wrappers.rhp_to_parent(index, res=1, dggs=rdggs) == "P1"
    assert rhp_wrappers.rhp_get_resolution(index, rdggs) == 3
    assert rhp_wrappers.rhp_get_base_cell(index, rdggs) == "P"


def test_center_child_uses_the_formatter():
    assert rhp_wrappers.rhp_to_center_child("P1", dggs=WGS84_003) == "P14"
    assert rhp_wrappers.rhp_to_center_child("P1", res=3, dggs=WGS84_003) == "P144"


def test_compaction_groups_by_parent():
    children = WGS84_002.index_children("P1")
    assert compact_cells(children, N_side=2) == {"P1"}
    assert compress_order_cells(children, N_side=2) == ["P1"]
    assert compact_cells(children[:-1], N_side=2) == set(children[:-1])


@pytest.mark.parametrize("rdggs", LARGE)
def test_large_n_side_has_no_index_strings(rdggs):
    assert not rdggs.has_index_strings
    cell = rdggs.cell(("N", 10, 15))
    assert str(cell) == "(N, 10, 15)"
    assert rdggs.cell(cell.suid[:-1]) == rdggs.cell(("N", 10))
    assert len(list(cell.subcells())) == rdggs.N_side**2


@pytest.mark.parametrize("rdggs", LARGE)
def test_every_string_entry_point_raises_for_large_n_side(rdggs):
    cell = rdggs.cell(("N", 10, 15))
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
        "polyfill_array": lambda: rhp_wrappers.polyfill_array(SQUARE, 1, dggs=rdggs),
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
        "compact_cells": lambda: compact_cells(["N1"], N_side=rdggs.N_side),
        "compress_order_cells": lambda: compress_order_cells(
            ["N1"], N_side=rdggs.N_side
        ),
    }
    for name, call in calls.items():
        with pytest.raises(ValueError, match="N_side 2 and 3"):
            call()
            pytest.fail(f"{name} did not raise")


def test_face_letters_are_the_resolution_0_index_strings():
    for rdggs in GRIDS:
        assert [str(c) for c in rdggs.grid(0)] == list(CELLS0)
        assert all(rdggs.parse_index(f) == (f,) for f in CELLS0)
        assert isinstance(Cell(rdggs, ["N"]), Cell)
