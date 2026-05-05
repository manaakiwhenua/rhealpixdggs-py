from shapely.geometry import Polygon, Point
from rhealpixdggs.dggs import RHEALPixDGGS, WGS84_003
from rhealpixdggs.cell import Cell
from itertools import compress
from collections.abc import Iterable


def get_finest_containing_cell(
    polygon: Polygon, rdggs: RHEALPixDGGS = WGS84_003
) -> Cell:
    """
    Finds the finest DGGS Cell containing a given cartesian polygon
    """

    def _cell_contains(cell: Cell, poly) -> bool:
        # Cap cells span all longitudes at the pole, so their boundary cannot
        # be represented as a simple polygon in geographic coordinates.
        # Project both sides to the rHEALPix plane and test containment there.
        if cell.ellipsoidal_shape() == "cap":
            plane_cell = Polygon(cell.vertices(plane=True))
            # MultiPolygon has no .exterior; check each part individually.
            parts = list(poly.geoms) if hasattr(poly, "geoms") else [poly]
            return all(
                plane_cell.contains(
                    Polygon([rdggs.rhealpix(x, y) for x, y in part.exterior.coords])
                )
                for part in parts
            )
        return Polygon(cell.vertices(plane=False)).contains(poly)

    def _get_finest_cell(polygon, suid, rdggs=rdggs):
        parent_cell = Cell(rdggs=rdggs, suid=suid)
        # get the children cells and polygons for these cells
        children_cells = [cell for cell in parent_cell.subcells()]
        # function and truth list for multipolygon / polygon (polygon) contained within multipolygon / polygon (cell)
        truth = [_cell_contains(cell, polygon) for cell in children_cells]
        # if we get something back, check the next level lower
        returned_cells = list(compress(children_cells, truth))
        if returned_cells:
            finest = _get_finest_cell(polygon, returned_cells[0].suid, rdggs)
        else:
            if _cell_contains(parent_cell, polygon):
                finest = parent_cell
            else:
                finest = None
        return finest

    for suid in [tuple(x) for x in ["N", "O", "P", "Q", "R", "S"]]:
        finest = _get_finest_cell(polygon, suid, rdggs)
        if finest is not None:
            return finest


# TODO class should be a general class for collections of cells (believe the term is 'zone'?)
class CellZoneFromPoly:
    def __init__(
        self,
        feature,
        res_limit,
        return_cells: bool,
        file=None,
        bounding_cell=None,
        rdggs=WGS84_003,
    ):
        self.label = feature[0]
        self.geometry = feature[1]
        self.res_limit = res_limit
        self.return_cells = return_cells
        if return_cells:
            self.cells_list = []
        # self.i = 0
        self.file = file
        if file:
            self.file.write(f"\n{self.label},")
        if bounding_cell is None:
            self._get_dggs_poly(get_finest_containing_cell(self.geometry, rdggs))
        else:
            self._get_dggs_poly(bounding_cell)

    def _get_dggs_poly(self, bounding_cell):
        bounding_poly = Polygon(bounding_cell.vertices(plane=False))
        if self.geometry.contains(
            bounding_poly
        ):  # edge case where the polygon is the same as the bounding cell
            self._write_cells(bounding_cell, bounding_poly, "bounding poly")
        else:
            if bounding_cell.resolution + 1 > self.res_limit:
                pass
            else:
                children_cells = [cell for cell in bounding_cell.subcells()]
                children_poly = [
                    Polygon(cell.vertices(plane=False)) for cell in children_cells
                ]
                together = list(zip(children_cells, children_poly))
                # print(f'processing chhildren for bounding cell {bounding_cell}')
                self._process_children(together)
        # print(f'*bounding cell* {bounding_cell}')
        # for i in self.cells_list:
        #     print(i)
        # print('****************')
        return self.cells_list

    def _process_children(self, together):
        for child_cell, child_poly in together:
            # 1: add contained cells
            if self.geometry.contains(child_poly):
                self._write_cells(child_cell, child_poly, "fully contained")
            # 2: check we're not at the limit, if we are, check centroids
            elif child_cell.resolution == self.res_limit:
                if self.geometry.contains(Point(child_cell.nucleus(plane=False))):
                    self._write_cells(child_cell, child_poly, "nucleus")
            # 3: check the children (call this same function on the children)
            else:
                if self.geometry.overlaps(child_poly):
                    self._get_dggs_poly(child_cell)

    def _write_cells(self, cell, poly, desc):
        """
        Writes cell / polygon details to either or both of a list or a file.
        """
        if self.file is not None:
            self.file.write(f"{str(cell)} {desc}\n ")
            # self.file.write(f"{poly.wkt}\n")
        if self.return_cells:
            # cell = self.add_int_suid(cell)
            self.cells_list.append(cell)


def compress_order_cells(cells: list[str], recursive: bool = False) -> list[str]:
    """
    Compress and sort a collection of cells.

    Each complete group of 9 siblings is replaced by their parent cell.
    By default only one level of compression is performed. Pass
    ``recursive=True`` to repeat until no further compression is possible.
    The result is always deduplicated and sorted alphanumerically.
    """
    import re

    def alphanum_sort(lst: list[str]) -> list[str]:
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
        return sorted(lst, key=alphanum_key)

    def _compress_once(cell_set: set[str]) -> set[str]:
        upper_cells: dict[str, list[str]] = {}
        for cell in cell_set:
            upper_cells.setdefault(cell[:-1], []).append(cell)
        result: set[str] = set()
        for k, v in upper_cells.items():
            if len(v) == 9:
                result.add(k)
            else:
                result.update(v)
        return result

    cell_set = set(cells)
    if recursive:
        while True:
            next_set = _compress_once(cell_set)
            if next_set == cell_set:
                break
            cell_set = next_set
    else:
        cell_set = _compress_once(cell_set)
    return alphanum_sort(list(cell_set))


def compact_cells(cells: Iterable[str]) -> set[str]:
    """
    Compact a set of rHEALPix DGGS cells by repeatedly merging complete groups
    of 9 siblings into their parent until no further merging is possible.

    All input cells must be at the same resolution (equal string length).
    Raises ValueError if the input contains cells at mixed resolutions.
    Returns a set; if the caller needs ordering, sort the result separately.
    """
    cell_set = set(cells)
    if not cell_set:
        return cell_set

    lengths = {len(c) for c in cell_set}
    if len(lengths) > 1:
        raise ValueError(
            f"All input cells must be at the same resolution; "
            f"got cells with string lengths {sorted(lengths)}"
        )

    while True:
        upper: dict[str, list[str]] = {}
        for cell in cell_set:
            upper.setdefault(cell[:-1], []).append(cell)

        next_set: set[str] = set()
        for parent, children in upper.items():
            if len(children) == 9:
                next_set.add(parent)
            else:
                next_set.update(children)

        if next_set == cell_set:
            return cell_set
        cell_set = next_set
