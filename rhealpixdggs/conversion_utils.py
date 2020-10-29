from shapely.geometry import Polygon, Point
from rhealpixdggs.dggs import Cell
from itertools import compress
from functools import reduce


def call_get_finest(polygon):
    """
    Returns the finest DGGS cell covering a cartesian polygon.
    Requires input as a shapely polygon.
    """
    for suid in [tuple(x) for x in ['N', 'O', 'P', 'Q', 'R', 'S']]:
        finest = _get_finest_cell(polygon, suid)
        if finest is not None:
            return finest


def _get_finest_cell(polygon, suid):
    parent_cell = Cell(suid=suid)
    # get the children cells and polygons for these cells
    children_cells = [cell for cell in parent_cell.subcells()]
    children_poly = [Polygon(cell.vertices(plane=False)) for cell in children_cells]
    # function and truth list for multipolygon / polygon (polygon) contained within multipolygon / polygon (cell)
    truth = [poly.contains(polygon) for poly in children_poly]
    # if we get something back, check the next level lower
    returned_cells = list(compress(children_cells, truth))
    if returned_cells:
        finest = _get_finest_cell(polygon, returned_cells[0].suid)
    else:
        parent_poly = Polygon(parent_cell.vertices(plane=False))
        if parent_poly.contains(polygon):
            finest = parent_cell
        else:
            finest = None
    return finest


# TODO class should be a general class for collections of cells (believe the term is 'zone'?)
class CellZoneFromPoly:
    def __init__(self, feature, res_limit, return_cells: bool, file=None, bounding_cell=None):
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
            self._get_dggs_poly(call_get_finest(self.geometry))
        else:
            self._get_dggs_poly(bounding_cell)

    def _get_dggs_poly(self, bounding_cell):
        bounding_poly = Polygon(bounding_cell.vertices(plane=False))
        if self.geometry.contains(bounding_poly):  # edge case where the polygon is the same as the bounding cell
            self._write_cells(bounding_cell, bounding_poly, 'bounding poly')
        else:
            if bounding_cell.resolution + 1 > self.res_limit:
                pass
            else:
                children_cells = [cell for cell in bounding_cell.subcells()]
                children_poly = [Polygon(cell.vertices(plane=False)) for cell in children_cells]
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
                self._write_cells(child_cell, child_poly, 'fully contained')
            # 2: check we're not at the limit, if we are, check centroids
            elif child_cell.resolution == self.res_limit:
                if self.geometry.contains(Point(child_cell.nucleus(plane=False))):
                    self._write_cells(child_cell, child_poly, 'nucleus')
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

def compress_order_cells(cells):
    """
    Compresses and order a set of cells
    """
    upper_cells = {}
    for cell in cells:
        upper_cells.setdefault(cell[:-1], []).append(cell)
    compressed_cells = []
    for k, v in upper_cells.items():
        if len(v) == 9:
            compressed_cells.append(k)
        else:
            compressed_cells.extend(v)
    ordered_compressed_cells = sorted(compressed_cells, key=len)
    return ordered_compressed_cells
