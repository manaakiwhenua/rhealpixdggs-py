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
    def __init__(self, feature, res_limit, return_cells: bool, file=None):
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

        self._get_dggs_poly(call_get_finest(self.geometry))

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
#
#
# # TODO extract compression function from below and benchmark against Nicks, also Nick's doesn't require integer suids
# # TODO remove dependence on dataframes; use lists - although dask + geopandas is worth in time exploring for performance ..
# with open(file_out_path, 'a') as file:
#     for row in act_geoms.itertuples():
#         reformatted_uri = row.Geom_ID.replace('http','https').replace('<','').replace('>','')
#         file.write(f"{reformatted_uri},")
#         cells_list = PolygonAsDGGS(polygon=(row.Geom_ID, row.Geometries),
#                                    res_limit=resolution,
#                                    return_cells=True,
#                                    file=None). \
#             get_dggs_poly(bounding_cell=row.finest_cell)
#         if cells_list:
#             # cell compression
#             suids_cells_list = list(zip([int(cell.integer_suid) for cell in cells_list], cells_list))
#             first_sort = sorted(suids_cells_list, key=lambda x: x[0])
#             for i, ele in enumerate(first_sort[0:-8]):
#                 if ele[0] % 10 == 0:
#                     try:
#                         if ele[0] + 8 == first_sort[i + 8][0]:
#                             # print(f'Full set of children cells found {ele[0]}:{first_sort[i + 8][0]}')
#                             first_sort.append((
#                                 int(str(ele[0])[:-1]),
#                                 Cell(TB16Pix, suid=ele[1].suid[:-1])))
#                             for n in range(i,i+9):
#                                 first_sort.pop(i)
#                     except IndexError:
#                         pass
#             # sort again to move compressed cells to the front of the list
#             second_sort = sorted(first_sort, key=lambda x: x[0])
#             for cell in second_sort:
#                 file.write(f"{str(cell[1])} ")
#         file.write('\n')
#
#
# def add_int_suid(self, cell):
#     """
#     Adds a 'numerical suid' to a cell. This is just the suid with the alphabetic character substituted for a digit.
#     Useful for:
#         1. Ordering cells (though not impossible to do with the alphanumeric suid)
#         2. Part of the cell compression process, to check whether all children cells are present in a sorted list.
#     Be careful with what integer suids are used for; they will only be functionally equivalent to regular suids in
#     specific cases!
#     """
#     replace_dict = {'N': 0, 'O': 1, 'P': 2, 'Q': 3, 'R': 4, 'S': 5}
#     tuple_int_suid = tuple([replace_dict[x] if type(x) is str else x for x in cell.suid])
#     cell.integer_suid = reduce(lambda sub, ele: sub * 10 + ele, tuple_int_suid)
#
#     # can this be changed to just replace the first character of a suid and append the other digits,
#     # something like: tuple(replace_dict[cell.suid[0]] + cell.suid[1:])
#     return cell