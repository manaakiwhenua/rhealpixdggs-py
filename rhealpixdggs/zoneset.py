"""
Sets of cells, and the zone-set operations of OGC Topic 21 v2.0 Table 53.

A ``ZoneSet`` holds cells of one grid, at any resolutions, exactly as
given: the spec's hierarchy operations return sets such as ``{40, 4}``
whose members nest. The set operations ``union``, ``intersection``,
``difference`` and ``sym_difference`` treat a set as the region its cells
cover and return disjoint cells: exact hierarchical algebra on the index
strings (the intersection of two nested cells is the finer one; a
difference refines a cell down to the resolution of what is removed from
it), then the spec's ``rangeRefine`` filter on resolution, with complete
sibling groups compacted no coarser than the range's lower bound.
``buffer`` grows a set by a distance; ``parent``, ``child`` and ``sibling``
walk the hierarchy by ``levels`` with the ``inheritID`` flag of clause
8.3.3; and ``geometry``, ``boundary``, ``convex_hull`` and
``boundary_type`` are the attributes of the query's combined geometry.
"""

import enum
from collections.abc import Iterable, Iterator
from typing import Any

import shapely

from rhealpixdggs import export
from rhealpixdggs.cell import Cell, _ring_step_neighbors
from rhealpixdggs.dggs import RHEALPixDGGS


class BoundaryType(enum.Enum):
    """
    The EA_BoundaryType code list of OGC Topic 21 v2.0 Table 76: the kinds
    of curve an equal-area cell boundary may be built from. Every edge of
    an rHEALPix cell is a straight line in the equal-area plane, which is
    PROJECTED_LINE; the edges of equatorial cells are also parallels
    (SMALL_CIRCLE) and meridians (GREAT_ELLIPSE).
    """

    GEODESIC = "geodesic"
    SMALL_CIRCLE = "smallCircle"
    GREAT_CIRCLE = "greatCircle"
    SMALL_ELLIPSE = "smallEllipse"
    GREAT_ELLIPSE = "greatEllipse"
    PROJECTED_LINE = "projectedLine"


def _covered(rdggs: RHEALPixDGGS, index: str, cells: frozenset[str] | set[str]) -> bool:
    """True if `index` is one of `cells` or a descendant of one."""
    return any(
        rdggs.index_ancestor(index, r) in cells
        for r in range(rdggs.index_resolution(index) + 1)
    )


def _within(rdggs: RHEALPixDGGS, index: str, ancestor: str) -> bool:
    """True if `index` is `ancestor` or a descendant of it."""
    return _covered(rdggs, index, {ancestor})


class ZoneSet:
    """
    A set of cells of one grid `rdggs`, given as index strings or ``Cell``
    objects and kept as index strings, deduplicated and otherwise exactly
    as given. Iterating yields ``Cell`` objects; ``cells`` is the sorted
    tuple of index strings. Raise ValueError for an invalid index or a
    cell of another grid.

    EXAMPLES::

        >>> from rhealpixdggs.dggs import WGS84_003
        >>> zs = ZoneSet(WGS84_003, ['P41', 'P40', 'P4'])
        >>> zs.cells, len(zs), 'P40' in zs, zs.resolutions
        (('P4', 'P40', 'P41'), 3, True, (1, 2))

    """

    def __init__(self, rdggs: RHEALPixDGGS, cells: Iterable[str | Cell] = ()) -> None:
        self.rdggs = rdggs
        indices = []
        for item in cells:
            if isinstance(item, Cell) and item.rdggs != rdggs:
                raise ValueError("Cells must belong to the ZoneSet's grid.")
            indices.append(str(item))
        invalid = [
            i for i, s in zip(indices, rdggs.shapes(indices), strict=True) if s == ""
        ]
        if invalid:
            raise ValueError(f"invalid cell indices: {invalid[:10]!r}")
        self._cells: frozenset[str] = frozenset(indices)

    @property
    def cells(self) -> tuple[str, ...]:
        """The index strings, sorted."""
        return tuple(sorted(self._cells))

    @property
    def resolutions(self) -> tuple[int, ...]:
        """The distinct resolutions present, ascending."""
        return tuple(sorted({self.rdggs.index_resolution(i) for i in self._cells}))

    def _cell(self, index: str) -> Cell:
        suid = self.rdggs.parse_index(index)
        if suid is None:
            raise ValueError(f"invalid cell index {index!r}")
        return self.rdggs.cell(suid)

    def _parent(self, index: str) -> str:
        """The parent's index string, the empty string for a resolution 0 cell."""
        resolution = self.rdggs.index_resolution(index)
        return (
            "" if resolution == 0 else self.rdggs.index_ancestor(index, resolution - 1)
        )

    def __len__(self) -> int:
        return len(self._cells)

    def __iter__(self) -> Iterator[Cell]:
        return (self._cell(index) for index in self.cells)

    def __contains__(self, item: object) -> bool:
        return str(item) in self._cells

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ZoneSet):
            return NotImplemented
        return self.rdggs == other.rdggs and self._cells == other._cells

    def __hash__(self) -> int:
        return hash(self._cells)

    def __repr__(self) -> str:
        return f"ZoneSet({list(self.cells)!r})"

    def _operand(self, other: "Cell | ZoneSet") -> "ZoneSet":
        if isinstance(other, Cell):
            return ZoneSet(self.rdggs, [other])
        if other.rdggs != self.rdggs:
            raise ValueError("ZoneSets must belong to the same grid.")
        return other

    def _disjoint(self) -> set[str]:
        """The cells with those nested in another member dropped."""
        return {
            i
            for i in self._cells
            if self.rdggs.index_resolution(i) == 0
            or not _covered(self.rdggs, self._parent(i), self._cells)
        }

    def _finish(
        self,
        cells: set[str],
        other: "ZoneSet",
        min_res: int | None,
        max_res: int | None,
    ) -> "ZoneSet":
        """
        Apply the spec's ``rangeRefine`` filter to a set operation's
        result and compact complete sibling groups no coarser than
        `min_res`. The default range spans both operands' resolutions.
        """
        resolutions = self.resolutions + other.resolutions or (
            0,
            self.rdggs.max_resolution,
        )
        low = min(resolutions) if min_res is None else min_res
        high = max(resolutions) if max_res is None else max_res
        resolution = self.rdggs.index_resolution
        kept = {i for i in cells if low <= resolution(i) <= high}
        full = self.rdggs.N_side**2
        while True:
            groups: dict[str, set[str]] = {}
            for index in kept:
                if resolution(index) - 1 >= low:
                    groups.setdefault(self._parent(index), set()).add(index)
            merged = {p for p, children in groups.items() if len(children) == full}
            if not merged:
                return ZoneSet(self.rdggs, kept)
            kept = {i for i in kept if self._parent(i) not in merged} | merged

    def _subtract(self, index: str, inside: list[str]) -> set[str]:
        """`index` minus the cells `inside` it (proper descendants)."""
        if not inside:
            return {index}
        result: set[str] = set()
        for child in self.rdggs.index_children(index):
            if child in inside:
                continue
            result |= self._subtract(
                child, [b for b in inside if _within(self.rdggs, b, child)]
            )
        return result

    def union(
        self,
        other: "Cell | ZoneSet",
        min_res: int | None = None,
        max_res: int | None = None,
    ) -> "ZoneSet":
        """
        The region covered by this set or `other`, as disjoint cells
        within the resolution range `min_res` to `max_res` (default: the
        operands' own range), compacted no coarser than `min_res`.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> eight = ZoneSet(WGS84_003, [f'P4{d}' for d in range(8)])
            >>> eight.union(ZoneSet(WGS84_003, ['P48']), min_res=1).cells
            ('P4',)
            >>> len(eight.union(ZoneSet(WGS84_003, ['P48'])))  # default range: resolution 2
            9

        """
        other = self._operand(other)
        return self._finish(
            ZoneSet(self.rdggs, self._cells | other._cells)._disjoint(),
            other,
            min_res,
            max_res,
        )

    def intersection(
        self,
        other: "Cell | ZoneSet",
        min_res: int | None = None,
        max_res: int | None = None,
    ) -> "ZoneSet":
        """
        The region covered by both this set and `other`: for two nested
        cells the finer one. Range and compaction as for ``union``.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> ZoneSet(WGS84_003, ['P4']).intersection(ZoneSet(WGS84_003, ['P40', 'P5'])).cells
            ('P40',)

        """
        other = self._operand(other)
        a, b = self._disjoint(), other._disjoint()
        cells = {i for i in a if _covered(self.rdggs, i, b)} | {
            i for i in b if _covered(self.rdggs, i, a)
        }
        return self._finish(cells, other, min_res, max_res)

    def difference(
        self,
        other: "Cell | ZoneSet",
        min_res: int | None = None,
        max_res: int | None = None,
    ) -> "ZoneSet":
        """
        The region covered by this set but not by `other`: a cell partly
        removed is refined to the resolution of the cells removed from
        it. Range and compaction as for ``union``.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> ZoneSet(WGS84_003, ['P4']).difference(ZoneSet(WGS84_003, ['P40'])).cells
            ('P41', 'P42', 'P43', 'P44', 'P45', 'P46', 'P47', 'P48')

        """
        other = self._operand(other)
        a, b = self._disjoint(), other._disjoint()
        cells: set[str] = set()
        for index in a:
            if _covered(self.rdggs, index, b):
                continue
            inside = [i for i in b if i != index and _within(self.rdggs, i, index)]
            cells |= self._subtract(index, inside)
        return self._finish(cells, other, min_res, max_res)

    def sym_difference(
        self,
        other: "Cell | ZoneSet",
        min_res: int | None = None,
        max_res: int | None = None,
    ) -> "ZoneSet":
        """
        The region covered by exactly one of this set and `other`. Range
        and compaction as for ``union``.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> ZoneSet(WGS84_003, ['P4']).sym_difference(ZoneSet(WGS84_003, ['P4', 'P5'])).cells
            ('P5',)

        """
        other = self._operand(other)
        return self._finish(
            set(self.difference(other, 0, self.rdggs.max_resolution)._cells)
            | set(other.difference(self, 0, self.rdggs.max_resolution)._cells),
            other,
            min_res,
            max_res,
        )

    def _shares_interior(self, other: "ZoneSet") -> bool:
        a, b = self._disjoint(), other._disjoint()
        return any(_covered(self.rdggs, i, b) for i in a) or any(
            _covered(self.rdggs, i, a) for i in b
        )

    def _only(self, other: "ZoneSet") -> bool:
        """True if this set covers ground `other` does not."""
        return bool(self.difference(other, 0, self.rdggs.max_resolution)._cells)

    def _cells_touch(self, other: "ZoneSet") -> bool:
        cells_a = [self._cell(i) for i in self._disjoint()]
        cells_b = [self._cell(i) for i in other._disjoint()]
        return any(a.touches(b) for a in cells_a for b in cells_b)

    def equals(self, other: "Cell | ZoneSet") -> bool:
        """
        DE-9IM `equals` between the regions this set and `other` cover:
        True if neither covers ground the other does not, however the
        cells are cut (nine children equal their parent).

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> children = [f'P4{d}' for d in range(9)]
            >>> ZoneSet(WGS84_003, children).equals(ZoneSet(WGS84_003, ['P4']))
            True

        """
        other = self._operand(other)
        return not self._only(other) and not other._only(self)

    def contains(self, other: "Cell | ZoneSet") -> bool:
        """
        DE-9IM `contains` between regions: True if the region `other`
        covers lies within the region this set covers. As for cells, this
        coincides with `covers`, since the regions are closed and tile.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> ZoneSet(WGS84_003, ['P4']).contains(ZoneSet(WGS84_003, ['P40', 'P41']))
            True

        """
        return not self._operand(other)._only(self)

    def within(self, other: "Cell | ZoneSet") -> bool:
        """DE-9IM `within` between regions: the converse of `contains`."""
        return self._operand(other).contains(self)

    def intersects(self, other: "Cell | ZoneSet") -> bool:
        """
        DE-9IM `intersects` between regions: True if they share a point,
        that is share interior or have a pair of cells that touch.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> ZoneSet(WGS84_003, ['P40']).intersects(ZoneSet(WGS84_003, ['P41']))
            True

        """
        other = self._operand(other)
        return self._shares_interior(other) or self._cells_touch(other)

    def disjoint(self, other: "Cell | ZoneSet") -> bool:
        """DE-9IM `disjoint` between regions: the negation of `intersects`."""
        return not self.intersects(other)

    def touches(self, other: "Cell | ZoneSet") -> bool:
        """
        DE-9IM `touches` between regions: True if they share boundary
        points but no interior.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> a = ZoneSet(WGS84_003, ['P40', 'P41'])
            >>> a.touches(ZoneSet(WGS84_003, ['P42'])), a.touches(ZoneSet(WGS84_003, ['P4']))
            (True, False)

        """
        other = self._operand(other)
        return not self._shares_interior(other) and self._cells_touch(other)

    def overlaps(self, other: "Cell | ZoneSet") -> bool:
        """
        DE-9IM `overlaps` between regions: True if they share interior and
        each covers ground the other does not. Unlike two single cells,
        two sets of cells can overlap; ``Cell.region_overlaps`` is the
        one-cell case.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> a = ZoneSet(WGS84_003, ['P40', 'P41'])
            >>> a.overlaps(ZoneSet(WGS84_003, ['P41', 'P42']))
            True
            >>> a.overlaps(ZoneSet(WGS84_003, ['P4']))
            False

        """
        other = self._operand(other)
        return self._shares_interior(other) and self._only(other) and other._only(self)

    def crosses(self, other: "Cell | ZoneSet") -> bool:
        """
        DE-9IM `crosses` between regions: always False, as both are
        two-dimensional; present because the interface mandates it.
        """
        self._operand(other)
        return False

    def buffer(self, dist: float, plane: bool = True, n: int = 8) -> "ZoneSet":
        """
        The cells at each member's resolution within `dist` of the set:
        every member, and every cell whose ``Cell.distance`` to a member is
        at most `dist` (planar distance with `plane` = True, geodesic with
        `plane` = False; see ``Cell.distance`` for units and `n`). With
        `dist` 0 that is the members and everything touching them.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> len(ZoneSet(WGS84_003, ['P44']).buffer(0))
            9

        """
        result: set[str] = set(self._cells)
        for cell in self:
            ring, seen = [cell], {cell.suid}
            while ring:
                candidates = []
                for member in ring:
                    for neighbor in _ring_step_neighbors(member):
                        if neighbor.suid not in seen:
                            seen.add(neighbor.suid)
                            candidates.append(neighbor)
                # Every cell beyond this ring is farther than the nearest
                # cell in it, so stop once none of a ring is within reach.
                ring = [
                    c for c in candidates if cell.distance(c, plane=plane, n=n) <= dist
                ]
                result.update(str(c) for c in ring)
        return ZoneSet(self.rdggs, result)

    def parent(self, levels: int = 1, inherit_id: bool = False) -> "ZoneSet":
        """
        The ancestors of the members from one to `levels` levels up, as
        many as exist: ``parent`` of Table 53. `inherit_id` is accepted for
        the spec's signature; in rHEALPix's nested hierarchy every ancestor
        shares its descendants' identifier prefix, so it changes nothing.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> ZoneSet(WGS84_003, ['P400']).parent(levels=2).cells
            ('P4', 'P40')

        """
        if levels < 1:
            raise ValueError("levels must be at least 1.")
        return ZoneSet(
            self.rdggs,
            {i[:-k] for i in self._cells for k in range(1, levels + 1) if len(i) > k},
        )

    def child(self, levels: int = 1, inherit_id: bool = False) -> "ZoneSet":
        """
        The descendants of the members from one to `levels` levels down:
        ``child`` of Table 53, including every intermediate level as the
        spec's examples do. `inherit_id` changes nothing in a nested
        hierarchy. Raise ValueError if a level would exceed the grid's
        ``max_resolution``.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> len(ZoneSet(WGS84_003, ['P4']).child(levels=2))
            90

        """
        if levels < 1:
            raise ValueError("levels must be at least 1.")
        result: set[str] = set()
        for cell in self:
            assert cell.resolution is not None
            for resolution in range(cell.resolution + 1, cell.resolution + levels + 1):
                result.update(str(c) for c in cell.subcells(resolution))
        return ZoneSet(self.rdggs, result)

    def sibling(self, levels: int = 1, inherit_id: bool = False) -> "ZoneSet":
        """
        The siblings of the members, themselves included: ``sibling`` of
        Table 53. With `inherit_id` True, the cells at a member's
        resolution under its ancestor `levels` levels up (its parent's
        children for `levels` 1); with `inherit_id` False, inheritance is
        ignored and the siblings are the cells within `levels` edge or
        corner steps at the same resolution, the k-ring of that radius, as
        the spec's example ``40.sibling()`` = {13, 22, 23, 31, 40, 41, 42,
        43} shows.

        EXAMPLES::

            >>> from rhealpixdggs.dggs import WGS84_003
            >>> len(ZoneSet(WGS84_003, ['P40']).sibling(inherit_id=True))
            9
            >>> len(ZoneSet(WGS84_003, ['P40']).sibling())
            9
            >>> len(ZoneSet(WGS84_003, ['P44']).sibling(levels=2))
            25

        """
        if levels < 1:
            raise ValueError("levels must be at least 1.")
        result: set[str] = set()
        for cell in self:
            assert cell.resolution is not None
            if inherit_id:
                ancestor = Cell(
                    self.rdggs, cell.suid[: max(1, len(cell.suid) - levels)]
                )
                result.update(str(c) for c in ancestor.subcells(cell.resolution))
                continue
            ring, seen = [cell], {cell.suid}
            for _ in range(levels):
                next_ring = []
                for member in ring:
                    for neighbor in _ring_step_neighbors(member):
                        if neighbor.suid not in seen:
                            seen.add(neighbor.suid)
                            next_ring.append(neighbor)
                ring = next_ring
            result.update(str(Cell(self.rdggs, suid)) for suid in seen)
        return ZoneSet(self.rdggs, result)

    @property
    def geometries(self) -> Any:
        """The members' geometries in longitude-latitude degrees (``export.geometries``)."""
        return export.geometries(self.cells, self.rdggs)

    @property
    def geometry(self) -> Any:
        """The union of the members' geometries: the query's combined geometry."""
        return shapely.union_all(self.geometries) if len(self) else shapely.Polygon()

    @property
    def boundary(self) -> Any:
        """The boundary of the combined geometry (Table 53's ``boundary``)."""
        return self.geometry.boundary

    @property
    def convex_hull(self) -> Any:
        """
        The convex hull of the combined geometry in longitude-latitude
        coordinates (Table 53's ``convexHull``); a set split by the
        antimeridian hulls across the map, so this is a planar notion.
        """
        return self.geometry.convex_hull

    @property
    def boundary_type(self) -> BoundaryType:
        """
        Table 53's ``boundaryType``: the kind of curve the boundary is made
        of. Every rHEALPix cell edge is a straight line in the equal-area
        plane, ``BoundaryType.PROJECTED_LINE``.
        """
        return BoundaryType.PROJECTED_LINE
