Zone Sets
=========

OGC Topic 21 v2.0 defines its zone queries (Table 53) between *zonal
identifiers*, and several of them return or compare sets of zones: the
set operations ``union``, ``intersection``, ``difference`` and
``symDifference``, ``buffer``, the hierarchy walks ``parent``, ``child``
and ``sibling``, and the attributes ``boundary``, ``boundaryType`` and
``convexHull`` of "the combined spatial geometries of the zones in the
query". :class:`~rhealpixdggs.zoneset.ZoneSet` is that set of zones, and
:class:`~rhealpixdggs.cell.Cell` has the same methods for the
zone-to-zone form, each returning a ``ZoneSet``. The DE-9IM predicates,
``distance`` and the one-dimensional relations are described on the
:doc:`cell` page.

A set of cells
--------------

A ``ZoneSet`` holds cells of one grid at any resolutions, exactly as
given: index strings, ``Cell`` objects or the array ``polyfill_array``
returns. Nested members are allowed, because the spec's own hierarchy
results contain them::

    >>> from rhealpixdggs.dggs import WGS84_003
    >>> from rhealpixdggs.zoneset import ZoneSet
    >>> zs = ZoneSet(WGS84_003, ['P40', 'P41', 'P43'])
    >>> zs.cells, len(zs), zs.resolutions
    (('P40', 'P41', 'P43'), 3, (2,))
    >>> [str(c) for c in zs][:2], 'P41' in zs
    (['P40', 'P41'], True)

Hierarchy
---------

``parent``, ``child`` and ``sibling`` take the spec's ``levels`` and
``inheritID`` arguments (clause 8.3.3). The worked example there, for a
zone ``40``, reads as follows for the cell P40::

    >>> p40 = WGS84_003.cell(['P', 4, 0])
    >>> p40.parent().cells
    ('P4',)
    >>> WGS84_003.cell(['P', 4, 0, 0]).parent(levels=2).cells
    ('P4', 'P40')
    >>> sorted(p40.sibling(inherit_id=True).cells)
    ['P40', 'P41', 'P42', 'P43', 'P44', 'P45', 'P46', 'P47', 'P48']
    >>> len(p40.sibling())   # inheritance ignored: the cell and its 8 neighbours
    9
    >>> len(WGS84_003.cell(['P', 4]).child(levels=2))   # children and grandchildren
    90

With ``inheritID`` True, siblings share a parent; with it False, the
default, inheritance is ignored and the siblings of a cell are the cells
at its resolution within ``levels`` edge or corner steps, its k-ring, as
the spec's ``40.sibling()`` = {13, 22, 23, 31, 40, 41, 42, 43} shows.
For ``parent`` and ``child`` the flag changes nothing, since every
rHEALPix cell's identifier extends its parent's. The Boolean forms are
:meth:`~rhealpixdggs.cell.Cell.parent_of`,
:meth:`~rhealpixdggs.cell.Cell.child_of` and
:meth:`~rhealpixdggs.cell.Cell.sibling_of`::

    >>> p32 = WGS84_003.cell(['P', 3, 2])   # left of P40, in the next parent
    >>> p40.sibling_of(p32), p40.sibling_of(p32, inherit_id=True)
    (True, False)

Set operations
--------------

The set operations treat a ``ZoneSet`` as the region its cells cover and
return disjoint cells. They are exact on the index strings: the
intersection of two nested cells is the finer one, and a difference
refines a cell down to the resolution of what is removed from it::

    >>> p4 = ZoneSet(WGS84_003, ['P4'])
    >>> p4.intersection(ZoneSet(WGS84_003, ['P40', 'P5'])).cells
    ('P40',)
    >>> p4.difference(ZoneSet(WGS84_003, ['P40'])).cells
    ('P41', 'P42', 'P43', 'P44', 'P45', 'P46', 'P47', 'P48')

The spec's ``rangeRefine`` argument, here ``min_res`` and ``max_res``,
filters the result to a range of resolutions, both ends included; its
default is the range spanned by the two operands. Within the range,
complete groups of siblings are compacted to their parent, but never
coarser than ``min_res``::

    >>> eight = p4.difference(ZoneSet(WGS84_003, ['P40']))
    >>> len(eight.union(ZoneSet(WGS84_003, ['P40'])))   # default range: resolution 2
    9
    >>> eight.union(ZoneSet(WGS84_003, ['P40']), min_res=1).cells
    ('P4',)

``buffer`` grows a set by a distance, using
:meth:`~rhealpixdggs.cell.Cell.distance`: every member, and every cell at
a member's resolution whose distance to it is at most ``dist``, in the
plane or on the ellipsoid::

    >>> len(WGS84_003.cell(['P', 4, 4]).buffer(0))   # touching cells are at distance 0
    9

Relations between sets
----------------------

The DE-9IM questions the :doc:`cell` page answers for two cells are also
answered between two sets, read as the regions their cells cover:
``equals``, ``contains``, ``within``, ``intersects``, ``disjoint``,
``touches``, ``overlaps`` and ``crosses``. Unlike two single cells, two
sets can genuinely overlap, sharing some ground while each has ground of
its own; ``Cell.region_overlaps`` is the one-cell case of that::

    >>> a = ZoneSet(WGS84_003, ['P40', 'P41'])
    >>> a.overlaps(ZoneSet(WGS84_003, ['P41', 'P42']))
    True
    >>> a.touches(ZoneSet(WGS84_003, ['P42'])), a.within(ZoneSet(WGS84_003, ['P4']))
    (True, True)
    >>> ZoneSet(WGS84_003, [f'P4{d}' for d in range(9)]).equals(ZoneSet(WGS84_003, ['P4']))
    True

Geometry of a query
-------------------

The combined geometry of the set, its boundary and its convex hull come
from the export geometries (:doc:`export`), so they are shapely objects in
longitude-latitude degrees with the antimeridian and pole conventions of
that page, and the convex hull is a planar hull of that map. The
``boundaryType`` attribute is a value of the EA_BoundaryType code list,
:class:`~rhealpixdggs.zoneset.BoundaryType`: every rHEALPix cell edge is a
straight line in the equal-area plane, ``PROJECTED_LINE``, and the edges
of equatorial cells are also parallels and meridians::

    >>> zs.geometry.geom_type, zs.boundary.geom_type
    ('Polygon', 'LineString')
    >>> zs.boundary_type.value
    'projectedLine'

Module reference
----------------

.. automodule:: rhealpixdggs.zoneset
    :members:
    :undoc-members:
    :show-inheritance:
