"""
Cell identifiers of an rHEALPix grid and the arithmetic behind them: index
strings, the orders of a cell's subcells, a cell's level, row and column in
the planar grid, and the uint64 and level-row-col text forms of that
position. All of it depends only on a grid's layout (``N_side``, the polar
squares and ``max_resolution``), not on its ellipsoid or projection. The
supported API is the ``RHEALPixDGGS`` and ``Cell`` methods built on it.
"""

import functools
import re
from collections.abc import Iterable, Sequence
from math import isqrt
from typing import Literal, Protocol, get_args

import numpy as np

#: The resolution 0 cells, the faces of the cube, which begin every identifier.
CELLS0 = ["N", "O", "P", "Q", "R", "S"]


class _GridLayout(Protocol):
    """What identifier arithmetic needs of a grid; ``RHEALPixDGGS`` has it."""

    N_side: int
    north_square: int
    south_square: int
    max_resolution: int


# Index strings ------------------------------------------------------------

#: The ``N_side`` values whose grids have index strings (issue #146).
INDEX_STRING_N_SIDES = (2, 3)


def _index_width(N_side: int) -> int:
    """
    Characters per digit in the index strings of a grid with `N_side`.
    Index strings, a face letter followed by one digit per resolution, are
    defined for ``N_side`` 2 and 3 (issue #146); raise ValueError otherwise.
    """
    if N_side not in INDEX_STRING_N_SIDES:
        raise ValueError(
            f"index strings are defined for N_side 2 and 3, not N_side={N_side}; "
            "address the cells of this grid by their suid tuples (issue #146)"
        )
    return 1


def _parse_index(index: str, N_side: int) -> tuple[str | int, ...] | None:
    """The suid for `index` in a grid with `N_side`; see ``RHEALPixDGGS.parse_index``."""
    _index_width(N_side)
    if not index or index[0] not in CELLS0:
        return None
    digits = index[1:]
    if not all(d in "0123456789" for d in digits):
        return None
    values = [int(d) for d in digits]
    if any(v >= N_side**2 for v in values):
        return None
    return (index[0], *values)


def _format_index(suid: Sequence[str | int], N_side: int) -> str:
    """The index string for `suid` in a grid with `N_side`; see ``RHEALPixDGGS.format_index``."""
    _index_width(N_side)
    return str(suid[0]) + "".join(str(d) for d in suid[1:])


def _parse_indices(
    indices: Iterable[str], N_side: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse cell index strings into arrays: a validity mask, the base cell
    as a code 0-5 (the position in ``CELLS0``), the digits as an
    ``(len(indices), max_resolution)`` int64 array (zero beyond a cell's
    resolution) and each index's resolution. Invalid indices -- empty,
    unknown base cell, or a character that is not a digit below
    ``N_side ** 2`` -- get resolution 0, base cell code -1 and zero digits.
    Raise ValueError if the grid has no index strings.
    """
    _index_width(N_side)
    strings = np.array(list(indices), dtype=str)
    count = len(strings)
    chars = max(strings.dtype.itemsize // 4, 2)
    codes = np.zeros((count, chars), dtype=np.uint32)
    if count:
        codes[:, : strings.dtype.itemsize // 4] = strings.view(np.uint32).reshape(
            count, -1
        )
    length = (codes != 0).sum(axis=1)
    face = np.full(count, -1, dtype=np.int64)
    for k, letter in enumerate(CELLS0):
        face[codes[:, 0] == ord(letter)] = k
    digits = codes[:, 1:].astype(np.int64) - ord("0")
    levels = np.arange(1, chars)
    active = levels[None, :] < length[:, None]
    digit_ok = ~active | ((digits >= 0) & (digits < N_side**2))
    valid = (length > 0) & (face >= 0) & digit_ok.all(axis=1)
    digits = np.where(active & valid[:, None], digits, 0)
    resolution = np.where(valid, length - 1, 0)
    face = np.where(valid, face, -1)
    return valid, face, digits, resolution


def _format_indices(
    face: np.ndarray, digits: np.ndarray, resolution: int, N_side: int
) -> np.ndarray:
    """
    The index strings of valid cells given as base cell codes and digit
    rows, all of resolution `resolution` (digit columns beyond it are
    ignored: ``_parse_planar_points`` always yields at least one).
    """
    _index_width(N_side)
    chars = np.empty((len(face), resolution + 1), dtype=np.uint32)
    chars[:, 0] = np.array([ord(c) for c in CELLS0])[face]
    # Column by column: a whole-array temporary would be int64 and, for
    # millions of cells, the largest allocation of the caller.
    for k in range(resolution):
        chars[:, k + 1] = digits[:, k] + ord("0")
    return chars.view(f"<U{resolution + 1}").ravel()


# Sub-zone orders ----------------------------------------------------------

# Each order maps a position to the (row, col) of a subcell in the n x n
# block of an ancestor's subcells, n = N**depth, rows counted from the top
# and columns from the left, and back.

#: The orders of a cell's subcells, by their names in OGC's DGGRS definition
#: schema; see ``Cell.subcells``.
SubZoneOrder = Literal["mortonCurve", "scanline", "hilbertCurve", "spiralFromCenter"]
#: The names in ``SubZoneOrder``.
SUB_ZONE_ORDERS: tuple[str, ...] = get_args(SubZoneOrder)


def _digits_to_rc(N: int, digits: Sequence[int]) -> tuple[int, int]:
    row = col = 0
    for d in digits:
        row, col = row * N + d // N, col * N + d % N
    return row, col


def _rc_to_digits(N: int, depth: int, row: int, col: int) -> list[int]:
    weights = [N ** (depth - 1 - i) for i in range(depth)]
    return [(row // w % N) * N + col // w % N for w in weights]


def _morton_at(N: int, depth: int, i: int) -> tuple[int, int]:
    digits = [i // (N * N) ** (depth - 1 - k) % (N * N) for k in range(depth)]
    return _digits_to_rc(N, digits)


def _morton_index(N: int, depth: int, row: int, col: int) -> int:
    index = 0
    for d in _rc_to_digits(N, depth, row, col):
        index = index * N * N + d
    return index


def _scanline_at(N: int, depth: int, i: int) -> tuple[int, int]:
    return divmod(i, N**depth)


def _scanline_index(N: int, depth: int, row: int, col: int) -> int:
    return row * N**depth + col


def _hilbert_at(N: int, depth: int, i: int) -> tuple[int, int]:
    # The iterative d-to-xy conversion, with x the column and y the row, so
    # the curve runs from the top left corner to the top right.
    n = N**depth
    x = y = 0
    s = 1
    while s < n:
        rx = 1 & (i // 2)
        ry = 1 & (i ^ rx)
        if ry == 0:
            if rx == 1:
                x, y = s - 1 - x, s - 1 - y
            x, y = y, x
        x += s * rx
        y += s * ry
        i //= 4
        s *= 2
    return y, x


def _hilbert_index(N: int, depth: int, row: int, col: int) -> int:
    n = N**depth
    x, y, index = col, row, 0
    s = n // 2
    while s > 0:
        rx = 1 if x & s else 0
        ry = 1 if y & s else 0
        index += s * s * ((3 * rx) ^ ry)
        if ry == 0:
            if rx == 1:
                x, y = n - 1 - x, n - 1 - y
            x, y = y, x
        s //= 2
    return index


# The spiral starts at the centre and goes right 1, down 1, left 2, up 2,
# right 3, ..., so ring k (cells k steps from the centre) starts at index
# (2k - 1)**2, k - 1 rows above the centre on its right side, and runs down
# its right side, left along its bottom, up its left side and right along
# its top, 2k cells each.
def _spiral_at(N: int, depth: int, i: int) -> tuple[int, int]:
    m = N**depth // 2
    if i == 0:
        return m, m
    k = (isqrt(i) + 1) // 2
    side, t = divmod(i - (2 * k - 1) ** 2, 2 * k)
    dr, dc = [
        (t - (k - 1), k),
        (k, k - 1 - t),
        (k - 1 - t, -k),
        (-k, t - (k - 1)),
    ][side]
    return m + dr, m + dc


def _spiral_index(N: int, depth: int, row: int, col: int) -> int:
    m = N**depth // 2
    dr, dc = row - m, col - m
    k = max(abs(dr), abs(dc))
    if k == 0:
        return 0
    if dc == k and dr > -k:
        offset = dr + k - 1
    elif dr == k:
        offset = 2 * k + k - 1 - dc
    elif dc == -k:
        offset = 4 * k + k - 1 - dr
    else:
        offset = 6 * k + dc + k - 1
    return (2 * k - 1) ** 2 + offset


_SUB_ZONE_AT = {
    "mortonCurve": _morton_at,
    "scanline": _scanline_at,
    "hilbertCurve": _hilbert_at,
    "spiralFromCenter": _spiral_at,
}
_SUB_ZONE_INDEX = {
    "mortonCurve": _morton_index,
    "scanline": _scanline_index,
    "hilbertCurve": _hilbert_index,
    "spiralFromCenter": _spiral_index,
}


def _check_sub_zone_order(order: str, N: int) -> None:
    if order not in SUB_ZONE_ORDERS:
        raise ValueError(
            f"order must be one of {', '.join(SUB_ZONE_ORDERS)}, not {order!r}"
        )
    if order == "hilbertCurve" and N & (N - 1):
        raise ValueError(f"hilbertCurve needs N_side a power of 2, not {N}")
    if order == "spiralFromCenter" and N % 2 == 0:
        raise ValueError(f"spiralFromCenter needs odd N_side, not {N}")


# Positions, uint64 and level-row-col text ---------------------------------

# The uint64 zone identifier of OGC's 2D tile matrix set form: level in the
# top 5 bits, row in the next 29, column in the low 30.
_UINT64_LEVEL_SHIFT = 59
_UINT64_ROW_SHIFT = 30
_UINT64_ROW_LIMIT = 1 << 29
_UINT64_COL_LIMIT = 1 << 30

# Cell.level_row_col_text: level, row and column in hexadecimal.
_LEVEL_ROW_COL_TEXT = re.compile(r"([0-9A-Fa-f]+)-([0-9A-Fa-f]+)-([0-9A-Fa-f]+)")

_Ints = Iterable[int] | np.ndarray


def _face_blocks(north_square: int, south_square: int) -> list[tuple[int, int]]:
    """(row, col) of each resolution 0 cell, in ``CELLS0`` order, in the 3 x 4
    frame of the planar grid."""
    return [(0, north_square), (1, 0), (1, 1), (1, 2), (1, 3), (2, south_square)]


def _level_row_col(
    grid: _GridLayout, suid: Sequence[str | int]
) -> tuple[int, int, int]:
    """The level, row and column of the cell with `suid`; see
    ``Cell.level_row_col``."""
    N = grid.N_side
    L = len(suid) - 1
    p = N**L
    block_row, block_col = _face_blocks(grid.north_square, grid.south_square)[
        CELLS0.index(str(suid[0]))
    ]
    row, col = _digits_to_rc(N, [int(d) for d in suid[1:]])
    return L, block_row * p + row, block_col * p + col


def _suid_from_level_row_col(
    grid: _GridLayout, level: int, row: int, col: int
) -> list[str | int]:
    """The suid of the cell at a position; see
    ``RHEALPixDGGS.cell_from_level_row_col``."""
    if not 0 <= level <= grid.max_resolution:
        raise ValueError(f"level must be from 0 to {grid.max_resolution}, not {level}")
    p = grid.N_side**level
    blocks = _face_blocks(grid.north_square, grid.south_square)
    if not (0 <= row < 3 * p and 0 <= col < 4 * p) or (
        (row // p, col // p) not in blocks
    ):
        raise ValueError(f"no level {level} cell at row {row}, column {col}")
    face = CELLS0[blocks.index((row // p, col // p))]
    return [face, *_rc_to_digits(grid.N_side, level, row % p, col % p)]


def _fits_uint64(row: int, col: int) -> bool:
    return row < _UINT64_ROW_LIMIT and col < _UINT64_COL_LIMIT


def _pack_uint64(level: int, row: int, col: int) -> int:
    return level << _UINT64_LEVEL_SHIFT | row << _UINT64_ROW_SHIFT | col


def _unpack_uint64(value: int) -> tuple[int, int, int]:
    if not 0 <= value < 1 << 64:
        raise ValueError(f"{value} is not a uint64")
    return (
        value >> _UINT64_LEVEL_SHIFT,
        value >> _UINT64_ROW_SHIFT & (_UINT64_ROW_LIMIT - 1),
        value & (_UINT64_COL_LIMIT - 1),
    )


def _level_row_col_text(level: int, row: int, col: int) -> str:
    return f"{level:X}-{row:X}-{col:X}"


def _parse_level_row_col_text(text: str) -> tuple[int, int, int]:
    match = _LEVEL_ROW_COL_TEXT.fullmatch(text)
    if match is None:
        raise ValueError(f"{text!r} is not level-row-col in hexadecimal")
    level, row, col = (int(v, 16) for v in match.groups())
    return level, row, col


def _level_row_col_from_indices(
    grid: _GridLayout, indices: Iterable[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """See ``RHEALPixDGGS.level_row_col_from_indices``."""
    strings = list(indices)
    valid, face, digits, resolution = _parse_indices(strings, grid.N_side)
    if not valid.all():
        raise ValueError(f"invalid index {strings[int(np.argmin(valid))]!r}")
    N = grid.N_side
    level = resolution.astype(np.int64)
    row = np.zeros(len(strings), dtype=np.int64)
    col = np.zeros(len(strings), dtype=np.int64)
    for k in range(digits.shape[1]):
        active = k < level
        row = np.where(active, row * N + digits[:, k] // N, row)
        col = np.where(active, col * N + digits[:, k] % N, col)
    blocks = np.array(_face_blocks(grid.north_square, grid.south_square))
    p = N**level
    return level, row + blocks[face, 0] * p, col + blocks[face, 1] * p


def _indices_from_level_row_col(
    grid: _GridLayout, level: _Ints, row: _Ints, col: _Ints
) -> list[str]:
    """See ``RHEALPixDGGS.indices_from_level_row_col``."""
    _index_width(grid.N_side)
    return _format_level_row_col(grid, *_checked_level_row_col(grid, level, row, col))


def _uint64_from_level_row_col(
    grid: _GridLayout, level: _Ints, row: _Ints, col: _Ints
) -> np.ndarray:
    """See ``RHEALPixDGGS.uint64_from_level_row_col``."""
    _, level_, row_, col_ = _checked_level_row_col(grid, level, row, col)
    if (row_ >= _UINT64_ROW_LIMIT).any() or (col_ >= _UINT64_COL_LIMIT).any():
        raise ValueError("cells too fine for a uint64 identifier")
    return (
        level_.astype(np.uint64) << np.uint64(_UINT64_LEVEL_SHIFT)
        | row_.astype(np.uint64) << np.uint64(_UINT64_ROW_SHIFT)
        | col_.astype(np.uint64)
    )


def _level_row_col_from_uint64(
    grid: _GridLayout, values: _Ints
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """See ``RHEALPixDGGS.level_row_col_from_uint64``."""
    z = np.asarray(values)
    if z.size and (z.dtype.kind not in "iu" or (z < 0).any()):
        raise ValueError("values must be non-negative integers")
    z = z.astype(np.uint64).ravel()
    level = (z >> np.uint64(_UINT64_LEVEL_SHIFT)).astype(np.int64)
    row = z >> np.uint64(_UINT64_ROW_SHIFT) & np.uint64(_UINT64_ROW_LIMIT - 1)
    col = z & np.uint64(_UINT64_COL_LIMIT - 1)
    row, col = row.astype(np.int64), col.astype(np.int64)
    face = _face_of_level_row_col(grid, level, row, col)
    if (face < 0).any():
        raise ValueError(f"{int(z[np.argmax(face < 0)])} names no cell of this grid")
    return level, row, col


def _level_row_col_text_from_indices(
    grid: _GridLayout, indices: Iterable[str]
) -> list[str]:
    """See ``RHEALPixDGGS.level_row_col_text_from_indices``."""
    level, row, col = _level_row_col_from_indices(grid, indices)
    return [
        _level_row_col_text(L, r, c)
        for L, r, c in zip(level.tolist(), row.tolist(), col.tolist(), strict=True)
    ]


def _indices_from_level_row_col_text(
    grid: _GridLayout, texts: Iterable[str]
) -> list[str]:
    """See ``RHEALPixDGGS.indices_from_level_row_col_text``."""
    _index_width(grid.N_side)
    strings = list(texts)
    positions = np.zeros((len(strings), 3), dtype=np.int64)
    for k, text in enumerate(strings):
        match = _LEVEL_ROW_COL_TEXT.fullmatch(text) if isinstance(text, str) else None
        values = [int(v, 16) for v in match.groups()] if match else []
        # Anything this large names no cell, and would overflow int64.
        if not values or max(values) >= 1 << 62:
            raise ValueError(f"{text!r} is not the level-row-col of a cell")
        positions[k] = values
    level, row, col = positions.T
    face = _face_of_level_row_col(grid, level, row, col)
    if (face < 0).any():
        text = strings[int(np.argmax(face < 0))]
        raise ValueError(f"{text!r} names no cell of this grid")
    return _format_level_row_col(grid, face, level, row, col)


def _checked_level_row_col(
    grid: _GridLayout, level: _Ints, row: _Ints, col: _Ints
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Face codes and flat int64 level, row and column arrays of positions
    broadcast together; raise ValueError if any is not integral or names no
    cell of the grid."""
    arrays = [np.asarray(a) for a in (level, row, col)]
    if any(a.size and a.dtype.kind not in "iu" for a in arrays):
        raise ValueError("level, row and column must be integers")
    level_, row_, col_ = (
        a.astype(np.int64).ravel() for a in np.broadcast_arrays(*arrays)
    )
    face = _face_of_level_row_col(grid, level_, row_, col_)
    if (face < 0).any():
        k = int(np.argmax(face < 0))
        raise ValueError(
            f"no level {level_[k]} cell at row {row_[k]}, column {col_[k]}"
        )
    return face, level_, row_, col_


def _face_of_level_row_col(
    grid: _GridLayout, level: np.ndarray, row: np.ndarray, col: np.ndarray
) -> np.ndarray:
    """The face code (position in ``CELLS0``) of each position, -1 where it
    names no cell: a level outside 0 to ``max_resolution``, or a row and
    column outside the frame or in one of its empty blocks."""
    ok = (level >= 0) & (level <= grid.max_resolution) & (row >= 0) & (col >= 0)
    p = grid.N_side ** np.where(ok, level, 0)
    ok &= (row < 3 * p) & (col < 4 * p)
    face_of_block = np.full((3, 4), -1)
    for k, (r, c) in enumerate(_face_blocks(grid.north_square, grid.south_square)):
        face_of_block[r, c] = k
    block = face_of_block[np.where(ok, row // p, 0), np.where(ok, col // p, 0)]
    return np.where(ok, block, -1)


def _format_level_row_col(
    grid: _GridLayout,
    face: np.ndarray,
    level: np.ndarray,
    row: np.ndarray,
    col: np.ndarray,
) -> list[str]:
    """The index strings of valid positions with face codes `face`."""
    N = grid.N_side
    result = np.empty(len(face), dtype=object)
    for L in np.unique(level):
        group = level == L
        p = N**L
        weights = N ** np.arange(L - 1, -1, -1, dtype=np.int64)
        r = (row[group] % p)[:, None] // weights % N
        c = (col[group] % p)[:, None] // weights % N
        result[group] = _format_indices(face[group], r * N + c, int(L), N)
    return result.tolist()


# Grid distance ------------------------------------------------------------

_STEPS = {"up": (-1, 0), "right": (0, 1), "down": (1, 0), "left": (0, -1)}

# Quarter turns of a polar face's rows and columns seen from the face across
# each of its edges: the net attaches it unturned to one band face, and each
# band face further round meets it a quarter turn on.
_POLAR_TURNS = {
    "N": {"down": 0, "right": 3, "up": 2, "left": 1},
    "S": {"up": 0, "left": 3, "down": 2, "right": 1},
}


def _face_neighbors(north_square: int, south_square: int) -> dict[str, dict[str, str]]:
    """The face across each edge of each face, as ``RHEALPixDGGS.atomic_neighbors``."""
    N, O, P, Q, R, S = CELLS0
    band = [O, P, Q, R]
    faces = {
        face: {"left": band[i - 1], "right": band[(i + 1) % 4], "up": N, "down": S}
        for i, face in enumerate(band)
    }
    for pole, square, around in (
        (N, north_square, ("down", "right", "up", "left")),
        (S, south_square, ("up", "right", "down", "left")),
    ):
        faces[pole] = {d: band[(square + k) % 4] for k, d in enumerate(around)}
    return faces


def _rotate_rc(turns: int, row: int, col: int, p: int) -> tuple[int, int]:
    """(row, col) of a p x p block after `turns` clockwise quarter turns."""
    for _ in range(turns % 4):
        row, col = col, p - 1 - row
    return row, col


@functools.cache
def _developments(
    north_square: int, south_square: int
) -> dict[tuple[str, str], tuple[tuple[int, int, int], ...]]:
    """
    For each ordered pair of faces, where the second lies in the plane of
    the first when the cube is unfolded along each simple path of at most
    three edges between them: the clockwise quarter turns of its rows and
    columns, and its row and column offsets in face widths.
    """
    faces = _face_neighbors(north_square, south_square)

    def crossing(face: str, direction: str) -> tuple[int, int, int]:
        across = faces[face][direction]
        if face in _POLAR_TURNS:
            turns = _POLAR_TURNS[face][direction]
        elif across in _POLAR_TURNS:
            back = next(d for d, f in faces[across].items() if f == face)
            turns = -_POLAR_TURNS[across][back] % 4
        else:
            turns = 0
        return turns, *_STEPS[direction]

    result: dict[tuple[str, str], list[tuple[int, int, int]]] = {
        (face, face): [(0, 0, 0)] for face in CELLS0
    }

    def extend(start: str, path: list[str], turns: int, dr: int, dc: int) -> None:
        face = path[-1]
        for direction in _STEPS:
            across = faces[face][direction]
            if across in path:
                continue
            t, ur, uc = crossing(face, direction)
            # Compose: the offset turns with the faces before it.
            for _ in range(turns):
                ur, uc = uc, -ur
            placed = ((turns + t) % 4, dr + ur, dc + uc)
            result.setdefault((start, across), []).append(placed)
            if len(path) < 3:
                extend(start, [*path, across], *placed)

    for face in CELLS0:
        extend(face, [face], 0, 0, 0)
    return {pair: tuple(placements) for pair, placements in result.items()}


def _grid_distance(
    grid: _GridLayout, a: Sequence[str | int], b: Sequence[str | int], diagonals: bool
) -> int:
    """
    Hops between the cells with suids `a` and `b` at the coarser one's
    resolution; see ``Cell.grid_distance``. A path unfolds into the plane
    along the faces it crosses, where it is no shorter than the planar hop
    count between its ends, and a shortest path crosses at most three edges.
    """
    N = grid.N_side
    level = min(len(a), len(b)) - 1
    p = N**level
    ra, ca = _digits_to_rc(N, [int(d) for d in a[1 : level + 1]])
    rb, cb = _digits_to_rc(N, [int(d) for d in b[1 : level + 1]])
    placements = _developments(grid.north_square, grid.south_square)
    best = None
    for turns, ur, uc in placements[str(a[0]), str(b[0])]:
        r, c = _rotate_rc(turns, rb, cb, p)
        dr, dc = abs(r + ur * p - ra), abs(c + uc * p - ca)
        hops = max(dr, dc) if diagonals else dr + dc
        best = hops if best is None else min(best, hops)
    assert best is not None
    return best
