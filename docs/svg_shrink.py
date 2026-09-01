"""
Shrink matplotlib's SVG output by dropping meaningless coordinate precision.

matplotlib writes path coordinates with six decimal places. The figures in
docs/source/images/ are a few hundred points across, so the sixth decimal of
a point is around a nanometre on the page -- far below anything any renderer
or printer can express, but it costs four bytes on every one of the hundreds
of thousands of numbers in a coastline or a grid of cell boundaries. Rounding
to two decimals (a hundredth of a point, ~1/25 of a pixel at 300 dpi) removes
roughly a third of the bytes with no visible change.

Only the ``d`` attribute of ``<path>`` elements drawn in the figure's own
coordinate system is touched. Paths inside ``<defs>`` are left alone: those
are text glyph outlines, which matplotlib emits at 64x and then scales down
with ``transform="scale(0.015625)"``, so a hundredth there is not a
hundredth on the page. Paths carrying their own transform, or sitting inside
a scaled group, are skipped for the same reason. Nothing else in the file --
styles, transforms, ids, metadata -- is rewritten.

Used by make_figures.py after it writes the figures, and runnable directly
to re-shrink files that are already committed::

    python docs/svg_shrink.py docs/source/images/*.svg

The transformation is idempotent: running it twice changes nothing the
second time.
"""

import pathlib
import re
import sys

# <defs>, <g ...>, </g> and <path .../>. Attribute values in matplotlib's
# output never contain '>', so matching up to the first '>' is safe here.
_TAG = re.compile(r"<(/?)(defs|g|path)\b([^>]*)>", re.S)
_D_ATTR = re.compile(r'(\sd=")([^"]*)(")', re.S)
_NUMBER = re.compile(r"-?\d+\.\d+")

DEFAULT_DECIMALS = 2


def _round_number(match: re.Match, decimals: int) -> str:
    """Round one decimal number, trimming the trailing zeros that creates."""
    text = f"{float(match.group(0)):.{decimals}f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    # Rounding can turn a small negative into "-0"; normalise it.
    return "0" if text in ("-0", "-") else text


def _round_path_data(data: str, decimals: int) -> str:
    return _NUMBER.sub(lambda m: _round_number(m, decimals), data)


def shrink(svg: str, decimals: int = DEFAULT_DECIMALS) -> str:
    """Return ``svg`` with unscaled path coordinates rounded to ``decimals``."""
    out = []
    end_of_previous = 0
    defs_depth = 0
    # One entry per open <g>: whether it, or anything containing it, scales.
    scale_stack = []

    for match in _TAG.finditer(svg):
        closing, tag, attributes = match.group(1), match.group(2), match.group(3)
        self_closing = attributes.rstrip().endswith("/")

        if tag == "defs":
            defs_depth += -1 if closing else 1
            continue

        if tag == "g":
            if closing:
                if scale_stack:
                    scale_stack.pop()
            elif not self_closing:
                inherited = scale_stack[-1] if scale_stack else False
                scale_stack.append(inherited or "scale(" in attributes)
            continue

        # tag == "path"
        if closing or defs_depth or (scale_stack and scale_stack[-1]):
            continue
        if "transform=" in attributes:
            continue

        rounded, count = _D_ATTR.subn(
            lambda m: m.group(1) + _round_path_data(m.group(2), decimals) + m.group(3),
            attributes,
        )
        if not count:
            continue
        out.append(svg[end_of_previous : match.start()])
        out.append(f"<{closing}{tag}{rounded}>")
        end_of_previous = match.end()

    out.append(svg[end_of_previous:])
    return "".join(out)


def shrink_file(path: pathlib.Path, decimals: int = DEFAULT_DECIMALS) -> tuple[int, int]:
    """Shrink ``path`` in place. Returns the byte count before and after."""
    original = path.read_text(encoding="utf-8")
    shrunk = shrink(original, decimals)
    if shrunk != original:
        path.write_text(shrunk, encoding="utf-8")
    return len(original.encode("utf-8")), len(shrunk.encode("utf-8"))


def main(argv: list[str]) -> int:
    paths = [pathlib.Path(a) for a in argv]
    if not paths:
        print(__doc__.strip().splitlines()[0], file=sys.stderr)
        print(f"usage: {sys.argv[0]} FILE.svg ...", file=sys.stderr)
        return 2

    before_total = after_total = 0
    for path in paths:
        before, after = shrink_file(path)
        before_total += before
        after_total += after
        saved = 100 * (before - after) / before if before else 0.0
        print(f"{path.name}: {before:>9,} -> {after:>9,} bytes ({saved:4.1f}% smaller)")

    if len(paths) > 1 and before_total:
        saved = 100 * (before_total - after_total) / before_total
        print(
            f"total: {before_total:,} -> {after_total:,} bytes ({saved:.1f}% smaller)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
