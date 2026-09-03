"""
Time the projection and the cell geometry built on it.

Usage, comparing a branch against master::

    git worktree add ../rhp-master master
    python ../rhp-master/scripts/benchmark.py master
    python scripts/benchmark.py branch
    git worktree remove ../rhp-master

The script measures the checkout it lives in, whatever package is
installed. Each run prints one Markdown table; paste the two side by
side. Timings
are the best of three repeats. ``--res 4`` adds the 39,366-cell grid
(tens of seconds per checkout); ``--count`` reports projected points
instead of time, which is deterministic and machine-independent.

Not part of the test suite: timings on shared CI runners are noise.
"""

import argparse
import sys
import timeit
from collections.abc import Callable
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rhealpixdggs.dggs import WGS84_003


class PointCounter:
    """Stand-in for RHEALPixDGGS.rhealpix that counts projected points."""

    def __init__(self, inner: Callable) -> None:
        self.inner = inner
        self.points = 0

    def __call__(self, *args, **kwargs):
        self.points += int(np.size(args[0]))
        return self.inner(*args, **kwargs)


def best_of(fn: Callable[[], object], number: int, repeat: int = 3) -> float:
    return min(timeit.repeat(fn, number=number, repeat=repeat)) / number


def workloads(rdggs, resolutions):
    quad, dart, skew = rdggs.cell(("P", 4)), rdggs.cell(("N", 0)), rdggs.cell(("N", 1))
    x_eq, y_eq = rdggs.rhealpix(20.0, 10.0)
    x_po, y_po = rdggs.rhealpix(20.0, 80.0)
    rows = [
        ("forward projection, equatorial", lambda: rdggs.rhealpix(20.0, 10.0), 2000),
        ("forward projection, polar", lambda: rdggs.rhealpix(20.0, 80.0), 2000),
        (
            "inverse projection, equatorial",
            lambda: rdggs.rhealpix(x_eq, y_eq, inverse=True),
            2000,
        ),
        (
            "inverse projection, polar",
            lambda: rdggs.rhealpix(x_po, y_po, inverse=True),
            2000,
        ),
        (
            "quad boundary(n=10, plane=False)",
            lambda: quad.boundary(n=10, plane=False),
            200,
        ),
        (
            "dart boundary(n=10, plane=False)",
            lambda: dart.boundary(n=10, plane=False),
            100,
        ),
        (
            "skew boundary(n=10, plane=False)",
            lambda: skew.boundary(n=10, plane=False),
            100,
        ),
        ("quad vertices(plane=False)", lambda: quad.vertices(plane=False), 500),
        (
            "dart interior(n=5, plane=False)",
            lambda: dart.interior(n=5, plane=False),
            100,
        ),
        ("quad centroid(plane=False)", lambda: quad.centroid(plane=False), 100),
        ("skew centroid(plane=False)", lambda: skew.centroid(plane=False), 5),
        ("dart centroid(plane=False)", lambda: dart.centroid(plane=False), 1),
        (
            "cell_from_point(plane=False)",
            lambda: rdggs.cell_from_point(3, (20.0, 80.0), plane=False),
            500,
        ),
    ]
    for res in resolutions:
        grid = list(rdggs.grid(res))
        n = 3 if res <= 2 else 1
        rows.append(
            (
                f"res {res}, {len(grid)} cells, per-cell boundary(n=10)",
                lambda grid=grid: [c.boundary(n=10, plane=False) for c in grid],
                n,
            )
        )
        rows.append(
            (
                f"res {res}, {len(grid)} cells, cell_boundaries(n=10)",
                lambda grid=grid: rdggs.cell_boundaries(grid, n=10, plane=False),
                n,
            )
        )
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("label", nargs="?", default="this checkout")
    parser.add_argument(
        "--res", type=int, nargs="*", default=[2, 3], help="grid resolutions"
    )
    parser.add_argument(
        "--count", action="store_true", help="count projected points, not time"
    )
    args = parser.parse_args(argv)
    rdggs = WGS84_003
    rows = workloads(rdggs, args.res)
    unit = "points" if args.count else "time"
    print(f"| workload | {args.label} ({unit}) |")
    print("|---|---|")
    for name, fn, number in rows:
        if args.count:
            counter = PointCounter(rdggs.rhealpix)
            rdggs.rhealpix = counter
            try:
                fn()
            finally:
                del rdggs.__dict__["rhealpix"]
            print(f"| {name} | {counter.points} |")
            continue
        seconds = best_of(fn, number)
        if seconds >= 0.1:
            text = f"{seconds:.2f} s"
        elif seconds >= 1e-3:
            text = f"{seconds * 1e3:.1f} ms"
        else:
            text = f"{seconds * 1e6:.1f} us"
        print(f"| {name} | {text} |", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
