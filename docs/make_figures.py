"""
Regenerate the static figures in docs/source/images/ (each written as
SVG for the HTML docs and PDF for the LaTeX/PDF docs; the figure
directives reference them with a wildcard so each Sphinx builder picks
the format it can use).

Run manually from the repository root whenever the figures need to change:

    python docs/make_figures.py

Requires matplotlib, which is deliberately NOT a dependency of the
package or of the docs build -- the generated SVGs are committed, and
regenerating them on every build would only add churn for figures that
only change when the grid geometry or the figure design does. It is
declared (and version-locked, keeping the SVG output consistent across
contributors) in the optional "figures" dependency group:

    poetry install --with figures

The SVGs are post-processed by svg_shrink (see that module), which strips
matplotlib's surplus coordinate precision; the committed figures are stored
in that form, so regenerating an unchanged figure still produces no diff.

All grid geometry is drawn from the library itself (Cell.ul_vertex,
Cell.boundary, Cell.nucleus). Coastlines are Natural Earth 1:110m data
(public domain), downloaded on first run (pinned to a specific upstream
commit for reproducibility) and cached under .cache/ at the repository
root (outside docs/, which is shipped in the sdist).
"""

import itertools
import json
import os
import pathlib
import sys
import urllib.request

# Make the SVG output deterministic, so regenerating unchanged figures
# produces no git diff: pin the embedded creation date and the salt
# matplotlib hashes into its clip-path/element ids.
os.environ.setdefault("SOURCE_DATE_EPOCH", "0")

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "rhealpixdggs"
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.patches import Polygon as PolygonPatch
from matplotlib.patches import Rectangle

# Draw from this repository's rhealpixdggs, not any installed copy, and pick
# up svg_shrink from alongside this file however the script was invoked.
sys.path.insert(0, str(pathlib.Path(__file__).parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).parent))


import svg_shrink

from rhealpixdggs.dggs import CELLS0, WGS84_003

OUT = pathlib.Path(__file__).parent / "source" / "images"

rdggs = WGS84_003
R = rdggs.ellipsoid.R_A

FACE_COLORS = {
    "N": "#e5735c",  # warm red
    "O": "#e8c34f",  # yellow
    "P": "#79bf6f",  # green
    "Q": "#5cc3c9",  # cyan
    "R": "#7c8fe0",  # blue-violet
    "S": "#d97cc0",  # magenta
}

COAST_COLOR = "#9a9a9a"

# Natural Earth 1:110m coastline (public domain), pinned to a specific
# upstream commit so regeneration stays reproducible.
COASTLINE_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
    "ca96624a56bd078437bca8184e78163e5039ad19/geojson/ne_110m_coastline.geojson"
)
COASTLINE_CACHE = (
    pathlib.Path(__file__).parents[1] / ".cache" / "ne_110m_coastline.geojson"
)


def coastline_segments():
    """
    Return the Natural Earth 1:110m coastlines as a list of segments,
    each a list of (longitude, latitude) pairs, split wherever a segment
    crosses the antimeridian.
    """
    if not COASTLINE_CACHE.exists():
        COASTLINE_CACHE.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(COASTLINE_URL, timeout=60) as response:
            COASTLINE_CACHE.write_bytes(response.read())
    collection = json.loads(COASTLINE_CACHE.read_text())
    segments = []
    for feature in collection["features"]:
        geometry = feature["geometry"]
        lines = (
            [geometry["coordinates"]]
            if geometry["type"] == "LineString"
            else geometry["coordinates"]
        )
        for line in lines:
            current = [line[0]]
            for prev, cur in itertools.pairwise(line):
                if abs(cur[0] - prev[0]) > 180:
                    segments.append(current)
                    current = []
                current.append(cur)
            segments.append(current)
    return [s for s in segments if len(s) > 1]


COASTLINES = coastline_segments()


# ---------------------------------------------------------------- figure 1
# Planar grid: the (0,0)-rHEALPix cube unfolding, resolution 0 faces
# labeled, resolution 1 sub-grid drawn, and one face's children labeled
# with their SUID digits.
fig, ax = plt.subplots(figsize=(8, 5))
for face in CELLS0:
    c0 = rdggs.cell([face])
    x, y = c0.ul_vertex()
    w = c0.width()
    ax.add_patch(
        Rectangle(
            (x / R, (y - w) / R),
            w / R,
            w / R,
            facecolor=FACE_COLORS[face],
            alpha=0.35,
            edgecolor="black",
            linewidth=1.8,
        )
    )
    ax.text(
        (x + 0.5 * w) / R,
        (y - 0.5 * w) / R,
        face,
        ha="center",
        va="center",
        fontsize=26,
        fontweight="bold",
        color="#333333",
        alpha=0.85,
    )
    # Resolution 1 sub-grid.
    for child in c0.subcells():
        cx, cy = child.ul_vertex()
        cw = child.width()
        ax.add_patch(
            Rectangle(
                (cx / R, (cy - cw) / R),
                cw / R,
                cw / R,
                facecolor="none",
                edgecolor="black",
                linewidth=0.4,
                alpha=0.6,
            )
        )
# Coastlines, projected through the DGGS's own rHEALPix projection.
# The polar-triangle rearrangement makes planar positions jump where a
# segment crosses between regions, so split on planar jumps far larger
# than any real 1:110m coastline step (~0.02 R at most).
for seg in COASTLINES:
    xy = [rdggs.rhealpix(lon, lat) for lon, lat in seg]
    run = [xy[0]]
    for prev, cur in itertools.pairwise(xy):
        if abs(cur[0] - prev[0]) > 0.15 * R or abs(cur[1] - prev[1]) > 0.15 * R:
            if len(run) > 1:
                ax.plot(
                    [p[0] / R for p in run],
                    [p[1] / R for p in run],
                    color=COAST_COLOR,
                    linewidth=0.5,
                    zorder=1,
                )
            run = []
        run.append(cur)
    if len(run) > 1:
        ax.plot(
            [p[0] / R for p in run],
            [p[1] / R for p in run],
            color=COAST_COLOR,
            linewidth=0.5,
            zorder=1,
        )
# Label P's children with their digits, teaching the SUID scheme.
for child in rdggs.cell(["P"]).subcells():
    cx, cy = child.ul_vertex()
    cw = child.width()
    ax.text(
        (cx + 0.15 * cw) / R,
        (cy - 0.22 * cw) / R,
        str(child.suid[1]),
        ha="center",
        va="center",
        fontsize=9,
        color="#222222",
    )
ax.set_xlim(-3.4, 3.4)
ax.set_ylim(-2.6, 2.6)
ax.set_aspect("equal")
ax.set_xlabel("x / authalic radius")
ax.set_ylabel("y / authalic radius")
ax.set_title("(0, 0)-rHEALPix planar grid: resolution 0 cells, resolution 1 sub-grid")
fig.tight_layout()
fig.savefig(OUT / "planar_grid.svg", bbox_inches="tight")
fig.savefig(OUT / "planar_grid.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- figure 2
# Ellipsoidal cells: every resolution 1 cell's boundary in longitude-
# latitude coordinates, colored by resolution 0 face; shows the quad,
# skew-quad, dart, and cap shapes.
fig, ax = plt.subplots(figsize=(9, 4.8))


def split_chart_discontinuities(points):
    """Split a closed boundary point list into segments wherever it
    jumps across the antimeridian or -- possible for grids recentred
    off the equator, whose cells can overflow past a pole -- wherever
    the latitude wraps from one pole to the other."""
    segs = [[points[0]]]
    for prev, cur in itertools.pairwise(points):
        if abs(cur[0] - prev[0]) > 180 or abs(cur[1] - prev[1]) > 90:
            segs.append([])
        segs[-1].append(cur)
    return segs


def fill_lonlat_polygon(ax, points, **kwargs):
    """Fill a closed lon-lat polygon on a plate carree axis, keeping the
    fill intact for polygons that straddle the antimeridian: unwrap the
    longitudes into one continuous ring, then draw it at every 360-degree
    offset that intersects the axis window and let clipping do the rest."""
    lons = [points[0][0]]
    for prev, cur in itertools.pairwise(points):
        d = cur[0] - prev[0]
        if d > 180:
            d -= 360
        elif d < -180:
            d += 360
        lons.append(lons[-1] + d)
    lats = [p[1] for p in points]
    offsets = {0.0}
    if max(lons) > 180:
        offsets.add(-360.0)
    if min(lons) < -180:
        offsets.add(360.0)
    for off in offsets:
        ax.fill([lon + off for lon in lons], lats, **kwargs)


def draw_coastlines_lonlat(ax, linewidth=0.5):
    for seg in COASTLINES:
        ax.plot(*zip(*seg), color=COAST_COLOR, linewidth=linewidth, zorder=1)


draw_coastlines_lonlat(ax)
for face in CELLS0:
    color = FACE_COLORS[face]
    for cell in rdggs.cell([face]).subcells():
        pts = cell.boundary(n=40, plane=False)
        pts = pts + [pts[0]]
        for seg in split_chart_discontinuities(pts):
            if len(seg) > 1:
                ax.plot(*zip(*seg), color=color, linewidth=0.9)
    # Face label at the resolution 0 nucleus.
    lon, lat = rdggs.cell([face]).nucleus(plane=False)
    lat = max(min(lat, 78), -78)  # nudge pole labels inward
    ax.text(
        lon,
        lat,
        face,
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold",
        color=color,
    )
ax.set_xlim(-180, 180)
ax.set_ylim(-90, 90)
ax.set_xticks(range(-180, 181, 60))
ax.set_yticks(range(-90, 91, 30))
ax.set_xlabel("longitude (degrees)")
ax.set_ylabel("latitude (degrees)")
ax.set_title("Resolution 1 ellipsoidal cells of the (0, 0)-rHEALPix DGGS (WGS84)")
ax.grid(True, linewidth=0.3, alpha=0.5)
fig.tight_layout()
fig.savefig(OUT / "ellipsoidal_cells.svg", bbox_inches="tight")
fig.savefig(OUT / "ellipsoidal_cells.pdf", bbox_inches="tight")
plt.close(fig)
print("figures written")


# ---------------------------------------------------------------- figure 3
# Globe views.
def ortho(lon_deg, lat_deg, lon0_deg, lat0_deg):
    """Orthographic projection centered on (lon0, lat0). Returns x, y and
    a visibility mask (front hemisphere)."""
    lon = np.radians(np.asarray(lon_deg))
    lat = np.radians(np.asarray(lat_deg))
    lon0 = np.radians(lon0_deg)
    lat0 = np.radians(lat0_deg)
    cosc = np.sin(lat0) * np.sin(lat) + np.cos(lat0) * np.cos(lat) * np.cos(lon - lon0)
    x = np.cos(lat) * np.sin(lon - lon0)
    y = np.cos(lat0) * np.sin(lat) - np.sin(lat0) * np.cos(lat) * np.cos(lon - lon0)
    return x, y, cosc >= 0


def draw_globe(ax, lon0, lat0, title):
    # Globe outline.
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(t), np.sin(t), color="#555555", linewidth=1.2)
    # Coastlines on the front hemisphere.
    for seg in COASTLINES:
        lons = [p[0] for p in seg]
        lats = [p[1] for p in seg]
        x, y, vis = ortho(lons, lats, lon0, lat0)
        x, y = np.where(vis, x, np.nan), np.where(vis, y, np.nan)
        ax.plot(x, y, color=COAST_COLOR, linewidth=0.5, zorder=1)
    # Light graticule for orientation.
    for glat in range(-60, 90, 30):
        lons = np.linspace(-180, 180, 361)
        x, y, vis = ortho(lons, np.full_like(lons, glat), lon0, lat0)
        x, y = np.where(vis, x, np.nan), np.where(vis, y, np.nan)
        ax.plot(x, y, color="#bbbbbb", linewidth=0.3)
    for glon in range(-180, 180, 30):
        lats = np.linspace(-90, 90, 181)
        x, y, vis = ortho(np.full_like(lats, glon), lats, lon0, lat0)
        x, y = np.where(vis, x, np.nan), np.where(vis, y, np.nan)
        ax.plot(x, y, color="#bbbbbb", linewidth=0.3)
    # Resolution 1 cells, colored by resolution 0 face: boundary lines
    # always; a subtle translucent fill for cells lying entirely on the
    # front hemisphere (filling horizon-clipped cells would need the
    # horizon arc stitched into the polygon, which isn't worth it for a
    # static illustration); a small SUID label for cells facing the
    # viewer squarely enough to read.
    for face in CELLS0:
        color = FACE_COLORS[face]
        for cell in rdggs.cell([face]).subcells():
            pts = cell.boundary(n=60, plane=False)
            pts = pts + [pts[0]]
            lons = [p[0] for p in pts]
            lats = [p[1] for p in pts]
            x, y, vis = ortho(lons, lats, lon0, lat0)
            if vis.all():
                ax.fill(x, y, color=color, alpha=0.18, linewidth=0)
            x, y = np.where(vis, x, np.nan), np.where(vis, y, np.nan)
            ax.plot(x, y, color=color, linewidth=1.0)
            nx, ny, nvis = ortho(*cell.nucleus(plane=False), lon0, lat0)
            ncos = np.sin(np.radians(lat0)) * np.sin(
                np.radians(cell.nucleus(plane=False)[1])
            ) + np.cos(np.radians(lat0)) * np.cos(
                np.radians(cell.nucleus(plane=False)[1])
            ) * np.cos(
                np.radians(cell.nucleus(plane=False)[0] - lon0)
            )
            if nvis and ncos > 0.45:
                ax.text(
                    nx,
                    ny,
                    str(cell),
                    ha="center",
                    va="center",
                    fontsize=5.5,
                    color="#333333",
                )
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=11)


fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.8))
draw_globe(axes[0], -45, 20, "P face (equatorial view)")
draw_globe(axes[1], 0, 90, "North polar view (cap cell N4)")
draw_globe(axes[2], 100, -35, "Oblique southern view")
fig.tight_layout()
fig.savefig(OUT / "globe_views.svg", bbox_inches="tight")
fig.savefig(OUT / "globe_views.pdf", bbox_inches="tight")
print("globe views written")


# ---------------------------------------------------------------- figure 4
# The H3-style wrappers in action over New Zealand: polyfill (cells whose
# centroids fall inside a polygon) and linetrace (cells touched by a
# linestring).

from shapely.geometry import LineString, Polygon

from rhealpixdggs import rhp_wrappers


def cell_from_address(address):
    return rdggs.cell([address[0]] + [int(d) for d in address[1:]])


def draw_cells_lonlat(ax, addresses, fill_alpha=0.35):
    for address in sorted(addresses):
        cell = cell_from_address(address)
        pts = cell.boundary(n=10, plane=False)
        pts = pts + [pts[0]]
        color = FACE_COLORS[address[0]]
        if all(abs(cur[0] - prev[0]) <= 180 for prev, cur in itertools.pairwise(pts)):
            ax.fill(
                [p[0] for p in pts],
                [p[1] for p in pts],
                color=color,
                alpha=fill_alpha,
                linewidth=0,
            )
        for seg in split_chart_discontinuities(pts):
            if len(seg) > 1:
                ax.plot(*zip(*seg), color=color, linewidth=0.8)


NZ_POLYGON = [
    (166.0, -46.8),
    (169.5, -47.5),
    (174.5, -42.0),
    (178.8, -37.8),
    (177.0, -35.5),
    (173.5, -34.0),
    (171.5, -36.5),
    (166.5, -43.5),
]
LINE = [(166.8, -46.4), (172.0, -41.0), (178.0, -38.2)]
RESOLUTION = 4

fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True)
for ax in axes:
    draw_coastlines_lonlat(ax, linewidth=0.7)
    ax.set_xlim(163, 181)
    ax.set_ylim(-49, -32)
    ax.set_aspect("equal")
    ax.grid(True, linewidth=0.3, alpha=0.5)
    ax.set_xlabel("longitude (degrees)")
axes[0].set_ylabel("latitude (degrees)")

filled = rhp_wrappers.polyfill(Polygon(NZ_POLYGON), RESOLUTION, plane=False, dggs=rdggs)
draw_cells_lonlat(axes[0], filled)
axes[0].plot(*zip(*(NZ_POLYGON + [NZ_POLYGON[0]])), color="#222222", linewidth=1.6)
axes[0].set_title(f"polyfill(polygon, res={RESOLUTION}, plane=False)")

traced = rhp_wrappers.linetrace(LineString(LINE), RESOLUTION, plane=False, dggs=rdggs)
draw_cells_lonlat(axes[1], traced)
axes[1].plot(*zip(*LINE), color="#222222", linewidth=1.6)
axes[1].set_title(f"linetrace(line, res={RESOLUTION}, plane=False)")

fig.tight_layout()
fig.savefig(OUT / "wrappers_nz.svg", bbox_inches="tight")
fig.savefig(OUT / "wrappers_nz.pdf", bbox_inches="tight")
plt.close(fig)
print("wrapper examples written")


# ---------------------------------------------------------------- figure 5
# linetrace across the north polar cap: the case where cells are not
# rectangles in longitude-latitude space, seen from above the pole.

from rhealpixdggs import rhp_wrappers

CAP_LINE = [
    (-134.998756, 86.549596),
    (-179.141527, 88.504030),
    (-44.874903, 86.549596),
    (-89.669615, 86.549596),
    (-134, 86),
]
CAP_RES = 3
traced = rhp_wrappers.linetrace(LineString(CAP_LINE), CAP_RES, plane=False)
traced_set = set(traced)

fig, ax = plt.subplots(figsize=(6.5, 6.5))
VIEW = (0.0, 90.0)  # straight down onto the north pole

# Arctic coastlines for context (northern Greenland, Svalbard, Franz
# Josef Land reach into this view).
for seg in COASTLINES:
    lons = [p[0] for p in seg]
    lats = [p[1] for p in seg]
    x, y, vis = ortho(lons, lats, *VIEW)
    ax.plot(
        np.where(vis, x, np.nan),
        np.where(vis, y, np.nan),
        color=COAST_COLOR,
        linewidth=0.7,
        zorder=1,
    )

# Graticule: parallels and meridians in the window.
for glat in (82, 84, 86, 88):
    lons = np.linspace(-180, 180, 721)
    x, y, vis = ortho(lons, np.full_like(lons, glat), *VIEW)
    ax.plot(
        np.where(vis, x, np.nan),
        np.where(vis, y, np.nan),
        color="#cccccc",
        linewidth=0.4,
    )
for glon in range(-180, 180, 30):
    lats = np.linspace(78, 90, 121)
    x, y, vis = ortho(np.full_like(lats, glon), lats, *VIEW)
    ax.plot(
        np.where(vis, x, np.nan),
        np.where(vis, y, np.nan),
        color="#cccccc",
        linewidth=0.4,
    )

# The resolution 4 sub-grid, as a faint outline.
for cell3 in rdggs.cell(["N", 4, 4]).subcells():
    for cell4 in cell3.subcells():
        pts = cell4.boundary(n=20, plane=False)
        pts = pts + [pts[0]]
        x, y, vis = ortho([p[0] for p in pts], [p[1] for p in pts], *VIEW)
        ax.plot(
            np.where(vis, x, np.nan),
            np.where(vis, y, np.nan),
            color=FACE_COLORS["N"],
            linewidth=0.35,
            alpha=0.45,
            zorder=2,
        )

# The resolution 3 cells around the pole (the children of N44), traced
# ones filled.
for cell in rdggs.cell(["N", 4, 4]).subcells():
    pts = cell.boundary(n=60, plane=False)
    pts = pts + [pts[0]]
    lons = [p[0] for p in pts]
    lats = [p[1] for p in pts]
    x, y, vis = ortho(lons, lats, *VIEW)
    name = str(cell)
    if name in traced_set:
        ax.fill(x, y, color=FACE_COLORS["N"], alpha=0.3, linewidth=0)
    ax.plot(
        np.where(vis, x, np.nan),
        np.where(vis, y, np.nan),
        color=FACE_COLORS["N"],
        linewidth=1.2,
        zorder=3,
    )
    nx, ny, _ = ortho(*cell.nucleus(plane=False), *VIEW)
    ax.text(
        nx,
        ny,
        name,
        ha="center",
        va="center",
        fontsize=9,
        color="#333333",
        bbox={"facecolor": "white", "alpha": 0.6, "linewidth": 0, "pad": 1},
        zorder=5,
    )

# The linestring, densified per leg in longitude-latitude space (its
# segments are straight in that space, so on the globe each leg curves
# around the pole).
for a, b in itertools.pairwise(CAP_LINE):
    ts = np.linspace(0, 1, 400)
    lons = a[0] + ts * (b[0] - a[0])
    lats = a[1] + ts * (b[1] - a[1])
    x, y, vis = ortho(lons, lats, *VIEW)
    ax.plot(
        np.where(vis, x, np.nan),
        np.where(vis, y, np.nan),
        color="#222222",
        linewidth=1.8,
    )
sx, sy, _ = ortho(*CAP_LINE[0], *VIEW)
ex, ey, _ = ortho(*CAP_LINE[-1], *VIEW)
ax.plot([sx], [sy], marker="o", color="#222222", markersize=5, zorder=4)
ax.plot([ex], [ey], marker="s", color="#222222", markersize=5, zorder=4)

lim = float(np.cos(np.radians(80.0)))
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title(
    f"linetrace across the north polar cap (resolution {CAP_RES}, view from above the pole)\n"
    f"traced: {' → '.join(traced)}",
    fontsize=10,
)
fig.tight_layout()
fig.savefig(OUT / "wrappers_cap_trace.svg", bbox_inches="tight")
fig.savefig(OUT / "wrappers_cap_trace.pdf", bbox_inches="tight")
plt.close(fig)
print("cap trace figure written")


# ---------------------------------------------------------------- figure 6
# Cell hierarchy: the SUID is the address. Zoom N -> N4 -> N44 -> N444 in
# the planar grid; each level's cell is the previous one's child, and each
# level appends one digit.
fig, ax = plt.subplots(figsize=(6.5, 6.5))
level_cells = [["N"], ["N", 4], ["N", 4, 4], ["N", 4, 4, 4]]
level_colors = ["#999999", "#7c8fe0", "#79bf6f", "#e5735c"]
n_face = rdggs.cell(["N"])
x0, y0 = n_face.ul_vertex()
w0 = n_face.width()


def to_fig(x, y):
    # Normalized coordinates within the N face: (0,0) bottom-left,
    # (1,1) top-right.
    return (x - x0) / w0, (y - (y0 - w0)) / w0


for suid, color, lw in zip(level_cells, level_colors, (1.2, 1.6, 2.0, 2.4)):
    cell = rdggs.cell(suid)
    x, y = cell.ul_vertex()
    w = cell.width()
    left, top = to_fig(x, y)
    ax.add_patch(
        Rectangle(
            (left, top - w / w0),
            w / w0,
            w / w0,
            facecolor=color if len(suid) == 4 else "none",
            alpha=0.5 if len(suid) == 4 else 1.0,
            edgecolor=color,
            linewidth=lw,
        )
    )
    # Sibling grid inside this cell, faintly.
    for child in cell.subcells():
        cx, cy = child.ul_vertex()
        cw = child.width()
        cleft, ctop = to_fig(cx, cy)
        ax.add_patch(
            Rectangle(
                (cleft, ctop - cw / w0),
                cw / w0,
                cw / w0,
                facecolor="none",
                edgecolor=color,
                linewidth=0.3,
                alpha=0.5,
            )
        )
    ax.annotate(
        str(cell),
        xy=(left, top),
        xytext=(left - 0.05, top + 0.04),
        fontsize=13,
        fontweight="bold",
        color=color,
        ha="right",
        va="bottom",
        arrowprops={"arrowstyle": "-", "color": color, "linewidth": 0.8},
    )
ax.set_xlim(-0.2, 1.02)
ax.set_ylim(-0.02, 1.15)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title(
    "The cell address is the location: each resolution appends one digit,\n"
    "and an address prefix is an ancestor (N contains N4 contains N44 ...)"
)
fig.tight_layout()
fig.savefig(OUT / "hierarchy.svg", bbox_inches="tight")
fig.savefig(OUT / "hierarchy.pdf", bbox_inches="tight")
plt.close(fig)
print("hierarchy figure written")

# ---------------------------------------------------------------- figure 7
# The four ellipsoidal cell shapes, each on an orthographic globe view
# centered on the cell.
SHAPE_EXAMPLES = [
    ("quad", ["Q", 4]),
    ("skew_quad", ["N", 7]),
    ("dart", ["N", 6]),
    ("cap", ["N", 4]),
]
fig, axes = plt.subplots(1, 4, figsize=(11.5, 3.1))
for ax, (shape, suid) in zip(axes, SHAPE_EXAMPLES):
    cell = rdggs.cell(suid)
    assert cell.ellipsoidal_shape == shape, (suid, cell.ellipsoidal_shape)
    lon_c, lat_c = cell.centroid(plane=False)
    view = (lon_c, min(lat_c, 55.0))  # keep some horizon context for the cap
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(t), np.sin(t), color="#555555", linewidth=1.0)
    for seg in COASTLINES:
        x, y, vis = ortho([p[0] for p in seg], [p[1] for p in seg], *view)
        ax.plot(
            np.where(vis, x, np.nan),
            np.where(vis, y, np.nan),
            color=COAST_COLOR,
            linewidth=0.4,
            zorder=1,
        )
    # Neighboring resolution 1 cells for context.
    for face in CELLS0:
        for other in rdggs.cell([face]).subcells():
            pts = other.boundary(n=40, plane=False)
            pts = pts + [pts[0]]
            x, y, vis = ortho([p[0] for p in pts], [p[1] for p in pts], *view)
            ax.plot(
                np.where(vis, x, np.nan),
                np.where(vis, y, np.nan),
                color="#bbbbbb",
                linewidth=0.4,
                zorder=2,
            )
    pts = cell.boundary(n=60, plane=False)
    pts = pts + [pts[0]]
    x, y, vis = ortho([p[0] for p in pts], [p[1] for p in pts], *view)
    if vis.all():
        ax.fill(x, y, color=FACE_COLORS[suid[0]], alpha=0.35, linewidth=0)
    ax.plot(
        np.where(vis, x, np.nan),
        np.where(vis, y, np.nan),
        color=FACE_COLORS[suid[0]],
        linewidth=1.6,
        zorder=3,
    )
    ax.set_xlim(-1.03, 1.03)
    ax.set_ylim(-1.03, 1.03)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"{shape}\n({cell})", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "cell_shapes.svg", bbox_inches="tight")
fig.savefig(OUT / "cell_shapes.pdf", bbox_inches="tight")
plt.close(fig)
print("cell shapes figure written")


# ---------------------------------------------------------------- figure 8
# DE-9IM predicates between cells, in the planar grid: edge touch, corner
# touch, containment, disjoint.
def draw_planar_cells(ax, parent, highlight, labels=True):
    px, py = parent.ul_vertex()
    pw = parent.width()
    for child in parent.subcells():
        x, y = child.ul_vertex()
        w = child.width()
        name = str(child)
        color, alpha = highlight.get(name, ("#ffffff", 0.0))
        left = (x - px) / pw
        bottom = (y - w - (py - pw)) / pw
        ax.add_patch(
            Rectangle(
                (left, bottom),
                w / pw,
                w / pw,
                facecolor=color,
                alpha=alpha if alpha else 1.0,
                edgecolor="#888888",
                linewidth=0.6,
            )
        )
        if labels:
            ax.text(
                left + 0.5 * w / pw,
                bottom + 0.5 * w / pw,
                name,
                ha="center",
                va="center",
                fontsize=8,
                color="#333333",
            )
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.axis("off")


P4 = rdggs.cell(["P", 4])
A_COLOR, B_COLOR = "#7c8fe0", "#e5735c"
cases = [
    ("touches (edge)", {"P41": (A_COLOR, 0.55), "P44": (B_COLOR, 0.55)}),
    ("touches (corner only)", {"P40": (A_COLOR, 0.55), "P44": (B_COLOR, 0.55)}),
    ("within / contains", {"P44": (B_COLOR, 0.55)}),
    ("disjoint", {"P40": (A_COLOR, 0.55), "P48": (B_COLOR, 0.55)}),
]
fig, axes = plt.subplots(1, 4, figsize=(11.5, 3.0))
for ax, (title, highlight) in zip(axes, cases):
    draw_planar_cells(ax, P4, highlight)
    if title == "within / contains":
        # The containing cell is P4 itself: outline the whole frame.
        ax.add_patch(
            Rectangle(
                (0, 0),
                1,
                1,
                facecolor=A_COLOR,
                alpha=0.18,
                edgecolor=A_COLOR,
                linewidth=2.0,
            )
        )
        ax.text(
            0.02,
            1.04,
            "P4",
            fontsize=11,
            fontweight="bold",
            color=A_COLOR,
            ha="left",
            va="bottom",
        )
    ax.set_title(title, fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "predicates.svg", bbox_inches="tight")
fig.savefig(OUT / "predicates.pdf", bbox_inches="tight")
plt.close(fig)
print("predicates figure written")

# --------------------------------------------------------------- figure 8b
# The same predicates on the north polar cap, where the players are caps,
# darts and skew quads rather than planar squares.


def draw_ortho_cells(ax, view, context, highlight, mark=None, label_at=None):
    """Draw cell boundaries in an orthographic view centered on
    ``view = (lon0, lat0)``; ``highlight`` maps cell names to fill
    colors, ``context`` lists extra cells drawn as thin outlines, and
    ``label_at`` optionally overrides a cell's label position (cap
    cells at the same pole share a centroid, so their labels collide)."""
    for name in list(highlight) + context:
        cell = rdggs.cell([name[0]] + [int(d) for d in name[1:]])
        pts = cell.boundary(n=40, plane=False)
        pts = pts + [pts[0]]
        x, y, vis = ortho([p[0] for p in pts], [p[1] for p in pts], *view)
        color = highlight.get(name)
        if color is not None and vis.all():
            ax.fill(x, y, color=color, alpha=0.45, linewidth=0, zorder=2)
        ax.plot(
            np.where(vis, x, np.nan),
            np.where(vis, y, np.nan),
            color=color or "#999999",
            linewidth=1.6 if color else 0.6,
            zorder=3 if color else 1,
        )
        lon_c, lat_c = (label_at or {}).get(name) or cell.centroid(plane=False)
        nx, ny, nvis = ortho(lon_c, lat_c, *view)
        if nvis:
            ax.text(
                nx,
                ny,
                name,
                ha="center",
                va="center",
                fontsize=9 if color else 7,
                fontweight="bold" if color else "normal",
                color="#333333",
                zorder=5,
            )
    if mark is not None:
        mx, my, mvis = ortho(mark[0], mark[1], *view)
        if mvis:
            ax.plot([mx], [my], marker="o", color="#222222", markersize=6, zorder=6)
    ax.set_aspect("equal")
    ax.axis("off")


POLE_VIEW = (0.0, 90.0)
N_CONTEXT = [f"N{i}" for i in range(9)]
CORNER = rdggs.cell(["N", 0]).ul_vertex(plane=False)
# Southernmost boundary latitudes of the two nested caps, used to slot
# the outer cap's label into the ring between their boundaries.
_lat_n4 = min(p[1] for p in rdggs.cell(["N", 4]).boundary(n=40, plane=False))
_lat_n44 = min(p[1] for p in rdggs.cell(["N", 4, 4]).boundary(n=40, plane=False))
polar_cases = [
    # (title, highlight, context, view, radius, mark, label_at)
    (
        "cap touches skew quad (edge)",
        {"N4": B_COLOR, "N1": A_COLOR},
        N_CONTEXT,
        POLE_VIEW,
        0.80,
        None,
        None,
    ),
    (
        "cap touches dart (corner only)",
        {"N4": B_COLOR, "N0": A_COLOR},
        N_CONTEXT,
        POLE_VIEW,
        0.80,
        rdggs.cell(["N", 4]).ul_vertex(plane=False),
        None,
    ),
    (
        "cap within cap",
        {"N44": B_COLOR, "N4": A_COLOR},
        N_CONTEXT,
        POLE_VIEW,
        0.80,
        None,
        {"N4": (0.0, 0.5 * (_lat_n4 + _lat_n44))},
    ),
    (
        "cube corner: all pairs touch (edges)",
        {"N0": A_COLOR, "Q2": B_COLOR, "R0": "#79b791"},
        ["N1", "N3", "Q1", "Q5", "R1", "R3"],
        CORNER,
        0.55,
        CORNER,
        None,
    ),
]
fig, axes = plt.subplots(1, 4, figsize=(12.5, 3.4))
for ax, (title, highlight, context, view, radius, mark, label_at) in zip(
    axes, polar_cases
):
    context = [name for name in context if name not in highlight]
    draw_ortho_cells(ax, view, context, highlight, mark=mark, label_at=label_at)
    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_title(title, fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "predicates_polar.svg", bbox_inches="tight")
fig.savefig(OUT / "predicates_polar.pdf", bbox_inches="tight")
plt.close(fig)
print("polar predicates figure written")

# ---------------------------------------------------------------- figure 9
# A cube corner is 3-valent: exactly three cells meet there, so one
# diagonal neighbor does not exist.
fig, ax = plt.subplots(figsize=(6.0, 6.0))
n0 = rdggs.cell(["N", 0])
corner_lon, corner_lat = n0.ul_vertex(plane=False)
view = (corner_lon, corner_lat)
t = np.linspace(0, 2 * np.pi, 400)
ax.plot(np.cos(t), np.sin(t), color="#555555", linewidth=1.0)
for seg in COASTLINES:
    x, y, vis = ortho([p[0] for p in seg], [p[1] for p in seg], *view)
    ax.plot(
        np.where(vis, x, np.nan),
        np.where(vis, y, np.nan),
        color=COAST_COLOR,
        linewidth=0.5,
        zorder=1,
    )
meeting = {"N0": "the cell", "Q2": "neighbor('up')", "R0": "neighbor('left')"}
context = ["N1", "N3", "N4", "Q1", "Q5", "R1", "R3"]
for name in list(meeting) + context:
    cell = rdggs.cell([name[0], int(name[1])])
    pts = cell.boundary(n=40, plane=False)
    pts = pts + [pts[0]]
    x, y, vis = ortho([p[0] for p in pts], [p[1] for p in pts], *view)
    color = FACE_COLORS[name[0]]
    bold = name in meeting
    if bold and vis.all():
        ax.fill(x, y, color=color, alpha=0.30, linewidth=0)
    ax.plot(
        np.where(vis, x, np.nan),
        np.where(vis, y, np.nan),
        color=color,
        linewidth=1.8 if bold else 0.6,
        zorder=3 if bold else 2,
    )
    lon_c, lat_c = cell.centroid(plane=False)
    nx, ny, nvis = ortho(lon_c, lat_c, *view)
    if nvis:
        label = name + (
            "\n" + meeting[name] if bold and meeting[name] != "the cell" else ""
        )
        ax.text(
            nx,
            ny,
            label,
            ha="center",
            va="center",
            fontsize=10 if bold else 8,
            fontweight="bold" if bold else "normal",
            color="#333333",
            zorder=5,
        )
cx, cy, _ = ortho(corner_lon, corner_lat, *view)
ax.plot([cx], [cy], marker="o", color="#222222", markersize=7, zorder=6)
ax.annotate(
    "cube corner: 3 cells meet here\nN0.diagonal_neighbor('up_left') is None",
    xy=(cx, cy),
    xytext=(cx - 0.55, cy + 0.62),
    fontsize=10,
    ha="center",
    color="#222222",
    arrowprops={"arrowstyle": "->", "color": "#222222"},
)
lim = 0.9
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("Three cells, not four, meet at a corner of the cube", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "cube_corner.svg", bbox_inches="tight")
fig.savefig(OUT / "cube_corner.pdf", bbox_inches="tight")
plt.close(fig)
print("cube corner figure written")

# --------------------------------------------------------------- figure 9b
# The same corner drawn on the cube itself: the three visible faces N, Q
# and R subdivided into their nine resolution 1 cells.


def iso(x, y, z):
    """Project cube coordinates in [0, 1]^3 onto the drawing plane."""
    return (x - y) * 0.866, z * 0.82 - (x + y) * 0.35


# Map each visible face's planar coordinates (p to the right, q downward,
# both in [0, 1]) onto the cube, oriented so shared edges line up the way
# the faces fold: the top (q=0) edges of Q and R meet N, N's top (q=0)
# edge lies along Q, N's left (p=0) edge lies along R, and Q's right
# (p=1) edge is R's left (p=0) edge. The corner vertex (1, 1, 1) is then
# shared by exactly the children N0, Q2 and R0.
CUBE_FACE_POINT = {
    "N": lambda p, q: (1 - p, 1 - q, 1.0),
    "Q": lambda p, q: (p, 1.0, 1 - q),
    "R": lambda p, q: (1.0, 1 - p, 1 - q),
}
CORNER_CELLS = ("N0", "Q2", "R0")

fig, ax = plt.subplots(figsize=(6.2, 6.8))
for face, to_cube in CUBE_FACE_POINT.items():
    color = FACE_COLORS[face]
    for row in range(3):
        for col in range(3):
            name = f"{face}{3 * row + col}"
            corners = [
                (col / 3, row / 3),
                ((col + 1) / 3, row / 3),
                ((col + 1) / 3, (row + 1) / 3),
                (col / 3, (row + 1) / 3),
            ]
            key = name in CORNER_CELLS
            ax.add_patch(
                PolygonPatch(
                    [iso(*to_cube(p, q)) for p, q in corners],
                    closed=True,
                    facecolor=to_rgba(color, 0.45 if key else 0.10),
                    edgecolor="#777777",
                    linewidth=0.6,
                    zorder=2,
                )
            )
            lx, ly = iso(*to_cube((col + 0.5) / 3, (row + 0.5) / 3))
            ax.text(
                lx,
                ly,
                name,
                ha="center",
                va="center",
                fontsize=8 if key else 7,
                fontweight="bold" if key else "normal",
                color="#333333",
                zorder=4,
            )

# Outer silhouette and the three edges that meet at the corner vertex.
silhouette = [(0, 1, 1), (0, 0, 1), (1, 0, 1), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
ax.plot(
    *zip(*[iso(*v) for v in silhouette + [silhouette[0]]]),
    color="#444444",
    linewidth=1.4,
    zorder=5,
)
for far in [(0, 1, 1), (1, 0, 1), (1, 1, 0)]:
    ax.plot(*zip(iso(1, 1, 1), iso(*far)), color="#444444", linewidth=1.0, zorder=5)

vx, vy = iso(1, 1, 1)
ax.plot([vx], [vy], marker="o", color="#222222", markersize=7, zorder=6)
ax.annotate(
    "cube corner: 3 cells meet here\nN0.diagonal_neighbor('up_left') is None",
    xy=(vx, vy),
    xytext=(0, -1.0),
    fontsize=10,
    ha="center",
    va="top",
    color="#222222",
    arrowprops={"arrowstyle": "->", "color": "#222222", "shrinkB": 6},
)
ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.45, 1.3)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("The same corner on the cube: faces N, Q and R", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "cube_corner_cube.svg", bbox_inches="tight")
fig.savefig(OUT / "cube_corner_cube.pdf", bbox_inches="tight")
plt.close(fig)
print("cube corner (cube view) figure written")

# --------------------------------------------------------------- figure 10
# The north_square/south_square parameters: where the polar squares sit.

from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID

fig, axes = plt.subplots(2, 1, figsize=(8, 6.4))
for ax, (ns, ss) in zip(axes, [(0, 0), (1, 2)]):
    rd = RHEALPixDGGS(
        ellipsoid=WGS84_ELLIPSOID, north_square=ns, south_square=ss, N_side=3
    )
    Rl = rd.ellipsoid.R_A
    for face in CELLS0:
        c0 = rd.cell([face])
        x, y = c0.ul_vertex()
        w = c0.width()
        ax.add_patch(
            Rectangle(
                (x / Rl, (y - w) / Rl),
                w / Rl,
                w / Rl,
                facecolor=FACE_COLORS[face],
                alpha=0.35,
                edgecolor="black",
                linewidth=1.2,
            )
        )
        ax.text(
            (x + 0.5 * w) / Rl,
            (y - 0.5 * w) / Rl,
            face,
            ha="center",
            va="center",
            fontsize=16,
            fontweight="bold",
            color="#333333",
            alpha=0.85,
        )
    for seg in COASTLINES:
        xy = [rd.rhealpix(lon, lat) for lon, lat in seg]
        run = [xy[0]]
        for prev, cur in itertools.pairwise(xy):
            if abs(cur[0] - prev[0]) > 0.15 * Rl or abs(cur[1] - prev[1]) > 0.15 * Rl:
                if len(run) > 1:
                    ax.plot(
                        [p[0] / Rl for p in run],
                        [p[1] / Rl for p in run],
                        color=COAST_COLOR,
                        linewidth=0.4,
                        zorder=1,
                    )
                run = []
            run.append(cur)
        if len(run) > 1:
            ax.plot(
                [p[0] / Rl for p in run],
                [p[1] / Rl for p in run],
                color=COAST_COLOR,
                linewidth=0.4,
                zorder=1,
            )
    ax.set_xlim(-3.4, 3.4)
    ax.set_ylim(-2.6, 2.6)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"north_square={ns}, south_square={ss}", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "polar_squares.svg", bbox_inches="tight")
fig.savefig(OUT / "polar_squares.pdf", bbox_inches="tight")
plt.close(fig)
print("polar squares figure written")


# --------------------------------------------------------------- figure 11
# Rotating the DGGS about the polar axis so that New Zealand sits in the
# middle of an equatorial face. Face edges and polar dart meridians lie on
# lon_0 + k*90, so lon_0 45 degrees west of Auckland keeps them clear of
# the country. Rotation is an isometry of the ellipsoid, so the rotated
# grid is congruent to the standard one. (Latitude recentring does not
# exist; Ellipsoid rejects a nonzero lat_0, see issue #93 and dggs.rst.)
from rhealpixdggs.ellipsoids import WGS84_A, WGS84_F, Ellipsoid

AKL = (174.0, -37.0)
LON_0_NZ = AKL[0] - 45.0  # 129
E_AKL = Ellipsoid(a=WGS84_A, f=WGS84_F, radians=False, lon_0=LON_0_NZ)
rd_akl = RHEALPixDGGS(E_AKL, N_side=3, north_square=0, south_square=0)
fig, ax = plt.subplots(figsize=(9, 4.8))
draw_coastlines_lonlat(ax, linewidth=0.5)
for face in CELLS0:
    color = FACE_COLORS[face]
    for cell in rd_akl.cell([face]).subcells():
        pts = cell.boundary(n=40, plane=False)
        pts = pts + [pts[0]]
        for seg in split_chart_discontinuities(pts):
            if len(seg) > 1:
                ax.plot(*zip(*seg), color=color, linewidth=0.8, zorder=2)
akl_cell = rd_akl.cell_from_point(1, AKL, plane=False)
pts = akl_cell.boundary(n=40, plane=False)
pts = pts + [pts[0]]
# The cell straddles the antimeridian, so fill it through the wrap.
fill_lonlat_polygon(
    ax, pts, color=FACE_COLORS[str(akl_cell)[0]], alpha=0.4, linewidth=0, zorder=3
)
ax.plot([AKL[0]], [AKL[1]], marker="*", color="#222222", markersize=12, zorder=4)
ax.annotate(
    f"Auckland, 45° east of lon_0 = {LON_0_NZ:.0f}:\n" f"mid-face, in cell {akl_cell}",
    xy=AKL,
    xytext=(AKL[0] - 120, AKL[1] - 30),
    fontsize=9,
    color="#222222",
    arrowprops={"arrowstyle": "->", "color": "#222222"},
)
ax.set_xlim(-180, 180)
ax.set_ylim(-90, 90)
ax.set_xticks(range(-180, 181, 60))
ax.set_yticks(range(-90, 91, 30))
ax.set_xlabel("longitude (degrees)")
ax.set_ylabel("latitude (degrees)")
ax.set_title(
    f"Resolution 1 cells of a DGGS rotated to put New Zealand mid-face "
    f"(lon_0={LON_0_NZ:.0f})"
)
ax.grid(True, linewidth=0.3, alpha=0.5)
fig.tight_layout()
fig.savefig(OUT / "recentred.svg", bbox_inches="tight")
fig.savefig(OUT / "recentred.pdf", bbox_inches="tight")
plt.close(fig)
print("recentred figure written")


# -------------------------------------------------------------- figure 11b
# Zoomed comparisons of the standard grid with a rotated one, cells
# coloured by ellipsoidal shape, so the effect of lon_0 on which shapes a
# region gets -- and which way its skew quads lean -- is visible. The
# equatorial/polar seam at +/-phi_0 cannot be moved by any rotation.
from matplotlib.patches import Patch

SHAPE_COLORS = {
    "quad": "#79bf6f",
    "skew_quad": "#7c8fe0",
    "dart": "#e5735c",
    "cap": "#e8c34f",
}
PHI_0 = rdggs.ellipsoid.phi_0  # degrees; the equatorial/polar seam


def cells_in_window(dggs, resolution, lon_range, lat_range, step):
    """Every resolution `resolution` cell of `dggs` that contains a point
    of the lon-lat lattice with spacing `step` over the window, in SUID
    order. `step` must be well under a cell width so no cell is missed."""
    cells = {}
    lons = np.arange(lon_range[0], lon_range[1] + step / 2, step)
    lats = np.arange(lat_range[0], lat_range[1] + step / 2, step)
    for lon in lons:
        lon = ((float(lon) + 180.0) % 360.0) - 180.0  # the window may pass 180
        for lat in lats:
            c = dggs.cell_from_point(resolution, (lon, float(lat)), plane=False)
            cells[str(c)] = c
    return [cells[k] for k in sorted(cells)]


def draw_cells_by_shape(ax, cells, lon_range, n=10):
    """Fill and outline each cell by its ellipsoidal shape. Longitudes are
    unwrapped into one continuous ring and the ring is drawn at whichever
    360-degree offsets meet the window, so a window that passes the
    antimeridian shows the cells on its far side whole."""
    for cell in cells:
        pts = cell.boundary(n=n, plane=False)
        pts = pts + [pts[0]]
        lons = [pts[0][0]]
        for prev, cur in itertools.pairwise(pts):
            d = cur[0] - prev[0]
            if d > 180:
                d -= 360
            elif d < -180:
                d += 360
            lons.append(lons[-1] + d)
        lats = [p[1] for p in pts]
        color = SHAPE_COLORS[cell.ellipsoidal_shape]
        for off in (-360.0, 0.0, 360.0):
            shifted = [lon + off for lon in lons]
            if max(shifted) < lon_range[0] or min(shifted) > lon_range[1]:
                continue
            ax.fill(shifted, lats, color=color, alpha=0.35, linewidth=0, zorder=2)
            ax.plot(shifted, lats, color=color, linewidth=0.7, zorder=3)


def draw_dart_meridians(ax, lon_0, lon_range):
    """Dashed verticals on the meridians lon_0 + k*90 inside the window:
    the polar dart meridians, which are also the equatorial face edges."""
    k0 = int(np.floor((lon_range[0] - lon_0) / 90))
    k1 = int(np.ceil((lon_range[1] - lon_0) / 90))
    for k in range(k0, k1 + 1):
        lon = lon_0 + 90 * k
        if lon_range[0] <= lon <= lon_range[1]:
            ax.axvline(lon, color="#e5735c", linestyle="--", linewidth=1.0, zorder=4)


def shape_legend(ax, shapes, loc):
    handles = [Patch(facecolor=SHAPE_COLORS[s], alpha=0.5, label=s) for s in shapes]
    handles.append(
        plt.Line2D(
            [], [], color="#444444", linestyle=":", label="equatorial/polar seam"
        )
    )
    handles.append(
        plt.Line2D([], [], color="#e5735c", linestyle="--", label="dart meridian")
    )
    ax.legend(handles=handles, loc=loc, fontsize=8, framealpha=0.9)


def rotated_comparison(
    lon_0, resolution, lon_range, lat_range, step, aspect, title, legend_loc
):
    dggs_rot = RHEALPixDGGS(
        Ellipsoid(a=WGS84_A, f=WGS84_F, radians=False, lon_0=lon_0),
        N_side=3,
        north_square=0,
        south_square=0,
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.6), sharey=True)
    shapes_seen = set()
    for ax, (dggs_i, label) in zip(
        axes,
        [
            (rdggs, "standard grid (lon_0=0)"),
            (dggs_rot, f"rotated grid (lon_0={lon_0:g})"),
        ],
    ):
        draw_coastlines_lonlat(ax, linewidth=0.7)
        cells = cells_in_window(dggs_i, resolution, lon_range, lat_range, step)
        shapes_seen.update(c.ellipsoidal_shape for c in cells)
        draw_cells_by_shape(ax, cells, lon_range)
        seam = PHI_0 if lat_range[0] >= 0 else -PHI_0
        ax.axhline(seam, color="#444444", linestyle=":", linewidth=1.0, zorder=4)
        draw_dart_meridians(ax, dggs_i.ellipsoid.lon_0, lon_range)
        ax.set_xlim(*lon_range)
        ax.set_ylim(*lat_range)
        ax.set_aspect(aspect)
        ax.grid(True, linewidth=0.3, alpha=0.5)
        ax.set_xlabel("longitude (degrees)")
        ax.set_title(f"{label}, resolution {resolution}", fontsize=11)
    axes[0].set_ylabel("latitude (degrees)")
    shape_legend(axes[1], [s for s in SHAPE_COLORS if s in shapes_seen], legend_loc)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    return fig


# New Zealand: the standard grid's dart meridian at 180 runs just east of
# the country and its face edge at 180 too; lon_0=129 moves both 45
# degrees away so NZ sits midway between dart meridians.
fig = rotated_comparison(
    LON_0_NZ,
    resolution=4,
    lon_range=(163, 181),
    lat_range=(-49, -32),
    step=0.25,
    aspect="equal",
    title="Cell shapes over New Zealand: standard grid vs rotated grid",
    legend_loc="upper right",
)
fig.savefig(OUT / "recentred_nz.svg", bbox_inches="tight")
fig.savefig(OUT / "recentred_nz.pdf", bbox_inches="tight")
plt.close(fig)
print("recentred NZ comparison written")

# Canada, after Bowater & Stefanakis (2018, section 4.2): the standard
# grid's dart meridian at 90 W runs through the middle of the country;
# lon_0=-50 moves the dart meridians to 50 W and 140 W, at its edges.
# Plate carree stretches longitude at these latitudes, so shorten the
# longitude axis to roughly restore local shape (cos 60 = 0.5).
fig = rotated_comparison(
    -50.0,
    resolution=3,
    lon_range=(-145, -50),
    lat_range=(40, 85),
    step=0.5,
    aspect=2.0,
    title="Cell shapes over Canada: standard grid vs rotated grid",
    legend_loc="lower left",
)
fig.savefig(OUT / "recentred_canada.svg", bbox_inches="tight")
fig.savefig(OUT / "recentred_canada.pdf", bbox_inches="tight")
plt.close(fig)
print("recentred Canada comparison written")

# --------------------------------------------------------------- figure 12
# cell_ring: rings by distance, and the 7-cell ring at a cube corner.
RING_COLORS = ["#e5735c", "#e8c34f", "#79bf6f", "#5cc3c9"]


def draw_ring_panel(ax, center_name, max_k, window_faces):
    rdggs.cell([c if c in CELLS0 else int(c) for c in center_name])
    rings = {0: [center_name]}
    for k in range(1, max_k + 1):
        rings[k] = rhp_wrappers.cell_ring(center_name, k)
    membership = {}
    for k, names in rings.items():
        for nm in names:
            membership[nm] = k
    for face in window_faces:
        for cell in rdggs.cell([face]).subcells():
            for sub in cell.subcells():
                x, y = sub.ul_vertex()
                w = sub.width()
                name = str(sub)
                k = membership.get(name)
                color = RING_COLORS[k] if k is not None else "#ffffff"
                ax.add_patch(
                    Rectangle(
                        (x, y - w),
                        w,
                        w,
                        facecolor=color,
                        alpha=0.55 if k is not None else 1.0,
                        edgecolor="#aaaaaa",
                        linewidth=0.3,
                    )
                )
    # face outlines
    for face in window_faces:
        c0 = rdggs.cell([face])
        x, y = c0.ul_vertex()
        w = c0.width()
        ax.add_patch(
            Rectangle(
                (x, y - w), w, w, facecolor="none", edgecolor="#444444", linewidth=1.2
            )
        )
        ax.text(
            x + 0.06 * w,
            y - 0.06 * w,
            face,
            ha="left",
            va="top",
            fontsize=12,
            fontweight="bold",
            color="#444444",
        )
    ax.set_aspect("equal")
    ax.axis("off")


fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4))
# Interior rings around Q44 (resolution 2), all within Q: planar view.
draw_ring_panel(axes[0], "Q44", 3, ["Q"])
q = rdggs.cell(["Q"])
qx, qy = q.ul_vertex()
qw = q.width()
axes[0].set_xlim(qx - 0.05 * qw, qx + 1.05 * qw)
axes[0].set_ylim(qy - 1.05 * qw, qy + 0.05 * qw)
axes[0].set_title("cell_ring('Q44', k=1..3): interior rings have 8k cells", fontsize=11)

# Rings around N00, whose corner digit chain puts it at a cube corner:
# drawn on the globe, because on the cube the ring wraps contiguously
# around the corner across the N, Q, and R faces.
ax = axes[1]
n00 = rdggs.cell(["N", 0, 0])
corner_lon, corner_lat = rdggs.cell(["N"]).ul_vertex(plane=False)
view = (corner_lon, corner_lat)
membership = {"N00": 0}
for k in (1, 2, 3):
    for nm in rhp_wrappers.cell_ring("N00", k):
        membership[nm] = k
for seg in COASTLINES:
    x, y, vis = ortho([p[0] for p in seg], [p[1] for p in seg], *view)
    ax.plot(
        np.where(vis, x, np.nan),
        np.where(vis, y, np.nan),
        color=COAST_COLOR,
        linewidth=0.5,
        zorder=1,
    )
for face in ("N", "Q", "R"):
    for cell1 in rdggs.cell([face]).subcells():
        for cell in cell1.subcells():
            pts = cell.boundary(n=25, plane=False)
            pts = pts + [pts[0]]
            x, y, vis = ortho([p[0] for p in pts], [p[1] for p in pts], *view)
            k = membership.get(str(cell))
            if k is not None and vis.all():
                ax.fill(x, y, color=RING_COLORS[k], alpha=0.55, linewidth=0)
            ax.plot(
                np.where(vis, x, np.nan),
                np.where(vis, y, np.nan),
                color="#aaaaaa",
                linewidth=0.3,
                zorder=2,
            )
cx, cy, _ = ortho(corner_lon, corner_lat, *view)
ax.plot([cx], [cy], marker="o", color="#222222", markersize=5, zorder=6)
lim = 0.42
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("around a cube corner ('N00'): the k=1 ring has 7 cells", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "rings.svg", bbox_inches="tight")
fig.savefig(OUT / "rings.pdf", bbox_inches="tight")
plt.close(fig)
print("rings figure written")


# --------------------------------------------------------------- figure 13
# The wrap_antimeridian flag: the same segment traced both ways.
SEG = ((179.0, 10.0), (-179.0, 10.0))
fig, axes = plt.subplots(2, 1, figsize=(9, 6.6))
for ax, wrap in zip(axes, (False, True)):
    draw_coastlines_lonlat(ax, linewidth=0.5)
    traced = rdggs.cells_from_line(
        1, SEG[0], SEG[1], plane=False, wrap_antimeridian=wrap
    )
    for cell in traced:
        pts = cell.boundary(n=20, plane=False)
        pts = pts + [pts[0]]
        color = FACE_COLORS[str(cell)[0]]
        fill_lonlat_polygon(ax, pts, color=color, alpha=0.4, linewidth=0)
        for seg in split_chart_discontinuities(pts):
            if len(seg) > 1:
                ax.plot(*zip(*seg), color=color, linewidth=0.9)
    if wrap:
        ax.plot([SEG[0][0], 180], [SEG[0][1], SEG[0][1]], color="#222222", linewidth=2)
        ax.plot([-180, SEG[1][0]], [SEG[1][1], SEG[1][1]], color="#222222", linewidth=2)
    else:
        ax.plot(
            [SEG[0][0], SEG[1][0]], [SEG[0][1], SEG[1][1]], color="#222222", linewidth=2
        )
    ax.plot([SEG[0][0]], [SEG[0][1]], marker="o", color="#222222", markersize=6)
    ax.plot([SEG[1][0]], [SEG[1][1]], marker="s", color="#222222", markersize=6)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-65, 65)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-60, 61, 30))
    ax.grid(True, linewidth=0.3, alpha=0.5)
    ax.set_title(
        "wrap_antimeridian={}: (179, 10) to (-179, 10) traced {}".format(
            wrap,
            (
                "the short way, across the antimeridian"
                if wrap
                else "the long way, through longitude 0 (the literal planar reading)"
            ),
        ),
        fontsize=10,
    )
axes[1].set_xlabel("longitude (degrees)")
fig.tight_layout()
fig.savefig(OUT / "wrap_antimeridian.svg", bbox_inches="tight")
fig.savefig(OUT / "wrap_antimeridian.pdf", bbox_inches="tight")
plt.close(fig)
print("wrap flag figure written")

# --------------------------------------------------------------- figure 14
# compact_cells: polyfill output before and after compaction.
from rhealpixdggs.conversion import compact_cells

FILL_RES = 5
filled = rhp_wrappers.polyfill(Polygon(NZ_POLYGON), FILL_RES, plane=False, dggs=rdggs)
compacted = compact_cells(filled, N_side=rdggs.N_side)
fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True)
for ax, cells, title in (
    (axes[0], filled, f"polyfill(polygon, res={FILL_RES}): {len(filled)} cells"),
    (
        axes[1],
        compacted,
        f"compact_cells(...): {len(compacted)} cells, mixed resolutions",
    ),
):
    draw_coastlines_lonlat(ax, linewidth=0.7)
    for address in sorted(cells):
        cell = cell_from_address(address)
        pts = cell.boundary(n=10, plane=False)
        pts = pts + [pts[0]]
        # Shade by resolution: coarser cells darker.
        depth = len(address) - 1
        alpha = 0.55 - 0.08 * (depth - 3)
        color = FACE_COLORS[address[0]]
        ax.fill(
            [p[0] for p in pts],
            [p[1] for p in pts],
            color=color,
            alpha=max(alpha, 0.2),
            linewidth=0,
        )
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color, linewidth=0.5)
    ax.plot(*zip(*(NZ_POLYGON + [NZ_POLYGON[0]])), color="#222222", linewidth=1.4)
    ax.set_xlim(163, 181)
    ax.set_ylim(-49, -32)
    ax.set_aspect("equal")
    ax.grid(True, linewidth=0.3, alpha=0.5)
    ax.set_xlabel("longitude (degrees)")
    ax.set_title(title, fontsize=11)
axes[0].set_ylabel("latitude (degrees)")
fig.tight_layout()
fig.savefig(OUT / "compaction.svg", bbox_inches="tight")
fig.savefig(OUT / "compaction.pdf", bbox_inches="tight")
plt.close(fig)
print("compaction figure written")

# Drop matplotlib's six decimal places of coordinate precision, which is
# around a nanometre on the page and about a third of every SVG's bytes.
# Committed figures are stored shrunk, so this keeps regeneration diff-free.
saved = 0
for svg in sorted(OUT.glob("*.svg")):
    before, after = svg_shrink.shrink_file(svg)
    saved += before - after
print(f"figures shrunk: {saved:,} bytes of surplus coordinate precision removed")
