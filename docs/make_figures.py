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

All grid geometry is drawn from the library itself (Cell.ul_vertex,
Cell.boundary, Cell.nucleus). Coastlines are Natural Earth 1:110m data
(public domain), downloaded on first run (pinned to a specific upstream
commit for reproducibility) and cached under .cache/ at the repository
root (outside docs/, which is shipped in the sdist).
"""
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
from matplotlib.patches import Rectangle

# Draw from this repository's rhealpixdggs, not any installed copy.
sys.path.insert(0, str(pathlib.Path(__file__).parents[1]))

from rhealpixdggs.dggs import WGS84_003, CELLS0

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
            for prev, cur in zip(line, line[1:]):
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
    for prev, cur in zip(xy, xy[1:]):
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
ax.set_title(
    "(0, 0)-rHEALPix planar grid: resolution 0 cells, resolution 1 sub-grid"
)
fig.tight_layout()
fig.savefig(OUT / "planar_grid.svg", bbox_inches="tight")
fig.savefig(OUT / "planar_grid.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- figure 2
# Ellipsoidal cells: every resolution 1 cell's boundary in longitude-
# latitude coordinates, colored by resolution 0 face; shows the quad,
# skew-quad, dart, and cap shapes.
fig, ax = plt.subplots(figsize=(9, 4.8))


def split_at_antimeridian(points):
    """Split a closed boundary point list into segments wherever the
    longitude jumps across the antimeridian."""
    segs = [[points[0]]]
    for prev, cur in zip(points, points[1:]):
        if abs(cur[0] - prev[0]) > 180:
            segs.append([])
        segs[-1].append(cur)
    return segs


def draw_coastlines_lonlat(ax, linewidth=0.5):
    for seg in COASTLINES:
        ax.plot(*zip(*seg), color=COAST_COLOR, linewidth=linewidth, zorder=1)


draw_coastlines_lonlat(ax)
for face in CELLS0:
    color = FACE_COLORS[face]
    for cell in rdggs.cell([face]).subcells():
        pts = cell.boundary(n=40, plane=False)
        pts = pts + [pts[0]]
        for seg in split_at_antimeridian(pts):
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
            ncos = (
                np.sin(np.radians(lat0)) * np.sin(np.radians(cell.nucleus(plane=False)[1]))
                + np.cos(np.radians(lat0))
                * np.cos(np.radians(cell.nucleus(plane=False)[1]))
                * np.cos(np.radians(cell.nucleus(plane=False)[0] - lon0))
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
        if all(abs(cur[0] - prev[0]) <= 180 for prev, cur in zip(pts, pts[1:])):
            ax.fill(
                [p[0] for p in pts],
                [p[1] for p in pts],
                color=color,
                alpha=fill_alpha,
                linewidth=0,
            )
        for seg in split_at_antimeridian(pts):
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

filled = rhp_wrappers.polyfill(
    Polygon(NZ_POLYGON), RESOLUTION, plane=False, dggs=rdggs
)
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
    ax.plot(np.where(vis, x, np.nan), np.where(vis, y, np.nan),
            color=COAST_COLOR, linewidth=0.7, zorder=1)

# Graticule: parallels and meridians in the window.
for glat in (82, 84, 86, 88):
    lons = np.linspace(-180, 180, 721)
    x, y, vis = ortho(lons, np.full_like(lons, glat), *VIEW)
    ax.plot(np.where(vis, x, np.nan), np.where(vis, y, np.nan),
            color="#cccccc", linewidth=0.4)
for glon in range(-180, 180, 30):
    lats = np.linspace(78, 90, 121)
    x, y, vis = ortho(np.full_like(lats, glon), lats, *VIEW)
    ax.plot(np.where(vis, x, np.nan), np.where(vis, y, np.nan),
            color="#cccccc", linewidth=0.4)

# The resolution 4 sub-grid, as a faint outline.
for cell3 in rdggs.cell(["N", 4, 4]).subcells():
    for cell4 in cell3.subcells():
        pts = cell4.boundary(n=20, plane=False)
        pts = pts + [pts[0]]
        x, y, vis = ortho([p[0] for p in pts], [p[1] for p in pts], *VIEW)
        ax.plot(np.where(vis, x, np.nan), np.where(vis, y, np.nan),
                color=FACE_COLORS["N"], linewidth=0.35, alpha=0.45, zorder=2)

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
    ax.plot(np.where(vis, x, np.nan), np.where(vis, y, np.nan),
            color=FACE_COLORS["N"], linewidth=1.2, zorder=3)
    nx, ny, _ = ortho(*cell.nucleus(plane=False), *VIEW)
    ax.text(
        nx,
        ny,
        name,
        ha="center",
        va="center",
        fontsize=9,
        color="#333333",
        bbox=dict(facecolor="white", alpha=0.6, linewidth=0, pad=1),
        zorder=5,
    )

# The linestring, densified per leg in longitude-latitude space (its
# segments are straight in that space, so on the globe each leg curves
# around the pole).
for a, b in zip(CAP_LINE, CAP_LINE[1:]):
    ts = np.linspace(0, 1, 400)
    lons = a[0] + ts * (b[0] - a[0])
    lats = a[1] + ts * (b[1] - a[1])
    x, y, vis = ortho(lons, lats, *VIEW)
    ax.plot(np.where(vis, x, np.nan), np.where(vis, y, np.nan),
            color="#222222", linewidth=1.8)
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
    "linetrace across the north polar cap (resolution %d, view from above the pole)\n"
    "traced: %s" % (CAP_RES, " → ".join(traced)),
    fontsize=10,
)
fig.tight_layout()
fig.savefig(OUT / "wrappers_cap_trace.svg", bbox_inches="tight")
fig.savefig(OUT / "wrappers_cap_trace.pdf", bbox_inches="tight")
plt.close(fig)
print("cap trace figure written")
