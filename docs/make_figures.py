"""
Regenerate the static figures in docs/source/images/.

Run manually from the repository root whenever the figures need to change:

    python docs/make_figures.py

Requires matplotlib, which is deliberately NOT a dependency of this
package or of the docs build -- the generated SVGs are committed, and
regenerating them on every build would only add churn (different
matplotlib versions emit slightly different SVGs) for figures that only
change when the grid geometry or the figure design does.

All geometry is drawn from the library itself (Cell.ul_vertex,
Cell.boundary, Cell.nucleus).
"""
import os
import pathlib
import sys

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
print("globe views written")
