# html/lib/galaxymap.py

"""
Galaxy-scale map: every galaxy-placed sector (`sectors.center_x/y/z_pc` not
NULL -- see `docs/design/galaxy-coordinate-system.md`) plotted as a bright,
star-like dot at its real position, split into four 90-degree **Quadrants**
(I-IV, azimuthal, this module's own `sector_quadrant`) and concentric
**Rings** (fixed-width bands of `shell_index`, `sector_ring`) so a person
can navigate from "the whole galaxy" down to one sector without scrolling a
flat list that could, in principle, run to hundreds of millions of rows.

Terminology note: `star_systems.quadrant` (Roman numerals I-VIII,
`spaceSector.classify_octant`) is a different, older concept -- an octant
classification of a *system's* position within its own *sector*. That
column/attribute name is left alone for schema stability, but every place
`html/` displays it now says "Octant" (`sector.py`, `system.py`,
`static/sectormap.js`) specifically so it doesn't collide with *this*
module's Quadrant, which is the real, galaxy-scale, 4-region azimuthal
concept the word ordinarily means (and which is a genuine 2D mathematical
standard, unlike the 8-region octant workaround -- see `spaceSector.py`'s
own module docstring).

Rendered as one flat 2D SVG rather than `starmap.py`'s rotatable 3D CSS
cube -- a Quadrant is inherently a flat, azimuthal split (blind to height),
so there's no "rotate to inspect" need the way a sector's cubic interior
has. Drill-down is plain server-side navigation (a different `?quadrant=`
query string), not client-side JS: `render_galaxy_map_panel` always draws
every placed sector into the same fixed coordinate space, and a Quadrant
view simply crops the same drawing to one quarter of it via the outer
`<svg>`'s own `viewBox` (SVG clips to its viewBox by default) -- so no new
static/*.js file is needed for this map at all.

Density is two genuinely different things here, deliberately drawn two
different ways:

- **Real, generated sectors** -- bright dots, sized/colored by their own
  actual `system_count`.
- **Un-generated space** -- a soft radial "cloud" gradient (brighter near
  the core, fading outward). This is purely illustrative shading, NOT a
  real population model -- no disk/bulge/spiral density envelope exists
  anywhere in this codebase yet (`galaxyGen.py` generates any requested
  shell/neighborhood uniformly; see
  `docs/design/galaxy-coordinate-system.md` section 7, question 2, still
  open). It exists only so an overwhelmingly-empty galaxy -- which is the
  normal, expected state, since the addressable volume runs into the
  hundreds of billions of sector slots -- doesn't render as a blank void.
"""

import math

from fmt import esc

try:
    from stellarObjects.program_constants import DEFAULT_SECTOR_EDGE_LY
    from stellarObjects.utils import pc_to_ly
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # duplicated fallback, matching starmap.py's own pattern (dbutil.py
    # itself no longer has one -- stellarObjects is a hard dependency
    # there now that its MySQL connection helpers are load-bearing).
    DEFAULT_SECTOR_EDGE_LY = 11.5
    pc_to_ly = None

QUADRANT_LABELS = ("I", "II", "III", "IV")
"""tuple[str]: The four galaxy-scale Quadrants, in azimuthal order starting
from `+X` -- `sector_quadrant` indexes into this. Deliberately the plain
2D I-IV convention (see this module's own docstring) rather than the
sector-internal octant's Roman-numeral scheme."""

RING_TARGET_LY = 100.0
"""float: The approximate light-year thickness a Ring should aim for.
`RING_SHELL_WIDTH` is derived from this so Ring boundaries land close to
round light-year milestones (~100 ly, ~200 ly, ...) while still being a
*fixed* shell-count width per Ring end to end (the preferred approach --
simple, uniform, and tracks the addressing scheme's own `shell_index`
directly, unlike a milestone-exact width, which would have to vary per
Ring)."""

RING_SHELL_WIDTH = max(1, round(RING_TARGET_LY / DEFAULT_SECTOR_EDGE_LY))
"""int: How many consecutive `shell_index` values make up one Ring."""

_SVG_SIZE = 560
_CENTER = _SVG_SIZE / 2
_MAP_RADIUS_PX = _SVG_SIZE / 2 - 40  # leaves margin for ring/quadrant labels
_MIN_RINGS_SHOWN = 3

_MIN_DOT_R = 2.6
_MAX_DOT_R = 12.0

# Bright, warm "star" coloring -- deliberately not spectral-accurate the
# way starmap.py's per-star dots are: one dot here represents a whole
# sector (a handful to a few dozen systems of every spectral type mixed
# together), so there's no single meaningful color to derive from real
# data. A halo (larger, low-opacity) behind a solid core is what reads as
# "bright" rather than a single flat circle -- both scale with the
# sector's own system_count via `_star_visual`.
_STAR_CORE_FILL = "#fff6df"
_STAR_CORE_STROKE = "#caa54d"
_STAR_HALO_FILL = "#ffd88a"

# theta=0/90/180/270 boundary -> which corner of the SVG's own coordinate
# square that Quadrant occupies, as (viewbox_x, viewbox_y) of a
# _SVG_SIZE/2-wide/tall crop -- see this module's docstring on how
# Quadrant zoom works (a plain viewBox crop, not a wedge clip-path, so a
# Quadrant view's corner nearest the core still shows a sliver of its two
# neighbors along the axis lines; accepted as a simplification, not worth
# a clip-path's extra complexity for a schematic navigation map).
_QUADRANT_VIEWBOX_ORIGIN = {
    "I": (_CENTER, 0.0),
    "II": (0.0, 0.0),
    "III": (0.0, _CENTER),
    "IV": (_CENTER, _CENTER),
}


def sector_quadrant(x_pc, y_pc):
    """
    Classifies a galaxy-frame `(x, y)` position into its Quadrant --
    `theta = atan2(y, x)` normalized to `[0, 2*pi)`, split into four
    90-degree bands starting at `+X` (matching
    `docs/design/galaxy-coordinate-system.md`'s own theta convention).
    Blind to `z`/height by design -- a Quadrant is purely azimuthal, the
    same way real astronomical galactic quadrants are.

    Args:
        x_pc (float): `sectors.center_x_pc`.
        y_pc (float): `sectors.center_y_pc`.

    Returns:
        str: One of `QUADRANT_LABELS` (`"I"`, `"II"`, `"III"`, `"IV"`).
    """
    theta = math.atan2(y_pc, x_pc) % (2 * math.pi)
    index = min(3, int(theta // (math.pi / 2)))
    return QUADRANT_LABELS[index]


def sector_ring(shell_index):
    """The Ring index (a group of `RING_SHELL_WIDTH` consecutive shells)
    that `shell_index` falls in."""
    return shell_index // RING_SHELL_WIDTH


def ring_bounds_ly(ring_index):
    """
    The `(inner_ly, outer_ly)` light-year bounds of Ring `ring_index` --
    exact conversion of its shell-index band's own boundary radius (design
    doc section 3: shell `k` spans `[k*edge_pc, (k+1)*edge_pc)`), not the
    Ring's nominal/average radius.

    Returns:
        tuple[float, float]: `(inner_ly, outer_ly)`.
    """
    inner_ly = ring_index * RING_SHELL_WIDTH * DEFAULT_SECTOR_EDGE_LY
    outer_ly = (ring_index + 1) * RING_SHELL_WIDTH * DEFAULT_SECTOR_EDGE_LY
    return inner_ly, outer_ly


def _rings_to_show(sectors):
    """At least `_MIN_RINGS_SHOWN`, or enough to cover every placed
    sector's own Ring plus one extra empty Ring of context beyond it."""
    max_shell = max((s["shell_index"] for s in sectors if s["shell_index"] is not None), default=None)
    if max_shell is None:
        return _MIN_RINGS_SHOWN
    return max(_MIN_RINGS_SHOWN, sector_ring(max_shell) + 2)


def _star_visual(system_count):
    """
    Maps a sector's `system_count` to `(radius_px, halo_opacity,
    core_opacity)` -- square-root scaled (matching `starmap.py`'s own
    `_star_dot_radius` reasoning: a linear scale would make a sparse
    sector's dot vanish next to a dense one) and clamped so the map stays
    legible at either extreme, including the `system_count == 0` case (a
    sector generated with zero placed systems is still real, generated
    content, and still gets a dim dot rather than disappearing).
    """
    count = system_count or 0
    radius_px = max(_MIN_DOT_R, min(_MAX_DOT_R, _MIN_DOT_R + 2.4 * math.sqrt(count)))
    brightness = max(0.35, min(1.0, 0.35 + 0.18 * math.log2(count + 1)))
    return radius_px, brightness * 0.3, brightness


def _project(sector, px_per_ly):
    """
    Converts one placed sector's stored galaxy-frame position into SVG
    pixel coordinates plus its Quadrant/Ring/display distance.

    Args:
        sector (dict): `x`, `y`, `galactic_radius_pc`, `shell_index`
                       (`sectors.center_x_pc`/`center_y_pc`/
                       `galactic_radius_pc`/`shell_index`).
        px_per_ly (float): Pixels per light-year at the current scale.

    Returns:
        dict: `svg_x`, `svg_y`, `radius_ly`, `quadrant`.
    """
    theta = math.atan2(sector["y"], sector["x"])
    radius_ly = pc_to_ly(sector["galactic_radius_pc"]) if pc_to_ly else sector["galactic_radius_pc"] * 3.2616
    radius_px = radius_ly * px_per_ly
    # +Y is "up" in galaxy-frame math, same flip starmap.py's own map
    # panel applies once, here, at the one place a normalized position
    # becomes a screen coordinate (SVG's own Y axis grows downward).
    svg_x = _CENTER + radius_px * math.cos(theta)
    svg_y = _CENTER - radius_px * math.sin(theta)
    return {
        "svg_x": svg_x, "svg_y": svg_y, "radius_ly": radius_ly,
        "quadrant": sector_quadrant(sector["x"], sector["y"]),
    }


def _cloud_defs():
    """The purely-illustrative "expected density" radial gradient -- see
    this module's docstring. `var(--accent)` keeps it theme-aware (this
    SVG is always inline in the page, so it shares the document's CSS
    cascade) without needing a light/dark-specific color of its own."""
    return (
        '<radialGradient id="galaxyCloud" cx="50%" cy="50%" r="50%">'
        '<stop offset="0%" stop-color="var(--accent)" stop-opacity="0.5"/>'
        '<stop offset="35%" stop-color="var(--accent)" stop-opacity="0.22"/>'
        '<stop offset="70%" stop-color="var(--accent)" stop-opacity="0.08"/>'
        '<stop offset="100%" stop-color="var(--accent)" stop-opacity="0"/>'
        "</radialGradient>"
    )


def _ring_elements(rings_to_show, px_per_ly):
    parts = []
    for ring_index in range(rings_to_show):
        _inner_ly, outer_ly = ring_bounds_ly(ring_index)
        radius_px = outer_ly * px_per_ly
        parts.append(
            f'<circle class="galaxymap-ring" cx="{_CENTER:.1f}" cy="{_CENTER:.1f}" r="{radius_px:.1f}"/>'
        )
        parts.append(
            f'<text class="galaxymap-ring-label" x="{_CENTER + 4:.1f}" y="{_CENTER - radius_px - 3:.1f}">'
            f"~{outer_ly:,.0f} ly</text>"
        )
    return "".join(parts)


def _quadrant_axis_elements():
    return (
        f'<line class="galaxymap-axis" x1="{_CENTER - _MAP_RADIUS_PX:.1f}" y1="{_CENTER:.1f}" '
        f'x2="{_CENTER + _MAP_RADIUS_PX:.1f}" y2="{_CENTER:.1f}"/>'
        f'<line class="galaxymap-axis" x1="{_CENTER:.1f}" y1="{_CENTER - _MAP_RADIUS_PX:.1f}" '
        f'x2="{_CENTER:.1f}" y2="{_CENTER + _MAP_RADIUS_PX:.1f}"/>'
    )


def _quadrant_label_elements(db_name, active_quadrant):
    """One clickable Roman-numeral label per Quadrant, near its outer
    bisector -- omitted for `active_quadrant` itself (already the whole
    view, per `render_galaxy_map_panel`'s own "Full galaxy view" link
    instead) since a Quadrant view's viewBox crop clips the other three
    away anyway."""
    parts = []
    label_radius = _MAP_RADIUS_PX + 18
    bisectors = {"I": 45, "II": 135, "III": 225, "IV": 315}
    for label, angle_deg in bisectors.items():
        if label == active_quadrant:
            continue
        angle = math.radians(angle_deg)
        x = _CENTER + label_radius * math.cos(angle)
        y = _CENTER - label_radius * math.sin(angle)
        href = f"galaxy.py?db={esc(db_name)}&amp;quadrant={label}"
        parts.append(
            f'<a href="{href}"><text class="galaxymap-quadrant-label" x="{x:.1f}" y="{y:.1f}" '
            f'text-anchor="middle" dominant-baseline="middle">{label}</text></a>'
        )
    return "".join(parts)


_PHENOMENON_DOT_R = 3.2
"""float: Fixed dot radius (px) for a galaxy-placed nebula/asteroid field
-- unlike a sector's own dot (`_star_visual`, scaled by `system_count`),
these carry no analogous "how much is here" quantity worth scaling by, and
at galaxy scale a nebula's own real `radius_ly` (which can itself span
several sectors) would be a wildly misleading dot size anyway -- see this
module's own docstring on the Sector Map being the right place to depict
one's actual physical extent instead."""

_PHENOMENON_COLORS = {"nebula": "#c9a8e0", "asteroid_field": "#b89a6e"}
_PHENOMENON_TYPE_LABELS = {"nebula": "Nebula", "asteroid_field": "Asteroid Field"}
_DEFAULT_PHENOMENON_COLOR = "#9aa0ac"


def _phenomenon_elements(phenomena, px_per_ly):
    """
    Plots every galaxy-placed nebula/asteroid field as a small, fixed-size
    dot -- the galaxy-scale counterpart to `_star_elements`, but with no
    click-through (a phenomenon has no detail page of its own, unlike a
    sector) -- just a hover tooltip naming it, its type, and its real size.

    Args:
        phenomena (list[dict]): `queryDb.galaxy_placed_phenomena`'s return
                                shape (`id`, `type`, `name`, `descriptor`,
                                `radius_ly`, `x`/`y`/`z`, `galactic_radius_pc`).
        px_per_ly (float): Pixels per light-year at the current scale.

    Returns:
        str: One `<g>` per phenomenon.
    """
    parts = []
    for phenomenon in phenomena:
        projected = _project(phenomenon, px_per_ly)
        x, y = projected["svg_x"], projected["svg_y"]
        color = _PHENOMENON_COLORS.get(phenomenon["type"], _DEFAULT_PHENOMENON_COLOR)
        type_label = _PHENOMENON_TYPE_LABELS.get(phenomenon["type"], phenomenon["type"])
        descriptor = (phenomenon["descriptor"] or "").capitalize()
        tooltip = (
            f'{esc(phenomenon["name"])} -- {esc(descriptor)} {esc(type_label)}, '
            f'~{phenomenon["radius_ly"]:,.1f} ly across, ~{projected["radius_ly"]:,.0f} ly from core'
        )
        parts.append(
            f'<g class="galaxymap-phenomenon"><title>{tooltip}</title>'
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{_PHENOMENON_DOT_R:.1f}" '
            f'fill="{color}" fill-opacity="0.85" stroke="#00000055" stroke-width="0.6"/>'
            "</g>"
        )
    return "".join(parts)


def _star_elements(db_name, sectors, px_per_ly):
    parts = []
    for sector in sectors:
        projected = _project(sector, px_per_ly)
        radius_px, halo_opacity, core_opacity = _star_visual(sector["system_count"])
        href = f'sector.py?db={esc(db_name)}&amp;id={sector["id"]}'
        tooltip = (
            f'{esc(sector["name"])} -- {sector["system_count"] or 0} system'
            f'{"s" if sector["system_count"] != 1 else ""}, '
            f'~{projected["radius_ly"]:,.0f} ly from core, Quadrant {projected["quadrant"]}'
        )
        x, y = projected["svg_x"], projected["svg_y"]
        parts.append(
            f'<a class="galaxymap-star" href="{href}"><title>{tooltip}</title>'
            f'<circle class="galaxymap-star-halo" cx="{x:.1f}" cy="{y:.1f}" r="{radius_px * 2.2:.1f}" '
            f'fill="{_STAR_HALO_FILL}" opacity="{halo_opacity:.2f}"/>'
            f'<circle class="galaxymap-star-core" cx="{x:.1f}" cy="{y:.1f}" r="{radius_px:.1f}" '
            f'fill="{_STAR_CORE_FILL}" stroke="{_STAR_CORE_STROKE}" opacity="{core_opacity:.2f}"/>'
            "</a>"
        )
    return "".join(parts)


def render_galaxy_map_panel(db_name, sectors, quadrant=None, phenomena=None):
    """
    Builds the "Galaxy Map" panel: a flat SVG plot of every galaxy-placed
    sector (a bright dot, sized/colored by `system_count`) inside four
    Quadrants and concentric Rings, over an illustrative "expected
    density" cloud gradient -- see this module's docstring.

    Args:
        db_name (str): The current `?db=` value, used to build every link.
        sectors (list[dict]): One entry per galaxy-placed sector (non-NULL
                              `center_x_pc`), each with `id`, `name`, `x`,
                              `y`, `galactic_radius_pc`, `shell_index`,
                              `system_count`.
        quadrant (str or None): One of `QUADRANT_LABELS` to zoom the map
                                into (a `viewBox` crop -- see this module's
                                docstring), or `None` for the full galaxy.
        phenomena (list[dict] or None): `queryDb.galaxy_placed_phenomena`'s
                                return shape -- every galaxy-placed nebula/
                                asteroid field, plotted as a small fixed-
                                size dot (`_phenomenon_elements`); its own
                                real physical size is instead depicted on
                                the Sector Map (`html/lib/starmap.py`) of
                                any sector its sphere reaches into. `None`/
                                empty plots none.

    Returns:
        str: A complete `<section class="panel">` block.
    """
    rings_to_show = _rings_to_show(sectors)
    _inner, outer_ly = ring_bounds_ly(rings_to_show - 1)
    px_per_ly = _MAP_RADIUS_PX / outer_ly if outer_ly else 1.0

    body = (
        f"<defs>{_cloud_defs()}</defs>"
        f'<circle cx="{_CENTER:.1f}" cy="{_CENTER:.1f}" r="{_MAP_RADIUS_PX:.1f}" fill="url(#galaxyCloud)"/>'
        f"{_ring_elements(rings_to_show, px_per_ly)}"
        f"{_quadrant_axis_elements()}"
        f"{_quadrant_label_elements(db_name, quadrant)}"
        f"{_phenomenon_elements(phenomena or [], px_per_ly)}"
        f"{_star_elements(db_name, sectors, px_per_ly)}"
    )

    if quadrant in QUADRANT_LABELS:
        vx, vy = _QUADRANT_VIEWBOX_ORIGIN[quadrant]
        view_box = f"{vx:.1f} {vy:.1f} {_CENTER:.1f} {_CENTER:.1f}"
        scope_label = f"Quadrant {quadrant}"
        nav_html = f'<p class="hint"><a href="galaxy.py?db={esc(db_name)}">&larr; Full galaxy view</a></p>'
    else:
        view_box = f"0 0 {_SVG_SIZE} {_SVG_SIZE}"
        scope_label = "Full galaxy"
        nav_html = ""

    svg = (
        f'<svg class="galaxymap-svg" viewBox="{view_box}" role="img" '
        f'aria-label="Galaxy map, {esc(scope_label)}. Click a bright dot for its sector, '
        f'a Roman-numeral label to jump to that Quadrant.">{body}</svg>'
    )

    if not sectors:
        legend_extra = '<p class="hint">No sectors have been placed in the galaxy yet -- see galaxyGen.py.</p>'
    else:
        legend_extra = ""

    return f"""
<section class="panel">
<div class="panel-header">
  <h2>Galaxy Map -- {esc(scope_label)}</h2>
  <span class="hint">Dot size/brightness &asymp; systems in that sector &middot; small purple/tan dots &asymp; nebulae/asteroid fields (hover for details) &middot; shading &asymp; illustrative expected density, not real data</span>
</div>
<div class="galaxymap-layout">
<div class="galaxymap-viewport">
{svg}
</div>
<div class="galaxymap-side">
{nav_html}
<p class="hint">Rings are fixed {RING_SHELL_WIDTH}-shell bands (~{RING_TARGET_LY:.0f} ly each); Quadrants I-IV split the galaxy into four 90-degree azimuthal arcs from the core.</p>
{legend_extra}
</div>
</div>
</section>
"""
