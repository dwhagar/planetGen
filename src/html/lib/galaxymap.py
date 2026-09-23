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
has. Drill-down is plain server-side navigation (a different `quadrant`
parameter), still no map-specific client-side JS of its own:
`render_galaxy_map_panel` always draws every placed sector into the same
fixed coordinate space, and a Quadrant view simply crops the same drawing
to one quarter of it via the outer `<svg>`'s own `viewBox` (SVG clips to
its viewBox by default). Each marker still needs `static/navform.js`
(loaded on every page, see `lib/page.py`'s `render`) to actually navigate,
though -- a `<form>` can't nest inside an SVG shape the way `page.
post_link` uses one everywhere else in `html/`, so a marker stays a real
`<a>` carrying `data-nav-target`/`data-nav-params` instead (`fmt.
data_nav_params`), intercepted by that shared script's click handler
rather than putting `sector.py?...`/`phenomenon.py?...` in the browser's
own address bar the way following its `href` directly would.

Density is two genuinely different things here, deliberately drawn two
different ways:

- **Real, generated sectors** -- bright dots, sized/colored by their own
  actual `system_count`.
- **Un-generated space** -- shaded by the galaxy's real predicted density,
  when it's known. `generate.py plan` computes and stores a singleton
  `galaxy_shape` row (`stellarObjects.galaxyDensity.GalaxyShape` -- the
  exponential-disk-plus-bulge-plus-spiral-arm model, see
  `docs/design/galaxy-disk-density.md`) that already gates/weights actual
  sector generation (`generate.py`'s `_BatchDensity`); `html/galaxy.py`
  passes that same shape down here (`galaxy_shape`, `queryDb.
  galaxy_density_shape`/`GET /api/galaxy/shape`) so `_expected_density_elements`
  can evaluate `galaxyDensity.relative_density` across the visible disk
  and shade a grid of small tiles by it -- un-filled space still reads as
  the actual predicted spiral, not a generic radial glow. If the skeleton
  has never been built (`galaxy_shape is None` -- `generate.py plan` was
  never run), this falls back to `_cloud_defs`'s original soft radial
  "cloud" gradient (brighter near the core, fading outward) -- purely
  illustrative shading, not a real population model, same as before this
  module read any actual density. Either way this exists so an
  overwhelmingly-empty galaxy -- the normal, expected state, since the
  addressable volume runs into the hundreds of billions of sector slots --
  doesn't render as a blank void.
"""

import math

from fmt import data_nav_params, esc, post_link

try:
    from stellarObjects.galaxyDensity import GalaxyShape, relative_density
    from stellarObjects.program_constants import DEFAULT_SECTOR_EDGE_LY
    from stellarObjects.utils import ly_to_pc, pc_to_ly
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # duplicated fallback, matching starmap.py's own pattern (dbutil.py
    # itself no longer has one -- stellarObjects is a hard dependency
    # there now that its MySQL connection helpers are load-bearing).
    # `GalaxyShape`/`relative_density` left `None` -- `_expected_density_elements`
    # treats that the same as "skeleton never built" and falls back to the
    # illustrative gradient, same as a deployment where it's just missing data.
    DEFAULT_SECTOR_EDGE_LY = 11.5
    pc_to_ly = None
    ly_to_pc = None
    GalaxyShape = None
    relative_density = None

QUADRANT_LABELS = ("I", "II", "III", "IV")
"""tuple[str]: The four galaxy-scale Quadrants, in azimuthal order starting
from `+X` -- `sector_quadrant` indexes into this. Deliberately the plain
2D I-IV convention (see this module's own docstring) rather than the
sector-internal octant's Roman-numeral scheme."""

_MIN_VIEW_SIZE_LY = 100.0
"""float: The deepest zoom-in this map allows, as a full viewBox-edge span
in light-years (not a radius) -- static/mapzoom.js clamps its own zoom
range to this, converted to SVG units via each render's own `px_per_ly`
(see `render_galaxy_map_panel`). Matches the ~100 ly scale
`RING_TARGET_LY` already uses for a Ring's own thickness -- deliberately
the same number (not a coincidence to reconcile): "zoom in to about one
Ring's width" is a natural, already-established scale reference on this
same map, not an arbitrary new one."""

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


_SPIRAL_VIEW_SCALE_LENGTHS = 2.0
"""float: How many `disk_scale_length_pc`s out the default (un-zoomed)
view reaches when the galaxy's real density model is known -- enough for
its spiral arms to complete a visible wind or two (see
`docs/design/galaxy-disk-density.md`'s own pitch-angle worked example),
so the map reads as a spiral by default even when almost nothing has been
generated yet. Without this, `_rings_to_show`'s sector-driven minimum
alone (`_MIN_RINGS_SHOWN`, ~300 ly) stays well inside the "trivially
solid" bulge core for any Milky-Way-scale shape (that design doc found it
extends to ~2,770 ly), so the density cloud would render as a uniform
bright disc -- accurate for that tiny a patch, but never showing the
actual spiral shape at all. Capped by the model's own real, stored outer
edge (`outer_shell_index`/`edge_pc`) when that's smaller, so a small toy
galaxy never gets zoomed out past its own true extent."""


def _rings_to_show(sectors, galaxy_shape=None):
    """At least `_MIN_RINGS_SHOWN`, or enough to cover every placed
    sector's own Ring plus one extra empty Ring of context beyond it --
    or, when `galaxy_shape` is known, enough to reach
    `_SPIRAL_VIEW_SCALE_LENGTHS` of its real disk scale length (capped at
    its own stored outer edge) if that's farther out still, so the
    "expected density" cloud actually has spiral structure to show by
    default -- see `_SPIRAL_VIEW_SCALE_LENGTHS`'s own docstring."""
    max_shell = max((s["shell_index"] for s in sectors if s["shell_index"] is not None), default=None)
    sector_rings = _MIN_RINGS_SHOWN if max_shell is None else max(_MIN_RINGS_SHOWN, sector_ring(max_shell) + 2)

    if not galaxy_shape:
        return sector_rings

    spiral_reach_pc = _SPIRAL_VIEW_SCALE_LENGTHS * galaxy_shape["disk_scale_length_pc"]
    outer_shell_index = galaxy_shape.get("outer_shell_index")
    edge_pc = galaxy_shape.get("edge_pc")
    if outer_shell_index is not None and edge_pc is not None:
        spiral_reach_pc = min(spiral_reach_pc, (outer_shell_index + 1) * edge_pc)

    spiral_reach_ly = pc_to_ly(spiral_reach_pc) if pc_to_ly else spiral_reach_pc * 3.2616
    spiral_rings = math.ceil(spiral_reach_ly / (RING_SHELL_WIDTH * DEFAULT_SECTOR_EDGE_LY))
    return max(sector_rings, spiral_rings)


_MAX_RINGS_DRAWN = 10
"""int: Never actually draw more than this many ring guide-circles/labels,
however many fixed-shell-width Rings a far-out placed sector's own real
distance implies (`_rings_to_show`, used unchanged to pick the map's scale
so that sector still fits) -- a sector placed deep in the galaxy can imply
hundreds of Rings (confirmed directly: shell_index ~1400 implied 157), which
only clutters the map into a dense, unreadable smear of overlapping guide
circles and labels instead of giving useful distance context. Past this
cap, `_ring_elements` switches from one guide per literal Ring to
`_MAX_RINGS_DRAWN` evenly-spaced distance markers covering the same span
instead -- still real, accurate distance labels, just no longer tied 1:1 to
`RING_SHELL_WIDTH`'s own fixed shell count once there'd be too many to
read."""


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


_DENSITY_GRID_CELLS = 64
"""int: Tiles per side of the "expected density" grid (see
`_expected_density_elements`) -- a compromise between a visibly-spiral
shape (needs enough resolution to show arm/inter-arm contrast rather than
smearing it into a uniform ring) and how many `<rect>` elements one page
render can afford (`_DENSITY_GRID_CELLS**2`, before the circular-crop skip
below removes the ~21% that fall in the square's corners -- a few thousand
either way, trivial for both server-side generation and the browser)."""

_DENSITY_MAX_OPACITY = 0.6
"""float: Cap on a single density tile's own `fill-opacity` -- keeps even
the single brightest tile in view from fully obscuring a "star" dot drawn
on top of it (same reasoning `_cloud_defs`'s own 0.5 center-stop cap
followed)."""

_DENSITY_RADIAL_GAMMA = 0.4
"""float: Exponent `_expected_density_elements` applies to each tile's
*azimuthally-averaged* density (`_radial_baseline_density` -- the same
`relative_density` formula with the spiral-arm term's own azimuthal
average substituted for its real value, i.e. `arm_factor -> 1`),
normalized against the single brightest tile currently in view, to get
that tile's smooth radial "how close to the core" glow. `relative_density`
falls off steeply (bulge/disk both exponential in radius, see
`docs/design/galaxy-disk-density.md`), so a linear (`gamma=1`) scale would
leave everything past the inner bulge reading as barely-there -- `< 1` is
a standard display-only contrast stretch (the same idea a telescope
image's own "curves" adjustment applies before its real structure is
visible by eye), same reasoning `_star_visual`'s own `log2` scale follows
for a sector dot's size. Kept separate from the spiral-arm contrast itself
(`_DENSITY_ARM_CONTRAST`) specifically so gamma-stretching one doesn't
also distort the other: the radius-spanning falloff and the arm/inter-arm
swing *at one radius* are different-sized effects (the model's own
worked example puts the swing at `~2.33x`, tiny next to the bulge-to-edge
falloff), so stretching both with one shared exponent leaves whichever
effect is smaller looking flat -- see `_expected_density_elements`'s own
docstring."""

_DENSITY_ARM_CONTRAST = 1.0
"""float: How strongly a tile's real spiral-arm modulation (`density /
_radial_baseline_density` at that same position -- exactly `arm_factor`,
`galaxyDensity`'s own `[1 - arm_amplitude, 1 + arm_amplitude]` range, once
the bulge-dilution near the core has been divided back out) multiplies
its radial glow (`_DENSITY_RADIAL_GAMMA`'s output) to get the final tile
opacity. `1.0` applies that real ratio unscaled -- already a genuine,
data-driven `~2.33x` arm/inter-arm swing for this model's own worked
example (`docs/design/galaxy-disk-density.md`), not fabricated contrast;
a value here would only ever be tuned up for a shape whose own
`arm_amplitude` is too subtle to read as a spiral at a glance, never used
to invent structure the model doesn't actually predict."""

_DENSITY_MIN_OPACITY = 0.01
"""float: Below this, `_expected_density_elements` skips emitting the tile
entirely rather than adding a practically-invisible `<rect>` -- keeps the
outer, empty-inter-arm-trough majority of the grid from bloating the SVG
with elements no one can see."""


def _radial_baseline_density(shape, r_pc):
    """
    `galaxyDensity._raw_density` at `(r_pc, 0, 0)` (disk plane, `z=0` --
    this map's own projection), but with the spiral-arm term's own
    azimuthal average (`arm_factor`'s mean over a full `theta` revolution
    is exactly `1`, since it's `1 + amplitude * cos(...)` and cosine
    averages to `0`) substituted for its real, `theta`-dependent value --
    i.e. this shape's smooth bulge+disk envelope *without* any spiral
    structure riding on it, at whatever radius a real tile's own
    `relative_density` (which *does* include the real arm term) can be
    compared against to isolate that tile's own arm/inter-arm contrast
    (`_expected_density_elements`: `density / _radial_baseline_density(...)
    == arm_factor` exactly, algebraically, once the shared bulge/disk
    terms cancel). Public formula, not a private `galaxyDensity` internal
    reused out of turn -- restated directly from
    `docs/design/galaxy-disk-density.md` section 1's own spec, the same
    one `galaxyDensity._raw_density`/`relative_density` implement.

    Args:
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        r_pc (float): In-plane (and, since `z=0`, also 3D) galactocentric
            radius, parsecs -- always `>= 0`.

    Returns:
        float: `>= 0`, normalized the same way `relative_density` is
            (`shape.k_norm` applied).
    """
    bulge = shape.bulge_amplitude * math.exp(-r_pc / shape.bulge_scale_radius_pc)
    disk_radial = math.exp(-r_pc / shape.disk_scale_length_pc)
    return shape.k_norm * (bulge + disk_radial)


def _expected_density_elements(galaxy_shape, px_per_ly):
    """
    The real "expected density" shading -- a `_DENSITY_GRID_CELLS`-square
    grid of small tiles, each shaded by `galaxyDensity.relative_density`
    at that tile's own galaxy-frame position (the disk plane, `z=0`: this
    map is a flat face-on projection, the same plane `sector_quadrant`
    already classifies by). Each tile's opacity is a smooth radial glow
    (`_radial_baseline_density`, peak-normalized and gamma-stretched --
    `_DENSITY_RADIAL_GAMMA`) multiplied by that tile's own real spiral-arm
    contrast (`_DENSITY_ARM_CONTRAST`) -- split this way, rather than
    gamma-stretching the raw `relative_density` directly, specifically
    because the radius-spanning falloff (bulge peak down to near-zero at
    the edge) is a far bigger effect than the arm/inter-arm swing *at any
    one radius*, so a single shared stretch flattens the arms into
    invisibility even though they're really there. Un-generated space this
    way still reads as the galaxy's actual
    predicted spiral shape -- bulge core, two arms, inter-arm troughs --
    rather than a generic radial glow, wherever that shape is actually
    known.

    Args:
        galaxy_shape (dict or None): `queryDb.galaxy_density_shape`'s
            return shape (every `galaxyDensity.GalaxyShape` field, plus a
            few this function ignores) -- `None` if `generate.py plan`
            has never been run against this database.
        px_per_ly (float): Pixels per light-year at the current scale
            (same value `_project`/`_ring_elements` already use), so a
            tile's screen position converts back to a real galaxy-frame
            position.

    Returns:
        str or None: A `<defs>` clip-path plus one clipped `<g>` of
            `<rect>` tiles, or `None` if `galaxy_shape` is `None` (no
            skeleton built yet), every in-view tile came back at a
            non-positive baseline density (never actually happens -- the
            model has no hard-zero region -- but guarded rather than
            assumed), or `stellarObjects` isn't importable in this
            deployment (see this module's own top-of-file try/except) --
            any of those, the caller falls back to `_cloud_defs`'s
            illustrative gradient instead.
    """
    if galaxy_shape is None or GalaxyShape is None:
        return None

    shape = GalaxyShape(**{field: galaxy_shape[field] for field in GalaxyShape._fields})
    cell_px = _SVG_SIZE / _DENSITY_GRID_CELLS
    max_r_px = _MAP_RADIUS_PX + cell_px  # generous margin so a tile centered just outside the
                                          # ring still gets clipped cleanly rather than cut mid-tile

    # First pass: every in-view tile's own real density and radial baseline, plus the
    # brightest baseline found -- the radial glow below needs that peak to normalize against.
    tiles = []
    peak_baseline = 0.0
    for row in range(_DENSITY_GRID_CELLS):
        cy = (row + 0.5) * cell_px
        dy_px = cy - _CENTER
        for col in range(_DENSITY_GRID_CELLS):
            cx = (col + 0.5) * cell_px
            dx_px = cx - _CENTER
            if dx_px * dx_px + dy_px * dy_px > max_r_px * max_r_px:
                continue

            x_ly, y_ly = dx_px / px_per_ly, -dy_px / px_per_ly  # +Y up in galaxy-frame, see _project
            x_pc = ly_to_pc(x_ly) if ly_to_pc else x_ly / 3.2616
            y_pc = ly_to_pc(y_ly) if ly_to_pc else y_ly / 3.2616

            density = relative_density((x_pc, y_pc, 0.0), shape)
            baseline = _radial_baseline_density(shape, math.hypot(x_pc, y_pc))
            tiles.append((cx, cy, density, baseline))
            if baseline > peak_baseline:
                peak_baseline = baseline

    if not tiles or peak_baseline <= 0:
        return None

    parts = [
        f'<defs><clipPath id="galaxyDensityClip">'
        f'<circle cx="{_CENTER:.1f}" cy="{_CENTER:.1f}" r="{_MAP_RADIUS_PX:.1f}"/>'
        f"</clipPath></defs>",
        '<g clip-path="url(#galaxyDensityClip)">',
    ]
    for cx, cy, density, baseline in tiles:
        radial_glow = max(0.0, baseline / peak_baseline) ** _DENSITY_RADIAL_GAMMA
        arm_contrast = (max(0.0, density) / baseline) ** _DENSITY_ARM_CONTRAST if baseline > 0 else 1.0
        opacity = min(_DENSITY_MAX_OPACITY, _DENSITY_MAX_OPACITY * radial_glow * arm_contrast)
        if opacity < _DENSITY_MIN_OPACITY:
            continue
        parts.append(
            f'<rect x="{cx - cell_px / 2:.2f}" y="{cy - cell_px / 2:.2f}" '
            f'width="{cell_px:.2f}" height="{cell_px:.2f}" fill="var(--accent)" '
            f'fill-opacity="{opacity:.3f}"/>'
        )
    parts.append("</g>")
    return "".join(parts)


def _ring_elements(rings_to_show, px_per_ly):
    if rings_to_show <= _MAX_RINGS_DRAWN:
        # The common case (browsing near the core): few enough real Rings
        # to show each one exactly, at its own fixed-shell-width boundary.
        outer_lys = [ring_bounds_ly(ring_index)[1] for ring_index in range(rings_to_show)]
    else:
        # Too many real Rings to draw individually -- fall back to
        # `_MAX_RINGS_DRAWN` evenly-spaced distance markers spanning the
        # same 0..outer_ly range instead (see `_MAX_RINGS_DRAWN`'s own
        # docstring), rather than one per literal Ring.
        _inner_ly, full_outer_ly = ring_bounds_ly(rings_to_show - 1)
        outer_lys = [full_outer_ly * (i + 1) / _MAX_RINGS_DRAWN for i in range(_MAX_RINGS_DRAWN)]

    parts = []
    for outer_ly in outer_lys:
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
        nav_params = data_nav_params({"db": db_name, "quadrant": label})
        parts.append(
            f'<a href="#" data-nav-target="galaxy.py" data-nav-params="{nav_params}">'
            f'<text class="galaxymap-quadrant-label" x="{x:.1f}" y="{y:.1f}" '
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

_PHENOMENON_COLORS = {
    "nebula": "#c9a8e0", "asteroid_field": "#b89a6e",
    # v21: black_holes/neutron_stars gained the same galaxy-frame placement
    # nebulae/asteroid_fields already had -- see queryDb._PHENOMENON_TABLES
    # and schema.sql's "v21" header note. Colors echo starmap.py's own
    # Sector Map markers for the same two types.
    "black_hole": "#1a1a1a", "neutron_star": "#cfe8ff",
}
_PHENOMENON_TYPE_LABELS = {
    "nebula": "Nebula", "asteroid_field": "Asteroid Field",
    "black_hole": "Black Hole", "neutron_star": "Neutron Star",
}
_DEFAULT_PHENOMENON_COLOR = "#9aa0ac"


def _phenomenon_elements(db_name, phenomena, px_per_ly):
    """
    Plots every galaxy-placed nebula/asteroid field/black hole/neutron
    star as a small, fixed-size dot -- the galaxy-scale counterpart to
    `_star_elements`, linking to `phenomenon.py` (this project's detail
    page for a standalone phenomenon) the same way a sector's own dot
    links to `sector.py` -- plus a hover tooltip naming it, its type, and
    its real size (a black hole/neutron star's own `radius_ly` is always
    0 -- point-like at this scale, see `queryDb._PHENOMENON_TABLES` -- so
    the tooltip omits the "ly across" figure for those two rather than
    showing a misleading "~0.0 ly across").

    Args:
        db_name (str): The current `?db=` value, used to build each dot's
                       `href`.
        phenomena (list[dict]): `queryDb.galaxy_placed_phenomena`'s return
                                shape (`id`, `type`, `name`, `descriptor`,
                                `radius_ly`, `x`/`y`/`z`, `galactic_radius_pc`).
        px_per_ly (float): Pixels per light-year at the current scale.

    Returns:
        str: One `<a>` per phenomenon.
    """
    parts = []
    for phenomenon in phenomena:
        projected = _project(phenomenon, px_per_ly)
        x, y = projected["svg_x"], projected["svg_y"]
        color = _PHENOMENON_COLORS.get(phenomenon["type"], _DEFAULT_PHENOMENON_COLOR)
        type_label = _PHENOMENON_TYPE_LABELS.get(phenomenon["type"], phenomenon["type"])
        descriptor = (phenomenon["descriptor"] or "").capitalize()
        size_bit = f'~{phenomenon["radius_ly"]:,.1f} ly across, ' if phenomenon["radius_ly"] else ""
        tooltip = (
            f'{esc(phenomenon["name"])} -- {esc(descriptor)} {esc(type_label)}, '
            f'{size_bit}~{projected["radius_ly"]:,.0f} ly from core'
        )
        nav_params = data_nav_params({"db": db_name, "type": phenomenon["type"], "id": phenomenon["id"]})
        parts.append(
            f'<a class="galaxymap-phenomenon" href="#" data-nav-target="phenomenon.py" '
            f'data-nav-params="{nav_params}"><title>{tooltip}</title>'
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{_PHENOMENON_DOT_R:.1f}" '
            f'fill="{color}" fill-opacity="0.85" stroke="#00000055" stroke-width="0.6"/>'
            "</a>"
        )
    return "".join(parts)


def _star_elements(db_name, sectors, px_per_ly):
    parts = []
    for sector in sectors:
        projected = _project(sector, px_per_ly)
        radius_px, halo_opacity, core_opacity = _star_visual(sector["system_count"])
        nav_params = data_nav_params({"db": db_name, "id": sector["id"]})
        tooltip = (
            f'{esc(sector["name"])} -- {sector["system_count"] or 0} system'
            f'{"s" if sector["system_count"] != 1 else ""}, '
            f'~{projected["radius_ly"]:,.0f} ly from core, Quadrant {projected["quadrant"]}'
        )
        x, y = projected["svg_x"], projected["svg_y"]
        parts.append(
            f'<a class="galaxymap-star" href="#" data-nav-target="sector.py" '
            f'data-nav-params="{nav_params}"><title>{tooltip}</title>'
            f'<circle class="galaxymap-star-halo" cx="{x:.1f}" cy="{y:.1f}" r="{radius_px * 2.2:.1f}" '
            f'fill="{_STAR_HALO_FILL}" opacity="{halo_opacity:.2f}"/>'
            f'<circle class="galaxymap-star-core" cx="{x:.1f}" cy="{y:.1f}" r="{radius_px:.1f}" '
            f'fill="{_STAR_CORE_FILL}" stroke="{_STAR_CORE_STROKE}" opacity="{core_opacity:.2f}"/>'
            "</a>"
        )
    return "".join(parts)


def render_galaxy_map_panel(db_name, sectors, quadrant=None, phenomena=None, galaxy_shape=None):
    """
    Builds the "Galaxy Map" panel: a flat SVG plot of every galaxy-placed
    sector (a bright dot, sized/colored by `system_count`) inside four
    Quadrants and concentric Rings, over an "expected density" cloud --
    see this module's docstring, and `galaxy_shape` below.

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
                                asteroid field/black hole/neutron star,
                                plotted as a small fixed-size dot
                                (`_phenomenon_elements`); a nebula/asteroid
                                field's own real physical size is instead
                                depicted on the Sector Map
                                (`html/lib/starmap.py`) of any sector its
                                sphere reaches into (a black hole/neutron
                                star is point-like even there). `None`/
                                empty plots none.
        galaxy_shape (dict or None): `queryDb.galaxy_density_shape`'s
                                return shape -- the galaxy's real, stored
                                disk/bulge/spiral-arm density model
                                (`generate.py plan`'s output). When given,
                                un-generated space is shaded by this
                                model's own `relative_density`
                                (`_expected_density_elements`) instead of
                                the generic illustrative gradient, so the
                                map still reads as the predicted spiral
                                even where nothing has been generated yet.
                                `None` (the skeleton was never built) falls
                                back to that illustrative gradient.

    Returns:
        str: A complete `<section class="panel">` block.
    """
    rings_to_show = _rings_to_show(sectors, galaxy_shape)
    _inner, outer_ly = ring_bounds_ly(rings_to_show - 1)
    px_per_ly = _MAP_RADIUS_PX / outer_ly if outer_ly else 1.0

    density_html = _expected_density_elements(galaxy_shape, px_per_ly)
    density_is_real = density_html is not None
    if not density_is_real:
        density_html = (
            f"<defs>{_cloud_defs()}</defs>"
            f'<circle cx="{_CENTER:.1f}" cy="{_CENTER:.1f}" r="{_MAP_RADIUS_PX:.1f}" fill="url(#galaxyCloud)"/>'
        )

    body = (
        f"{density_html}"
        f"{_ring_elements(rings_to_show, px_per_ly)}"
        f"{_quadrant_axis_elements()}"
        f"{_quadrant_label_elements(db_name, quadrant)}"
        f"{_phenomenon_elements(db_name, phenomena or [], px_per_ly)}"
        f"{_star_elements(db_name, sectors, px_per_ly)}"
    )

    if quadrant in QUADRANT_LABELS:
        vx, vy = _QUADRANT_VIEWBOX_ORIGIN[quadrant]
        view_box = f"{vx:.1f} {vy:.1f} {_CENTER:.1f} {_CENTER:.1f}"
        scope_label = f"Quadrant {quadrant}"
        nav_html = f'<p class="hint">{post_link("galaxy.py", {"db": db_name}, "&larr; Full galaxy view")}</p>'
    else:
        view_box = f"0 0 {_SVG_SIZE} {_SVG_SIZE}"
        scope_label = "Full galaxy"
        nav_html = ""

    # Real interactive zoom/pan (static/mapzoom.js), not just the quadrant
    # crop's own fixed viewBox above -- that crop only decides where this
    # map STARTS; from there, scroll/drag/the +/-/Reset buttons let a
    # visitor zoom in far past it. `_MIN_VIEW_SIZE_LY` (~100 ly across, via
    # this module's own `RING_TARGET_LY`) converted to this render's own
    # `px_per_ly` is the deepest zoom-in this map allows -- independent of
    # `quadrant`, so the same 100 ly floor applies whether zoom starts from
    # the full galaxy or a Quadrant crop. Floored/capped so a sparse galaxy
    # (whose own full extent is already under 100 ly) can't compute a
    # min-zoom LARGER than the starting view, which would make "zoom in"
    # zoom out instead.
    min_view_size = min(_SVG_SIZE, max(20.0, _MIN_VIEW_SIZE_LY * 2 * px_per_ly))

    svg = (
        f'<svg class="galaxymap-svg" id="galaxymap-svg" viewBox="{view_box}" role="img" '
        f'data-min-view-size="{min_view_size:.2f}" data-max-view-size="{_SVG_SIZE:.2f}" '
        f'data-px-per-ly="{px_per_ly:.6f}" '
        f'aria-label="Galaxy map, {esc(scope_label)}. Scroll or use the +/- buttons to zoom, '
        f'drag to pan. Click a bright dot for its sector, a Roman-numeral label to jump to that '
        f'Quadrant.">{body}</svg>'
    )

    hint_parts = []
    if not sectors:
        hint_parts.append('<p class="hint">No sectors have been placed in the galaxy yet -- see galaxyGen.py.</p>')
    if not density_is_real:
        hint_parts.append(
            '<p class="hint">The galaxy\'s density skeleton hasn\'t been built yet '
            "(<code>generate.py plan</code>) -- shading below is illustrative only, not the galaxy's "
            "real predicted density.</p>"
        )
    legend_extra = "".join(hint_parts)

    density_hint = (
        "shading &asymp; the galaxy's real predicted stellar density (generate.py plan), including "
        "space not generated yet"
        if density_is_real
        else "shading &asymp; illustrative expected density, not real data"
    )

    return f"""
<section class="panel">
<div class="panel-header">
  <h2>Galaxy Map -- {esc(scope_label)}</h2>
  <span class="hint">Dot size/brightness &asymp; systems in that sector &middot; small purple/tan/dark/blue dots &asymp; nebulae/asteroid fields/black holes/neutron stars (hover for details) &middot; {density_hint}</span>
</div>
<div class="galaxymap-layout">
<div class="galaxymap-viewport">
{svg}
<div class="starmap-scale" id="galaxymap-scale"></div>
</div>
<div class="galaxymap-side">
<div class="starmap-controls" id="galaxymap-controls">
  <button type="button" class="starmap-btn" data-action="zoom-out" aria-label="Zoom out">&minus;</button>
  <button type="button" class="starmap-btn" data-action="zoom-in" aria-label="Zoom in">+</button>
  <button type="button" class="starmap-btn" data-action="reset">Reset view</button>
</div>
{nav_html}
<p class="hint">Rings are fixed {RING_SHELL_WIDTH}-shell bands (~{RING_TARGET_LY:.0f} ly each); Quadrants I-IV split the galaxy into four 90-degree azimuthal arcs from the core. Scroll/drag/+/- to zoom in, down to about {_MIN_VIEW_SIZE_LY:.0f} ly wide.</p>
{legend_extra}
</div>
</div>
</section>
<script src="static/mapzoom.js"></script>
<script src="static/galaxymap.js"></script>
"""
