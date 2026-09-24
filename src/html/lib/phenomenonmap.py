# html/lib/phenomenonmap.py

"""
Stellar phenomenon diagram: a flat, zoomable SVG showing one standalone
phenomenon's own real physical extent, drawn directly to astronomical-unit
(AU) scale -- `phenomenon.py`'s counterpart to the Galaxy Map/Sector Map,
for the one thing neither of those can show: what a nebula/asteroid field/
supernova remnant's own real size actually looks like next to a familiar
AU-scale yardstick, rather than as a same-size dot on a galaxy-wide plot
(`html/lib/galaxymap.py`'s own docstring explains why that map deliberately
never scales a phenomenon's dot by its real `radius_ly`).

Same `viewBox`-is-the-camera approach as `html/lib/galaxymap.py` (see
`static/mapzoom.js`'s own module docstring), but simpler: this diagram's
own SVG user units ARE astronomical units directly, one-to-one -- no
separate px-per-ly scale factor is needed the way the Galaxy Map's much
larger, dynamically-sized canvas requires, since this diagram's dynamic
range (1 AU to 1 ly, i.e. ~63,241x -- the user-facing zoom bounds, matching
the ~1 AU to ~1 ly range explicitly asked for) is already anchored to one
fixed, real physical unit end to end.

A black hole/neutron star has no `radius_ly` column of its own (its real
size -- an event horizon a few km across, a neutron star ~10-20 km -- is
utterly negligible at this diagram's own AU scale, exactly as at the
Galaxy/Sector Map's much larger scales); it renders as a small fixed dot
rather than a to-scale circle, same convention `galaxymap.py`'s own
`_PHENOMENON_TABLES` comment documents.

A supernova remnant has no galaxy-frame placement at all (see
`queryDb._SUPERNOVA_REMNANT_TABLE`'s own docstring) and so never appears
on the Galaxy Map or as a NAV endpoint -- but it does have its own real
`radius_ly`, which is all this diagram needs, so it renders here exactly
like a nebula/asteroid field.
"""

from fmt import esc

try:
    from stellarObjects.utils import ly_to_au
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # duplicated fallback, matching galaxy.py's/phenomenon.py's own
    # identical pattern for pc_to_ly.
    ly_to_au = None

_LY_TO_AU_FALLBACK = 63241.1


def _to_au(radius_ly):
    return ly_to_au(radius_ly) if ly_to_au else radius_ly * _LY_TO_AU_FALLBACK


_MIN_VIEW_SIZE_AU = 1.0
"""float: The deepest zoom-in this diagram allows, as a full viewBox-edge
span in AU -- the "max zoom in as 1 AU" bound."""

_MAX_VIEW_SIZE_AU = _to_au(1.0)
"""float: The furthest zoom-out this diagram allows, as a full viewBox-edge
span in AU (one light-year, converted) -- the "max zoom out of 1 ly"
bound. A phenomenon whose own real radius is larger than this (a big
nebula can be dozens of ly across) is simply not shown in full at the
default view -- this diagram exists to give an AU-scale sense of the
object, not to replace the Galaxy/Sector Map's own whole-object view."""

_POINT_OBJECT_RADIUS_AU = 0.15
"""float: Fixed, illustrative dot radius (AU) for a point-like phenomenon
(`radius_ly` 0 or `None` -- a black hole/neutron star). Not a real
physical size (see this module's own docstring); chosen purely to read as
a small, clearly-marked point at this diagram's own default zoom."""

_POINT_OBJECT_DEFAULT_VIEW_AU = 12.0
"""float: The starting view span (AU) for a point-like phenomenon, since
(unlike a nebula/asteroid field/supernova remnant) there's no real
physical extent to size a sensible default view from -- a modest
"local neighborhood" view, well within `_MIN_VIEW_SIZE_AU`/
`_MAX_VIEW_SIZE_AU`, that the visitor can freely zoom in/out from."""

_DEFAULT_VIEW_MARGIN = 2.4
"""float: How much wider than the phenomenon's own diameter the default
starting view is, for a real-radius (non-point-like) phenomenon -- enough
margin that the whole circle reads clearly inside the viewport with room
to spare, without being so wide it looks like a speck."""

_TYPE_COLORS = {
    "nebula": "#c9a8e0", "asteroid_field": "#b89a6e",
    "black_hole": "#1a1a1a", "neutron_star": "#cfe8ff",
    "supernova_remnant": "#e08a5c",
    "rogue_planet": "#7a8ba0", "interstellar_comet": "#a8d0e0",
}
_DEFAULT_COLOR = "#9aa0ac"


def render_phenomenon_map_panel(phenomenon_type, name, radius_ly):
    """
    Builds the "Diagram" panel: a flat, zoomable SVG showing one
    phenomenon's own real extent at AU scale (or a fixed illustrative dot
    for a point-like compact remnant) -- see this module's docstring.

    Args:
        phenomenon_type (str): One of `queryDb._PHENOMENON_TYPE_TO_TABLE`'s
            keys (`"nebula"`, `"asteroid_field"`, `"black_hole"`,
            `"neutron_star"`, `"supernova_remnant"`, `"rogue_planet"`, or
            `"interstellar_comet"`).
        name (str): The phenomenon's own name, for the SVG's `aria-label`.
        radius_ly (float): The phenomenon's own real radius, in
            light-years -- `0`/`None` for a point-like compact remnant
            (`html/phenomenon.py`'s caller passes `detail.get("radius_ly")
            or 0`, since `black_holes`/`neutron_stars` have no such column
            of their own at all).

    Returns:
        str: A complete `<section class="panel">` block, including its own
            `static/mapzoom.js`/`static/phenomenonmap.js` `<script>` tags
            (same self-contained convention `galaxymap.py`'s own panel
            uses).
    """
    color = _TYPE_COLORS.get(phenomenon_type, _DEFAULT_COLOR)
    radius_au = _to_au(radius_ly) if radius_ly else 0.0
    is_point_like = radius_au <= 0

    if is_point_like:
        default_view_au = _POINT_OBJECT_DEFAULT_VIEW_AU
        stroke_width = _POINT_OBJECT_RADIUS_AU * 0.35
        marker = (
            f'<circle cx="0" cy="0" r="{_POINT_OBJECT_RADIUS_AU:.4g}" '
            f'fill="{color}" stroke="#ffffff" stroke-width="{stroke_width:.4g}"/>'
        )
        hint_extra = " -- shown as a fixed illustrative point; its own real size is negligible at this scale"
    else:
        default_view_au = min(_MAX_VIEW_SIZE_AU, max(_MIN_VIEW_SIZE_AU, radius_au * _DEFAULT_VIEW_MARGIN))
        stroke_width = default_view_au * 0.004
        marker = (
            f'<circle cx="0" cy="0" r="{radius_au:.6g}" '
            f'fill="{color}" fill-opacity="0.45" stroke="{color}" stroke-width="{stroke_width:.4g}"/>'
        )
        hint_extra = ""

    half = default_view_au / 2
    view_box = f"{-half:.6g} {-half:.6g} {default_view_au:.6g} {default_view_au:.6g}"

    svg = (
        f'<svg class="phenomenonmap-svg" id="phenomenonmap-svg" viewBox="{view_box}" role="img" '
        f'data-min-view-size="{_MIN_VIEW_SIZE_AU:.6g}" data-max-view-size="{_MAX_VIEW_SIZE_AU:.6g}" '
        f'aria-label="{esc(name)} diagram, drawn to astronomical-unit scale. Scroll or use the '
        f'+/- buttons to zoom, drag to pan.">{marker}</svg>'
    )

    return f"""
<section class="panel">
<div class="panel-header">
  <h2>Diagram</h2>
  <span class="hint">Drawn to real astronomical-unit (AU) scale{hint_extra}</span>
</div>
<div class="phenomenonmap-layout">
<div class="phenomenonmap-viewport">
{svg}
<div class="starmap-scale" id="phenomenonmap-scale"></div>
</div>
<div class="phenomenonmap-side">
<div class="starmap-controls" id="phenomenonmap-controls">
  <button type="button" class="starmap-btn" data-action="zoom-out" aria-label="Zoom out">&minus;</button>
  <button type="button" class="starmap-btn" data-action="zoom-in" aria-label="Zoom in">+</button>
  <button type="button" class="starmap-btn" data-action="reset">Reset view</button>
</div>
<p class="hint">Scroll/drag/+/- to zoom, from about 1 AU up to 1 ly across.</p>
</div>
</div>
</section>
<script src="static/mapzoom.js"></script>
<script src="static/phenomenonmap.js"></script>
"""
