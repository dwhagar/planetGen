# html/lib/phenomenonrender.py

"""
The "View" panel on a phenomenon page: a rendered picture of the object
itself, chosen per type (UX.17), instead of the flat AU-scale diagram
`phenomenonmap.py` draws.

- Neutron star: the star spinning with its two radio beams (none for a
  non-pulsing one), in rough time with its real spin period: slowed to a
  visible rate on a log scale, so a faster pulsar still turns faster, and
  the caption states both periods.
- Black hole: a 3D accretion disk around the black shadow (a bare shadow
  when it has no disk). A stellar-mass hole's disk is hotter and bluer
  than an intermediate one's, as a smaller hole's inner disk is.
- Quasar: a bigger, brighter disk inside a dusty torus, with jets when it
  is radio-loud.
- Rogue planet: a dark, starless body, glowing faintly when it keeps
  internal heat.
- Interstellar comet: a nucleus, with a coma and tail when it is active.
- Asteroid field: no view at all (`view_kind` returns `None`).
- Nebula and supernova remnant: the AU diagram until a render is designed
  for them.

The panel holds a static SVG still that reads without JavaScript; the
script (`static/phenomenonrender.js`) swaps in a three.js canvas and
draws one still frame under `prefers-reduced-motion`.
"""

import json
import math

from fmt import esc, format_duration_seconds, format_number

_RENDERED = {"neutron_star", "black_hole", "quasar", "rogue_planet", "interstellar_comet"}
_NO_VIEW = {"asteroid_field"}

_SHOWN_PERIOD_MIN_S = 0.6
"""float: How long the fastest pulsar (1 ms) takes per turn on screen."""

_SHOWN_PERIOD_LOG_DIVISOR = 4.0
"""float: Each factor of 10 in the real period stretches the shown one by
10**(1/4) (about 1.8x), so 1 ms shows as 0.6 s and 10 s as 6 s."""


def view_kind(phenomenon_type):
    """`"render"`, `"map"` (the AU diagram) or `None` (no view at all)."""
    if phenomenon_type in _NO_VIEW:
        return None
    return "render" if phenomenon_type in _RENDERED else "map"


def shown_spin_period_s(spin_period_ms):
    """Seconds per turn on screen for a real spin period in ms: a log
    mapping that keeps faster pulsars faster but always visible."""
    period_ms = max(1.0, float(spin_period_ms or 1000.0))
    return _SHOWN_PERIOD_MIN_S * 10 ** (math.log10(period_ms) / _SHOWN_PERIOD_LOG_DIVISOR)


def _params(phenomenon_type, detail):
    """The numbers the script needs, plus the caption under the view."""
    if phenomenon_type == "neutron_star":
        pulsing = detail.get("pulsar_type") != "non-pulsing"
        real_s = (detail.get("spin_period_ms") or 1000.0) / 1000
        shown_s = shown_spin_period_s(detail.get("spin_period_ms"))
        caption = (f"One turn every {format_duration_seconds(real_s)} in reality, slowed here to "
                   f"one every {shown_s:.1f} s (about {format_number(shown_s / real_s)} times slower).")
        if pulsing:
            caption += " Each flash is a radio beam sweeping past you."
        else:
            caption += " It no longer pulses, so it shows no beams."
        return {"pulsing": pulsing, "period_s": round(shown_s, 3),
                "magnetar": (detail.get("magnetic_field_gauss") or 0) >= 1e14}, caption
    if phenomenon_type == "black_hole":
        disk = bool(detail.get("has_accretion_disk"))
        stellar = (detail.get("mass_class") or "stellar") == "stellar"
        caption = ("The black shadow is the event horizon and the light it bends; the disk "
                   "is gas spiralling in, hottest at the inner edge." if disk else
                   "With no accretion disk, the black hole shows only as a shadow that "
                   "bends the starlight behind it.")
        return {"disk": disk, "stellar": stellar,
                "spin": detail.get("spin") or 0}, caption
    if phenomenon_type == "quasar":
        jets = bool(detail.get("is_radio_loud"))
        caption = ("A supermassive black hole's blazing accretion disk inside a ring of dust"
                   + (", with jets of plasma along its spin axis." if jets else
                      "; it is radio-quiet, so it shows no jets."))
        return {"jets": jets}, caption
    if phenomenon_type == "rogue_planet":
        gas = detail.get("planet_type") == "g"
        heat = bool(detail.get("has_internal_heat"))
        caption = ("A planet with no star, lit only by distant starlight"
                   + (", glowing faintly with its own internal heat." if heat else "."))
        return {"gas": gas, "heat": heat}, caption
    active = bool(detail.get("is_active"))
    caption = ("The nucleus is shedding gas and dust into a coma and tail." if active else
               "Far from any star it is frozen and inactive: a bare nucleus with no tail.")
    return {"active": active}, caption


_STILL_COLORS = {
    "neutron_star": "#cfe8ff", "black_hole": "#ffb070", "quasar": "#fff3d6",
    "rogue_planet": "#5a6878", "interstellar_comet": "#a8d0e0",
}


def _still_svg(phenomenon_type, params, name):
    """A simple picture for when the script or WebGL can't run."""
    color = _STILL_COLORS[phenomenon_type]
    parts = []
    if phenomenon_type == "neutron_star" and params["pulsing"]:
        parts.append('<path d="M0 0 L-18 -95 L18 -95 Z M0 0 L-18 95 L18 95 Z" '
                     'fill="#9fd0ff" fill-opacity="0.35" transform="rotate(20)"/>')
    if phenomenon_type in ("black_hole", "quasar") and (params.get("disk", True)):
        parts.append(f'<ellipse cx="0" cy="0" rx="85" ry="22" fill="none" stroke="{color}" '
                     'stroke-width="14" stroke-opacity="0.8"/>')
    if phenomenon_type == "quasar" and params["jets"]:
        parts.append('<path d="M-4 0 L-10 -98 L10 -98 L4 0 Z M-4 0 L-10 98 L10 98 L4 0 Z" '
                     'fill="#a8c8ff" fill-opacity="0.5"/>')
    if phenomenon_type == "interstellar_comet" and params["active"]:
        parts.append('<path d="M0 -8 L95 -30 L95 30 L0 8 Z" fill="#a8d0e0" fill-opacity="0.3"/>')
    core = {"black_hole": ("22", "#000"), "quasar": ("16", "#fffbe8"),
            "neutron_star": ("12", color), "rogue_planet": ("40", color),
            "interstellar_comet": ("7", "#8a8f96")}[phenomenon_type]
    parts.append(f'<circle cx="0" cy="0" r="{core[0]}" fill="{core[1]}" stroke="{color}" '
                 'stroke-width="2"/>')
    return (f'<svg class="phenomrender-still" viewBox="-100 -100 200 200" role="img" '
            f'aria-label="{esc(name)}, a still drawing">'
            '<rect x="-100" y="-100" width="200" height="200" fill="#05070c"/>'
            + "".join(parts) + "</svg>")


def render_phenomenon_view_panel(phenomenon_type, detail):
    """
    The phenomenon page's "View" panel for a rendered type (see
    `view_kind`).

    Args:
        phenomenon_type (str): A `queryDb._PHENOMENON_TYPE_TO_TABLE` key.
        detail (dict): The phenomenon's API detail (name and the columns
            `_params` reads).

    Returns:
        str: A complete `<section class="panel">`; the page loads
            `static/phenomenonrender.js` itself.
    """
    params, caption = _params(phenomenon_type, detail)
    data = json.dumps({"kind": phenomenon_type, **params})
    name = detail.get("name") or "Phenomenon"
    return f"""
<section class="panel" aria-labelledby="phenomrender-heading">
<div class="panel-header">
  <h2 id="phenomrender-heading">View</h2>
  <span class="hint">An artist's impression from its properties; drag to turn it.</span>
</div>
<div class="phenomrender-viewport" id="phenomrender" data-view="{esc(data)}">
{_still_svg(phenomenon_type, params, name)}
<canvas class="phenomrender-canvas" id="phenomrender-canvas" hidden
  role="img" aria-label="{esc(name)}, rendered in 3D. Drag to turn it."></canvas>
</div>
<p class="hint phenomrender-caption">{esc(caption)}</p>
</section>
"""
