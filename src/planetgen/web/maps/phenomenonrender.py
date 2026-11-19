# planetgen/web/maps/phenomenonrender.py

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
- Rogue planet: a dark, starless body painted by the same painters as a
  planet (`bodylook.rogue_planet_look`, MAP.278), glowing faintly when it
  keeps internal heat.
- Interstellar comet: a nucleus (grey-white when active, dark red-brown
  when not), with a coma and tail when it is active.
- Asteroid field (GEN.361): its own bodies as points in a cloud, each at the place `field_dispersal.member_position_m` gives it at the field's age (a
  ball for a collision, a thin ring for a belt a supernova unbound), the heaviest ranks drawn larger. The caption states how many bodies the field
  really holds above 1 km, that each point is drawn far larger than any body, and the field's column density, how soon after it began to spread its
  pieces stopped colliding and how far galactic shear has drawn it out.
- Nebula (MAP.105, MAP.142): the nebula drawn from its own points in 3D, with the brightest stars of the
  galaxy round it dimmed for reference (`render_nebula_view_panel`,
  `static/nebulaview.js`).
- Supernova remnant: the AU diagram until a render is designed for it.

The panel holds a static SVG still that reads without JavaScript; the
script (`static/phenomenonrender.js`) swaps in a three.js canvas and
draws one still frame under `prefers-reduced-motion`.
"""

import json
import math

from planetgen.physics import belt_ledger, constants, field_dispersal
from planetgen.web.lib.fmt import esc, format_duration_seconds, format_number
from planetgen.web.maps import bodylook
from planetgen.web.maps.starmap import NEBULA_SURROUNDINGS_PATH, NEBULA_VOLUME_PATH, nebula_look

_RENDERED = {"neutron_star", "black_hole", "quasar", "rogue_planet", "interstellar_comet", "asteroid_field"}
_NO_VIEW = set()

_SHOWN_PERIOD_MIN_S = 0.6
"""float: How long the fastest pulsar (1 ms) takes per turn on screen."""

_SHOWN_PERIOD_LOG_DIVISOR = 4.0
"""float: Each factor of 10 in the real period stretches the shown one by
10**(1/4) (about 1.8x), so 1 ms shows as 0.6 s and 10 s as 6 s."""


def view_kind(phenomenon_type):
    """`"render"`, `"nebula"` (the 3D shape view), `"map"` (the AU diagram)
    or `None` (no view at all)."""
    if phenomenon_type in _NO_VIEW:
        return None
    if phenomenon_type == "nebula":
        return "nebula"
    return "render" if phenomenon_type in _RENDERED else "map"


FIELD_VIEW_POINTS = 480
"""int: How many bodies of a field the view draws (the first ranks of the field's seed): enough to show the shape of the cloud, few enough that the record stays about 11 KB."""
FIELD_VIEW_HEROES = 12
"""int: How many of those, the heaviest ranks, the view draws larger."""
FIELD_REFERENCE_DIAMETER_M = 1000.0
"""float: The size the caption counts bodies above, and the pieces whose freeze-out it states: 1 km."""
FIELD_RING_RADIUS_M = constants.AU_M
"""float: The ring a belt that a supernova unbound left from, for the view only: the cloud is a ring of light-years, so the ring it started on (a belt of tens of AU, not stored) is a point
beside it; a stand-in for the drawing, not a stored number."""


def _field_bodies_text(count):
    """`count` as a rounded phrase, "about 2 billion", "about 3.4 million"; the power of ten from a trillion."""
    for scale, word in ((1e12, "trillion"), (1e9, "billion"), (1e6, "million"), (1e3, "thousand")):
        if count >= scale:
            return f"about {format_number(count / scale, ',.2g')} {word}"
    return f"about {format_number(count, ',.2g')}"


def field_view(detail):
    """
    The numbers of an asteroid field's view and the line under it (GEN.361): `(parameters, caption)`.

    The bodies are the field's own: body `rank` (1 to `FIELD_VIEW_POINTS`) of the seed made from the field's Object ID is at `field_dispersal.member_position_m` at the field's age and rms
    speed, the same cloud after a database rebuild. The positions are scaled so that their rms distance is 1 (the page draws the cloud at its own scale and states the size in the caption). A field
    stored without its age or speed (an old store) is drawn as a ball of the radius it has, the speed the age implies.
    """
    uid = str(detail.get("uid") or detail.get("id"))
    radius_m = max(float(detail.get("radius_ly") or 0.0), 1e-9) * constants.LIGHTYEAR_M
    age_yr = float(detail.get("birth_age_yr") or 0.0)
    speed_ms = float(detail.get("rms_speed_kms") or 0.0) * 1000.0
    if not (age_yr > 0.0 and speed_ms > 0.0):
        age_yr = 1.0e6
        speed_ms = radius_m / (age_yr * constants.SECONDS_PER_YEAR)
    shape = detail.get("shape") if detail.get("shape") in field_dispersal.SHAPES else "ball"
    positions = [field_dispersal.member_position_m(uid, rank, age_yr, speed_ms, shape, FIELD_RING_RADIUS_M) for rank in range(1, FIELD_VIEW_POINTS + 1)]
    scale = math.sqrt(sum(x * x + y * y + z * z for x, y, z in positions) / len(positions)) or 1.0
    points = [round(value / scale, 3) for position in positions for value in position]
    mass_kg = float(detail.get("mass_kg") or 0.0)
    law = belt_ledger.ICY_LAW if detail.get("composition_family") == "icy" else belt_ledger.ROCKY_LAW
    caption = (f"{format_number(len(positions))} of its bodies, each at the place its field's seed gives it after "
               f"{format_duration_seconds(age_yr * constants.SECONDS_PER_YEAR)} of drifting apart; the heaviest are drawn larger. ")
    if mass_kg > 0.0:
        reference_kg = math.pi / 6.0 * FIELD_REFERENCE_DIAMETER_M ** 3 * law.density_kg_m3
        caption += (f"It really holds {_field_bodies_text(belt_ledger.count_above(law, mass_kg, reference_kg))} bodies larger than 1 km, each point here "
                    "drawn far larger than any of them: seen close, the field is almost empty space. "
                    f"Seen end on, its mass comes to {format_number(field_dispersal.column_density_kg_m2(mass_kg, radius_m), '.2g')} kg per square metre; "
                    f"pieces of 1 km stopped colliding {format_duration_seconds(field_dispersal.freeze_out_time_s(mass_kg, FIELD_REFERENCE_DIAMETER_M, speed_ms))} after it began to spread. ")
    if age_yr <= field_dispersal.SHEAR_VALID_YEARS:
        caption += (f"Galactic shear has drawn its far side "
                    f"{format_number(field_dispersal.shear_displacement_m(radius_m, age_yr) / constants.LIGHTYEAR_M, ',.2g')} ly along its orbit.")
    else:
        caption += "It is older than the galaxy's shear is followed for (26 Myr), so the shear is left out of this account."
    return {"points": points, "heroes": FIELD_VIEW_HEROES, "ring": shape == "annulus", "seed": bodylook.seed_for("asteroid_field", uid)}, caption


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
                   f"one every {format_duration_seconds(shown_s)} (about {format_number(shown_s / real_s)} times slower).")
        if pulsing:
            caption += " Each flash is a radio beam sweeping past you."
        else:
            caption += " It no longer pulses, so it shows no beams."
        return {"pulsing": pulsing, "period_s": round(shown_s, 3),
                "magnetar": (detail.get("magnetic_field_gauss") or 0) >= 1e14,
                "surface": bodylook.blackbody_hex(detail.get("surface_temperature_k")),
                "seed": bodylook.seed_for("neutron_star", detail.get("id"))}, caption
    if phenomenon_type == "black_hole":
        disk = bool(detail.get("has_accretion_disk"))
        stellar = (detail.get("mass_class") or "stellar") == "stellar"
        caption = ("The black shadow is the event horizon and the light it bends; the disk "
                   "is gas spiralling in, hottest at the inner edge." if disk else
                   "With no accretion disk, the black hole shows only as a shadow that "
                   "bends the starlight behind it.")
        return {"disk": disk, "stellar": stellar, "spin": detail.get("spin") or 0,
                "diskColor": bodylook.blackbody_hex(detail.get("temperature_k")) if disk else None,
                "seed": bodylook.seed_for("black_hole", detail.get("id"))}, caption
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
        if detail.get("effective_temperature_k"):
            caption += f" It radiates at about {format_number(detail['effective_temperature_k'], ',.0f')} K."
        look = bodylook.rogue_planet_look(detail)
        return {"gas": gas, "heat": heat, "look": look, "seed": look["seed"]}, caption
    if phenomenon_type == "asteroid_field":
        return field_view(detail)
    active = bool(detail.get("is_active"))
    caption = ("The nucleus is shedding gas and dust into a coma and tail." if active else
               "Far from any star it is frozen and inactive: a bare nucleus with no tail.")
    return {"active": active, "nucleus": bodylook.comet_nucleus_color(detail),
            "seed": bodylook.seed_for("interstellar_comet", detail.get("id"))}, caption


def view_params(phenomenon_type, detail):
    """`(parameters, caption)` of a rendered type (`view_kind` is `"render"`), the record `static/phenomenonrender.js` paints and the line under
    it; the Sensor registry (`web/lib/sensor.py`, UX.275) reads the same pair. `KeyError`-free for the five rendered types, `ValueError` for any other."""
    if phenomenon_type not in _RENDERED:
        raise ValueError(f"{phenomenon_type!r} has no rendered view")
    return _params(phenomenon_type, detail)


_STILL_COLORS = {
    "neutron_star": "#cfe8ff", "black_hole": "#ffb070", "quasar": "#fff3d6",
    "rogue_planet": "#5a6878", "interstellar_comet": "#a8d0e0", "asteroid_field": "#b9a88c",
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
    if phenomenon_type == "asteroid_field":
        for number in range(60):      # a still scatter from the field's own points: x and y of the first 60 bodies, drawn at 70 units across
            x, y = params["points"][3 * number], params["points"][3 * number + 1]
            parts.append(f'<circle cx="{max(-95.0, min(95.0, x * 40)):.1f}" cy="{max(-95.0, min(95.0, y * 40)):.1f}" r="{2.6 if number < 8 else 1.2}" fill="{color}" fill-opacity="0.8"/>')
        return (f'<svg class="phenomrender-still" viewBox="-100 -100 200 200" role="img" '
                f'aria-label="{esc(name)}, a still drawing"><rect x="-100" y="-100" width="200" height="200" fill="#000000"/>' + "".join(parts) + "</svg>")
    core = {"black_hole": ("22", "#000"), "quasar": ("16", "#fffbe8"),
            "neutron_star": ("12", color), "rogue_planet": ("40", color),
            "interstellar_comet": ("7", "#8a8f96")}[phenomenon_type]
    parts.append(f'<circle cx="0" cy="0" r="{core[0]}" fill="{core[1]}" stroke="{color}" '
                 'stroke-width="2"/>')
    return (f'<svg class="phenomrender-still" viewBox="-100 -100 200 200" role="img" '
            f'aria-label="{esc(name)}, a still drawing">'
            '<rect x="-100" y="-100" width="200" height="200" fill="#000000"/>'
            + "".join(parts) + "</svg>")


def render_phenomenon_view_panel(phenomenon_type, detail):
    """
    The phenomenon page's "View" panel for a rendered type (see
    `view_kind`).

    Args:
        phenomenon_type (str): A `planetgen.db.query._PHENOMENON_TYPE_TO_TABLE` key.
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
  <h2 id="phenomrender-heading">Sensor view</h2>
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


def render_nebula_view_panel(detail):
    """
    A nebula page's "View" panel (MAP.105, MAP.142): the nebula drawn from its own
    points in 3D with the brightest stars of the galaxy round it, dimmed. The script
    (`static/nebulaview.js`) draws it from the volume and surroundings
    endpoints; the still drawing reads without it.

    Args:
        detail (dict): The nebula's API detail (`id`, `name`,
            `nebula_type`, `radius_ly`).

    Returns:
        str: A complete `<section class="panel">`; the page loads
            `static/nebulaview.js` itself.
    """
    color, core, edge = nebula_look(detail.get("nebula_type"))
    name = detail.get("name") or "Nebula"
    data = json.dumps({
        "nebulaId": str(detail["id"]),
        "volumePath": NEBULA_VOLUME_PATH.replace("{id}", str(detail["id"])),
        "surroundingsPath": NEBULA_SURROUNDINGS_PATH.replace("{id}", str(detail["id"])),
        "color": color, "coreOpacity": core, "edgeOpacity": edge,
    })
    still = (f'<svg class="phenomrender-still" viewBox="-100 -100 200 200" role="img" '
             f'aria-label="{esc(name)}, a still drawing"><rect x="-100" y="-100" width="200" height="200" '
             f'fill="#000000"/><ellipse cx="0" cy="0" rx="62" ry="40" fill="{color}" fill-opacity="{edge:.2f}" '
             f'stroke="{color}" stroke-opacity="0.6" transform="rotate(-18)"/></svg>')
    return f"""
<section class="panel" aria-labelledby="nebulaview-heading">
<div class="panel-header">
  <h2 id="nebulaview-heading">Sensor view</h2>
  <span class="hint">Its cloud, with the galaxy's brightest stars round it; drag to turn, scroll to zoom.</span>
</div>
<div class="phenomrender-viewport" id="nebulaview" data-view="{esc(data)}">
{still}
<canvas class="phenomrender-canvas" id="nebulaview-canvas" hidden
  role="img" aria-label="{esc(name)}, its cloud in 3D with the brightest stars round it. Drag to turn it."></canvas>
</div>
<p><button type="button" class="starmap-btn" id="nebulaview-recenter" hidden
  title="Back to the opening angle and zoom, about the nebula's center (C)">Re-center</button></p>
<p class="hint phenomrender-caption">An irregular cloud drawn from its properties. The dimmed points are the
brightest stars of the galaxy near it, there for reference.</p>
</section>
"""
