# planetgen/web/maps/ringlook.py

"""
The look of a planet's ring system (MAP.275): what `static/ringpaint.js` needs to draw the rings a planet stores.

`ring_look(row)` turns the stored ring (`ring_kind`, `ring_inner_km`, `ring_outer_km`, `ring_optical_depth`, `ring_ice_fraction`, GEN.732) and the planet's
moons into one record of edges in planet radii, a colour, a peak opacity and the gaps the moons dig at their resonances. Colour and brightness follow the
ice share (GEN.670 research, `reports/gen670-body-rendering-research-2026-11-04.md` section 5): a reddish slope below 0.55 micrometers, steepest for an
icy ring; brightness set by the ice share, a dusty ring dark. Nothing here is stored or feeds generation; it is a picture of what is stored.

Every constant is a stand-in registered in `docs/design/research/gap-filling-plan.md` section 6.170 (rows RNG1 to RNG8).
"""

import math

ICY_FRACTION = 0.985
"""float: An ice share from this up (Saturn's A and B rings hold 99.5 percent or more by volume) reddens at the full measured slope [R5090] [R5091]."""
BLUE_SLOPE_ICY = 3.5
"""float: The reflectance slope between 0.35 and 0.55 micrometers (per micrometer, normalized to 1 at 0.55) of an icy ring: 3 to 4 in the B ring [R5091]."""
BLUE_SLOPE_DUST = 1.5
"""float: The same slope for a ring of dust and no ice (a dark red-brown ring); between the two ice shares it is interpolated. A stand-in."""
RED_SLOPE_ICY = 0.05
"""float: The slope between 0.55 and 0.95 micrometers of a thick icy ring: 0.0 to 0.1 in the A and B rings [R5091]."""
RED_SLOPE_THIN = 0.25
"""float: The same for a ring of low optical depth: about 0.3 in the C ring, 0.2 in the Cassini Division [R5091]."""
THIN_TAU, THICK_TAU = 0.1, 0.5
"""float: The optical depths below and above which the red slope is the thin and the thick ring's: the C ring's 0.05 to 0.10 and the A ring's 0.5 [R5090]."""
WAVELENGTH_UM = {"r": 0.65, "g": 0.55, "b": 0.45}
"""dict: The wavelength, micrometers, taken for each display channel."""

ALBEDO_ICE = 0.55
"""float: The Bond albedo of a ring with an ice share of 0.95 or more: the A and B rings' thermal-model value [R5090]. Bolometric, used for the visible look."""
ALBEDO_HALF = 0.1
"""float: The albedo at an ice share of 0.5 (the low end of the C ring's 0.4 to 0.1 [R5090])."""
ALBEDO_DUST = 0.04
"""float: The albedo of a ring with no ice, dark gray-brown: an educated guess, no source gives it."""
DISPLAY_GAIN = 1.6
"""float: How much brighter the picture is than the albedo (so the albedo 0.55 of an icy ring reads as a bright ring; the visible I/F of Saturn's rings is about 0.5 [R5090])."""

DUST_ALPHA_FLOOR = 0.06
"""float: The least opacity a stored ring is drawn with, so a ring of very low optical depth (Jupiter's dust ring) is still seen as a faint haze."""
MAX_ALPHA = 0.95
"""float: The most opaque a ring is drawn."""

RESONANCES = ((2, 1, 1.0), (3, 2, 0.6), (3, 1, 0.35), (4, 3, 0.35))
"""tuple: The resonances with a moon that clear a gap, `(p, q, weight)`: a ring particle goes p times round for each q turns of the moon; the radius is
`a_moon * (q / p) ** (2/3)`. The 2:1 resonance of Mimas is the Cassini Division; the weight is the strength against it (first order strongest). Stand-in."""
GAP_COEFFICIENT = 2500.0
"""float: A gap's width over its radius is `weight * GAP_COEFFICIENT * (m_moon / M_planet) ** (2/3)`: fitted to the Cassini Division (0.04 of the radius for Mimas, 6.6e-8). Stand-in."""
GAP_MAX_SHARE = 0.25
"""float: No gap is wider than this share of the ring's width (a big moon clears more of the ring than a gap shows)."""
GAP_MIN_SHARE = 0.006
"""float: A gap narrower than this share of the ring's width is not drawn."""


def _clamp(x, low, high):
    return low if x < low else high if x > high else x


def albedo(ice_fraction):
    """The Bond albedo of a ring that is `ice_fraction` ice: 0.55 from 0.95, falling to 0.1 at 0.5 and to 0.04 at none."""
    ice = _clamp(ice_fraction, 0.0, 1.0)
    if ice >= 0.95:
        return ALBEDO_ICE
    if ice >= 0.5:
        return ALBEDO_HALF + (ALBEDO_ICE - ALBEDO_HALF) * (ice - 0.5) / 0.45
    return ALBEDO_DUST + (ALBEDO_HALF - ALBEDO_DUST) * ice / 0.5


def slopes(ice_fraction, optical_depth):
    """`(blue, red)` reflectance slopes (per micrometer from 0.55): the blue one is the icy ring's flat 3.5 from an ice share of 0.985 and falls to the dust
    ring's below; the red one is the thin ring's 0.25 and falls to the icy ring's 0.05 as the ring thickens and the ice share rises."""
    ice = _clamp((ice_fraction - 0.5) / (ICY_FRACTION - 0.5), 0.0, 1.0)
    blue = BLUE_SLOPE_DUST + (BLUE_SLOPE_ICY - BLUE_SLOPE_DUST) * ice
    thick = _clamp((optical_depth - THIN_TAU) / (THICK_TAU - THIN_TAU), 0.0, 1.0)
    red = RED_SLOPE_THIN + (RED_SLOPE_ICY - RED_SLOPE_THIN) * ice * thick
    return blue, red


def rgb(ice_fraction, optical_depth):
    """The ring's colour, three 0 to 255 integers: the albedo times the reflectance `1 + slope * (wavelength - 0.55)` of each channel."""
    blue, red = slopes(ice_fraction, optical_depth)
    base = albedo(ice_fraction) * DISPLAY_GAIN
    out = []
    for channel in ("r", "g", "b"):
        step = WAVELENGTH_UM[channel] - WAVELENGTH_UM["g"]
        slope = red if step > 0 else blue
        out.append(int(round(255 * _clamp(base * (1.0 + slope * step), 0.0, 1.0))))
    return out


def alpha(optical_depth):
    """The peak opacity of the picture: `1 - exp(-tau)`, no less than `DUST_ALPHA_FLOOR` and no more than `MAX_ALPHA`."""
    tau = optical_depth if optical_depth and optical_depth > 0 else 0.0
    return round(_clamp(1.0 - math.exp(-tau), DUST_ALPHA_FLOOR, MAX_ALPHA), 3)


def resonance_gaps(inner_km, outer_km, planet_mass_kg, moons):
    """
    The gaps the moons dig in a ring between `inner_km` and `outer_km`, each `[center, width, depth]` with the center and width as shares of the ring's width
    (0 at the inner edge) and the depth 0 to 1 (1 empty). A moon's resonance falls in the ring at `a * (q / p) ** (2/3)` for each of `RESONANCES`.

    Args:
        inner_km, outer_km (float): The ring's edges, km from the planet's center.
        planet_mass_kg (float): The planet's mass.
        moons (iterable): Dicts with `distance_km` and `mass_kg`.
    """
    width_km = outer_km - inner_km
    if not (width_km > 0 and planet_mass_kg and planet_mass_kg > 0):
        return []
    gaps = []
    for moon in moons or ():
        distance, mass = moon.get("distance_km"), moon.get("mass_kg")
        if not (distance and distance > 0 and mass and mass > 0):
            continue
        for p, q, weight in RESONANCES:
            radius = distance * (q / p) ** (2.0 / 3.0)
            if not inner_km < radius < outer_km:
                continue
            share = weight * GAP_COEFFICIENT * (mass / planet_mass_kg) ** (2.0 / 3.0) * radius / width_km
            share = min(share, GAP_MAX_SHARE)
            if share < GAP_MIN_SHARE:
                continue
            gaps.append([round((radius - inner_km) / width_km, 4), round(share, 4), round(_clamp(share / 0.04, 0.35, 0.95), 3)])
    gaps.sort()
    return gaps


def ring_look(row):
    """
    The `ring` entry of a body's look record, or `None` for a body that stores no ring: `{kind, origin, in, out, tau, ice, alpha, rgb, gaps}` with the edges
    `in` and `out` in the body's radii.

    Args:
        row (dict): A planet or moon row with the `ring_*` columns, `radius_km`, `mass_kg` and (for the gaps) `moons`.
    """
    inner, outer, radius = row.get("ring_inner_km"), row.get("ring_outer_km"), row.get("radius_km")
    kind = row.get("ring_kind")
    if not kind or not (inner and outer and radius and radius > 0 and outer > inner):
        return None
    ice = row.get("ring_ice_fraction")
    ice = 0.0 if ice is None else float(ice)
    tau = row.get("ring_optical_depth")
    tau = 0.0 if tau is None else float(tau)
    return {
        "kind": kind, "origin": row.get("ring_origin"),
        "in": round(inner / radius, 4), "out": round(outer / radius, 4),
        "tau": float(f"{tau:.3g}"), "ice": round(ice, 3), "alpha": alpha(tau), "rgb": rgb(ice, tau),
        "gaps": resonance_gaps(inner, outer, row.get("mass_kg"), row.get("moons")),
    }
