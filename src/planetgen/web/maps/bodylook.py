# planetgen/web/maps/bodylook.py

"""
The look record (MAP.271): one small JSON object per rendered body that holds every number the page needs to
paint it, built here from the stored row so the page never parses a label.

`systemmap.py` puts it on each marker as `data-look`; `static/bodytexture.js` and `static/solidpaint.js` turn it
plus `seed` into a texture. A field with no stored value is left out (not zero), so the page can tell "unknown"
from "none"; the painters fall back to the plain tinted disc for anything missing.

Fields (all optional but `v`, `kind`, `body`, `seed`):

    v           record version (1)
    kind        "star", "planet" or "moon"
    body        "star", "solid" or "gas"
    seed        32 bits from the body's id and kind (`seed_for`); display noise only, never generation
    cls         planet class letter
    color       the base colour as `#rrggbb`: a star's blackbody colour, a solid world's land colour, a giant's cloud top
    quality     "estimate" while `color` comes from temperature and body family alone (MAP.274 replaces it with
                the composition colour), absent for a star
    glow        atmosphere limb colour (`#rrggbb`), present only with an atmosphere
    atm         true when the body has an atmosphere
    pkpa        surface pressure, kPa
    tempK       mean surface temperature, K (a star's effective temperature)
    radiusKm    radius
    hydro       `hydrosphere` (dry, vapour, ice, ice-covered ocean, surface ocean, hycean)
    ocean       fraction of the surface under ocean, 0 to 1
    land        fraction of the surface above sea level, 0 to 1
    oceanClass  ice-sealed, chloride brine, acid sulfate, soda or neutral
    oceanColor  the sea colour for that class, `#rrggbb`
    oceanKm     ocean depth, km
    iceKm       ice shell thickness, km
    lumSolar    a star's luminosity in solar units
    spec        a star's `star_type`
    yerkes      a star's Yerkes luminosity class (V dwarf, III giant, VII white dwarf, ...)

The colour tables are stand-ins for hue only (the gap plan registry lists them): the colours are picked to read as
the thing they are, and they are deliberately not photographs of any one world.
"""

import hashlib
import json
import math

from planetgen.physics import photometry
from planetgen.web.maps.starmap import star_color

LOOK_VERSION = 1

_SOLAR_LUMINOSITY_W = 3.828e26

# Ocean colour by `ocean_class`; the neutral blue is the unlabelled default.
OCEAN_COLORS = {
    "neutral": "#1f4f8f",
    "chloride brine": "#2a6a7a",
    "acid sulfate": "#8a9a3a",
    "soda": "#3a8a8a",
    "ice-sealed": "#cfe0ec",
}
_DEFAULT_OCEAN_COLOR = OCEAN_COLORS["neutral"]

# Bare land by mean surface temperature (K): the upper bound of each band and the colour below it.
_LAND_BANDS = (
    (150.0, "#b8c4cc"),   # frozen, pale grey-blue ice-dusted rock
    (250.0, "#8f8a84"),   # cold grey-brown regolith
    (330.0, "#a39070"),   # temperate tan
    (700.0, "#a65a3a"),   # hot, rust
    (math.inf, "#6a2c20"),  # molten-dark
)
_DEFAULT_LAND_COLOR = "#8a8f9c"

# A giant's cloud top by temperature: methane blue, ammonia cream, water-cloud white, clear blue, alkali dark,
# silicate-cloud hot. Bands only; never read from the class description text.
_GIANT_BANDS = (
    (130.0, "#9fd0dc"),
    (250.0, "#d9c9a3"),
    (800.0, "#e8e2d6"),
    (1400.0, "#3b4a7a"),
    (2200.0, "#3a3340"),
    (math.inf, "#a0522d"),
)


def seed_for(kind, body_id):
    """32 bits of a stable hash of the body's kind and id: the same body gets the same picture on every machine."""
    digest = hashlib.blake2b(f"{kind}:{body_id}".encode(), digest_size=4).digest()
    return int.from_bytes(digest, "big")


def _band(bands, value):
    for upper, color in bands:
        if value < upper:
            return color
    return bands[-1][1]


def land_color(temp_k):
    """Bare-land colour estimate from the mean surface temperature alone."""
    return _DEFAULT_LAND_COLOR if temp_k is None else _band(_LAND_BANDS, temp_k)


def giant_color(temp_k):
    """Cloud-top colour estimate for a gas or ice giant from temperature alone."""
    return _DEFAULT_LAND_COLOR if temp_k is None else _band(_GIANT_BANDS, temp_k)


def glow_color(temp_k):
    """The atmosphere rim colour: orange-hazy when hot, pale blue when cold, sky blue between."""
    if temp_k is None:
        return "#bcdfff"
    if temp_k >= 320.0:
        return "#ffb066"
    if temp_k <= 200.0:
        return "#bcd7ff"
    return "#bfe3ff"


def _put(look, key, value):
    if value is not None:
        look[key] = value


def _rounded(value, digits):
    return None if value is None else round(float(value), digits)


def planet_look(planet, kind="planet"):
    """The look record of a planet or moon row."""
    atmosphere = planet.get("atmosphere")
    has_atmosphere = bool(atmosphere) and atmosphere != "None"
    temp_k = planet.get("surface_temperature_k")
    gas = planet.get("body_type") == "g"
    look = {
        "v": LOOK_VERSION, "kind": kind, "body": "gas" if gas else "solid",
        "seed": seed_for(kind, planet["id"]),
        "color": giant_color(temp_k) if gas else land_color(temp_k),
        "quality": "estimate",
    }
    _put(look, "cls", (planet.get("planet_class") or "").upper() or None)
    if has_atmosphere:
        look["atm"] = True
        look["glow"] = glow_color(temp_k)
    pressure_pa = planet.get("atmospheric_pressure_pa")
    _put(look, "pkpa", _rounded(pressure_pa / 1000.0 if pressure_pa is not None else None, 3))
    _put(look, "tempK", _rounded(temp_k, 1))
    _put(look, "radiusKm", _rounded(planet.get("radius_km"), 1))
    if not gas:
        _put(look, "hydro", planet.get("hydrosphere"))
        _put(look, "ocean", _rounded(planet.get("ocean_fraction"), 4))
        _put(look, "land", _rounded(planet.get("land_fraction"), 4))
        _put(look, "oceanClass", planet.get("ocean_class"))
        look["oceanColor"] = ocean_color(planet.get("ocean_class"))
        _put(look, "oceanKm", _rounded(planet.get("ocean_depth_km"), 3))
        _put(look, "iceKm", _rounded(planet.get("ice_shell_km"), 3))
    return look


def star_surface_color(star):
    """The colour of a star's disc: the blackbody-based `photometry.star_colour` of its effective temperature (the
    same function the sky colour assumes for the star, so the two agree), else the spectral-letter colour."""
    temperature_k = star.get("temperature_k")
    if temperature_k and temperature_k > 0:
        try:
            return photometry.star_colour(float(temperature_k)).hex
        except ValueError:
            pass
    return star_color(star["star_type"], temperature_k, star["luminosity_w"])[0]


def star_look(star):
    """The look record of a star row."""
    look = {"v": LOOK_VERSION, "kind": "star", "body": "star", "seed": seed_for("star", star["id"]),
            "color": star_surface_color(star)}
    _put(look, "yerkes", star.get("yerkes_class"))
    _put(look, "tempK", _rounded(star.get("temperature_k"), 1))
    luminosity_w = star.get("luminosity_w")
    _put(look, "lumSolar", None if luminosity_w is None else float(f"{luminosity_w / _SOLAR_LUMINOSITY_W:.4g}"))
    _put(look, "radiusKm", _rounded(star.get("radius_km"), 1))
    _put(look, "spec", star.get("star_type"))
    return look


def ocean_color(ocean_class):
    return OCEAN_COLORS.get(ocean_class or "", _DEFAULT_OCEAN_COLOR)


def look_json(look):
    """The record as compact, key-sorted JSON (the `data-look` value)."""
    return json.dumps(look, separators=(",", ":"), sort_keys=True)


def look_note(look):
    """The info panel's line about how the picture was worked out (UX.158 replaces it with the full account)."""
    if look.get("quality") == "estimate":
        return "Estimate from temperature and stored water fractions; not yet from composition"
    return None


# --- Phenomena (MAP.278) ----------------------------------------------------

# Grain colour of an asteroid field by composition family (stand-ins for hue only): carbonaceous dark, stony tan,
# metallic grey, icy pale blue-grey, basaltic dark brown, mixed and collisional between, dust beige.
ASTEROID_GRAIN_COLORS = {
    "carbonaceous": "#4a4540",
    "stony": "#a89878",
    "metallic": "#8c9096",
    "icy": "#b8c8d4",
    "basaltic": "#6b5a4e",
    "mixed": "#8a7a66",
    "dust": "#9c8f7c",
    "collisional": "#8a8480",
}
_DEFAULT_GRAIN_COLOR = "#b89a6e"


def asteroid_grain_color(class_code):
    """The grain colour of an asteroid field from its class code (`B3`: the letter names the composition family)."""
    from planetgen.generation.phenomena.asteroid_field import composition_family_for_letter
    letter = (class_code or "")[:1].upper()
    try:
        return ASTEROID_GRAIN_COLORS.get(composition_family_for_letter(letter), _DEFAULT_GRAIN_COLOR)
    except ValueError:
        return _DEFAULT_GRAIN_COLOR


def blackbody_hex(temperature_k):
    """The `#rrggbb` blackbody colour of `temperature_k`, `None` when it is not a usable temperature."""
    if not temperature_k or temperature_k <= 0:
        return None
    try:
        return photometry.star_colour(float(temperature_k)).hex
    except ValueError:
        return None


def rogue_planet_look(rogue):
    """The look record of a rogue planet row: the same painters as a planet, lit only by distant starlight, so
    no sea unless it is stored as liquid under ice, and ice where the surface is cold."""
    temp_k = rogue.get("surface_temperature_k")
    if temp_k is None:
        temp_k = rogue.get("effective_temperature_k")
    gas = rogue.get("planet_type") == "g"
    pressure_pa = rogue.get("surface_pressure_pa")
    row = {
        "id": rogue.get("id"), "planet_class": rogue.get("planet_class"), "body_type": "g" if gas else "t",
        "surface_temperature_k": temp_k, "radius_km": rogue.get("radius_km"),
        "atmosphere": "yes" if pressure_pa else None, "atmospheric_pressure_pa": pressure_pa,
    }
    if not gas:
        if rogue.get("has_liquid_water") and rogue.get("ice_shell_thickness_km"):
            row["hydrosphere"] = "ice-covered ocean"
        elif temp_k is not None and temp_k < 200.0:
            row["hydrosphere"] = "ice"
        else:
            row["hydrosphere"] = "dry"
    look = planet_look(row, "rogue_planet")
    look.pop("glow", None)
    return look


def comet_nucleus_color(comet):
    """The nucleus colour of an interstellar comet: grey-white when it is active (volatile-rich, shedding ice),
    dark red-brown when it is inactive (a crust of dark organics)."""
    return "#9ea4ac" if comet.get("is_active") else "#4a3a32"
