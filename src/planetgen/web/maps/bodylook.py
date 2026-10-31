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
    quality     "composition" when a solid world's `color` is the lab colour of its stored ground types (MAP.274),
                "estimate" while it comes from temperature alone (no stored ground, no land, or too little of the land
                with a measured colour), absent for a star
    landKnown   with quality "composition": the share of the land whose ground types have a measured colour, 0.5 to 1
    glow        atmosphere limb colour (`#rrggbb`): the atmosphere colour of the body's air under its star (MAP.244,
                `physics/atmosphere_colour.py`), present only when the body has air and the row carries its star's
                temperature (`_star_temperature_k`, put there by `systemmap.py`)
    atm         true when the body has an atmosphere
    glowStrength  how strongly the limb glows, 0 to `GLOW_STRENGTH_MAX` (MAP.277): from the column of the air (`air_glow`); a
                thin air has almost none, a thick one the most, and a body with no air has no glow fields at all
    glowWidth   how far the glow reaches past the limb, in radii: the air's scale height over the radius, exaggerated
                (`GLOW_EXAGGERATION`) and held to `GLOW_WIDTH_RANGE`; left out when the air's composition is not stored
    veil        (solid world with air, MAP.350) how far the surface is tinted toward the sky colour `glow`, 0 to 1: `1 - exp(-tau_550)`,
                the vertical optical depth of the air at 550 nm (gas and aerosol, `AtmosphereColour.tau_550`); a thin air leaves the ground as
                it is, a thick one drowns it in its own colour
    haze        the washing-out of the ground's contrast by haze and dust, 0 to 1: `1 - exp(-(tau_haze + tau_dust))`; left out when neither is stored
    cloudCover  the share of the sky the cloud covers (the stored `cloud_cover`), left out when none is stored
    cloudTau    the vertical optical depth of the cloud at 550 nm (the stored `tau_cloud`), left out when none is stored; the painter
                draws a cloud layer over the ground, thicker and whiter with the depth
    twilight    how far past the day and night line the light of the star is carried by the air, as a cosine of the light
                on the ground (0 for no air to `TWILIGHT_MAX` for a thick one): the terminator of a world with air is softer
                and lies a little onto the night side
    pkpa        surface pressure, kPa
    tempK       mean surface temperature, K (a star's effective temperature)
    radiusKm    radius
    roundness   the roundness ratio of a solid body (GEN.1158, GEN.1159): its radius over the radius at which its weight beats
                its strength, 1 or more round; under 1 the body is drawn as a lump (MAP.279, `static/lumpymesh.js`); left out
                when none is stored
    hydro       `hydrosphere` (dry, vapour, ice, ice-covered ocean, surface ocean, brine seep, hycean)
    ocean       fraction of the surface under ocean, 0 to 1
    land        fraction of the surface above sea level, 0 to 1
    oceanClass  ice-sealed, chloride brine, acid sulfate, soda or neutral
    oceanColor  the sea colour for that class, `#rrggbb`
    oceanKm     ocean depth, km
    iceKm       ice shell thickness, km
    lumSolar    a star's luminosity in solar units
    spec        a star's `star_type`
    yerkes      a star's Yerkes luminosity class (V dwarf, III giant, VII white dwarf, ...)
    spinH       rotation period, hours (MAP.328); left out when none is stored, and the body is drawn still
    axis        the unit spin axis [x, y, z] in the system's frame, right-handed (a retrograde body's points against its
                orbit's normal)
    tiltDeg     axial tilt, 0 to 180 degrees (past 90 the body turns backwards)
    locked      true when the body is tidally locked (its rotation period is its orbital period, UX.244); the views then keep
                its day side toward what it orbits, and the object page says it is locked

The colour tables are stand-ins for hue only (the gap plan registry lists them): the colours are picked to read as
the thing they are, and they are deliberately not photographs of any one world.
"""

import hashlib
import json
import math
from types import SimpleNamespace

from planetgen.physics import (
    atmosphere, atmosphere_colour, gases, giant_colour, photometry, regional_temperature, stellar_library, surface_composition)
from planetgen.web.maps.starmap import star_color
from planetgen.physics import unit_registry as units

_pa_to_kpa = units.converter("Pa", "kPa")

LOOK_VERSION = 1

# MAP.277: the glow of an atmosphere from its column. The column of the air is its mass over an area, `P / g`; the glow is the share of the light it
# scatters out, taken as `1 - exp(-column / EARTH_COLUMN / GLOW_COLUMN_SCALE)`. Display constants (a picture, not physics): the numbers are picked so
# that the Earth's air glows clearly, Mars's (about 1.5 percent of the Earth's column) hardly at all, and Venus's or Titan's most.
EARTH_COLUMN_KG_M2 = 101325.0 / 9.80665
GLOW_COLUMN_SCALE = 0.6
GLOW_STRENGTH_MAX = 0.85
GLOW_STRENGTH_MIN = 0.01
"""float: A glow weaker than this is not drawn (and is not on the look record): the air is too thin to show."""
GIANT_GLOW_STRENGTH = 0.7
"""float: A giant's air is always deep (the 1 bar level is far above its cloud tops' column), so its glow is a fixed strong one."""
GLOW_EXAGGERATION = 50.0
"""float: The scale height over the radius is about 0.001 for the Earth, too thin to see; the glow is drawn this many times wider."""
GLOW_WIDTH_RANGE = (0.03, 0.2)
"""tuple: The least and most the glow reaches past the limb, in radii."""
TWILIGHT_MAX = 0.2
"""float: The twilight (a cosine of the light on the ground) of the thickest air."""
VEIL_TAU_SCALE = 1.0
"""float: The air's vertical optical depth at 550 nm over which the surface is half-way to the sky colour: `veil = 1 - exp(-tau / VEIL_TAU_SCALE)`."""
HAZE_TAU_SCALE = 1.0
"""float: The same for the haze and dust depths, which wash out the ground's contrast: `haze = 1 - exp(-(tau_haze + tau_dust) / HAZE_TAU_SCALE)`."""
_BOLTZMANN = 1.380649e-23

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


# The share of the land that must have a measured colour for the composition colour to stand in for the estimate.
MIN_KNOWN_LAND_SHARE = 0.5


def _srgb_byte(linear):
    """The 8-bit sRGB value of a linear reflectance (clipped to 0 to 1)."""
    linear = min(1.0, max(0.0, linear))
    encoded = 12.92 * linear if linear <= 0.0031308 else 1.055 * linear ** (1.0 / 2.4) - 0.055
    return round(255.0 * encoded)


def composition_land_color(row):
    """
    The colour of a solid world's bare land from its stored ground types (MAP.274): the lab reflectance in the 450, 550 and 650 nm bands
    (`surface_composition.band_albedo`) of the areal mix of the land types, read as blue, green and red under white light, linear to sRGB.
    The ocean and the ice are not land (the painter draws them on their own), and a land type with no measured visible spectrum (the metal, the
    intermediate mixes) is left out of the average.

    Returns:
        tuple: `(#rrggbb, known share of the land)`, or `None` when the row stores no ground, has no land, or fewer than
        `MIN_KNOWN_LAND_SHARE` of its land has a measured colour.
    """
    mix = surface_composition.mix_of(planet_row_namespace(row))
    if not mix:
        return None
    land = {name: share for name, share in mix.items() if name not in (surface_composition.OCEAN, surface_composition.ICE)}
    total = math.fsum(land.values())
    if not total > 0.0:
        return None
    composition = surface_composition.SurfaceComposition(
        ground_type=row.get("ground_type") or max(land, key=land.get), mix=land, components={},
        grain=row.get("ground_grain") or surface_composition.POWDER)
    bands, covered = surface_composition.band_albedo(composition)
    known = covered / total
    values = [bands[band] for band in surface_composition.BANDS]
    if known < MIN_KNOWN_LAND_SHARE or any(value is None for value in values):
        return None
    blue, green, red = values
    return "#{:02x}{:02x}{:02x}".format(_srgb_byte(red), _srgb_byte(green), _srgb_byte(blue)), min(1.0, known)


def planet_row_namespace(row):
    """The row's ground columns as attributes, the shape `surface_composition.mix_of` reads."""
    return SimpleNamespace(**{field: row.get(field) for field in surface_composition.FRACTION_FIELDS})


def giant_color(temp_k):
    """Cloud-top colour estimate for a gas or ice giant from temperature alone."""
    return _DEFAULT_LAND_COLOR if temp_k is None else _band(_GIANT_BANDS, temp_k)


def air_colour(planet):
    """
    The `AtmosphereColour` (a solid world) or `GiantColour` (a giant) of a planet or moon row: its stored mixture and ozone column
    under its star (`_star_temperature_k`). `None` when it has no air, no star temperature or no stored mixture.
    """
    atmosphere_text = planet.get("atmosphere")
    temperature = planet.get("_star_temperature_k")
    if not atmosphere_text or atmosphere_text == "None" or not temperature:
        return None
    body = SimpleNamespace(
        body_type=planet.get("body_type"), gravity=planet.get("gravity_g"), atmospheric_pressure=planet.get("atmospheric_pressure_pa"),
        planet_class=planet.get("planet_class"), ozone_column_du=planet.get("ozone_column_du"), star=SimpleNamespace(temperature=temperature),
        **{atmosphere.GAS_FIELD: atmosphere.gas_partials_pa(planet)})
    if planet.get("body_type") == "g":
        return giant_colour.body_giant_colour(body)
    return atmosphere_colour.body_atmosphere_colour(body)


def _mean_molar_mass_kg_mol(partials_pa):
    """The mean molar mass of an air from its partial pressures (formula -> Pa), kg/mol; `None` for none stored or a gas the table lacks."""
    total = sum(partials_pa.values())
    if total <= 0.0:
        return None
    try:
        return sum(pressure * gases.gas(formula).molar_mass for formula, pressure in partials_pa.items()) / total
    except KeyError:
        return None


def air_glow(planet):
    """
    The glow of a body's air (MAP.277): `{"glowStrength", "glowWidth", "twilight"}` from its surface pressure, gravity and the gases it stores;
    `{}` for a body with no air or air too thin to show (`GLOW_STRENGTH_MIN`). A giant's glow is a fixed strong one (`GIANT_GLOW_STRENGTH`).

    The strength is `GLOW_STRENGTH_MAX (1 - exp(-column / EARTH_COLUMN_KG_M2 / GLOW_COLUMN_SCALE))`, the column being `P / g`: the Earth's air
    gets 0.69 and Mars's 0.02, Venus's and Titan's the most. The width is the scale height `k T / (m g)` over the radius times `GLOW_EXAGGERATION`
    (needs the stored gases for the mean molar mass, the temperature and the radius). The twilight grows with the strength.
    """
    pressure_pa, gravity_g = planet.get("atmospheric_pressure_pa"), planet.get("gravity_g")
    if planet.get("body_type") == "g":
        strength = GIANT_GLOW_STRENGTH
    elif pressure_pa and gravity_g and pressure_pa > 0 and gravity_g > 0:
        column = pressure_pa / (gravity_g * 9.80665) / EARTH_COLUMN_KG_M2
        strength = GLOW_STRENGTH_MAX * -math.expm1(-column / GLOW_COLUMN_SCALE)
    else:
        return {}
    if strength < GLOW_STRENGTH_MIN:
        return {}
    glow = {"glowStrength": round(strength, 3), "twilight": round(TWILIGHT_MAX * strength / GLOW_STRENGTH_MAX, 3)}
    mass = _mean_molar_mass_kg_mol(atmosphere.gas_partials_pa(planet))
    temperature, radius_km = planet.get("surface_temperature_k"), planet.get("radius_km")
    if mass and temperature and radius_km and gravity_g:
        height_m = _BOLTZMANN * temperature * 6.02214076e23 / (mass * gravity_g * 9.80665)
        width = GLOW_EXAGGERATION * height_m / (radius_km * 1000.0)
        glow["glowWidth"] = round(min(GLOW_WIDTH_RANGE[1], max(GLOW_WIDTH_RANGE[0], width)), 3)
    return glow


def _put(look, key, value):
    if value is not None:
        look[key] = value


def _rounded(value, digits):
    return None if value is None else round(float(value), digits)


def _put_veil(look, row, colour):
    """
    What the air does to the ground seen through it (MAP.350): `veil` from the air's optical depth, `haze` from the stored haze and dust depths, and
    `cloudCover` and `cloudTau` as stored. Each is left out when its input is not stored, and the painter then leaves the ground as it was.
    """
    if colour is not None and colour.tau_550 is not None and math.isfinite(colour.tau_550):
        look["veil"] = round(-math.expm1(-max(0.0, colour.tau_550) / VEIL_TAU_SCALE), 3)
    depths = [row.get(name) for name in ("tau_haze", "tau_dust")]
    depths = [d for d in depths if d is not None and math.isfinite(d)]
    if depths:
        look["haze"] = round(-math.expm1(-max(0.0, math.fsum(depths)) / HAZE_TAU_SCALE), 3)
    cover, depth = row.get("cloud_cover"), row.get("tau_cloud")
    if cover is not None and depth is not None and math.isfinite(cover) and math.isfinite(depth) and cover > 0.0 and depth > 0.0:
        look["cloudCover"] = round(min(1.0, cover), 3)
        look["cloudTau"] = round(depth, 3)


def _put_spin(look, row):
    """The spin the row stores, for `static/spin.js` (MAP.328): `spinH` (rotation period, hours), `axis` (the unit spin
    axis, right-handed, in the system's frame) and `tiltDeg` (axial tilt, 0 to 180; past 90 the body turns backwards).
    A key is left out when the row stores nothing for it; a body with no `spinH` is drawn still, never with a made-up
    period."""
    period = row.get("rotation_period_hours")
    if period is not None and math.isfinite(period) and period > 0:
        look["spinH"] = float(f"{period:.4g}")
    axis = [row.get(f"spin_axis_{name}") for name in "xyz"]
    if all(v is not None and math.isfinite(v) for v in axis):
        length = math.sqrt(sum(v * v for v in axis))
        if length > 0:
            look["axis"] = [round(v / length, 4) for v in axis]
    _put(look, "tiltDeg", _rounded(row.get("axial_tilt_deg"), 2))
    if regional_temperature.is_locked(period, row.get("period_years")):
        look["locked"] = True


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
    if not gas:
        found = composition_land_color(planet)
        if found is not None:
            look["color"], known_land = found
            look["quality"] = "composition"
            look["landKnown"] = round(known_land, 3)
    _put(look, "cls", (planet.get("planet_class") or "").upper() or None)
    _put_spin(look, planet)
    if has_atmosphere:
        look["atm"] = True
        glow = air_glow(planet)
        colour = air_colour(planet) if glow else None       # an air too thin to show has no glow at all
        if colour is not None:
            look["glow"] = (colour.atmosphere if gas else colour).hex
        look.update(glow)
        if not gas:
            _put_veil(look, planet, colour)
    pressure_pa = planet.get("atmospheric_pressure_pa")
    _put(look, "pkpa", _rounded(_pa_to_kpa(pressure_pa) if pressure_pa is not None else None, 3))
    _put(look, "tempK", _rounded(temp_k, 1))
    _put(look, "radiusKm", _rounded(planet.get("radius_km"), 1))
    if not gas:
        _put(look, "roundness", _rounded(planet.get("roundness_ratio"), 4))
        _put(look, "hydro", planet.get("hydrosphere"))
        _put(look, "ocean", _rounded(planet.get("ocean_fraction"), 4))
        _put(look, "land", _rounded(planet.get("land_fraction"), 4))
        _put(look, "oceanClass", planet.get("ocean_class"))
        look["oceanColor"] = ocean_color(planet.get("ocean_class"))
        _put(look, "oceanKm", _rounded(planet.get("ocean_depth_km"), 3))
        _put(look, "iceKm", _rounded(planet.get("ice_shell_km"), 3))
    return look


def star_surface_color(star):
    """The colour of a star's disc: that of its model atmosphere (`stellar_library.star_colour`, GEN.971: the spectrum the sky colour
    assumes for the star, so the two agree) at its effective temperature and surface gravity, else the spectral-letter colour."""
    temperature_k = star.get("temperature_k")
    if temperature_k and temperature_k > 0:
        mass_kg, radius_km = star.get("mass_kg"), star.get("radius_km")
        logg = stellar_library.surface_logg(mass_kg, radius_km) if mass_kg and radius_km else stellar_library.DEFAULT_LOGG
        try:
            return stellar_library.star_colour(float(temperature_k), logg).hex
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
    _put_spin(look, star)
    return look


def ocean_color(ocean_class):
    return OCEAN_COLORS.get(ocean_class or "", _DEFAULT_OCEAN_COLOR)


def look_json(look):
    """The record as compact, key-sorted JSON (the `data-look` value)."""
    return json.dumps(look, separators=(",", ":"), sort_keys=True)


def colour_account(look):
    """
    The info panel's account of how the colour was worked out (UX.158): which stored inputs it used and whether it is a partial estimate
    or the full account. A pure function of the look record, which carries exactly the inputs the painters had, so it cannot say more than
    the picture used. `None` for a record that is not a body.

    A star is full when its colour is a model atmosphere's at its temperature, partial when only its spectral class gave it. A solid world
    or a giant is a partial estimate for as long as the record says `quality` "estimate" (the colour from temperature and body family, not
    from the surface composition, cloud and haze) or "composition" (a solid world whose land colour is the lab reflectance of its stored ground
    types, MAP.274; the cloud and haze are still not stored).
    """
    body = look.get("body")
    temperature = look.get("tempK")
    if body == "star":
        if temperature is None:
            return "Partial: colour from its spectral class alone; no temperature is stored for a model atmosphere"
        return f"Full: colour of a model atmosphere at its effective temperature ({temperature:,.0f} K) and gravity"
    if body not in ("solid", "gas"):
        return None
    kind = look.get("kind")
    layer = "Cloud-top" if body == "gas" else "Ground"
    if temperature is not None:
        parts = [f"{layer} colour from its mean {'cloud-top' if body == 'gas' else 'surface'} temperature ({temperature:,.0f} K) alone"]
    else:
        parts = [f"{layer} colour is a default: no temperature is stored"]
    if body == "solid" and look.get("quality") == "composition":
        known = look.get("landKnown")
        parts[0] = "Ground colour from the stored ground types (lab reflectance)" + (
            "" if (known or 1.0) >= 0.995 else "; %d percent of the land has no measured colour and is left out" % round(100.0 * (1.0 - known)))
    if body == "solid" and look.get("ocean"):
        if look.get("oceanClass"):
            parts.append(f"sea colour from its ocean class ({look['oceanClass']})")
        else:
            parts.append("sea colour is a default: no ocean class is stored")
    if look.get("atm") and kind != "rogue_planet":
        if look.get("glow"):
            parts.append("limb glow from the composition of its air under its star")
        else:
            parts.append("no limb glow: the colour of its air could not be worked out")
    if look.get("quality") == "estimate":
        missing = "surface composition, clouds and haze" if body == "solid" else "chemistry bands and haze"
        head = "Partial estimate: "
        tail = f". Not used yet: {missing} (not stored)"
    elif look.get("quality") == "composition":
        # MAP.274: the land colour is the stored ground types' lab reflectance; the clouds and haze are still not stored
        head = "Partial estimate: "
        tail = ". Not used yet: clouds and haze (not stored)"
    else:
        head, tail = "Full: ", ""
    return head + "; ".join(parts) + tail


def look_note(look):
    """The info panel's line about how the picture was worked out (`colour_account`), then how it is turned or locked."""
    notes = []
    account = colour_account(look)
    if account:
        notes.append(account)
    if look.get("locked") and look.get("body") != "star":
        notes.append("Tidally locked: its day side stays toward " + ("its planet" if look.get("kind") == "moon" else "its star"))
    elif look.get("spinH") and look.get("body") == "solid":
        notes.append("Seen from the side of its orbit, with the day and night line down the middle, turning much faster than its real rotation")
    elif look.get("spinH"):
        notes.append("Turning much faster than its real rotation")
    return "; ".join(notes) or None


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
    for key in ("glow", "glowStrength", "glowWidth", "twilight"):     # no star: nothing lights the air
        look.pop(key, None)
    return look


def comet_nucleus_color(comet):
    """The nucleus colour of an interstellar comet: grey-white when it is active (volatile-rich, shedding ice),
    dark red-brown when it is inactive (a crust of dark organics)."""
    return "#9ea4ac" if comet.get("is_active") else "#4a3a32"
