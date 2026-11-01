# planetgen/physics/mixture_subject.py

"""
The description of a planet or moon that the atmosphere mixture draw reads (GEN.756; docs/design/composition/mixture-draw.md)
==============================================================================================================

`mixture_draw.draw_mixture` takes a `Subject`: the body's class and the quantities the recipe tables and the chemistry read
(`mixture_draw.required_inputs(class)`). This module builds that subject from a generated body and its star, so that the generator
(GEN.411) calls the draw without writing the translation itself:

    subject = subject_of_body(planet, ts_k=288.0, pressure_bar=1.0, ocean_present=True, long_lived_geology=True)
    mixture_draw.draw_mixture(subject)

What the builder derives, from the body's stored fields and its star: the mass, radius and escape speed, the age of the star, the
starlight and the equilibrium temperature at the body's Bond albedo, the XUV exposure index, the position against the snow line, the
runaway-greenhouse limit and the habitable zone of the host, and for a giant the temperature of its 1 bar layer, its dayside and its
effective temperature. What it does not derive, and takes as arguments: the surface temperature and the surface pressure (the climate
loop's, which the mixture feeds back into), whether an ocean can condense and whether the body keeps a long-lived geology (the surface
pressure code's liquid-water range and the interior's), and whether the water was lost (the water-loss clock's).

The subject holds exactly `required_inputs(class)`; a quantity that cannot be built (a rogue body with no star, a moon whose distance from
its star is not given) is refused by name (`MixtureError`), never defaulted.
"""

import math

from planetgen import tuning
from planetgen.generation.possibility_table import Subject
from planetgen.physics import activity, bulk_composition, bulk_structure, climate_draws, constants, giant_draw, giant_interior, ice_species, mixture_draw, orbit_elements, orbits, unit_registry as units
from planetgen.physics.mixture_base import MixtureError

_kg_to_mearth = units.converter("kg", "Mearth")
_kg_to_msun = units.converter("kg", "Msun")

PHOTOSPHERE_BAR = 0.1
"""float: The pressure of the layer a giant's regime is read at, bar: the dayside emission of a hot giant comes from about 0.1 bar (the
recipe quantity's own definition, `atmosphere_recipes.INPUTS`). A stand-in (AM27)."""

DAYSIDE_OVER_EQUILIBRIUM = (8.0 / 3.0) ** 0.25
"""float: The dayside temperature of a body with no heat redistribution over its equilibrium temperature: the starlight is spread over
the day hemisphere alone (`f = 2/3` against the `f = 1/4` of the whole surface), so `T` rises by `(f / (1/4))^(1/4) = (8/3)^(1/4)`.
Derived; the generator has no redistribution model yet, so every body gets the hot end."""

FILE_HOST_TYPES = {"F": "F", "G": "G", "K": "K", "M": "M"}
"""dict: The spectral letters of the recipe tables' hosts. A hotter star (O, B, A, or a Wolf-Rayet W) takes the F rows and a brown dwarf (L, T, Y) the M rows,
the nearest the tables have; a white dwarf is `WD`, a pulsar or a black hole `PSR`."""


def host_type_of(star):
    """
    The recipe tables' host type (`F`, `G`, `K`, `M`, `WD`, `PSR`) of `star`.

    Raises:
        MixtureError: The star carries no spectral type.
    """
    kind = getattr(star, "type", None)
    if not kind:
        raise MixtureError("the star has no spectral type, so the host type of the mixture draw cannot be built")
    letter = str(kind).strip()[:1].upper()
    name = type(star).__name__.lower()
    if "pulsar" in name or "neutron" in name or "blackhole" in name:  # a black hole is a bare remnant like a pulsar: no photosphere
        return "PSR"
    if letter == "D" or getattr(star, "yerkes_class", None) in ("VII", "D"):
        return "WD"
    if letter in FILE_HOST_TYPES:
        return FILE_HOST_TYPES[letter]
    if letter in "OBAW":
        return "F"
    if letter in "LTY":
        return "M"
    raise MixtureError(f"the star's spectral type {kind!r} has no host type in the recipe tables")


def runaway_ratio(insolation, star_temperature_k):
    """The starlight over the runaway-greenhouse limit of the host (the inner edge of its habitable zone, [R2447]), both
    against the Earth's."""
    low, high = orbits.KOPPARAPU_TEFF_RANGE_K
    limit = orbits.kopparapu_flux(orbits.KOPPARAPU_RUNAWAY_GREENHOUSE, min(max(star_temperature_k, low), high))
    return insolation / limit


def quantities_of_body(body, *, ts_k, pressure_bar, ocean_present=False, long_lived_geology=False, water_lost=False,
                       star_distance_au=None):
    """
    Every quantity of the mixture draw that can be built for `body`, as a dict (a quantity that cannot is absent).

    Args:
        body: A generated planet or moon (`generation.planet.Planet`): `mass` (kg), `radius` (km), `distance` (AU, from its primary),
            `star`, `bond_albedo`, `is_moon`; a giant's `core_mass_earth` or a rocky planet's `envelope_mass_fraction`.
        ts_k (float): The surface temperature, K.
        pressure_bar (float): The surface pressure, bar.
        ocean_present (bool): True when an ocean can condense (the climate loop's test against `surface_pressure.LIQUID_WATER_K`).
        long_lived_geology (bool): True for a tectonic or long-lived body.
        water_lost (bool): True when the water was lost.
        star_distance_au (float or None): The distance of the body (a moon: of its parent) from the star, as the circle that gives the same
            mean light as its orbit (`orbit_elements.flux_distance_au`, GEN.727); for a planet the same circle for its own `distance` and
            `eccentricity` when not given (GEN.842).
    """
    star = getattr(body, "star", None)
    mass_earth = _kg_to_mearth(body.mass)
    quantities = {
        "mass_earth": mass_earth,
        "mass_jupiter": mass_earth / giant_interior.EARTH_MASS_JUPITER,
        "radius_km": body.radius,
        "v_esc_kms": math.sqrt(2.0 * constants.G * body.mass / (body.radius * constants.KM_TO_M_FACTOR)) / 1000.0,
        "ts_k": ts_k, "pressure_bar": pressure_bar,
        "ocean_present": bool(ocean_present), "long_lived_geology": bool(long_lived_geology), "water_lost": bool(water_lost),
        "is_moon": bool(getattr(body, "is_moon", False)),
        "photosphere_bar": PHOTOSPHERE_BAR,
    }
    fractions = bulk_composition.bulk_fractions(body)
    if fractions is not None and fractions.rock > 0.0:   # GEN.821: the ammonia ice a Titan-like air takes its nitrogen from
        shares = fractions.ice_shares or ice_species.median_shares()
        quantities["nh3_mass_fraction"] = fractions.ice * shares.get("NH3", 0.0)
        quantities["water_mass_fraction"] = fractions.water
        quantities["rock_mass_fraction"] = fractions.rock
    core = getattr(body, "core_mass_earth", None)
    envelope = getattr(body, "envelope_mass_fraction", None)
    if not core and envelope is not None:  # a rocky planet's core is the mass its envelope (none, for a stripped one) leaves
        core = mass_earth * (1.0 - envelope)
    if core:
        quantities["core_radius_earth"] = bulk_structure.core_radius_earth(core)
    if star is None:
        return quantities
    quantities["age_gy"] = getattr(star, "age", None) or tuning.ROGUE_EARTH_REFERENCE_AGE_GY
    quantities["host_type"] = host_type_of(star)
    quantities["host_mass_solar"] = _kg_to_msun(star.mass)
    quantities["luminosity_solar"] = star.luminosity / constants.SOLAR_LUMINOSITY
    distance_au = star_distance_au if star_distance_au is not None else (
        None if quantities["is_moon"] else orbit_elements.flux_distance_au(body.distance, getattr(body, "eccentricity", None)))
    if not distance_au:
        return quantities
    insolation = quantities["luminosity_solar"] / distance_au ** 2
    quantities["insolation"] = insolation
    quantities["runaway_ratio"] = runaway_ratio(insolation, star.temperature)
    snow_line = getattr(star, "snow_line_au", None)
    if snow_line:
        quantities["snow_ratio"] = distance_au / snow_line
    zone_edges = getattr(star, "habitable_zone", None)
    if zone_edges:
        inner, outer = zone_edges
        quantities["zone"] = "h" if distance_au < inner else ("c" if distance_au > outer else "e")
    # GEN.806: the index `magnetism.update_exposure` stored for the distance the body sits at; the one implementation (`activity.exposure_index`) answers
    # for a distance it was not stored for (the mantle's redox is read at another)
    stored = getattr(body, "xuv_exposure_index", None)
    same_distance = getattr(body, "_star_distance_au", None) == distance_au
    index = stored if stored is not None and same_distance else activity.exposure_index(star, distance_au)
    if index is not None:
        quantities["xuv_index"] = index
    if getattr(star, "log_lx_lbol", None) is not None:
        quantities["log_lx_lbol"] = star.log_lx_lbol
    flux = climate_draws.starlight_w_m2(star.luminosity, distance_au * constants.AU_TO_KM)
    albedo = getattr(body, "bond_albedo", None)
    if albedo is not None:
        teq = climate_draws.equilibrium_temperature_k(albedo, flux)
        quantities["teq_k"] = teq
        quantities["tday_k"] = teq * DAYSIDE_OVER_EQUILIBRIUM
        quantities["tp_k"] = teq * DAYSIDE_OVER_EQUILIBRIUM
    effective = giant_draw.effective_temperature_k(flux, giant_draw.own_flux_w_m2(body))
    quantities["teff_k"] = effective
    quantities["t_layer_k"] = giant_draw.ONE_BAR_OVER_EFFECTIVE * effective
    return quantities


OPTIONAL_INPUTS = ("log_lx_lbol", "host_mass_solar", "nh3_mass_fraction", "water_mass_fraction", "rock_mass_fraction")
"""tuple: The quantities the escape filter reads when the body has them (GEN.804): the star's own activity at its age, which the XUV history follows
(the table's nominal decay without it). The magnetopause (`mixture_escape` reads `magnetopause_rp` when a subject carries it, 0 for no field) is not built here: the
field is drawn after the air until GEN.967 orders it first, and a pass that ran after it would close the air to the wind for some bodies and not others.
`nh3_mass_fraction`, `water_mass_fraction` and `rock_mass_fraction` (GEN.821) are the ammonia ice, the water and the rock of the stored bulk composition, as shares
of the whole body's mass: the nitrogen of a Titan-like air (`mixture_evolved.titan_like`)."""


KIND_OF_BODY_TYPE = {"t": "terrestrial", "g": "giant"}
"""dict: The recipe tables' kind of body (`atmosphere_recipes.BODY_KINDS`) by a planet's or moon's stored body type."""


def kind_of_body(body):
    """The recipe tables' kind of body of a planet or moon: `terrestrial` or `giant`, by its stored body type.

    Raises:
        MixtureError: The body has no body type of `t` or `g`.
    """
    try:
        return KIND_OF_BODY_TYPE[getattr(body, "body_type", None)]
    except KeyError:
        raise MixtureError(f"the body type {getattr(body, 'body_type', None)!r} is neither 't' nor 'g', so the recipe kind cannot be built") from None


def subject_of_body(body, *, ts_k, pressure_bar, ocean_present=False, long_lived_geology=False, water_lost=False,
                    star_distance_au=None):
    """
    The `Subject` of `body` for `mixture_draw.draw_mixture`: its kind and exactly `mixture_draw.required_inputs(kind)`. The arguments
    are those of `quantities_of_body`.

    Raises:
        MixtureError: A quantity the class's rows read cannot be built for this body (named).
    """
    kind = kind_of_body(body)
    built = quantities_of_body(body, ts_k=ts_k, pressure_bar=pressure_bar, ocean_present=ocean_present,
                               long_lived_geology=long_lived_geology, water_lost=water_lost, star_distance_au=star_distance_au)
    needed = mixture_draw.required_inputs(kind)
    missing = [name for name in needed if name not in built]
    if missing:
        raise MixtureError(f"the mixture draw of a {kind} reads {', '.join(missing)}, which cannot be built for this body"
                           f" (no star, no distance from the star, no Bond albedo or no stored interior)")
    optional = {name: built[name] for name in OPTIONAL_INPUTS if built.get(name) is not None}
    return Subject(kind, {**{name: built[name] for name in needed}, **optional})
