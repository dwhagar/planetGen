# planetgen/physics/planets.py

"""
Planet Physical and Orbital Generation
=======================================

This module determines what a planet or moon physically *is*: its orbital
zone, kind, mass, radius, density, composition, atmosphere, surface
gravity, and atmospheric conditions, plus the generation of any moons it
has. No class letter is read: the mass comes from the formation law
(`planet_kind.draw_formation`) or a moon's size law, the composition from the
star and where the body formed, and the class is derived from the finished
body (`apply_class`, `planet_classes.classify`). It operates on `Planet` instances (from `planet`) passed in as
arguments, rather than as methods on the class, so that generation logic is
kept separate from the `Planet` object's own definition and presentation.

Life chemistry and evolutionary data are intentionally NOT determined here —
see `life.apply_life_data`, which `StarSystem` applies to every planet
and moon in a system after all of them have been generated.
"""

import contextlib
import math

from planetgen.physics import moon_history, body_state, bulk_draw, clearance_ledger, composition_text, constants, envelope_capture, giant_draw, giant_interior, giant_tilt, kepler, metallicity, moon_composition, moon_disk_zone, orbit_elements, planet_classes, planet_kind, planet_mass, rings, oddities
from planetgen.physics import aerosols, atmosphere, body_colour, body_spin, body_strength, class_code, daylight, envelope_share, planet_type, spin_torques, climate_draws, habitability_world, hosting_gate, hydrosphere, interior_heat, magnetism, ocean_tides, phosphorus_budget, radiation, radiogenic, relations, spin, surface_air, surface_composition, tidal_spin, air_branches, venus_spin
from planetgen.physics import disk_truncation, eclipses, isotope_record, ore_deposits, radiation_field, sky_suns
from planetgen import tuning
from planetgen.generation import irregular_moons
from planetgen.util import draw
from planetgen.util import log
from planetgen.physics import orbits
from planetgen.physics.orbits import calculate_hill_sphere, calculate_reflex_offset, minimum_update_interval_years
from planetgen.physics.state_vectors import state_from_elements
from planetgen.util.checks import finite_domain
from planetgen.util.random import sample_bounded_bell, sample_power_law
from planetgen.physics import unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")
_kg_to_msun = units.converter("kg", "Msun")
_mearth_to_kg = units.converter("Mearth", "kg")

_km_to_m = units.converter("km", "m")
_m_to_km = units.converter("m", "km")
_au_to_m = units.converter("au", "m")
_au_to_km = units.converter("au", "km")
_km_to_au = units.converter("km", "au")
_gcm3_to_kgm3 = units.converter("g_cm3", "kg_m3")
_h_to_s = units.converter("h", "s")
_yr_to_s = units.converter("yr", "s")
_gyr_to_yr = units.converter("Gyr", "yr")
_yr_to_h = units.converter("yr", "h")
_au_yr_to_kms = units.converter("au_yr", "km_s")
_TWO_PI_RAD = 2.0 * math.pi

CAP_TOLERANCE = 1e-9
"""float: How far above the stored temperature the air's uncapped warming must reach to count as held at `surface_air.MAX_SURFACE_TEMPERATURE_K` (GEN.1236)."""


def uses_giant_mass_radius(planet):
    """Whether `planet` gets its mass and radius from the giant-planet mass-radius relation (GEN.34): a giant."""
    return planet.body_type == "g"


def _giant_transition_earth():
    """`(mass, radius)` in Earth units where the Neptunian and Jovian
    branches of the giant mass-radius relation meet."""
    m_n, r_n, s_n = constants.GIANT_NEPTUNIAN_MASS_RADIUS
    m_j, r_j, s_j = constants.GIANT_JOVIAN_MASS_RADIUS
    # r_n * (M/m_n)^s_n == r_j * (M/m_j)^s_j, solved in log space.
    log_m = (math.log(r_j / r_n) + s_n * math.log(m_n) - s_j * math.log(m_j)) / (s_n - s_j)
    mass = math.exp(log_m)
    return mass, r_n * (mass / m_n) ** s_n


GIANT_TRANSITION_MASS_EARTH, GIANT_TRANSITION_RADIUS_EARTH = _giant_transition_earth()


def giant_regime(mass_kg):
    """`"neptunian"` below the relation's transition mass, else `"jovian"`."""
    mass_earth = _kg_to_mearth(mass_kg)
    return "neptunian" if mass_earth < GIANT_TRANSITION_MASS_EARTH else "jovian"


def giant_radius_km(mass_kg):
    """The giant mass-radius relation's median radius (km) for `mass_kg`
    (`constants.GIANT_NEPTUNIAN_MASS_RADIUS`), with no scatter."""
    mass_earth = _kg_to_mearth(mass_kg)
    if giant_regime(mass_kg) == "neptunian":
        m0, r0, slope = constants.GIANT_NEPTUNIAN_MASS_RADIUS
    else:
        m0, r0, slope = constants.GIANT_JOVIAN_MASS_RADIUS
    return r0 * (mass_earth / m0) ** slope * constants.EARTH_RADIUS_KM


def _sample_giant_mass_kg(low_kg, high_kg):
    """A giant's mass in [low_kg, high_kg], from dN/dlogM ~
    M^-GIANT_MASS_FUNCTION_SLOPE."""
    return sample_power_law(low_kg, high_kg, tuning.GIANT_MASS_FUNCTION_SLOPE)


def sample_giant_radius_km(mass_kg):
    """`giant_radius_km(mass_kg)` with the relation's own scatter
    (`GIANT_RADIUS_SCATTER`, cut at 3 sigma). Bound giants
    (`_sample_giant_radius_km`) and rogue gas giants (GEN.60) both draw
    their radius here."""
    sigma = constants.GIANT_RADIUS_SCATTER[giant_regime(mass_kg)]
    factor = 1 + max(-3.0, min(3.0, draw.gauss(0.0, 1.0))) * sigma
    return giant_radius_km(mass_kg) * factor


def _giant_mass_for_radius_kg(radius_km, star_mass_kg):
    """
    A mass (kg) for a giant whose radius was fixed first (a caller's radius): the Neptunian branch's inverse below the transition
    radius, never under the core-critical mass, else a draw from the runaway mass function (the Jovian branch is nearly flat, so
    radius barely constrains mass there), never under the transition mass.
    """
    transition_kg = _mearth_to_kg(GIANT_TRANSITION_MASS_EARTH)
    radius_earth = radius_km / constants.EARTH_RADIUS_KM
    if radius_earth < GIANT_TRANSITION_RADIUS_EARTH:
        m0, r0, slope = constants.GIANT_NEPTUNIAN_MASS_RADIUS
        mass_kg = m0 * (radius_earth / r0) ** (1 / slope) * constants.EARTH_MASS_TO_KG
        return max(mass_kg, _mearth_to_kg(tuning.GIANT_CORE_CRITICAL_MASS_EARTH))
    return max(transition_kg, planet_mass.runaway_mass_kg(star_mass_kg))


def apply_class(planet):
    """
    Sets the class of a finished planet or moon from its stored conditions (GEN.689, `planet_classes.classify`; the composition text is
    written from the stored fractions, GEN.430). Called when the body is built (so a request can test it) and again once the system is
    final (the seasonal swing of its starlight, `flux_variation`, is known only then). The words that describe a class are read at view
    time (`physics/planet_class_text.py`); a body keeps none (GEN.1212).
    """
    planet.planet_class = None
    planet.class_variant = planet.class_match = planet.runner_up_class = planet.runner_up_match = None
    planet_type.fill(planet)            # GEN.1219: terrestrial, transitional, super-puff or gas, from the envelope share (GEN.1226) and the super-puff gate (GEN.1343)
    body_strength.generate(planet)      # GEN.1159: the roundness ratio the class reads (a giant body type has none)
    detail = planet_classes.classify_detail(planet)
    if detail is None:
        return      # GEN.1160: a body whose roundness ratio is under 1 is an Asteroid, which is an object type and has no class letter
    planet.planet_class = detail.letter
    # GEN.1207: the best-fit class, its match, the runner-up and the digit of the class code (0, or 9 when the gravity or the air amount is the limit)
    planet.class_variant, planet.class_match = detail.digit, detail.match
    planet.runner_up_class, planet.runner_up_match = detail.runner_up, detail.runner_up_match


@finite_domain()
def calculate_orbital_period_years(distance_au, primary_mass_kg):
    """
    Kepler's third law: T(years) = sqrt(a(AU)^3 / M_primary(Msun)).

    "Primary" is whatever body this orbit is actually around -- the host
    star for an ordinary planet, but the parent planet for a moon (see
    `Planet.__init__`'s `primary_mass_kg` parameter). A single shared
    helper so this formula is computed the same way wherever an orbital
    distance is set or changed -- also `StarSystem.validate_system`, which
    adjusts `distance` after generation to resolve overlapping orbits and
    must keep `period` in sync with it.

    Args:
        distance_au (float): Orbital distance from the primary, in AU.
        primary_mass_kg (float): The primary's mass, in kg.

    Returns:
        float: Orbital period in years.

    Raises:
        ValueError: If `distance_au` or `primary_mass_kg` isn't positive --
            every real caller always has a positive orbital distance and
            primary mass, so this guards against a `ZeroDivisionError`/
            `math domain error` escaping from a coding mistake upstream
            (e.g. an orbital-overlap correction pushing `distance` to zero
            or negative) with the same clean-`ValueError` treatment
            `calculate_surface_gravity` already gives a non-positive
            gravity, instead of a confusing low-level exception.
    """
    if distance_au <= 0 or primary_mass_kg <= 0:
        raise ValueError(
            f"calculate_orbital_period_years: distance_au ({distance_au}) and primary_mass_kg "
            f"({primary_mass_kg}) must both be positive."
        )
    primary_mass_sol = _kg_to_msun(primary_mass_kg)
    return math.sqrt(distance_au ** 3 / primary_mass_sol)


def formation_distance_au(planet):
    """Where the body formed, in AU: its own distance, or its parent's for a moon (`Planet.formation_distance_au`)."""
    distance = getattr(planet, "formation_distance_au", None)
    return planet.distance if distance is None else distance


def formed_beyond_snow_line(planet):
    """Whether `planet` formed at or beyond its star's snow line (False for a star with none)."""
    line = getattr(getattr(planet, "star", None), "snow_line_au", None)
    return line is not None and formation_distance_au(planet) >= line


def generate_planet_properties(planet, zone_override=None):
    """
    Generates a planet's physical and orbital properties: zone, kind, mass, radius, density, composition, and the interior
    (GEN.646, no class letter is read).

    A planet draws its kind and mass from the formation law at its distance (`planet_kind.draw_formation`); a moon's mass is
    given by its parent (`generate_moons`, the moon size law); a caller may give a mass, a radius or both (kept). A terrestrial
    body then draws its bulk composition from its star and its origin (`bulk_draw.draw_terrestrial`) and a giant its heavy
    elements and radius (`giant_draw.draw_giant`). The air, the water and the class follow (`calculate_atmospheric_conditions`,
    `apply_class`).

    Args:
        planet (Planet): The planet to generate properties for. Any of `planet.radius`/`planet.mass` already set are kept; unset
                         ones are filled in.
        zone_override (str, optional): A character ('h', 'c', 'e') to manually set the planet's zone, overriding the
                                       calculation based on distance.
    """
    # Determine the planet's zone (hot, ecosphere, or cold)
    zone = zone_for(planet.habitable_zone, planet.distance, planet.star)

    if zone_override and zone_override.lower() in "hce":
        zone = zone_override.lower()
    planet.zone = zone

    radius_given = planet.radius is not None
    mass_given = planet.mass is not None
    beyond = formed_beyond_snow_line(planet)
    _check_given_bounds(planet.mass if mass_given else None, planet.radius if radius_given else None)

    if planet.is_moon:
        kind = "terrestrial"
        if not mass_given and not radius_given:
            raise ValueError("A moon is built from a mass or a radius its parent gives (generate_moons).")
    elif mass_given or radius_given:
        kind = planet_kind.kind_of_given(planet.mass if mass_given else None, planet.radius if radius_given else None, beyond)
    else:
        formation = planet_kind.draw_formation(planet.distance, planet.star.mass, planet.star.snow_line_au,
                                               _star_fe_h(planet.star), bulk_draw.metal_ratio_to_solar(planet.star),
                                               gas_factor=disk_truncation.gas_survival_factor(planet.star))
        kind = formation.kind
        planet.mass = max(formation.mass_kg, tuning.PLANET_MIN_MASS_KG)
        mass_given = True
        log.choice("Planet kind", kind, f"formation at {planet.distance:.4g} AU gave {planet.mass:.4g} kg "
                   f"(core {formation.core_kg:.4g} kg, runaway {formation.runaway}, beyond the snow line {formation.beyond_snow_line})")

    if mass_given and radius_given:
        _check_given_density(planet.mass, planet.radius, kind)

    planet.body_type = "g" if kind == "giant" else "t"
    rings.clear(planet)  # a new draw draws its ring again (GEN.732)
    planet.oddities = ()  # ... and its oddities (GEN.1047)
    planet._venus_spin = None
    planet._resonance_e = None

    if uses_giant_mass_radius(planet):
        apply_giant(planet, radius_given, mass_given)
    else:
        # The composition comes first and the radius, density and gravity follow from it (GEN.428, GEN.429; `bulk_draw`).
        bulk_draw.draw_terrestrial(planet, radius_given, mass_given)
        giant_interior.clear(planet)

    oddities.apply(planet, levers=("L3",))     # GEN.1047: the bulk inputs (core, water, envelope) are settled; a row of lever L3 acts here

    # The air is made later, with the temperature it gives (`calculate_atmospheric_conditions`, GEN.411 and GEN.422): whether the
    # body holds one is the solver's call. Clear what a previous draw left behind -- `reconcile_zone` regenerates a body
    # in place -- and the class, which is derived once the body is finished (`apply_class`).
    planet._bare_core = None             # a new draw has no envelope yet (`apply_envelope`, GEN.1323)
    planet._captured_envelope = 0.0
    planet.atmosphere = None
    planet.atm_density = None
    planet.atm_molar_density = None
    planet.scale_height = None
    planet.planet_class = None

    planet.volume = (4 / 3) * math.pi * planet.radius ** 3
    # GEN.459: the potassium, thorium, uranium and aluminium-26 the interior heat and the ground's dose both read (a stored draw of the
    # body, so the later recalculation of its air and water, which moves nothing in its rock, does not redraw it).
    # GEN.428: Th and U go with the rock's Al2O3 (part E 2.8b), the body's over the bulk silicate Earth's.
    radiogenic.generate_content(planet, refractory_enrichment=radiogenic.rock_of(planet).aluminium / radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
    # GEN.699: the rheology of the mantle that the heat it was born with leaves through (the heat itself: `internal_heat_flux_w_m2`).
    interior_heat.generate_interior(planet)

    update_hill_sphere(planet)


ECCENTRICITY_SETTLE_TOLERANCE = 1e-4
"""float: A planet or captured moon whose eccentricity at the star's final age is within this of the one drawn keeps its orbit (`settle_orbit`);
below it the orbit cannot be told apart. A stand-in, registry row OE14."""


def settle_orbit(planet):
    """
    Takes the tides of a planet's (or captured moon's) orbit again at the star's final age (GEN.985): `generate_orbit_shape` wore the
    birth eccentricity down for the age the star had when the planet was drawn, and the age rises when the planets are in. The birth
    eccentricity is kept with the body for the generation (`_orbit_birth`, not stored), and the circularisation time is read again from the
    body's radius as it is now. Ages only rise, so the eccentricity only falls (it is never raised above what the body has, which the spacing of the orbits may have cut), and a lower one clears the star and the neighbours better.
    The position and velocity are made again on the new ellipse (`update_orbital_position`, the body keeps its place in its plane) and
    a planet's moons are told (`host_eccentricity`). The air, the tide's heat and the spin of a rocky body are read again by the
    re-finish of GEN.992, not here; the star's XUV flux, exposure and flares at the new mean light are set here, so the stored flux follows the stored orbit.

    Returns:
        bool: Whether the eccentricity moved by more than `ECCENTRICITY_SETTLE_TOLERANCE`.
    """
    birth = getattr(planet, "_orbit_birth", None)
    if birth is None:
        return False
    # never above what the planet has: the spacing of the orbits may have cut it (`validation`) and ages only rise
    eccentricity = min(planet.eccentricity, orbit_elements.settled_eccentricity(planet, birth, _gyr_to_yr(planet.star.age), planet.primary_mass_kg, planet.distance))
    if planet.eccentricity - eccentricity <= ECCENTRICITY_SETTLE_TOLERANCE:
        return False
    planet.eccentricity = eccentricity
    for moon in getattr(planet, "moons", ()):
        moon.host_eccentricity = eccentricity
    update_orbital_position(planet)
    # the star's XUV and flares at the mean light of the new orbit (no draws), so the stored flux follows the stored ellipse; the air waits for the re-finish
    light_au = orbit_elements.flux_distance_au(planet.distance, eccentricity)
    for body in (planet, *getattr(planet, "moons", ())):
        if getattr(body, "_star_distance_au", None):
            magnetism.update_exposure(body, light_au)
    return True


AGE_SETTLE_TOLERANCE = 0.02
"""float: A rocky planet or moon made when the star's age was within this share of its final age keeps its air, class and dose
(`refinish_rocky`); the age the decay clock and the heat history read is then within 2 percent of the right one. A stand-in, registry
row OE15."""


def age_moved(body):
    """Whether the star's age has risen by more than `AGE_SETTLE_TOLERANCE` of the age `body` was made at."""
    drawn = getattr(body, "_age_drawn", None)
    return drawn is not None and body.star.age - drawn > AGE_SETTLE_TOLERANCE * max(drawn, 0.1)


def orbit_moved(body):
    """Whether the eccentricity the body's light goes by (a moon's host's) has moved by more than `ECCENTRICITY_SETTLE_TOLERANCE` since its air was made (GEN.1402): `settle_orbit` wore the
    orbit down after the air was drawn, so the air, temperature, dose and tide's heat were made on another orbit. A body with no air made on an orbit is not moved."""
    made_on = getattr(body, "_air_eccentricity", None)
    if made_on is None:
        return False
    now = body.host_eccentricity if body.is_moon else body.eccentricity
    return abs(float(now) - made_on) > ECCENTRICITY_SETTLE_TOLERANCE


def needs_remake(body):
    """Whether a rocky body's end stages are to be made again: the star's age rose past `AGE_SETTLE_TOLERANCE` (`age_moved`) or the orbit settled past `ECCENTRICITY_SETTLE_TOLERANCE` (`orbit_moved`)."""
    return age_moved(body) or orbit_moved(body)


def remake_end_stages(body, planet):
    """
    Makes a rocky planet's or moon's air, water, surface, class and spin again from a stream of the body's own (`refinish_rocky`). `planet` is the body
    itself for a planet and the parent for a moon.
    """
    stream = draw.Stream(f"refinish-rocky:{float(body.mass)!r}:{float(body.distance)!r}:{float(planet.distance)!r}:"
                         f"{float(body.star.mass)!r}:{float(body.star.age)!r}:{float(body.radius)!r}:{int(body.is_moon)}")
    with draw.bound(stream):
        if body.is_moon:
            calculate_atmospheric_conditions(body, planet.distance, planet.eccentricity)
            apply_class(body)
            generate_spin(body, planet.mass, birth=getattr(body, "_birth_spin", None))
        else:
            calculate_atmospheric_conditions(body)
            apply_class(body)
            generate_spin(body, body.primary_mass_kg, birth=getattr(body, "_birth_spin", None))


def refinish_rocky(planet):
    """
    Makes a rocky planet's end stages again at the star's final age, and its moons' (GEN.992, found by GEN.986): the air, the water and the
    surface (`calculate_atmospheric_conditions`), the class, the spin at that age (its own birth spin worn by the tide again, `generate_spin` with `birth`, GEN.961) and what ends it (the field, the dose, the colour and the
    habitability score; `generate_spin`), the same sequence a moved planet gets (`validation.reconcile_moved_planet`), but from a stream of
    the body's own (mass, distance, star and the age it settles to), so the run's stream does not move. A body keeps its mass, radius,
    composition and orbit. The forced requests are judged on the classes this makes (`StarSystem`), and the moons' induced fields, rings and
    life follow in later passes. A planet or moon whose star's age rose by under `AGE_SETTLE_TOLERANCE` of the age it was made at, and whose eccentricity moved by no more than `ECCENTRICITY_SETTLE_TOLERANCE` since its air was made, is left (GEN.1402).

    Args:
        planet (Planet): A rocky planet (`body_type` "t", not a moon) with its moons; a giant's moons are made again too.

    Returns:
        int: How many bodies were made again (the planet and its moons).
    """
    made = 0
    # a rocky planet on a moved age or orbit; a giant only on a moved orbit (its temperature reads the light of the orbit, GEN.1402; its age does not move its end stages)
    with_planet = not planet.is_moon and ((planet.body_type == "t" and needs_remake(planet)) or (planet.body_type == "g" and orbit_moved(planet)))

    def bodies():
        if with_planet:
            yield planet
        yield from list(planet.moons)      # read after the planet is made again: its moons' history (GEN.959) can replace a moon with a new body

    for body in bodies():
        if not needs_remake(body):
            continue
        remake_end_stages(body, planet)
        body._age_drawn = body.star.age
        made += 1
    if made:  # GEN.767: the spin and period changed, the moons' day and night shares follow (no draws)
        for moon in planet.moons:
            daylight.generate(moon, planet)
    if made:
        # GEN.1034: the sea and the spin are new, and the tide worked out from the old ones follows them (no draws)
        ocean_tides.generate(planet)
        for moon in planet.moons:
            ocean_tides.generate(moon, planet)
    return made


ENVELOPE_SETTLE_TOLERANCE = 0.02
"""float: A planet whose envelope at the star's final age and its final orbit is within this share of the one it was drawn with keeps it and its air, class and
spin (`settle_envelope`); the radius moves by less than a hundredth. A stand-in, registry row EC6."""


def apply_envelope(planet):
    """
    Gives a rocky planet the hydrogen-helium envelope it captured and kept (GEN.1323, `envelope_capture`), from the star's age and the planet's
    orbit as they are now: it draws the envelope the planet started with from a stream of its own (the run's stream does not move; the key is
    the planet's mass, formation distance and star, not the age or the orbit, so every pass draws the same start), works out what the loss to
    the star's X-rays and heat leaves, and when something is, stores the envelope as a fraction of the mass (the solid fractions give way to
    it), sets the radius of the core under the envelope (the law of Lopez and Fortney, as the valley survey does) and the volume, density and
    gravity that follow, and drops the moons that no longer fit. A planet that lost it all, a planet out of the mass range, a moon and a giant
    are left as they were. Called again it starts from the planet as it was drawn (`_bare_core`), so it never compounds, and a planet whose
    envelope is within `ENVELOPE_SETTLE_TOLERANCE` of the one it has is left. It does not make the air, class or spin again: the caller does
    (`Planet`, before the air is made the first time; `settle_envelope`, after).

    Returns:
        bool: Whether the planet's envelope (and so its radius) changed.
    """
    if planet.body_type != "t" or planet.is_moon:
        return False
    mass_earth = _kg_to_mearth(planet.mass)
    if not envelope_capture.in_range(mass_earth):
        return False
    if getattr(planet, "_late_core", False):
        start = 0.0                     # a core that finished after the gas disk went captured nothing
    else:
        stream = draw.Stream(f"envelope-start:{float(planet.mass)!r}:{float(formation_distance_au(planet))!r}:{float(planet.star.mass)!r}")
        with draw.bound(stream):
            start = envelope_capture.draw_start(mass_earth)
            puff = envelope_capture.draw_puff_start(mass_earth, envelope_capture.period_days(planet))     # GEN.1343 (W15): one sub-Neptune core in 20 starts as a super-puff
            start = puff or start
    result = envelope_capture.kept_with_tide(planet, start, tides=True)       # GEN.1347: the obliquity tide of the planet's own tilt, once it has one (the settle pass)
    held = getattr(planet, "_captured_envelope", 0.0)
    if result is None and not held:
        return False                    # a bare core, as drawn: nothing is touched
    if result is not None and held and abs(result[0] - held) <= ENVELOPE_SETTLE_TOLERANCE * held:
        inflated = getattr(planet, "_tidal_radius_km", None)
        new_radius = result[1] * constants.EARTH_RADIUS_KM
        if inflated is None and new_radius <= planet.radius * (1.0 + ENVELOPE_SETTLE_TOLERANCE):
            return False
        if inflated is not None and abs(new_radius - inflated) <= ENVELOPE_SETTLE_TOLERANCE * inflated:
            return False
    bare = getattr(planet, "_bare_core", None)
    if bare is not None:                # it had an envelope: back to the drawn core first
        for name, value in bare["fractions"].items():
            setattr(planet, name, value)
        planet.radius = bare["radius"]
        planet._bare_core = None
        planet._captured_envelope = 0.0
        if hasattr(planet, "_formation_fractions"):
            del planet._formation_fractions
    if result is None:                  # it kept an envelope before and keeps none now: the drawn core again
        planet.envelope_helium_fraction = planet.envelope_z_fraction = None
        planet._tidal_radius_km = None
    else:
        envelope_share.restore(planet)  # the drawn fractions, before the envelope takes its share
        fraction, radius_earth, tide_ratio = result
        planet._bare_core = {"radius": planet.radius, "fractions": {name: getattr(planet, name) for name in (*envelope_share.SOLIDS, "envelope_mass_fraction")}}
        solids = sum(getattr(planet, name, None) or 0.0 for name in envelope_share.SOLIDS)
        keep = (1.0 - fraction) / solids if solids > 0.0 else 1.0       # the solids as drawn may already have made room for a hydrogen air's share: they give the envelope its share of the whole
        for name in envelope_share.SOLIDS:
            value = getattr(planet, name, None)
            if value is not None:
                setattr(planet, name, value * keep)
        planet.envelope_mass_fraction = fraction
        planet._captured_envelope = fraction
        envelope_share.set_gas_detail(planet)       # the captured gas is the disk's: the star's helium and heavy elements
        planet.radius = radius_earth * constants.EARTH_RADIUS_KM
        planet._tidal_radius_km = planet.radius if tide_ratio > 1.0 + 1e-9 else None     # the radius the obliquity tide gave it, or None for a planet it did not inflate (GEN.1480: the ratio itself; the loss follows the tide, so an uninflated envelope is not the one to compare with)
    planet.volume = (4 / 3) * math.pi * planet.radius ** 3
    planet.density = planet.mass / (planet.volume * 1.0e9) / 1000.0     # g/cm3: kg per km3 over 1e9 m3 per km3, then to g/cm3
    calculate_surface_gravity(planet)
    planet.composition = composition_text.text(planet)
    if getattr(planet, "moons", None) and drop_unstable_moons(planet):
        planet.reflex_offset_x, planet.reflex_offset_y, planet.reflex_offset_z = calculate_reflex_offset(
            planet.mass, [(m.mass, m.position_x, m.position_y, m.position_z) for m in planet.moons])
    return True


def settle_envelope(planet):
    """
    The envelope of a rocky planet once the star's age and the planet's orbit are final (GEN.1323; the loss needs both, so the planet was given
    the one the provisional age and orbit leave when it was drawn, `apply_envelope`). When the final one differs, the radius, density and
    gravity follow and the air, class and spin are made again (the gravity and radius changed), like a giant's settled radius (`settle_giant`).

    Returns:
        bool: Whether the planet's envelope changed.
    """
    if not apply_envelope(planet):
        return False
    remake_end_stages(planet, planet)
    planet._age_drawn = planet.star.age
    return True


def settle_giant(planet):
    """
    Settles a giant sized from the radius tables once the star's age and the orbit are final (GEN.961, `giant_draw.settle`): when its
    radius moves, the volume and the gravity follow, and the air is made again from a stream of its own (the gravity sets its scale
    height; the run's stream is not touched), then the class and the spin (which ends with the field, the dose, the colour and the habitability score).

    Returns:
        bool: Whether the giant's radius moved.
    """
    moved = giant_draw.settle(planet)
    stale = bool(getattr(planet, "_air_stale", False))       # the secular exchange moved the orbit and with it the mean light (GEN.1505)
    planet._air_stale = False
    if not moved and not stale:
        # The radius stayed (the tables moved it by under the tolerance, or the giant was sized another way), but the star's age rose after the flux was
        # stored and the giant's own heat falls with the age: the flux it radiates is made again at the final age.
        if not planet.is_moon and getattr(planet, "energy_flux_w_m2", None) is not None:
            planet.energy_flux_w_m2 = habitability_world.giant_energy_flux_w_m2(planet)
        return False
    if moved:
        planet.volume = (4 / 3) * math.pi * planet.radius ** 3
        calculate_surface_gravity(planet)
    stream = draw.Stream(f"settle-giant-air:{float(planet.mass)!r}:{float(planet.distance)!r}:{float(planet.star.mass)!r}:"
                         f"{float(planet.star.age)!r}:{float(planet.radius)!r}")
    with draw.bound(stream):
        calculate_atmospheric_conditions(planet)
        apply_class(planet)
        generate_spin(planet, planet.primary_mass_kg, birth=getattr(planet, "_birth_spin", None))  # the birth spin again worn to the final age, then the field, the dose, the colour and the score
    if drop_unstable_moons(planet):  # a smaller giant holds a smaller moon and a smaller range of orbits
        planet.reflex_offset_x, planet.reflex_offset_y, planet.reflex_offset_z = calculate_reflex_offset(
            planet.mass, [(m.mass, m.position_x, m.position_y, m.position_z) for m in planet.moons])
    return moved


def _check_given_bounds(mass_kg, radius_km):
    """The mass and radius a caller gave must be inside `tuning.GIVEN_MASS_KG_BOUNDS` and `tuning.GIVEN_RADIUS_KM_BOUNDS`.

    Raises:
        ValueError: One is outside its bounds (or not a number).
    """
    for name, value, unit, (low, high) in (("mass", mass_kg, "kg", tuning.GIVEN_MASS_KG_BOUNDS),
                                           ("radius", radius_km, "km", tuning.GIVEN_RADIUS_KM_BOUNDS)):
        if value is not None and not low <= value <= high:
            raise ValueError(f"The given {name} ({value:.4g} {unit}) is outside {low:.3g} to {high:.3g} {unit}.")


def _check_given_density(mass_kg, radius_km, kind):
    """A mass and a radius a caller gave together must make a body of a density that kind of body can have.

    Raises:
        ValueError: The mean density is outside the bounds.
    """
    density = mass_kg / ((4 / 3) * math.pi * (radius_km * constants.KM_TO_M_FACTOR) ** 3) / 1000.0
    low, high = tuning.GIANT_DENSITY_BOUNDS_G_CM3 if kind == "giant" else bulk_draw.DENSITY_BOUNDS_G_CM3
    if not low <= density <= high:
        raise ValueError(f"A {kind} body of {mass_kg:.4g} kg and {radius_km:.4g} km has a density of {density:.4g} g/cm3; "
                         f"a {kind} body has {low:g} to {high:g}.")


def _star_fe_h(star):
    """The star's [Fe/H] for the formation law (0 for a star without one, a remnant)."""
    fe_h = getattr(star, "metallicity_fe_h", None)
    return 0.0 if fe_h is None else fe_h


def apply_giant(planet, radius_given, mass_given):
    """
    Sets a giant's mass, radius, density and composition: the mass first (`_apply_giant_mass`), then its heavy elements, radius
    and interior (`giant_draw.draw_giant`, GEN.460 and GEN.642), with the anchored mass-radius relation (GEN.34) for a body the
    radius tables do not hold.
    """
    _apply_giant_mass(planet, radius_given, mass_given)
    giant_draw.draw_giant(planet, radius_given, sample_giant_radius_km)


def _apply_giant_mass(planet, radius_given, mass_given):
    """
    Sets a giant's mass: the formation law's or the caller's, kept; with only a radius, the anchored relation's mass for it
    (`_giant_mass_for_radius_kg`).
    """
    if not mass_given:
        planet.mass = _giant_mass_for_radius_kg(planet.radius, planet.star.mass)


def update_hill_sphere(planet):
    """
    Sets `planet.hill_radius` (km) and `planet.min_orbit_distance` (5 Hill
    radii, in AU: the clearance the next body out must keep) from the
    planet's current distance and mass about `planet.primary_mass_kg` (the
    star for a planet, the parent planet for a moon). Call again whenever
    either changes, e.g. after `StarSystem.validate_system` moves the planet.
    """
    distance_m = _au_to_m(planet.distance)
    # A moon's Hill sphere is within its planet's pull, not the star's (GEN.138).
    primary_mass_kg = getattr(planet, "primary_mass_kg", None) or planet.star.mass
    planet.hill_radius = _m_to_km(calculate_hill_sphere(distance_m, planet.mass, primary_mass_kg))
    planet.min_orbit_distance = _km_to_au(5 * planet.hill_radius)


def calculate_surface_gravity(planet):
    """
    Calculates the surface gravity of the planet in g's (multiples of Earth's gravity).

    This computes the surface gravity based on the planet's mass and
    radius. The result is normalized to Earth's gravity (g's). It includes
    special adjustments for certain planet classes to ensure realistic values.

    Args:
        planet (Planet): The planet to calculate surface gravity for. Its
                         `gravity` attribute is set in place.

    Raises:
        ValueError: If the computed gravity is zero or negative.
    """
    radius_meters = planet.radius * 1000
    surface_gravity = relations.surface_gravity(planet.mass, radius_meters)
    surface_gravity_g = surface_gravity / constants.EARTH_GRAVITY
    if surface_gravity_g <= 0:
        raise ValueError('Invalid value for gravity.')
    planet.gravity = surface_gravity_g


def calculate_atmospheric_conditions(planet, distance_override=None, host_eccentricity=None):
    """
    Makes the air of the planet or moon and the surface temperature it gives (GEN.411, GEN.422; `surface_air`).

    The Bond albedo is drawn from the range of the body's regime (`climate_draws.albedo_range`); the heat the interior has carried out (`interior_heat.update_heat`) says how much
    of the body's volatiles has reached the surface (`surface_pressure`); the water it holds is the reservoir its air's water vapor
    stands in contact with; and the mixture, the pressure with the freeze-out, the density, the scale height, the greenhouse and the
    temperature are solved together (`surface_air.solve`). A body that holds no air, and a body whose air is lost, have the temperature the
    starlight alone gives.

    Args:
        planet (Planet): The planet to calculate atmospheric conditions for. Its `surface_temperature`, `atmospheric_pressure`,
                         `atmosphere`, `atm_density`, `atm_molar_density`, `scale_height`, gases, hydrosphere and aerosols are set
                         in place.
        distance_override (float, optional): An override for the planet's distance, used for special cases like moons.
        host_eccentricity (float, optional): For a moon, its parent planet's eccentricity, which the light it gets goes by (GEN.727);
                                             kept on the moon (`Planet.host_eccentricity`) and left as it was when not given.
    """
    if not planet.is_moon:
        planet._air_replay = (draw.getstate(), distance_override, host_eccentricity, None)    # GEN.1493: kept, not stored; where this air's draws began, for `remake_air` once the system's belts are known
    envelope_share.restore(planet)      # GEN.1226: the drawn fractions, before the air is solved again (the air share is put in again at the end)
    if host_eccentricity is not None:
        planet.host_eccentricity = host_eccentricity
    distance = float(distance_override) if distance_override is not None else float(planet.distance)
    # GEN.727: the light of an eccentric orbit, averaged over it, is that of a circle of radius a (1 - e^2)^(1/4)
    distance = orbit_elements.flux_distance_au(distance, planet.host_eccentricity if planet.is_moon else planet.eccentricity)
    planet._air_eccentricity = float(planet.host_eccentricity if planet.is_moon else planet.eccentricity)      # kept, not stored: the orbit this air was made on (`orbit_moved`, GEN.1402)
    planet._air_light_au = distance     # kept, not stored: the light this air, temperature, dose and spectrum were made on (a later settled orbit moves `_star_distance_au` and leaves them)
    planet.light_factor = eclipses.light_factor_of(planet)       # GEN.1484: the share of its light a circumbinary planet keeps after its two stars eclipse each other; None for every other body
    solar_output_at_orbit = radiation_field.starlight_w_m2(planet.star, distance) * (planet.light_factor or 1.0)     # GEN.1441: the host's light and every other star's, by position
    # The range of the body's regime (a frozen or a wet world reflects more than dry rock): it reads only what the body is made of and
    # the starlight, so it can be drawn before the air it will be solved with.
    albedo = draw.uniform(*climate_draws.albedo_range(planet.body_type, planet.water_mass_fraction, solar_output_at_orbit,
                                                      getattr(planet.star, "temperature", None)))
    # Floored at the cosmic microwave background: starlight alone can put a
    # body around a very dim star (or far out) below ~2.7 K, but nothing in
    # space is colder than the background it sits in.
    surface_temperature_no_atmosphere = climate_draws.equilibrium_temperature_k(albedo, solar_output_at_orbit)

    # GEN.699: the heat the interior has carried out (no draws), which the degassing of the air follows; GEN.88: the surface water, drawn
    # before the air because the ocean is the reservoir the air's water vapor and its carbonate cycle stand in contact with.
    magnetism.update_exposure(planet, distance)     # GEN.86: the star's XUV and flares at this distance (no draws); the water loss reads it
    heat = interior_heat.update_heat(planet)
    giant = not planet_type.has_ground(planet)      # GEN.1334: a giant, or a rocky planet under a deep envelope, has no surface; the air is solved at its 1 bar level
    own_heat_w_m2 = (getattr(planet, "internal_heat_flux_w_m2", None) or 0.0) + (getattr(planet, "tidal_flux_w_m2", None) or 0.0)
    if getattr(planet.star, "starless", False) and not giant and own_heat_w_m2 > 0.0:
        # GEN.1111: a body with no star (an ejected one) lives on its own heat, which is a flux through its surface; under a star it is below a part in
        # a thousand of the absorbed flux for an ordinary rocky planet and is left out, so a lit planet is as it was (a giant's is added by its setting)
        surface_temperature_no_atmosphere = climate_draws.equilibrium_temperature_k(albedo, solar_output_at_orbit, internal_flux_w_m2=own_heat_w_m2)
    surface_w = None if giant else hydrosphere.draw_surface_water(planet)
    planet._surface_w = surface_w                   # GEN.839: kept for `refresh_at_final_age`, not stored
    water_lost = (not giant) and hydrosphere.water_loss_kg(planet) > 0.0

    planet.surface_temperature = surface_temperature_no_atmosphere
    planet.atmospheric_pressure = 0.0
    climate_draws.record(planet, albedo, 0.0)       # the mixture subject reads the albedo
    setting = surface_air.setting_of(planet, albedo=albedo, starlight_w_m2=solar_output_at_orbit, star_distance_au=distance,
                                     heat=heat, surface_w=surface_w, water_lost=water_lost)
    start_k = surface_temperature_no_atmosphere
    if _past_runaway_with_water(planet, solar_output_at_orbit, surface_w):
        start_k = max(start_k, tuning.RUNAWAY_START_TEMPERATURE_K)    # GEN.1238: the steam branch, not the cool one (see `_past_runaway_with_water`)
    air = surface_air.solve(setting, start_k)
    if not planet.is_moon and not giant and air is not None and air.partials_pa:       # GEN.1493: the late accretion this air was made with and what it is small or large against, for `remake_air`
        planet._air_replay = (*planet._air_replay[:3], (setting, air.late_state, air.late, air.degassed))
    planet._retention_cap_pa = None if air is None else air.retention_cap_pa     # GEN.1340: kept, not stored; the life stage keeps its oxygen inside it
    air = _air_oddities(planet, setting, air, start_k, solar_output_at_orbit, distance, albedo, giant)
    planet._co2_sink = None if air is None else air.co2_sink                      # GEN.1309 A: kept, not stored; the life stage balances the CO2 again at the temperature the oxygen rise leaves
    if air is None or not air.partials_pa:
        if setting.albedo_law is not None:
            # GEN.972: with no air the ice the cold leaves is the only feedback; the temperature and the albedo that keep each other
            planet.surface_temperature, albedo = surface_air.airless_state(setting.albedo_law, solar_output_at_orbit, surface_temperature_no_atmosphere)
            climate_draws.record(planet, albedo, 0.0)
        planet.atm_density = planet.atm_molar_density = planet.scale_height = None
        atmosphere.store_air(planet, {}, None if giant else _mantle_redox(planet, air, distance, water_lost))
    else:
        planet.surface_temperature = air.temperature_k
        planet.atmospheric_pressure = air.pressure_pa
        planet.atm_molar_density = air.state.mean_molar_mass_kg_mol
        planet.atm_density = air.state.density_kg_m3
        planet.scale_height = _m_to_km(air.state.scale_height_m)
        factor = air.greenhouse_factor
        albedo_kept = albedo if air.albedo is None else air.albedo      # a body with ice has the albedo of its last step (GEN.972)
        # GEN.1236: the fixed point holds the surface at `surface_air.MAX_SURFACE_TEMPERATURE_K` (a stand-in guard, the equilibrium of the gas does
        # not converge far above it), but the warming it keeps is the uncapped one, so a body at the cap kept a greenhouse factor that gave 2370 K
        # or 2600 K; its stored G is taken from the stored temperature too, so the three numbers agree.
        capped = air.warming is not None and air.warming.surface_temperature_k > air.temperature_k * (1.0 + CAP_TOLERANCE)
        if giant or capped or air.knee:        # (GEN.1344: the knee of an unsettled air is not the temperature its last step found)
            # GEN.412 reads the temperature as starlight alone times (1 + G)^(1/4); a giant also radiates its own heat, which the air's factor
            # leaves out of T_e, so the stored G is taken against the starlight alone (GEN.1030) and the pair reproduces the temperature
            factor = climate_draws.greenhouse_from_temperature(air.temperature_k, albedo_kept, solar_output_at_orbit)
            factor = air.greenhouse_factor if factor is None else max(factor, 0.0)
        climate_draws.record(planet, albedo_kept, factor)     # GEN.412: the two numbers the temperature is made of
        atmosphere.store_air(planet, air.partials_pa, None if giant else _mantle_redox(planet, air, distance, water_lost))

    # GEN.88: water, land, ice and the ocean's chemistry.
    hydrosphere.generate_hydrosphere(planet, 0.0 if heat is None else heat.flux_w_m2, surface_w)
    # GEN.721: stagnant or mobile lid, which the magnetic field, the heat and the ground below all read (the hydrosphere above
    # read the heat of a stagnant lid, the regime not being drawn yet; it needs the hydrosphere). The heat is made again with it.
    interior_heat.generate_lid(planet)
    interior_heat.update_heat(planet)
    # GEN.441: the minerals, ices and ocean on the ground, from the bulk fractions, the temperature, the air and the water.
    surface_composition.generate_surface(planet)
    phosphorus_budget.generate(planet)         # GEN.1367: the crust's P2O5 and hosts, the ocean's phosphate
    # GEN.1380: the crust's thickness, craton and sulfur, and the magmatic ore types the body admits (own streams, no draw of the run's)
    ore_deposits.generate(planet)
    # GEN.1377, GEN.1420: the isotope record of the elements the body holds (own streams, no draw of the run's)
    isotope_record.generate(planet)
    # GEN.414: cloud, haze and dust depths, on the body's own stream (no draw of the run's moves).
    aerosols.generate_aerosols(planet)
    # GEN.1226: the true share of the hydrogen-rich air, in the envelope share of every rocky body (after the water and the ground, which read the drawn fractions)
    envelope_share.fill(planet)


def _air_oddities(planet, setting, air, start_k, starlight_w_m2, distance_au, albedo, giant):
    """
    The rows of the air lever (L7, GEN.1051: W19 the nightside collapse of a locked planet's air, W20 the other branch of a bistable climate) on the air just solved: they may replace it
    with the airless state (`None`) or with the air of the other branch (`physics/air_branches.py`). The rows are rolled again each time the air is made, from their own streams.
    """
    air_branches.clear(planet)
    if giant or planet.is_moon or air is None:
        return air
    planet._air_context = air_branches.Context(setting, air, start_k, starlight_w_m2, distance_au, albedo)
    oddities.drop(planet, ("L7",))
    oddities.apply(planet, levers=("L7",))
    choice = getattr(planet, air_branches.CHOICE_ATTRIBUTE, None)
    air_branches.clear(planet)
    if choice == "lost":
        return None
    return air if choice is None else choice


def _past_runaway_with_water(planet, starlight_w_m2, surface_w):
    """
    Whether the body holds surface water under starlight past its runaway greenhouse limit (GEN.1238): an effective flux (in the Sun's at the
    Earth, `starlight_w_m2` over `RUNAWAY_FLUX_REFERENCE_W_M2`) above `orbits.runaway_greenhouse_flux` of the star's temperature and the body's mass,
    the limit the water-loss history already follows (`runaway_history`, GEN.639). Past it no ocean can stand, but the air's fixed point has two
    answers there (a cool air with a few bar of vapor, whose greenhouse the table gives no limit to, and the steam air that holds the water),
    and which one it finds depends on where it starts, so a body past the limit starts hot and finds the steam one; a body below it, or with no
    water, starts from the airless temperature as before. No draw.
    """
    if not surface_w or planet.body_type == "g" or not getattr(planet.star, "temperature", None) or not planet.mass:
        return False
    effective = starlight_w_m2 / tuning.RUNAWAY_FLUX_REFERENCE_W_M2
    return effective > orbits.runaway_greenhouse_flux(planet.star.temperature, units.converter("kg", "Mearth")(planet.mass))


LATE_ACCRETION_GATE = 0.5
"""float: The air is not made again for a new late accretion when it moves the carbon and the nitrogen the late accretion left by less than this share of what the volcanoes gave. The share is how far the
inventory moves, not how far the air does (a carbonate sink holds the carbon of a wet world whatever is added): in 246 planets of 240 seeds 201 are identical to the full replay, 13 differ by more than 1 percent
(the worst by a branch of a bistable climate, 99 percent in pressure), and the pass costs 8.8 percent of generation against 18 for every planet; at gate 1.0 a 30-seed trial had a worst of 29 percent. [Stand-in LA3, tier Marginal educated guess from physics, tweak L3.]"""


def remake_air(planet, late_accreted_percent, gate=LATE_ACCRETION_GATE):
    """
    Makes the air of a planet again with the late accretion the system's own belts give it (GEN.1493, `generation/late_accretion.py`): `calculate_atmospheric_conditions` is run again on the
    draws it began with (`planet._air_replay`, taken at its last call), in a stream of its own, so every other draw of the air, the water and the ground is the one it was, the late accretion is
    `late_accreted_percent` of the planet's mass instead of the law's draw, and no draw of the run's stream moves. A planet whose air was never made here (a moon, a loaded one) or has no air is left
    as it is, and so is one whose new late accretion moves the carbon and the nitrogen it left by less than `gate` of what the volcanoes gave (it reaches the air only through them).

    Args:
        planet (Planet): The planet.
        late_accreted_percent (float): The late accretion, percent of the planet's mass.
        gate (float): The share of the degassed carbon or nitrogen the new late accretion moves it by, under which the air is left as it is; 0 makes it for every planet (the full replay the census compares with).

    Returns:
        bool: Whether the air was made again.
    """
    replay = getattr(planet, "_air_replay", None)
    if replay is None or planet.is_moon:
        return False
    state, distance_override, host_eccentricity, settled = replay
    planet._air_replay = None       # the replay is spent; the memory it held goes
    planet._late_accreted_percent = late_accreted_percent
    if settled is None:      # no air was made, so no late accretion reached one
        return False
    if gate > 0.0 and surface_air.late_accretion_moves(*settled, late_accreted_percent) < gate:
        return False
    stream = draw.Stream(0)
    stream.setstate(state)
    with draw.bound(stream):
        calculate_atmospheric_conditions(planet, distance_override, host_eccentricity)
    planet._air_replay = None       # the air made again kept a replay of its own; the memory it held goes
    _remake_what_reads_the_air(planet)
    return True


def _remake_what_reads_the_air(planet):
    """
    The class, the spin and what ends it (the field, the dose, the ultraviolet index, the colour, the equipment and habitability scores; `generate_spin`) made again on the air that
    `remake_air` just made, the way `remake_end_stages` follows the air (GEN.992), from a stream of the planet's own so the run's stream does not move, and the moons' day and night shares (no
    draws) after the planet's new spin. The planet keeps its birth spin, so a planet whose air is made again with the same pressure and clouds keeps its day.
    """
    stream = draw.Stream(f"late-accretion-end:{float(planet.mass)!r}:{float(planet.distance)!r}:{float(planet.star.mass)!r}:{float(planet.star.age)!r}:{float(planet.radius)!r}")
    with draw.bound(stream):
        apply_class(planet)
        generate_spin(planet, planet.primary_mass_kg, birth=getattr(planet, "_birth_spin", None))
    for moon in planet.moons:
        daylight.generate(moon, planet)


def refresh_at_final_age(body):
    """
    GEN.839: makes again, at the heat of the star's final age, what read the heat of the provisional one: the ice shell and the rest of
    the hydrosphere (`hydrosphere.generate_hydrosphere`), the ground with its volcanism gate (`surface_composition.generate_surface`),
    the ground dose when the ground changed (`radiation.update_ground_dose`), and the habitability score that reads them all. Called in the life pass, after
    `interior_heat.update_heat` has run again at the final age. The draws these make (the brine, the ocean chemistry, the ground and its
    gates) come from the body's own stream (`util.draw.Stream`, seeded from what is final about it), so no draw of the run's stream moves;
    a body with no water drawn (a giant, a belt, or one loaded from the database) is left as it is.
    """
    surface_w = getattr(body, "_surface_w", None)
    if body.body_type in ("g", "a") or surface_w is None:
        return
    seed = f"final-age:{body.is_moon}:{body.mass!r}:{body.radius!r}:{body.distance!r}:{getattr(body.star, 'age', None)!r}"
    flux = internal_heat_flux_w_m2(body)            # the heat at the star's age now (no draws; the life pass has run it already)
    ground = surface_composition.mix_of(body)
    with draw.bound(draw.Stream(seed)):
        hydrosphere.generate_hydrosphere(body, flux, surface_w)
        surface_composition.generate_surface(body)
        phosphorus_budget.generate(body)        # GEN.1367: the heat, the water and the ground moved
        ore_deposits.generate(body)         # GEN.1380: the heat, the water and the ground moved
        isotope_record.generate(body)       # GEN.1377: the air and the interior moved
        if surface_composition.mix_of(body) != ground:
            radiation.update_ground_dose(body)
    envelope_share.fill(body)       # GEN.1226: the air may have changed with the water and the ground; the share follows the final air
    habitability_world.generate(body, flux)


def _mantle_redox(planet, air, distance, water_lost):
    """The mantle redox (log10 above iron-wuestite) of a rocky body: the air mixture's when it has one, else (no air, or an envelope's
    mixture, which has none) the recipe's redox rows at the body's own mass and age (`surface_air.mantle_redox_diw`)."""
    if air is not None and air.mixture.redox_diw is not None:
        return air.mixture.redox_diw
    return surface_air.mantle_redox_diw(planet, distance, water_lost)


def internal_heat_flux_w_m2(planet):
    """A rocky body's own heat flow, W/m^2, at its star's age now: radiogenic plus primordial heat (`interior_heat.update_heat`, GEN.699),
    stored with the body (its Urey ratio is an output), and with it the flow of its tide from its orbit (GEN.737, `tidal_flux_w_m2`, not in
    this figure: the habitability energy term adds the two). 0 for a body with no stored content (a giant, or a row from before GEN.699)."""
    found = interior_heat.update_heat(planet)
    return found.own_flux_w_m2 if found is not None else 0.0


def generate_orbit_shape(planet, primary_mass_kg, planned_objects=None, parent_distance_au=None):
    """
    Draws `planet.eccentricity` (GEN.727, `physics.orbit_elements`): a planet's from the law of the system it is in (`planned_objects`
    is the number of planets and belts planned round its star, one for a single planet; `None` counts as several), cut where its
    periapsis would fall inside the star's Roche limit; a moon's from how it came to be there (`planet.moon_origin`), cut where its
    apoapsis would leave its parent's stable Hill sphere (`parent_distance_au` is the parent's distance from the star). The tides of the
    orbit then take what they can in the age of the star, from a planet's and a captured moon's eccentricity, which nothing renews; a
    regular moon's is forced (by its neighbours and the resonances among them, as Io's is by Europa and Ganymede) and stays, which is what
    makes its tide a lasting heat source (GEN.737). Draws the eccentricity before the climate, which goes by the mean light of
    the orbit; the plane and the angles come with the rest of the motion (`generate_orbital_motion_properties`).

    Args:
        planet (Planet): The planet or moon; needs `distance`, `mass`, `radius` and `star`.
        primary_mass_kg (float): The mass (kg) of what it orbits.
        planned_objects (int, optional): As above.
        parent_distance_au (float, optional): For a moon, the parent's distance from the star, AU.
    """
    if planet.is_moon:
        hill_au = None
        if parent_distance_au is not None and primary_mass_kg and planet.star.mass:
            hill_au = _km_to_au(_m_to_km(calculate_hill_sphere(_au_to_m(parent_distance_au), primary_mass_kg, planet.star.mass)))
        e_max = orbit_elements.max_moon_eccentricity(planet.distance, hill_au, getattr(planet, "host_eccentricity", None) or 0.0)
        birth = orbit_elements.draw_moon_eccentricity(planet.moon_origin, e_max)
    else:
        e_max = orbit_elements.max_planet_eccentricity(planet.distance, planet.star, planet)
        birth = orbit_elements.draw_planet_eccentricity(None if planned_objects is None else planned_objects > 1, e_max, giant=planet.body_type == 'g')
    if planet.is_moon and planet.moon_origin != "captured":
        planet.eccentricity = birth
        return
    planet.eccentricity = orbit_elements.settled_eccentricity(planet, birth, _gyr_to_yr(planet.star.age), primary_mass_kg, planet.distance)
    planet._orbit_birth = birth  # not stored: `settle_orbit` takes the tides again once the star's age is final (GEN.985)


def generate_orbital_motion_properties(planet, primary_mass_kg):
    """
    Sets a planet's (or moon's) 3D orbital orientation, position, orbital
    speed, and rotation period.

    Together with `planet.distance` (its semi-major axis) and
    `planet.eccentricity` (`generate_orbit_shape`), `orbital_inclination_deg`
    (the tilt of the orbital plane) and `orbital_ascending_node_deg` (where
    that plane crosses its primary's reference plane) fully orient the
    orbit in 3D, with a periapsis drawn here (uniform round the plane); the
    mean anomaly drawn here is where the body starts around it (uniform
    in time, so a body is found nearer apoapsis most often), and neither
    is kept: the body is its state vector (`position_x/y/z`, `velocity_*_kms`,
    see `update_orbital_position`), which `planetgen.cli.orbits` advances
    over time on `planet.period` with `kepler.universal_step` (GEN.289).

    `position_x/y/z` (AU, via `update_orbital_position`) are this body's
    Cartesian position *relative to its orbital anchor* -- the star (or,
    for a binary system, the `BinaryStarProxy` standing in for the system's
    combined center) for a planet, the parent planet for a moon -- the same
    "each body positioned relative to its immediate primary, not some
    absolute frame" convention `docs/design/galaxy/galaxy-coordinate-system.md`
    already uses one level up for sectors/systems relative to the galactic
    center.

    `rotation_period_hours` is a separate, purely descriptive "day length"
    stat (this generator doesn't track rotational phase -- nothing consumes
    "which side currently faces the primary"). `generate_spin` draws the
    body's birth spin and wears it down by the tide of what it orbits for the
    age of its star (`body_spin`, GEN.731 and GEN.957): the moon ends up locked
    (`rotation_period_hours` == its own orbital period, converted to hours)
    or in a spin-orbit resonance when the tide has had time to bring it there,
    which is the norm for large or close-in moons, not a rare special case,
    but genuinely not universal: a large moon far from a low-mass primary can
    easily have a locking time longer than the system itself, which is
    exactly why not every moon comes out tidally locked here either.

    Args:
        planet (Planet): The planet or moon to set these properties on.
                         Must already have `distance`, `period`, `mass`,
                         `radius`, and `body_type` set (see
                         `Planet.__init__`).
        primary_mass_kg (float): The mass (kg) of the body this one
                                 actually orbits -- the star for an
                                 ordinary planet, the parent planet for a
                                 moon (see `Planet.__init__`'s parameter of
                                 the same name, which this is threaded
                                 straight through from).
    """
    if planet.is_moon and getattr(planet, "moon_origin", None) is not None:
        # GEN.714: the orbit follows from how the moon came to be there (a captured one can be retrograde)
        planet.orbital_inclination_deg = moon_composition.draw_inclination_deg(planet.moon_origin)
    elif planet.is_moon:
        planet.orbital_inclination_deg = draw.uniform(0, constants.MOON_ORBITAL_INCLINATION_MAX_DEG)
    else:
        planet.orbital_inclination_deg = orbit_elements.draw_planet_inclination_deg()  # GEN.727: Rayleigh, as the Kepler multiples are
    if not planet.is_moon:
        apply_orbit_oddities(planet)      # GEN.1049: a close giant round a hot star can be given a polar orbit before the spin reads the plane
    planet.orbital_ascending_node_deg = draw.uniform(0, 360)
    arg_periapsis_deg = orbit_elements.draw_argument_of_periapsis_deg()
    mean_anomaly_rad = draw.uniform(0, _TWO_PI_RAD)

    generate_spin(planet, primary_mass_kg)
    true_anomaly_rad, _ = kepler.true_anomaly_and_distance_elliptical(mean_anomaly_rad, planet.eccentricity, planet.distance)
    update_orbital_position(planet, arg_periapsis_deg + math.degrees(true_anomaly_rad), arg_periapsis_deg)


def apply_orbit_oddities(planet):
    """
    Rolls the rows of the orbit and history levers (L2: W3 a highly eccentric orbit, W4 a polar orbit, W5 a resonant chain's forced eccentricity, GEN.1049; L5: W14 a disintegrating planet, GEN.1050) again on the planet as it is
    now, taking off the ones it had: a row whose gate still holds is given again from its own stream, so the ids always match the orbit the planet has.
    """
    if planet.body_type not in ("t", "g") or planet.is_moon:
        return
    oddities.drop(planet, ("L2", "L5"))
    oddities.apply(planet, levers=("L2", "L5"))


def apply_ring_oddities(planet):
    """Rolls the rows of the ring lever (L6: W9 a ring from a fallen moon, GEN.1050) on a planet whose ring (`rings.apply`) is made, taking off the ones it had."""
    if planet.body_type not in ("t", "g") or planet.is_moon:
        return
    oddities.drop(planet, ("L6",))
    oddities.apply(planet, levers=("L6",))


def finish_orbit_oddities(planet):
    """
    The orbit rows again once the orbit is final (after `settle_orbit`, GEN.1049): the eccentricity has worn down for the star's final age and the planets' places in a
    chain are known (`resonant_chains.give_forced_eccentricities`). A chain planet held on a forced eccentricity has its position, its moons' host eccentricity and the light of
    the orbit made again, and a rocky one its air, class and spin (the tide of the orbit is a heat source, `remake_end_stages`).

    Returns:
        bool: Whether the eccentricity moved.
    """
    before = planet.eccentricity
    apply_orbit_oddities(planet)
    if abs(planet.eccentricity - before) <= ECCENTRICITY_SETTLE_TOLERANCE:
        return False
    remake_after_orbit_change(planet)
    return True


def follow_orbit_change(planet):
    """What follows a planet's eccentricity having moved once the orbit is final, short of its air: its moons' host eccentricity, its place on the new ellipse and the exposure of it and its moons
    at the new mean light."""
    for moon in getattr(planet, "moons", ()):
        moon.host_eccentricity = planet.eccentricity
    update_orbital_position(planet)
    light_au = orbit_elements.flux_distance_au(planet.distance, planet.eccentricity)
    for body in (planet, *getattr(planet, "moons", ())):
        if getattr(body, "_star_distance_au", None):
            magnetism.update_exposure(body, light_au)


def remake_after_orbit_change(planet):
    """`follow_orbit_change` and, for a rocky planet, its air, class and spin made again at once (the tide of the orbit is a heat source, `remake_end_stages`). Used by `finish_orbit_oddities`,
    after which nothing moves the orbit again."""
    follow_orbit_change(planet)
    if planet.body_type == "t":
        remake_end_stages(planet, planet)


def stale_air(planet):
    """`follow_orbit_change`, then the planet's and its moons' air marked as made on an orbit that has since moved: `age_moved` says so (the age stamp `_age_drawn` is zeroed), so the final-age refresh,
    `refinish_rocky`, makes the air, class and spin again once, on the orbit the settle pass leaves (GEN.1450: the exchange comes first, the tide wears the eccentricity down after it,
    and a planet made again at the exchange's eccentricity would be left on the wrong one)."""
    follow_orbit_change(planet)
    for body in (*([planet] if planet.body_type == "t" else []), *getattr(planet, "moons", ())):        # a giant is not made again by the refresh (its moons are)
        if getattr(body, "_age_drawn", None) is not None:
            body._age_drawn = 0.0
    if planet.body_type == "g":
        planet._air_stale = True                 # `settle_giant` makes a giant's air again on the orbit it has now, whether or not its radius moves (GEN.1505)


def carry_spin_with_orbit(planet, old_inclination_deg, old_node_deg):
    """
    GEN.1504: the spin axis of a planet whose orbital plane the secular exchange has moved (the elements are already the new ones) follows the plane: the axis is turned with the orbit normal,
    so the obliquity is kept (the planet's spin precesses about its orbit's normal far faster than the plane's nodes regress, so the obliquity is the adiabatic invariant; the angle round the normal
    is kept too), and a locked or resonant body, whose obliquity is its Cassini state's share of the inclination (`tidal_spin.cassini_tilt_deg`), keeps that share of the new inclination.
    The sky of the other suns follows the new normal, the day and night shares the new tilt, and the moons are set in the new equator again. No draws.
    """
    axis = tuple(getattr(planet, name, None) for name in ("spin_axis_x", "spin_axis_y", "spin_axis_z"))
    if any(c is None for c in axis):
        return
    new_inclination = planet.orbital_inclination_deg or 0.0
    tilt = None
    if spin.spin_state(planet.rotation_period_hours, planet.period)[0] in ("locked", "resonant"):
        fold = lambda inc: min(abs(inc), 180.0 - abs(inc))
        if fold(old_inclination_deg) > 1.0e-12:
            tilt = planet.axial_tilt_deg * fold(new_inclination) / fold(old_inclination_deg)
    old_normal = spin.orbit_normal(old_inclination_deg, old_node_deg)
    new_normal = spin.orbit_normal(new_inclination, planet.orbital_ascending_node_deg or 0.0)
    planet.spin_axis_x, planet.spin_axis_y, planet.spin_axis_z = spin.carry_axis(axis, old_normal, new_normal, tilt)
    if tilt is not None:
        planet.axial_tilt_deg = tilt
        daylight.generate(planet)
    planet.sky_suns = sky_suns.encode(sky_suns.suns_of(
        planet.star, (planet.spin_axis_x, planet.spin_axis_y, planet.spin_axis_z), new_normal,
        orbit_elements.flux_distance_au(planet.distance, getattr(planet, "eccentricity", 0.0))))
    if getattr(planet, "moons", None):
        orient_moons(planet)


def settle_obliquity(planets, host_mass_kg):
    """
    GEN.1454: a free-spinning planet whose tilt lies in a chaotic zone of the spin-orbit secular resonances (`obliquity_variation`: the torques of the star and the planet's moons on its
    bulge against the nodal modes of the other planets' pull on its orbit) is not at the tilt it was made with: over the age of the star its axis has wandered over the zone, so it is
    given a tilt drawn from the zone (even in the cosine) from a stream of the body's own, the axis kept at its angle round the orbit's normal. Called with the orbits final
    (`StarSystem._settle_to_final_age`) for one star's planets and belts in orbital order. A planet whose tilt is not chaotic keeps it (the regular wander of an Earth or a Mars
    is a few degrees to a few tens, and the stored tilt is one place in it); a locked or resonant body takes the Cassini state instead; a lone planet has no modes. The sky of the other
    suns follows the new tilt, the day and night shares it, and the moons are set in the new equator, as `carry_spin_with_orbit` does. Returns how many planets changed tilt.
    """
    from planetgen.physics import obliquity_variation, tidal_spin

    bodies = [p for p in planets if getattr(p, "body_type", None) in ("t", "g") and not getattr(p, "is_moon", False) and getattr(p, "period", None) and p.mass and p.distance]
    if len(bodies) < 2 or not host_mass_kg > 0.0:
        return 0
    a_au = [b.distance for b in bodies]
    masses = [b.mass for b in bodies]
    inclinations = [math.radians(b.orbital_inclination_deg or 0.0) for b in bodies]
    nodes = [math.radians(b.orbital_ascending_node_deg or 0.0) for b in bodies]
    changed = 0
    all_modes = None
    for index, body in enumerate(bodies):
        tilt = getattr(body, "axial_tilt_deg", None)
        rotation = getattr(body, "rotation_period_hours", None)
        axis = tuple(getattr(body, name, None) for name in ("spin_axis_x", "spin_axis_y", "spin_axis_z"))
        if tilt is None or not rotation or not body.radius or any(c is None for c in axis):
            continue
        if spin.spin_state(rotation, body.period)[0] != "free":
            continue
        omega = 2.0 * math.pi / (abs(rotation) * 3600.0)
        torquers = [obliquity_variation.Torquer(host_mass_kg, body.distance * constants.AU_TO_M, getattr(body, "eccentricity", 0.0) or 0.0)]
        torquers.extend(obliquity_variation.Torquer(moon.mass, moon.distance * constants.AU_TO_M, getattr(moon, "eccentricity", 0.0) or 0.0)
                        for moon in getattr(body, "moons", ()) if moon.mass and moon.distance)
        ellipticity = obliquity_variation.dynamical_ellipticity(body.mass, body.radius * 1000.0, omega, tidal_spin.moment_factor(body))
        alpha = obliquity_variation.precession_constant_rad_yr(ellipticity, omega, torquers)
        if all_modes is None:                    # the matrix is made once for the star's planets, and only if a planet is free to wander
            all_modes = obliquity_variation.all_sibling_modes(a_au, masses, host_mass_kg, inclinations, nodes)
        variation = obliquity_variation.analyse(alpha, tilt, all_modes[index])
        if not variation.chaotic:
            continue
        stream = draw.Stream(f"chaotic-obliquity:{float(body.mass)!r}:{float(body.distance)!r}:{float(host_mass_kg)!r}:{float(body.radius)!r}")
        new_tilt = obliquity_variation.chaotic_tilt_deg(variation, stream.random())
        normal = spin.orbit_normal(body.orbital_inclination_deg or 0.0, body.orbital_ascending_node_deg or 0.0)
        body.spin_axis_x, body.spin_axis_y, body.spin_axis_z = spin.carry_axis(axis, normal, normal, new_tilt)
        body.axial_tilt_deg = new_tilt
        daylight.generate(body)
        body.sky_suns = sky_suns.encode(sky_suns.suns_of(
            body.star, (body.spin_axis_x, body.spin_axis_y, body.spin_axis_z), normal,
            orbit_elements.flux_distance_au(body.distance, getattr(body, "eccentricity", 0.0))))
        if getattr(body, "moons", None):
            orient_moons(body)
        changed += 1
    return changed


SECULAR_OUTCOMES = ("exchanged", "held by the spacing", "unchanged", "outside the theory")
"""tuple: What `secular_exchange` reports for a star's planets: the orbits took the exchange; some eccentricity was held short of what the exchange gave because the spacing would not hold;
nothing moved by more than `ECCENTRICITY_SETTLE_TOLERANCE`; or the planets or the result were beyond the first-order theory's limits."""
SECULAR_BISECTION_STEPS = 12
"""int: The steps of the bisection for the largest eccentricity the spacing allows a planet (the interval narrows to 2.4e-4 of the rise)."""


def secular_exchange(planets, host_mass_kg, age_yr, keeps_spacing):
    """
    GEN.1450: gives the planets of one star the orbits the Laplace-Lagrange secular exchange (`secular_chain`) leaves them at the star's age, from the eccentricities and inclinations the
    generator made them with (after `generate_orbit_shape`, the spacing pass and the moons; the exchange comes before the tide wears an eccentricity down again, so `settle_orbit` takes the
    lower of the tide's and the exchange's). Semi-major axes stay; the eccentricity, the longitude of periapsis, the inclination and the node move. Draws nothing.

    The spacing the generator keeps between neighbours (mutual Hill radii between the inner one's apoapsis and the outer one's periapsis) is its screen for a chain that stays stable, and it
    is not broken: the eccentricities the exchange lowers are lowered, then those it raises are raised planet by planet, each to what the exchange gave or, if the spacing will not hold,
    to the most it will (a bisection). A chain packed to the screen's limit therefore keeps its eccentricities and takes the rest of the exchange (the inclinations, the periapses, the nodes).

    Args:
        planets (list): One star's planets and belts in orbital order (belts and moons do not take part).
        host_mass_kg (float): The star's mass.
        age_yr (float): The star's age, years.
        keeps_spacing (callable): Given `planets`, whether their spacing and the moons' stable ranges hold (no problem beyond those the system already had), asked with the new
            eccentricities set.

    Returns:
        str: One of `SECULAR_OUTCOMES`.
    """
    from planetgen.physics import secular_chain

    bodies = [p for p in planets if getattr(p, "body_type", None) in ("t", "g") and not getattr(p, "is_moon", False) and getattr(p, "period", None) and p.mass and p.distance]
    if len(bodies) < 2 or not host_mass_kg > 0.0 or not age_yr > 0.0:
        return "unchanged"
    start = []
    for body in bodies:
        orbit = body.orbit_from_vector()
        start.append((orbit["eccentricity"], (orbit["ascending_node"] + orbit["arg_periapsis"]) % (2.0 * math.pi), orbit["inclination"], orbit["ascending_node"]))
    result = secular_chain.evolve([b.distance for b in bodies], [b.mass for b in bodies], host_mass_kg,
                                  [x[0] for x in start], [x[1] for x in start], [x[2] for x in start], [x[3] for x in start], age_yr)
    if not result.linear:
        return "outside the theory"
    old = [body.eccentricity for body in bodies]
    ceilings = [orbit_elements.max_planet_eccentricity(body.distance, body.star, body) for body in bodies]
    targets = [min(e, ceiling) for e, ceiling in zip(result.eccentricity, ceilings)]
    moved = [i for i, body in enumerate(bodies) if abs(targets[i] - old[i]) > ECCENTRICITY_SETTLE_TOLERANCE
             or abs(math.degrees(result.inclination[i]) - (body.orbital_inclination_deg or 0.0)) > 1e-3]
    if not moved:
        return "unchanged"
    for i, body in enumerate(bodies):
        if targets[i] < old[i]:
            body.eccentricity = targets[i]
    held = False
    for i in sorted((i for i in range(len(bodies)) if targets[i] > old[i]), key=lambda i: old[i] - targets[i]):
        body = bodies[i]
        body.eccentricity = targets[i]
        if keeps_spacing(planets):
            continue
        held = True
        low, high = old[i], targets[i]
        for _ in range(SECULAR_BISECTION_STEPS):
            body.eccentricity = 0.5 * (low + high)
            low, high = (body.eccentricity, high) if keeps_spacing(planets) else (low, body.eccentricity)
        body.eccentricity = low
    for i in sorted(set(moved) | {i for i, body in enumerate(bodies) if body.eccentricity != old[i]}):      # any change of the eccentricity moves the mean light
        body = bodies[i]
        phase = body_state.argument_of_latitude_deg((body.position_x, body.position_y, body.position_z), body.orbital_inclination_deg or 0.0, body.orbital_ascending_node_deg or 0.0)
        old_inclination, old_node = body.orbital_inclination_deg or 0.0, body.orbital_ascending_node_deg or 0.0
        if getattr(body, "moons", None):
            release_moons(body)                  # the moons lie in the planet's equator: written back in it, turned onto the new axis below
        body.orbital_inclination_deg = math.degrees(result.inclination[i])
        body.orbital_ascending_node_deg = math.degrees(result.node_longitude[i])
        carry_spin_with_orbit(body, old_inclination, old_node)
        update_orbital_position(body, phase, math.degrees((result.periapsis_longitude[i] - result.node_longitude[i]) % (2.0 * math.pi)))
        for moon in getattr(body, "moons", ()):
            moon.host_eccentricity = body.eccentricity
        if body.eccentricity != old[i] or i in moved:
            stale_air(body)
    return "held by the spacing" if held else "exchanged"


def draw_birth_spin(planet):
    """The spin a planet or moon is born with, `(rate rad/s, tilt degrees, capture units)`: a rocky planet's from the giant impacts that built it
    (`spin.accretion_spin`), a moon's or giant's at a rate from `spin.planet_day_hours` with a Rayleigh tilt (GEN.104, GEN.979). The capture units are one uniform
    for each resonance the spin can pass without being caught (`spin_torques.CAPTURE_RATIOS`, GEN.1172), drawn here so that running the history again decides the same."""
    if planet.body_type == "t" and not planet.is_moon:
        rate, tilt = spin.accretion_spin(planet.mass, planet.radius)
    else:
        rate, tilt = 2.0 * math.pi / _h_to_s(spin.planet_day_hours(planet.body_type)), spin.rayleigh_tilt_deg()
        if planet.body_type == "g" and not planet.is_moon:
            tilt = giant_tilt.tilt_after_merger_deg(planet.mass, planet.radius, tilt)    # GEN.1165: an ice giant's merger can lay it on its side
    return rate, tilt, tuple(draw.random() for _ in spin_torques.CAPTURE_RATIOS)


def generate_spin(planet, primary_mass_kg, birth=None, found=None):
    """
    Draws a planet's or moon's `rotation_period_hours` and spin axis (GEN.104, GEN.731; `physics.spin.SPIN_FIELDS`). It is born spinning
    (`draw_birth_spin`, kept on the body as `_birth_spin` so that a later pass can run its history again without drawing it) and the tide
    of what it orbits then wears the spin toward the orbit's rate for the age of its star (`body_spin.evolve`: the gravitational tide with
    the body's own frequency-dependent rheology and, for a rocky planet with an air, the thermal tide of the air, GEN.957; the orbit's
    eccentricity from `generate_orbit_shape`). A planet with moons runs the history of the moons with its spin
    (`moon_history.apply`, GEN.959), which calls this again with the `found` state it reached. A body
    the tide has brought all the way down turns once an orbit (locked) or in a spin-orbit resonance such as Mercury's 3:2
    (`spin.spin_state`), with its axis in the first Cassini state (a small share of the orbit's inclination off the orbit normal, `tidal_spin.cassini_tilt_deg`); one the thermal and gravitational tides hold in balance away from synchronism turns at that
    rate, upright or upside down; any other keeps its tilt, a rocky planet's from the impacts and a moon's or giant's a Rayleigh one. Called again when a move changes the body's period, since that can change the state.

    Args:
        planet (Planet): The planet or moon.
        primary_mass_kg (float): What it orbits.
        birth (tuple, optional): `(rate, tilt, capture units)` to start from instead of drawing it.
        found (spin_torques.SpinResult, optional): The state to set instead of running the history (`moon_history.apply`'s).
    """
    if birth is None:
        birth = draw_birth_spin(planet)
    planet._birth_spin = birth
    spin_rate, birth_tilt, capture_units = birth
    if found is None:
        if not planet.is_moon and (getattr(planet, "moons", None) or getattr(planet, "_born_moons", None)):
            release_moons(planet)      # GEN.1352: the history moves them against the equator they were drawn in
            return moon_history.apply(planet, primary_mass_kg)
        found = body_spin.evolve(planet, primary_mass_kg, _au_to_m(planet.distance), spin_rate, _yr_to_s(_gyr_to_yr(planet.star.age)),
                                 getattr(planet, "eccentricity", 0.0) or 0.0, capture_units=capture_units)
    period_hours = planet.period * _yr_to_h(1.0)
    if found.ratio is not None:
        planet.rotation_period_hours = period_hours / found.ratio
        tilt = tidal_spin.cassini_tilt_deg(planet.orbital_inclination_deg, draw.random())    # GEN.1165: the first Cassini state, not exactly 0
    else:
        planet.rotation_period_hours = 2.0 * math.pi / abs(found.rate) / 3600.0
        if found.kind == "asynchronous":
            tilt = 0.0 if found.rate > 0.0 else 180.0      # the tides have brought it to its balance: upright, or upside down (Venus)
        else:
            tilt = birth_tilt
            if found.rate < 0.0:
                tilt = 180.0 - tilt
    # GEN.1048 (oddity W2): the rows of the spin lever are rolled again on the body as it is now (its air, its moons and the age change between calls) and a
    # thick-air planet the tides leave turning slowly backwards takes the Venus-type state in place of the tide history's
    planet._venus_spin = None
    if planet.body_type in ("t", "g") and not planet.is_moon:
        planet.axial_tilt_deg = tilt             # the rows read the tilt the planet is about to have (set below)
        oddities.drop(planet, ("L1",))
        oddities.apply(planet, levers=("L1",))
        venus_tilt = venus_spin.set_on(planet)
        if venus_tilt is not None:
            tilt = venus_tilt
    spin.set_spin(planet, spin.orbit_normal(planet.orbital_inclination_deg,
                                            planet.orbital_ascending_node_deg), tilt)
    frame = getattr(planet, "_frame_axis", None) if planet.is_moon else None
    if frame is not None:      # GEN.1352: a moon already set in the system's frame has its new axis turned the way its position was
        planet.spin_axis_x, planet.spin_axis_y, planet.spin_axis_z = spin.turn_to_axis((planet.spin_axis_x, planet.spin_axis_y, planet.spin_axis_z), frame)
    if not planet.is_moon:
        planet.sky_suns = sky_suns.encode(sky_suns.suns_of(
            planet.star, (planet.spin_axis_x, planet.spin_axis_y, planet.spin_axis_z),
            spin.orbit_normal(planet.orbital_inclination_deg, planet.orbital_ascending_node_deg),
            orbit_elements.flux_distance_au(planet.distance, getattr(planet, "eccentricity", 0.0))))     # GEN.1458: the other suns of its sky
    if not planet.is_moon and getattr(planet, "moons", None):
        orient_moons(planet)      # GEN.1352: its moons lie in the equator this axis gives
    # GEN.767: the day and night shares follow the rotation and tilt just set (a planet's; a moon's need its parent's axis, `refinish_rocky`)
    if not planet.is_moon:
        daylight.generate(planet)
    # GEN.86: the field depends on the rotation just set.
    magnetism.generate_field(planet, None if planet.body_type == "g" else interior_heat.core_of(planet))
    # GEN.87: the surface dose, from the air, water and field above.
    radiation.generate_dose(planet)
    # DB.96: the atmosphere colour, from the air, the ozone and the aerosols above.
    body_colour.generate(planet)
    # GEN.89: the habitability score, from all of the above.
    habitability_world.generate(planet, internal_heat_flux_w_m2(planet))



def _turn_moon(moon, axis, turn):
    """Turns a moon's position, velocity and spin axis by `turn` (`spin.turn_to_axis` or `spin.turn_from_axis`) about the planet's pole `axis`."""
    moon.set_position_au(*turn((moon.position_x, moon.position_y, moon.position_z), axis))
    moon.set_velocity_kms(*turn((moon.velocity_x_kms, moon.velocity_y_kms, moon.velocity_z_kms), axis))
    if moon.spin_axis_x is not None:
        moon.spin_axis_x, moon.spin_axis_y, moon.spin_axis_z = turn((moon.spin_axis_x, moon.spin_axis_y, moon.spin_axis_z), axis)


@contextlib.contextmanager
def _in_equator_frame(body):
    """
    GEN.1352: holds the position and velocity of a moon `orient_moons` has set in the system's frame in the frame of its planet's equator while the block runs
    (its elements, which `move_orbit` and `update_orbital_position` read, are measured there) and turns them back after. A body not so set is left as it is.
    """
    axis = getattr(body, "_frame_axis", None)
    if axis is None:
        yield
        return
    body.set_position_au(*spin.turn_from_axis((body.position_x, body.position_y, body.position_z), axis))
    body.set_velocity_kms(*spin.turn_from_axis((body.velocity_x_kms, body.velocity_y_kms, body.velocity_z_kms), axis))
    try:
        yield
    finally:
        body.set_position_au(*spin.turn_to_axis((body.position_x, body.position_y, body.position_z), axis))
        body.set_velocity_kms(*spin.turn_to_axis((body.velocity_x_kms, body.velocity_y_kms, body.velocity_z_kms), axis))


def release_moons(planet):
    """
    Writes the moons of `planet` that `orient_moons` set in the system's frame back in the frame of the planet's equator, the one their
    elements (`orbital_inclination_deg`, `orbital_ascending_node_deg`) are measured in and the history of the moons moves them in (no draws).
    """
    for moon in (*planet.moons, *(getattr(planet, "_born_moons", None) or ())):      # the born ones too: the history draws on them again, and one dropped since may be back
        axis = getattr(moon, "_frame_axis", None)
        if axis is not None:
            _turn_moon(moon, axis, spin.turn_from_axis)
            moon._frame_axis = None


def orient_moons(planet):
    """
    GEN.1352: sets the moons of `planet` in the system's frame. A moon's inclination and node are drawn and kept against the equator of its planet
    (`generate_orbital_motion_properties`, `moon_history.retrograde`), so its position, velocity and spin axis are turned from that frame onto the planet's
    spin axis: a moon of a planet tipped by 98 degrees goes round in its tipped equator, and one of a planet with no tilt goes round in the plane of its
    planet's orbit. A turn only, no draws; the planet's wobble from its moons (`reflex_offset_*`) follows the new positions. The irregular moons keep
    their own frame (`generation/irregular_moons.py`).
    """
    axis = tuple(getattr(planet, name, None) for name in ("spin_axis_x", "spin_axis_y", "spin_axis_z"))
    if any(c is None for c in axis):
        return
    release_moons(planet)
    for moon in planet.moons:
        _turn_moon(moon, axis, spin.turn_to_axis)
        moon._frame_axis = axis
    if planet.moons:
        planet.reflex_offset_x, planet.reflex_offset_y, planet.reflex_offset_z = calculate_reflex_offset(
            planet.mass, [(m.mass, m.position_x, m.position_y, m.position_z) for m in planet.moons])


def move_orbit(body, distance_au, primary_mass_kg):
    """
    Sets the semi-major axis of a planet's or moon's orbit and what follows from it (period, Hill sphere, position and speed), keeping the
    plane, the shape and where the body is along the orbit (no draws): the move of a moon by the tide of its planet (`moon_history`, GEN.959).
    """
    with _in_equator_frame(body):
        orbit = body.orbit_from_vector()
        arg_periapsis_deg = math.degrees(orbit["arg_periapsis"])
        phase_deg = body_state.argument_of_latitude_deg((body.position_x, body.position_y, body.position_z),
                                                        body.orbital_inclination_deg, body.orbital_ascending_node_deg)
        body.distance = distance_au
        body.period = calculate_orbital_period_years(distance_au, primary_mass_kg)
        update_hill_sphere(body)
        _update_orbital_position(body, phase_deg, arg_periapsis_deg)


def update_orbital_position(planet, phase_deg=None, arg_periapsis_deg=0.0):
    """
    Recomputes `planet.position_x/y/z` (AU, relative to its orbital
    anchor), its velocity, `planet.orbital_speed_kms`, and
    `planet.min_update_interval_years` from its current `distance` (the
    semi-major axis), `eccentricity`, `period`, and orbital plane
    (`orbital_inclination_deg`/`orbital_ascending_node_deg`), as the
    ellipse `state_vectors.state_from_elements` makes of them (a circle
    when the eccentricity is 0).

    Called both at initial generation (`generate_orbital_motion_properties`,
    which passes the angles it rolled) and again by
    `StarSystem.validate_system` whenever it corrects a top-level planet's
    `distance` post-hoc to resolve an orbital overlap -- the same
    "recompute anything that depends on distance" treatment `period`
    already gets there (see that method's docstring). There the body keeps
    where it is along its plane and the orbit keeps its periapsis: no
    angle is passed and the ones the current state has are used.
    Position/speed/`min_update_interval_years` would
    otherwise silently go stale relative to the corrected `distance`/
    `period` the same way `period` itself used to before that fix.

    `planetgen.cli.orbits`'s periodic time-based advancement is a separate
    path (`planetgen.db.store.advance_orbital_phases`) that moves the
    stored state vector in the database rather than through this function
    -- this generator has no live "simulation loop" over in-memory
    objects, only one-shot generation (this function) followed by periodic
    database updates (see that module's own docstring).
    `min_update_interval_years` depends only on `period`, fixed once
    generated, so the update never touches it -- it only *reads* the
    column, to decide whether a given call's `elapsed_years` is even worth
    writing for this row.

    Args:
        planet (Planet): The planet or moon to update in place. Must
                         already have `distance`, `eccentricity`, `period`,
                         `orbital_inclination_deg` and
                         `orbital_ascending_node_deg` set, and, without
                         `phase_deg`, its position and velocity.
        phase_deg (float, optional): Where along its orbit plane the body
                         is (the argument of latitude, from the ascending
                         node, which is the argument of periapsis plus the
                         true anomaly); the one its current position has
                         when not given.
        arg_periapsis_deg (float, optional): The angle from the ascending
                         node to periapsis (meaningless for a circle);
                         the one its current state has when `phase_deg` is
                         not given.
    """
    with _in_equator_frame(planet):
        _update_orbital_position(planet, phase_deg, arg_periapsis_deg)


def _update_orbital_position(planet, phase_deg, arg_periapsis_deg):
    """`update_orbital_position` for a body whose position and velocity are in the frame its elements are measured in."""
    if phase_deg is None:
        orbit = planet.orbit_from_vector()
        arg_periapsis_deg = math.degrees(orbit["arg_periapsis"])
        phase_deg = body_state.argument_of_latitude_deg(
            (planet.position_x, planet.position_y, planet.position_z),
            planet.orbital_inclination_deg, planet.orbital_ascending_node_deg)
    position, velocity = state_from_elements(
        (2.0 * math.pi / planet.period) ** 2 * planet.distance ** 3, orbit_elements.periapsis_au(planet.distance, planet.eccentricity), planet.eccentricity,
        math.radians(planet.orbital_inclination_deg), math.radians(planet.orbital_ascending_node_deg), math.radians(arg_periapsis_deg),
        math.radians(phase_deg - arg_periapsis_deg))
    planet.set_position_au(*position)
    to_kms = _au_yr_to_kms(1.0)
    planet.set_velocity_kms(*(v * to_kms for v in velocity))
    planet.orbital_speed_kms = math.sqrt(sum(v * v for v in velocity)) * to_kms
    planet.min_update_interval_years = minimum_update_interval_years(planet.period)


def reconcile_zone(planet, primary_mass_kg, distance_override=None):
    """
    Re-derives `planet`'s zone from its *current* orbital distance and records the mass of the body it orbits.

    Meant to be called after something external to normal generation -- `StarSystem.validate_system`'s orbital-overlap
    correction -- has moved `planet.distance` (for an ordinary planet) or its parent's distance (for a moon, via
    `distance_override`). The mass, radius and composition of the body do not change with where it sits; its climate does, so the
    caller makes the air again (`calculate_atmospheric_conditions`) and the class follows (`apply_class`).

    Args:
        planet (Planet): The planet or moon to reconcile in place. Must already have `distance`, `zone` and `habitable_zone` set.
        primary_mass_kg (float): The mass (kg) of the body this one actually orbits -- the star for an ordinary planet, the parent
            planet for a moon.
        distance_override (float, optional): The distance (AU) to determine the zone from, when it differs from `planet.distance`
            -- a moon's own `distance` is its orbit around its *parent planet*, so a moon passes its parent's `distance` here.

    Returns:
        bool: True if the zone changed.
    """
    distance = distance_override if distance_override is not None else planet.distance
    new_zone = zone_for(planet.habitable_zone, distance, planet.star)
    zone_changed = new_zone != planet.zone
    planet.zone = new_zone
    planet.primary_mass_kg = primary_mass_kg  # what its Hill sphere is measured against (GEN.138)
    return zone_changed


def zone_for(habitable_zone, distance_au, star=None):
    """The zone letter ('h' hot, 'e' ecosphere, 'c' cold) of `distance_au` against `habitable_zone` (inner, outer) in AU. With the host `star` the other
    stars' light counts (GEN.1441): the body is placed at the distance of a body of a lone star that gets the same total light."""
    inner_bound, outer_bound = habitable_zone
    if star is not None:
        distance_au = radiation_field.lit_distance_au(star, distance_au)
    if distance_au < inner_bound:
        return 'h'
    if distance_au > outer_bound:
        return 'c'
    return 'e'


def hosting_population(planet):
    """
    What the hosting gate (`physics.hosting_gate`) calls `planet`: a small body (under `tuning.HOSTING_SMALL_BODY_MAX_RADIUS_KM` in
    radius) is a `trans-neptunian-cold` or `trans-neptunian-hot` body past `tuning.HOSTING_TRANS_NEPTUNIAN_MIN_SNOW_LINES`
    snow lines of its star (cold under `tuning.HOSTING_COLD_MAX_INCLINATION_DEG` of inclination, [R4207]), a
    `belt-small-body` in the cold zone short of that, and an `inner-small-body` in the hot or habitable zone; every other
    planet is a `planet`, which the gate's share table has no rows for (its host-class table decides).
    """
    if planet.radius >= tuning.HOSTING_SMALL_BODY_MAX_RADIUS_KM:
        return "planet"
    if planet.zone != "c":
        return "inner-small-body"
    if planet.distance < tuning.HOSTING_TRANS_NEPTUNIAN_MIN_SNOW_LINES * planet.star.snow_line_au:
        return "belt-small-body"
    inclination = planet.orbital_inclination_deg or 0.0
    return "trans-neptunian-cold" if inclination < tuning.HOSTING_COLD_MAX_INCLINATION_DEG else "trans-neptunian-hot"


def hosting_host(planet):
    """`planet` as the hosting gate reads it (`hosting_gate.Host`): its population, size, mass and orbit about what it
    orbits (the planets of this generator are on circular orbits); a gas or ice giant is fluid."""
    return hosting_gate.Host(
        hosting_population(planet), 2.0 * _km_to_m(planet.radius), planet.mass, parent_mass_kg=planet.primary_mass_kg,
        orbit_radius_m=_au_to_m(planet.distance), fluid=planet.body_type != "t")


def satellite_share(planet):
    """The chance that `planet` has satellites from the gate's `satellite-share` table, or `None` for a planet the table has
    no rows for (the host-class chance of the draw decides)."""
    return hosting_gate.satellite_share(hosting_host(planet))


def moon_orbit_bounds_km(planet):
    """
    The range of orbits, in km from `planet`'s center, where a moon of it
    can sit.

    Innermost: clear of the planet's body and the largest moon it could
    have (`planet.radius / 10**(1/3)`, so no moon touches it; this is also
    past the rigid-body Roche limit of ~1.26 planet radii for similar
    densities), plus 15 atmospheric scale heights (or 100 km with no
    atmosphere) of drag-free margin. Outermost: the prograde stability
    limit of [R1138], a fraction of the Hill radius that the planet's eccentricity lowers
    (0.4895 on a circular orbit, `clearance_ledger.moon_bounds_km`); past it the star strips the moon.

    Returns:
        tuple: `(low_km, high_km)`; `low_km >= high_km` means no room.
    """
    max_moon_radius = planet.radius / (10 ** (1 / 3))
    atmosphere_margin_km = planet.scale_height * 15 if planet.scale_height else 100
    low_km = planet.radius + max_moon_radius + atmosphere_margin_km
    high_km = clearance_ledger.moon_bounds_km(planet.radius, planet.hill_radius, getattr(planet, "eccentricity", None) or 0.0)[1]
    return low_km, high_km


def moon_stable_limit_km(planet, moon_eccentricity=0.0):
    """The farthest semi-major axis, km from `planet`'s center, a moon of this eccentricity keeps ([R1138]: the
    fit is for the semi-major axis, and the moon's own eccentricity lowers it)."""
    return clearance_ledger.moon_bounds_km(planet.radius, planet.hill_radius, getattr(planet, "eccentricity", None) or 0.0,
                                           moon_eccentricity or 0.0)[1]


def drop_unstable_moons(planet):
    """
    Removes the moons of `planet` that no longer fit after its class (and
    so its radius and mass) was regenerated: one orbiting outside
    `moon_orbit_bounds_km`, or one larger than the planet could hold
    (`moon_size_limits`: `planet.radius / 10**(1/3)` and a tenth of its
    mass). A real planet that lost mass this way
    would lose those moons to the star or to a collision.

    Returns:
        int: How many moons were dropped.
    """
    low_km, _high_km = moon_orbit_bounds_km(planet)
    max_moon_radius, max_moon_mass = moon_size_limits(planet)
    kept = [moon for moon in planet.moons
            if low_km <= _au_to_km(moon.distance) <= moon_stable_limit_km(planet, moon.eccentricity)
            and moon.radius <= max_moon_radius and moon.mass <= max_moon_mass]
    dropped = len(planet.moons) - len(kept)
    planet.moons[:] = kept
    irregular_moons.keep_stable(planet)      # the irregular moons are held to the same zone (GEN.1084); they do not count in the return, no reflex offset comes of them
    return dropped


def moon_size_limits(parent):
    """The largest moon `parent` can hold: `(radius_km, mass_kg)`, at most
    `radius / 10**(1/3)` and a tenth of its mass."""
    return parent.radius / (10 ** (1 / 3)), parent.mass / 10


MOON_REJECTION_LIMIT = 5
"""int: How many moons in a row `generate_moons` may discard for being larger than their parent can hold before it stops."""


def generate_moons(planet, moon_count=None):
    """
    Generates a system of moons for the given planet.

    This procedurally generates moons for the planet. It determines the
    number, size, and orbital distance of the moons based on the planet's
    properties, such as its mass and Hill radius. The generated moons are
    themselves `Planet` instances, with the `is_moon` flag set, appended to
    `planet.moons`.

    Args:
        planet (Planet): The planet to generate moons for.
        moon_count (int, optional): An exact number of moons to attempt to
            generate. If None (the default), moons are generated until no
            more orbital room is available, as before. If given, generation
            stops once `moon_count` moons exist, even if more room remains;
            if the planet doesn't have room for `moon_count` moons, fewer
            than requested may be generated.
    """
    if moon_count == 0:
        return
    # GEN.761: a body asks the hosting gate first (a satellite must fit, a stable orbit must exist, and the body's
    # population and size must have a chance of satellites at all); the share of the gated table is the draw's, in
    # `Planet._draw_moons`
    if not hosting_gate.can_host(hosting_host(planet)).allowed:
        return
    # The lightest moon must fit under the largest the parent can hold (`moon_size_limits`).
    limit_radius_km, limit_mass_kg = moon_size_limits(planet)
    if limit_mass_kg <= tuning.MOON_MIN_MASS_KG:
        return
    # GEN.1334: a rocky planet under a deep envelope (a sub-Neptune, 1 percent or more) has a disk of its own, so it is a giant parent for its moons' masses and origins alike (it
    # drew massive rocky impact moons and then disk moons of ice that were too large for it to hold)
    parent_kind = moon_composition.parent_kind(planet.envelope_mass_fraction or 0.0, planet.mass)
    parent_is_giant = planet.body_type == "g" or parent_kind == moon_composition.GIANT
    # GEN.962: a rocky planet (no row in the `satellite-share` table, not a giant) holds massive moons from giant impacts at a rate of its own; its other moons
    # are the captured ones of `captured_moons`
    impact_moons = not parent_is_giant and hosting_gate.satellite_share(hosting_host(planet)) is None
    impact_floor_kg = tuning.ROCKY_MASSIVE_MOON_MASS_RATIO * planet.mass if impact_moons else None
    # a massive moon is often denser or lighter than its parent can hold (the radius limit), so its draws are repeated further before the slot is given up
    rejection_limit = MOON_REJECTION_LIMIT * 4 if impact_moons else MOON_REJECTION_LIMIT

    low_orbit, high_orbit = moon_orbit_bounds_km(planet)
    total_orbit_distance = low_orbit
    parent_kind = moon_composition.parent_kind(planet.envelope_mass_fraction or 0.0, planet.mass)
    # GEN.466: a giant's regular moons share one circumplanetary disk (the lifetime, the ice line, the ice ramp)
    disk_system = moon_composition.draw_disk_system(planet.mass) if parent_kind == moon_composition.GIANT else None
    if disk_system is not None:
        moon_composition.store_disk(planet, disk_system)

    # Deferred import: planet imports this module at load time, so Planet
    # can't be imported here at module level without a circular import.
    from planetgen.generation.planet import Planet

    disk_zone = (moon_disk_zone.zone_km(planet.mass, planet.radius, (getattr(planet, "_birth_spin", None) or (None,))[0], low_orbit, high_orbit)
                 if disk_system is not None else None)
    rejected = 0
    while total_orbit_distance < high_orbit and total_orbit_distance < _au_to_km(planet.distance):
        if moon_count is not None and len(planet.moons) >= moon_count:
            break

        moon_mass = planet_kind.draw_moon_mass_kg(planet.mass, parent_is_giant, limit_mass_kg, impact_floor_kg)
        if moon_mass is None:
            break

        # Log-uniform, not linear-uniform: [total_orbit_distance, high_orbit]
        # can span many orders of magnitude (high_orbit reaches out to about
        # half the planet's own Hill radius, which for a large planet is tens
        # of millions of km -- far beyond where any real large moon
        # actually orbits, e.g. our Moon at ~384,400 km), and a plain
        # random.uniform over that range spends almost all its density in the
        # single largest order of magnitude, so nearly every moon landed
        # implausibly far out. Real moon systems are much closer to
        # log-spaced -- e.g. the Galilean moons run 421,700 / 671,100 /
        # 1,070,400 / 1,882,700 km, each roughly 1.5-1.6x the last, not a
        # near-flat distribution across the whole possible range. Sampling
        # log-uniformly gives every order of magnitude equal weight instead,
        # which both matches that real spacing pattern better and -- since
        # tidal-locking timescale scales with distance^6
        # (`body_spin.evolve`) -- stops real tidal-locking
        # physics from calling almost every moon unlocked purely because the
        # old distribution pushed it implausibly far from its primary.
        # GEN.714, GEN.1244: the origin comes first and the orbit follows from it; a moon of the disk forms inside the disk's zone (`moon_disk_zone`), the rest anywhere a moon can stay
        origin = moon_composition.draw_origin(parent_kind, moon_mass)
        span_low, span_high = total_orbit_distance, high_orbit
        if origin == "disk" and disk_zone is not None:
            span_low, span_high = max(total_orbit_distance, disk_zone[0]), min(high_orbit, disk_zone[1])
            if span_low >= span_high:    # the disk has no room left past the last moon
                draw.redrawn("moon.disk_full")
                rejected += 1
                if rejected >= MOON_REJECTION_LIMIT:
                    break
                continue
        moon_distance_km = math.exp(draw.uniform(math.log(span_low), math.log(span_high)))
        moon_distance = _km_to_au(moon_distance_km)
        new_moon = Planet(planet.system_config, planet.star, planet.habitable_zone, moon_distance,
                          mass=moon_mass, zone_override=planet.zone,
                          distance_override=planet.distance, is_moon=True, primary_mass_kg=planet.mass,
                          moon_origin=origin, moon_disk_system=disk_system, host_eccentricity=planet.eccentricity)
        if new_moon.radius > limit_radius_km:    # denser or lighter than the parent could hold; its slot stays open
            draw.redrawn("moon.too_large")
            rejected += 1
            if rejected >= rejection_limit:
                break
            continue
        rejected = 0
        planet.moons.append(new_moon)
        total_orbit_distance = _au_to_km(new_moon.distance) + _au_to_km(new_moon.min_orbit_distance)
        if impact_moons and moon_count is None and draw.random() >= tuning.ROCKY_IMPACT_MOON_CHANCE:    # GEN.962: a further impact, or no more
            break

    if planet.moons:
        # This planet's own reflex-offset "wobble" (schema v20) from the
        # combined pull of its own moons -- a proper two-body treatment
        # alongside each moon's own unchanged position_x/y/z (relative to
        # this planet); see Planet.reflex_offset_x's own docstring and
        # orbits.calculate_reflex_offset. Left at its 0.0 default when this
        # planet ends up with no moons (every early-return path above).
        planet.reflex_offset_x, planet.reflex_offset_y, planet.reflex_offset_z = calculate_reflex_offset(
            planet.mass, [(m.mass, m.position_x, m.position_y, m.position_z) for m in planet.moons]
        )
