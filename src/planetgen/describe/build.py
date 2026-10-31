# planetgen/describe/build.py

"""
`describe(system)`: the Document of a star system (docs/design/architecture/system-description.md).

A pure function of the `StarSystem` object graph, whether it was generated a
moment ago or rebuilt from the database (`db.store.load_star_system`). Every
number goes in as a typed quantity; nothing here formats one. The words are
the words the old per-class `to_paragraph_list()` methods used (wording is
reworked by GEN.221, not here); what changed is that the facts are table rows
and the sentences carry quantities, not finished strings.
"""

import math

from planetgen import tuning
from planetgen.describe import life_facts, model
from planetgen.describe.hierarchy import hierarchy_sentences
from planetgen.describe.model import (
    fact,
    facts,
    inline,
    measured,
    prose,
    qty,
    section,
    swatch,
)
from planetgen.generation.belt import format_composition_summary
from planetgen.generation.comet import PERIOD_CLASS_LABELS
from planetgen.generation.star import compressed_heliosphere_radius
from planetgen.generation.stellar_object import supported_evolutionary_scales
from planetgen.names.bodies import close_pair_label
from planetgen.physics import atmosphere as atmosphere_physics
from planetgen.physics import (
    atmosphere_colour,
    equipment,
    gases,
    giant_interior,
    habitability,
    habitability_explain,
    meteorite_groups,
    outside_hours,
    pre_main_sequence,
    protection,
    regional_temperature,
    unit_registry,
)
from planetgen.physics.quantity import Quantity
from planetgen.util.format import three_figures_or_tiny

ECCENTRICITY_NOTE = "(0 is a circle, 1 is the escape orbit)"
"""str: What a bare eccentricity reads beside its number (GEN.772)."""

HELIOSPHERE_COMPRESSION_NOTED = 0.95
"""float: A compressed heliosphere is mentioned when it is under this share of the open-space one."""

# --- quantities in the units the model stores -------------------------------------------------------

def _held(unit, value, shown=None, **extra):
    """The inline quantity of `value` in the registry unit `unit`, shown as `shown` (the registry kind's name when omitted)."""
    registry_kind = unit_registry.unit(unit).kind
    return measured(shown or registry_kind, None if value is None else Quantity(registry_kind, value, unit), **extra)


def _au(value):
    return _held("au", value)


def _km(value):
    return _held("km", value)


def _radius_km(value):
    return _held("km", value, shown="radius")


def _years(value):
    return _held("yr", value, shown="time")


def _gy(value):
    return _held("Gyr", value, shown="time")


def _kms(value):
    return _held("km_s", value)


def _kelvin(value, companions=True):
    return qty("temperature", value) if companions else qty("temperature", value, companions=False)


def _pascal(value):
    return qty("pressure", value)


def _gravity(value_g):
    return _held("g0", value_g, shown="gravity")


def _mass(kilograms, solar=False):
    quantity = None if kilograms is None else Quantity("mass", kilograms, "kg")
    return measured("mass", quantity, **({"sol": quantity.as_("Msun")} if solar else {}))


def _mass_solar(solar_masses):
    """A mass held in solar masses, shown with its ratio to the Sun."""
    quantity = Quantity("mass", solar_masses, "Msun")
    return measured("mass", quantity, sol=quantity.Msun)


def _luminosity(watts):
    quantity = Quantity("power", watts, "W")
    return measured("luminosity", quantity, sol=quantity.as_("Lsun"))


def _galactic_orbit(speed_kms, period_gy):
    return inline(_kms(speed_kms), " (", _gy(period_gy), " per orbit)")


def _offset(*components_au):
    """The length of a vector given in AU, as a length in AU."""
    return math.sqrt(sum(component ** 2 for component in components_au))


def _ref(kind, obj):
    """`kind:id` for a body that has a database row, else `None`."""
    db_id = getattr(obj, "db_id", None)
    return None if db_id is None else f"{kind}:{db_id}"


def _noun(planet):
    return "moon" if planet.is_moon else "planet"


# --- stars ---------------------------------------------------------------------------------------------

def _star_facts(star):
    rows = [
        fact("type", "Type", star.type),
        fact("radius", "Radius", _radius_km(star.radius)),
        fact("mass", "Mass", _mass(star.mass, solar=True)),
        fact("temp", "Temperature", _kelvin(star.temperature, companions=False)),
        fact("lum", "Luminosity", _luminosity(star.luminosity)),
        fact("hab", "Habitable zone", "Between ", _au(star.habitable_zone[0]), " and ", _au(star.habitable_zone[1])),
        fact("orbit", "Galactic orbit", _galactic_orbit(star.galactic_orbital_speed_kms, star.galactic_orbital_period_gy)),
    ]
    if star.reflex_offset_x or star.reflex_offset_y or star.reflex_offset_z:
        offset = _offset(star.reflex_offset_x, star.reflex_offset_y, star.reflex_offset_z)
        rows.append(fact("wobble", "Planetary wobble", _au(offset), " from its nominal position, pulled by its own planets"))
    if star.disk_kind:
        rows.append(fact("disk", "Disk", _star_disk_text(star)))
    rows.append(fact("loc", "Location", star.name))
    return facts("Star", rows)


_DISK_LABELS = {"post-agb": "Circumbinary post-AGB"}
"""dict: The words for a `disk_kind` whose name does not read as a sentence."""


def _star_disk_text(star):
    if star.disk_mass_solar is None:    # a debris disk (GEN.787): its dust light is measured, its mass is not
        return inline(_DISK_LABELS.get(star.disk_kind) or star.disk_kind.replace("-", " ").capitalize(), " disk, from ",
                      _au(star.disk_inner_au), " to ", _au(star.disk_outer_au), ", with no bodies in it")
    return inline(_DISK_LABELS.get(star.disk_kind) or star.disk_kind.replace("-", " ").capitalize(), " disk of ",
                  _mass_solar(star.disk_mass_solar), ", from ",
                  _au(star.disk_inner_au), " to ", _au(star.disk_outer_au), ", with no bodies in it")


LIFESPAN_SENTENCES = {
    "fast": "Its lifespan is too short for anything slower than the fastest evolution, so only the simplest light-driven chemistry could arise.",
    "normal": "That is long enough for life that evolves as fast as Earth's did, but not for the slower paces of evolution.",
    "slow": "That is long enough for even the slowest paces of evolution to reach a technological civilization.",
}
"""dict: What a main-sequence star's lifespan means for life, by the slowest pace it allows (`stellar_object.supported_evolutionary_scales`, GEN.942)."""


def _star_age_sentence(star):
    """Two short sentences: the star's age and lifespan, then what that means for life (the slowest pace its lifespan allows)."""
    age = _gy(star.age)
    if star.yerkes_class in pre_main_sequence.STAGE_CLASSES:
        # GEN.242: still contracting.
        sentence = inline("The star is about ", age, " old and is still contracting towards the main sequence. ",
                          "It should last about ", _gy(star.lifespan), ".")
    elif star.lifespan == float("inf"):
        sentence = inline("The star is about ", age, " old and is now a white dwarf, which will cool for trillions of years.")
    else:
        sentence = inline("The star is about ", age, " old and should last about ", _gy(star.lifespan), ".")
    # What the lifespan allows is told of a main-sequence star only: an evolved star's letter and lifespan are not those of its main sequence.
    if star.yerkes_class in ("V", "VI"):
        scales = supported_evolutionary_scales(star)
        if scales:
            sentence = inline(sentence, " ", LIFESPAN_SENTENCES[scales[-1]])
    return sentence


def _galactic_orbit_or_centre(star):
    if getattr(star, "ejecta_velocity_kms", None) is not None:
        speed = math.sqrt(sum(component * component for component in star.ejecta_velocity_kms))
        return inline("Ejected from the galaxy, leaving at ", _kms(speed), " on a straight line")
    if star.galactic_orbital_period_gy is None:
        return inline("None (the galaxy's centre)")
    return _galactic_orbit(star.galactic_orbital_speed_kms, star.galactic_orbital_period_gy)


def _wobble_row(star):
    if star.reflex_offset_x or star.reflex_offset_y or star.reflex_offset_z:
        offset = _offset(star.reflex_offset_x, star.reflex_offset_y, star.reflex_offset_z)
        return [fact("wobble", "Planetary wobble", _au(offset), " from its nominal position, pulled by its own planets")]
    return []


def _black_hole_blocks(hole):
    rows = [
        fact("type", "Type", hole.type),
        fact("mass", "Mass", _mass(hole.mass, solar=True)),
        fact("event_horizon", "Event horizon radius", _radius_km(hole.radius)),
        fact("spin", "Spin parameter", qty("ratio", hole.spin, note="of the maximum (0 is not spinning)")),
        fact("disk", "Accretion disk", "Active accretion disk detected" if hole.has_accretion_disk else "None detected"),
        fact("orbit", "Galactic orbit", _galactic_orbit_or_centre(hole)),
        *_wobble_row(hole),
        fact("loc", "Location", hole.name),
    ]
    glow = _kelvin(hole.temperature, companions=False)
    if hole.mass_class == "supermassive":
        sentences = [inline("It is the supermassive black hole at the centre of the galaxy, and every star in the galaxy "
                            "orbits it in the end."),
                     inline("Its nucleus is quiet: a thin, hot accretion flow feeds it at a tiny fraction of its Eddington "
                            "limit and glows at about ", glow, ".")]
    else:
        if hole.mass_class == "intermediate":
            origin = inline("It is an intermediate-mass black hole, heavier than any one star could leave behind, probably "
                            "grown in a dense star cluster ", _gy(hole.age), " ago.")
        else:
            origin = inline("It is the collapsed core of a massive star, formed in a supernova ", _gy(hole.age), " ago.")
        disk = (inline("A faint accretion disk of infalling matter surrounds it, its inner edge glowing at about ", glow, ".")
                if hole.has_accretion_disk else
                inline("No companion star or debris disk remains, so only its gravity shows it is there."))
        sentences = [origin, inline("Nothing, not even light, escapes from inside its event horizon."), disk]
    return [facts("Black hole", rows), prose("summary", *sentences)]


def _neutron_star_blocks(star):
    spin = _held("ms", star.spin_period_ms, shown="time")
    rows = [
        fact("type", "Type", star.type),
        fact("mass", "Mass", _mass(star.mass, solar=True)),
        fact("radius", "Radius", _radius_km(star.radius)),
        fact("spin_period", "Spin period", spin),
        fact("magnetic_field", "Magnetic field", _held("G", star.magnetic_field_gauss, shown="magnetic field")),
        fact("surface_temp", "Surface temperature", _kelvin(star.surface_temperature_k, companions=False)),
        fact("orbit", "Galactic orbit", _galactic_orbit_or_centre(star)),
        *_wobble_row(star),
        fact("loc", "Location", star.name),
    ]
    sentences = [inline("It is a neutron star, the collapsed core of a massive star left by a supernova ", _gy(star.age),
                        " ago. It holds more mass than the Sun in a sphere ", _radius_km(star.radius), " in radius.")]
    if star.pulsar_type == "magnetar":
        sentences.append(inline("It is a magnetar, one of the rarest objects in the galaxy. Its magnetic field of ",
                                _held("G", star.magnetic_field_gauss, shown="magnetic field"), " is among the strongest known, it "
                                "turns once every ", spin, ", and it is young enough to erupt in giant flares."))
    elif star.pulsar_type != "non-pulsing":
        sentences.append(inline(f"It is a {star.pulsar_type} pulsar and sweeps a beam of radiation past any observer once "
                                "every ", spin, "."))
    else:
        sentences.append(inline("Its beam no longer sweeps across this vantage point."))
    return [facts("Neutron star", rows), prose("summary", *sentences)]


def _star_section(star, children=()):
    if star.yerkes_class == "BH":
        return section("black-hole", _ref("star", star), star.name, _black_hole_blocks(star), children)
    if star.yerkes_class == "NS":
        return section("neutron-star", _ref("star", star), star.name, _neutron_star_blocks(star), children)
    return section("star", _ref("star", star), star.name, [_star_facts(star), prose("summary", _star_age_sentence(star))],
                   children)


def _close_pair_section(pair):
    primary, secondary = pair._primary, pair._secondary
    offsets = []
    for member, parts in ((primary, (pair.binary_primary_position_x, pair.binary_primary_position_y,
                                     pair.binary_primary_position_z)),
                          (secondary, (pair.binary_secondary_position_x, pair.binary_secondary_position_y,
                                       pair.binary_secondary_position_z))):
        offsets.append(inline(f"{member.name}: ", _au(_offset(*parts))))
    rows = [
        fact("type", "Type", pair.type),
        fact("mass", "Mass", _mass(pair.mass)),
        fact("lum", "Luminosity", _luminosity(pair.luminosity)),
        fact("hab", "Habitable zone", "Between ", _au(pair.habitable_zone[0]), " and ", _au(pair.habitable_zone[1])),
        fact("separation", "Stellar separation", _au(pair.binary_separation_au)),
        fact("mutual_orbit", "Mutual orbit", _kms(pair.binary_mutual_orbital_speed_kms), " (",
             _years(pair.binary_mutual_orbital_period_years), " per orbit)"),
        fact("wobble", "Barycenter offset", offsets[0], ", ", offsets[1], " from the barycenter"),
        fact("orbit", "Galactic orbit", _galactic_orbit(pair.galactic_orbital_speed_kms, pair.galactic_orbital_period_gy)),
        fact("loc", "Location", f"{close_pair_label(primary.name, secondary.name)} Binary System"),
    ]
    age = inline("The binary system is about ", _gy(pair.age), " old.")
    return section("pair", None, f"{close_pair_label(primary.name, secondary.name)} Binary System",
                   [facts("Binary system", rows), prose("summary", age)])


def _wide_pair_section(pair):
    primary, secondary = pair.primary, pair.secondary
    offsets = [inline(f"{member.name}: ", _au(_offset(*parts))) for member, parts in (
        (primary, (pair.primary_position_x_au, pair.primary_position_y_au, pair.primary_position_z_au)),
        (secondary, (pair.secondary_position_x_au, pair.secondary_position_y_au, pair.secondary_position_z_au)))]
    rows = [
        fact("separation", "Stellar separation", _au(pair.separation_au), ", ranging from ", _au(pair.periapsis_au),
             " at periapsis to ", _au(pair.apoapsis_au), " at apoapsis"),
        fact("eccentricity", "Orbital eccentricity", qty("ratio", pair.eccentricity, note=ECCENTRICITY_NOTE)),
        fact("mutual_orbit", "Mutual orbit", _kms(pair.speed_kms), " (", _years(pair.period_years), " per orbit)"),
        fact("wobble", "Barycenter offset", offsets[0], ", ", offsets[1], " from the barycenter"),
        fact("primary_limit", "Planetary limit of the primary", f"{primary.name}'s planetary system is stable out to ",
             _au(pair.primary_a_crit_au), f" due to {secondary.name}'s gravity."),
        fact("secondary_limit", "Planetary limit of the secondary", f"{secondary.name}'s planetary system is stable out to ",
             _au(pair.secondary_a_crit_au), f" due to {primary.name}'s gravity."),
    ]
    return section("pair", None, f"{primary.name} and {secondary.name}", [facts("Wide binary system", rows)])


# --- planets, moons, belts, comets ----------------------------------------------------------------------

def _gas_giant_depth_sentences(planet):
    depth = planet.gas_giant_depth()
    if depth is None:
        return []
    km = depth["depth_km"]
    if km > 0:
        first = inline("That is ", _km(km), " below the 1 atmosphere level, where pressure roughly doubles every ",
                       _km(depth["doubling_km"]), " of depth.")
    else:
        first = inline("That is ", _km(abs(km)), " above the 1 atmosphere level, where pressure roughly doubles ",
                       "every ", _km(depth["doubling_km"]), " of depth.")
    return [
        first,
        inline("Around ", _km(depth["transition_depth_km"]), " down, the gas is squeezed into liquid metal at about ",
               _kelvin(depth["transition_temperature_k"]), "."),
        inline("The pressure at the core is of the order of ", _pascal(depth["core_pressure_pa"]), "."),
    ]


def _atmosphere_composition(planet):
    """The air as shares of its partial pressures, with the stored phrase; `None` for an airless body."""
    if planet.atmosphere in (None, "None"):
        return None
    air = atmosphere_physics.gas_partials_pa(planet)
    total = sum(air.values())
    parts = []
    if total > 0:
        shares = [(gases.gas(formula).name, value) for formula, value in air.items() if value]
        for name, value in sorted(shares, key=lambda item: -item[1]):
            parts.append({"name": name, "share": value / total})
    return model.composition("Atmosphere", parts, planet.atmosphere.lower())


def _timeline_sentences(life):
    """The system's age, from the `life_facts.LifeFacts` alone (the "timeline" block)."""
    return [inline("The system is about ", _gy(life.timeline.system_age_gyr), " old.")]


def _milestone_sentences(life):
    """The latest milestone the world's life has reached (the "milestone" block; a generated world states what is there)."""
    timeline = life.timeline
    if timeline.stage is None:
        return [inline("No milestone has been reached yet at the system's age.")]
    return [inline("Its latest milestone is ", timeline.stage_name, ", reached when the system was ", _gy(timeline.milestone_age_gyr),
                   " old. ", timeline.stage_text)]


def _colours(names):
    """Spectrum names as a plain list: "None (Charcoal-black)" is just "charcoal-black", and an acronym like UV keeps its case."""
    words = []
    for name in names:
        if name.lower().startswith("none (") and name.endswith(")"):
            name = name[6:-1]
        words.append(name[0].lower() + name[1:] if not name[:2].isupper() else name)
    return ", ".join(words[:-1]) + " and " + words[-1] if len(words) > 1 else words[0]


def _life_class_sentences(life):
    """The first sentence of the life paragraph, from the `life_facts.LifeFacts` alone, in the present tense: life that has reached a
    milestone lives on the world; a world that suits a life class and has no life yet says so."""
    if life is None:
        return []
    name = life.class_name[0].lower() + life.class_name[1:]
    if life.timeline is None or life.timeline.stage is not None:
        pace = "" if life.pace is None else f", evolving at a {life.pace} pace"
        return [f"Class {life.class_id} life, the {name}, lives on this {life.noun}{pace}."]
    pace = "" if life.pace is None else f", which evolves at a {life.pace} pace"
    return [f"Conditions on this {life.noun} suit Class {life.class_id} life, the {name}{pace}, but none has appeared yet."]


def _life_sentences(life):
    """What the life looks like, from the `life_facts.LifeFacts` alone: its colors in visible light and outside it."""
    if life is None:
        return []
    lives = life.timeline is None or life.timeline.stage is not None
    subject = "that life" if lives else "such life"
    sentences = []
    if life.visible:
        sentences.append(f"In visible light {subject} looks {_colours(life.visible)}.")
    if life.non_visible:
        sentences.append(f"Outside visible light {'it' if lives else 'such life'} shows {_colours(life.non_visible)}.")
    return sentences


def _planet_facts(planet):
    rows = [
        fact("class", "Class", planet.planet_class),
        fact("distance", "Distance", _au(planet.distance)),
        fact("period", "Period", _years(planet.period)),
        fact("speed", "Speed", _kms(planet.orbital_speed_kms)),
        fact("eccentricity", "Eccentricity", qty("ratio", planet.eccentricity, note=ECCENTRICITY_NOTE)),
        fact("radius", "Radius", _radius_km(planet.radius)),
        fact("gravity", "Gravity", _gravity(planet.gravity)),
        fact("mass", "Mass", _mass(planet.mass)),
    ]
    if planet.reflex_offset_x or planet.reflex_offset_y or planet.reflex_offset_z:
        rows.append(fact("moon_wobble", "Moon Wobble",
                         _au(_offset(planet.reflex_offset_x, planet.reflex_offset_y, planet.reflex_offset_z)),
                         " from its nominal position, pulled by its own moons"))
    if not planet.is_moon:
        rows.append(fact("moon_count", "Moons", qty("count", len(planet.moons))))
    if planet.body_type == "g":
        rows.append(fact("temperature", "Temperature", _kelvin(planet.surface_temperature)))
        rows.append(fact("internal_pressure", "Internal Pressure", _pascal(planet.atmospheric_pressure)))
    else:
        rows.append(fact("temperature", "Surface Temperature", _kelvin(planet.surface_temperature)))
        airless = planet.atmosphere in (None, "None")
        rows.append(fact("pressure", "Surface Pressure", _pascal(0.0 if airless else planet.atmospheric_pressure)))
        colour = None if airless else atmosphere_colour.body_atmosphere_colour(planet)
        if colour is not None:
            rows.append(fact("atmosphere_colour", "Atmosphere Colour", swatch(colour.hex, colour.label)))
    return facts("Moon" if planet.is_moon else "Planet", rows)


def _ratings_block(planet):
    """The six habitability categories of a scored rocky body (from the same world the stored scores come from)."""
    scored = planet.habitability_world()
    if scored is None:
        return None
    regions = {field: getattr(planet, field, None) for field in regional_temperature.FIELDS}
    result = protection.rating(*scored, regions=regions)
    rows = [{"domain": row["domain"], "tier": row["tier"], "score": row["score"], "text": inline([row["text"]])}
            for row in result["rows"]]
    return model.ratings({"name": "PHI-6", "score": result["overall"], "tier": habitability.tier(result["overall"])},
                         rows, model.chem_text(result["gear"]))


def _radiation_block(planet):
    """GEN.445: the terms behind the Radiation rating of a scored rocky body, each with the composition that sets it: the dose from
    cosmic rays (the air above and the magnetic field), stellar particles, the ground (its potassium, thorium and uranium) and
    the ultraviolet (the ozone). None for a gas giant or a body not scored."""
    if planet.body_type == "g" or getattr(planet, "surface_dose_msv_yr", None) is None:
        return None
    terms = habitability_explain.radiation_terms(planet)

    def dose(msv_yr):
        return qty("dose", msv_yr / 1000.0, "Sv/yr")

    air = ", ".join(f"{name} {three_figures_or_tiny(share * 100.0)} percent" for name, share in terms.air)
    shield = (f" under {three_figures_or_tiny(terms.column_g_cm2)} g/cm2 of air ({air})" if terms.column_g_cm2 > 0.0 and air
              else f" under {three_figures_or_tiny(terms.column_g_cm2)} g/cm2 of air" if terms.column_g_cm2 > 0.0
              else " with no air to stop them")
    field = (f" and a {three_figures_or_tiny(terms.cutoff_gv)} GV magnetic cutoff" if terms.cutoff_gv > 0.0
             else " and no magnetic field")
    ground = ""
    if terms.potassium_ppm is not None:
        ground = (f": the gamma rays and radon of {three_figures_or_tiny(terms.potassium_ppm)} ppm potassium, "
                  f"{three_figures_or_tiny(terms.thorium_ppm)} ppm thorium and {three_figures_or_tiny(terms.uranium_ppm)} ppm uranium "
                  f"in the rock at formation")
    ozone = ("" if terms.ozone_du is None else f", under {three_figures_or_tiny(terms.ozone_du)} Dobson units of ozone"
             if terms.ozone_du > 0.0 else ", with no ozone")
    rows = [
        fact("radiation-dose", "Total dose", dose(terms.total_msv_yr)),
        fact("radiation-limit", "Sets the rating", f"The {terms.limit}, {terms.tier}"),
        fact("radiation-cosmic", "Cosmic rays", dose(terms.gcr_msv_yr), shield + field),
        fact("radiation-stellar", "Stellar particles", dose(terms.sep_msv_yr), ", from the star's flares, cut by the same air and field"),
        fact("radiation-ground", "Ground", dose(terms.ground_msv_yr), ground),
        *([fact("radiation-belts", "Parent planet's radiation belts", dose(terms.belt_msv_yr),
                ", trapped particles of the parent's field at this moon's orbit, cut by its air")] if terms.belt_msv_yr > 0.0 else []),
        fact("radiation-uv", "Ultraviolet", qty("ratio", terms.uv_index, note="times Earth's DNA-damaging dose"), ozone),
    ]
    return facts("Radiation", rows)


SLOT_LABELS = {"body": "Body protection", "filter": "Air filter", "supply": "Air supply", "radiation": "Radiation shelter"}
"""dict: The label of each slot of the survival equipment in a body's details."""


def _equipment_block(planet):
    """UX.163: the survival equipment of a scored rocky body as a facts table, one row per item with how badly it is
    needed, the PHI-6 category it comes from and the tech level to make it, then the shelter's shielding and the hours
    outside a year when there are any. None for a gas giant or a body not scored."""
    scored = planet.habitability_world()
    if scored is None:
        return None
    levels = equipment.compute(scored[0])
    rows = []
    for index, item in enumerate(equipment.items(levels)):
        tech = equipment.item_tech(item.name, tuning.SPECIES_EQUIPMENT_TECH)
        made = "no tech level makes it" if tech >= 99 else f"made at tech level {tech:g}"
        rows.append(fact(f"equipment-{index}", SLOT_LABELS[item.slot],
                         model.chem_text(item.name[:1].upper() + item.name[1:]),
                         f": {item.severity}; {' and '.join(domain.lower() for domain in item.domains)}; {made}"))
    if levels.get("shelter_column_g_cm2"):
        rows.append(fact("shelter-shielding", "Shelter shielding", f"{round(levels['shelter_column_g_cm2']):,} g/cm2"))
    if levels.get("outside_hours") is not None:
        rows.append(fact("outside-hours", "Hours outside a year", f"{round(levels['outside_hours']):,} hours"))
    day, limit = getattr(planet, "outside_hours_day", None), getattr(planet, "outside_hours_limit", None)
    if day is not None:
        cause = outside_hours.LIMIT_WORDS.get(limit, limit)
        rows.append(fact("outside-hours-day", "Hours outside in a day", f"{day:.1f} hours", f", set by {cause}"))
    return facts("Survival equipment", rows) if rows else None


def _planet_section(planet):
    blocks = [_planet_facts(planet)]
    summary = [f"It is {planet.description.lower()}."]
    blocks.append(prose("description", *summary))
    if planet.body_type == "g":
        blocks.append(prose("depth", *_gas_giant_depth_sentences(planet)))
        interior = giant_interior.interior_of(planet)
        if interior is not None:
            blocks.append(prose("interior", giant_interior.summarize_interior(interior), *giant_interior.describe_interior(interior)))
    else:
        blocks.append(_atmosphere_composition(planet))
    blocks.append(model.composition("Composition", summary=planet.composition.lower()))
    blocks.append(_ratings_block(planet))
    blocks.append(_radiation_block(planet))
    blocks.append(_equipment_block(planet))
    life = life_facts.resolve(planet)
    blocks.append(prose("life-class", *_life_class_sentences(life)))
    blocks.append(prose("life", *_life_sentences(life)))
    if life is not None and life.timeline is not None:
        blocks.append(prose("timeline", *_timeline_sentences(life)))
        blocks.append(prose("milestone", *_milestone_sentences(life)))
    if planet.flavor_text:
        blocks.append(prose("flavour", f"Sensors show {planet.flavor_text}"))
    kind = _noun(planet)
    children = [_planet_section(moon) for moon in planet.moons]
    return section(kind, _ref(kind, planet), planet.name, [b for b in blocks if b and (b["type"] != "prose" or b["text"])],
                   children)


def meteorite_group_fact(body):
    """The row naming the chondrite group a belt's or field's mean composition is nearest (`small_body_composition.chondrite_group`), or `None` for one stored without it."""
    group = body.chondrite_group
    if group is None:
        return None
    return fact("meteorite-group", "Nearest Meteorite Group", "None fits (not a chondrite mix)" if group == meteorite_groups.NONE else f"{group} chondrites")


def _belt_section(belt):
    rows = [
        fact("inner", "Inner Edge", _au(belt.lower_limit)),
        fact("outer", "Outer Edge", _au(belt.upper_limit)),
        fact("density", "Density", belt.density),
        fact("mass", "Mass", _mass(belt.mass_kg)),
    ]
    rows.append(meteorite_group_fact(belt))
    parts = [{"name": component, "level": level} for component, level in belt.composition]
    return section("belt", _ref("belt", belt), "Asteroid Belt",
                   [facts("Asteroid belt", rows),
                    model.composition("Composition", parts, format_composition_summary(belt.composition))])


def _comet_section(comet):
    parabolic = comet.orbit_type == "parabolic"
    rows = [
        fact("type", "Type", "Parabolic Comet" if parabolic else "Comet"),
        fact("nucleus", "Nucleus Diameter", _km(comet.nucleus_diameter_km)),
        fact("perihelion", "Perihelion", _au(comet.perihelion_distance_au)),
        fact("eccentricity", "Eccentricity", qty("ratio", comet.eccentricity, note=ECCENTRICITY_NOTE)),
    ]
    if not parabolic:
        rows.append(fact("period", "Period", _years(comet.orbital_period_years)))
    rows += [
        fact("distance", "Distance from the Star", _au(comet.distance_au)),
        fact("speed", "Speed", _kms(comet.orbital_speed_kms)),
        fact("activity", "Activity", "Active" if comet.is_active else "Dormant"),
    ]
    if parabolic:
        orbit = ("It is on a near-parabolic orbit that is only just unbound. It began in this system's outer reaches "
                 "and will not return.")
    else:
        orbit = (f"It follows a {PERIOD_CLASS_LABELS.get(comet.period_class, 'periodic')} elliptical orbit around the "
                 "star and returns to the inner system every orbit.")
    activity = ("It has an active coma and tail as its surface ices sublimate." if comet.is_active
                else "It is inactive, its surface ices dormant.")
    return section("comet", _ref("comet", comet), comet.name,
                   [facts("Comet", rows), prose("description", orbit, activity),
                    model.composition("Composition", summary=comet.get_composition_summary())])


def _body_section(body):
    return _belt_section(body) if body.body_type == "a" else _planet_section(body)


# --- the system ------------------------------------------------------------------------------------------

def _habitable_count(system):
    """How many rocky planets and moons are habitable: an overall PHI-6 rating that `habitability.is_habitable` accepts
    (GEN.643)."""
    count = 0
    for planet in system.all_planets:
        if planet.body_type == "a":
            continue
        for body in (planet, *planet.moons):
            scored = body.habitability_world()
            if scored is not None and habitability.is_habitable(protection.rating(*scored)["overall"]):
                count += 1
    return count


def _hierarchy_sentences(system):
    """What the system's stars are to one another (a triple's close pair and the star that orbits it), as short sentences, or `[]` for a single
    star or a system stored without its tree of orbits."""
    slots = system.star_slots
    if not system.orbit_nodes or not slots or len(slots) != len(system.stars):
        return []
    return hierarchy_sentences(
        [{"name": star.name, "type": star.type, "node_no": node_no, "slot": slot} for star, (node_no, slot) in zip(system.stars, slots)],
        [{"node_no": node.node_no, "parent_node_no": node.parent_node_no, "slot": node.slot, "kind": node.kind}
         for node in system.orbit_nodes]) or []


def _overview_sentences(system):
    sentences = []
    sentences += _hierarchy_sentences(system)
    segments = []
    if system.planet_count > 0:
        segments.append(f"{system.planet_count} planet{'s' if system.planet_count > 1 else ''}")
    if system.belt_count > 0:
        segments.append(f"{system.belt_count} asteroid belt{'s' if system.belt_count > 1 else ''}")
    if system.comet_count > 0:
        segments.append(f"{system.comet_count} comet{'s' if system.comet_count > 1 else ''}")
    if system.moon_count > 0:
        segments.append(f"{system.moon_count} moon{'s' if system.moon_count > 1 else ''}")
    if segments:
        if len(segments) == 1:
            joined = segments[0]
        elif len(segments) == 2:
            joined = f"{segments[0]} and {segments[1]}"
        else:
            joined = ", ".join(segments[:-1]) + f", and {segments[-1]}"
        sentences.append(f"This system contains {joined}.")
    elif system.planet_count == 0 and system.belt_count == 0 and system.comet_count == 0 and system.moon_count == 0:
        sentences.append("There are no stellar objects in this system.")

    habitable = _habitable_count(system)
    if habitable == 1:
        sentences.append("There is 1 habitable world in the system.")
    elif habitable > 1:
        sentences.append(f"There are {habitable} habitable worlds in the system.")
    else:
        sentences.append("There are no habitable worlds in this system.")

    star = system.star
    cloud = getattr(system, "surrounding_cloud", None)
    wind = star.heliosphere_radius > 0
    compressed = (compressed_heliosphere_radius(star.heliosphere_radius, cloud["density_cm3"], cloud["temperature_k"])
                  if cloud and wind else None)
    if compressed is not None and compressed < star.heliosphere_radius * HELIOSPHERE_COMPRESSION_NOTED:
        sentences.append(inline(
            "The star's wind blows a bubble called the heliosphere. In open space it would reach about ",
            _au(star.heliosphere_radius), f", but the gas of {cloud['name'] or 'the cloud'} around the system presses it in to about ",
            _au(compressed), "."))
    elif wind:
        sentences.append(inline("The star's wind blows a bubble called the heliosphere, which reaches about ",
                                _au(star.heliosphere_radius), "."))
    if wind:
        sentences.append(inline("Beyond the bubble, the star's gravity still holds on out to ",
                                _au(star.system_perimeter), ", the outer edge of the system."))
    else:
        # A black hole or neutron star blows no stellar wind: only its gravity bounds the system.
        sentences.append(inline("The star's gravity holds on out to ", _au(star.system_perimeter), ", the edge of the system."))
    return sentences


def _member_sections(system):
    """The stars beyond the pair, each with its own planets and belts as children (GEN.332)."""
    return [_star_section(star, [_body_section(body) for body in bodies]) for star, bodies in system.member_hosts()]


def describe(system):
    """
    The Document of `system` (a `StarSystem`): one `system` section holding the overview, whose children are
    the stars (or the pair and its stars), planets with their moons, belts and comets in the order the page
    lists them. A wide binary's planets and comets are children of the star they orbit.
    """
    blocks = [prose("overview", *_overview_sentences(system))]
    if system.system_flavor_text:
        blocks.append(prose("flavour", f"Sensors show {system.system_flavor_text}"))

    if system.binary_type == "close":
        children = [_close_pair_section(system.star)]
        children += [_star_section(star) for star in system.stars[:2]] + _member_sections(system)
    elif system.binary_type == "wide":
        children = [_wide_pair_section(system.wide_binary)]
        for star, planets, comets in ((system.primary_star, system.planets, system.comets),
                                      (system.secondary_star, system.secondary_planets, system.secondary_comets)):
            children.append(_star_section(star, [_body_section(body) for body in planets]
                                          + [_comet_section(comet) for comet in comets]))
        children += _member_sections(system)    # GEN.631: the stars beyond the pair
    else:
        children = [_star_section(system.star)]
    if system.binary_type != "wide":
        children += [_body_section(body) for body in system.planets]
        children += [_comet_section(comet) for comet in system.comets]

    return model.document("system", _ref("system", system), system.name,
                          [section("system", _ref("system", system), system.name, blocks, children)])
