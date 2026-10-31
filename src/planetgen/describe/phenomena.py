# planetgen/describe/phenomena.py

"""
The Document of a phenomenon (docs/design/architecture/system-description.md, UX.108): a nebula, a quasar, a rogue planet, an
interstellar comet, an asteroid field, a supernova remnant, or a black hole or neutron star on its own. One section
holds a facts table and a few short sentences; a supernova remnant's compact remnant is its child section.

Like `build.describe`, this holds no Markdown or HTML and no raw number formatting: every number is a
typed quantity the renderers write.
"""

from planetgen import tuning
from planetgen.describe import build, model
from planetgen.describe.model import fact, facts, inline, prose, qty, section
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.generation.belt import format_composition_summary


def _ly(value):
    return build._held("ly", value)


def _years(value):
    return build._years(value)


def _density(per_cm3):
    return build._held("cm-3", per_cm3, shown="number density")


def _article(word):
    return "an" if word[:1].lower() in "aeiou" else "a"


def _orbit(obj):
    """The galactic-orbit sentence every free-floating object shares."""
    return inline("It is bound to no star, and it orbits the galactic centre at ",
                  build._galactic_orbit(obj.galactic_orbital_speed_kms, obj.galactic_orbital_period_gy), ".")


def _orbit_fact(obj):
    return fact("orbit", "Galactic orbit", build._galactic_orbit(obj.galactic_orbital_speed_kms,
                                                                obj.galactic_orbital_period_gy))


def _gas_facts(obj):
    """The class, gas and light rows a nebula and a supernova remnant share."""
    return [
        fact("species", "Dominant Gas", model.chem_text(obj.dominant_species)),
        fact("density", "Density", _density(obj.density_cm3)),
        fact("temperature", "Temperature", build._kelvin(obj.temperature_k, companions=False)),
        fact("extinction", "Extinction", qty("ratio", obj.extinction_av, note="magnitudes")),
    ]


def _nebula(nebula):
    class_name = tuning.NEBULA_CLASSES[nebula.nebula_class]["name"]
    rows = [
        fact("class", "Class", f"{nebula.nebula_class} ({class_name.lower()})"),
        fact("type", "Type", nebula.nebula_type.capitalize()),
        fact("radius", "Radius", _ly(nebula.radius_ly)),
        fact("composition", "Composition", nebula.composition),
        *_gas_facts(nebula),
    ]
    if nebula.hosts:
        rows.append(fact("hosts", "Stars", nebula.hosts))
        rows.append(fact("age", "Age of the Stars", build._held("Myr", nebula.age_myr, shown="time")))
        if nebula.q_total:
            rows.append(fact("ionizing", "Ionizing Photons", qty("count", nebula.q_total), " a second"))
    rows.append(_orbit_fact(nebula))
    sentences = [inline(f"{nebula.name} is {_article(nebula.nebula_type)} {nebula.nebula_type} nebula about ",
                        _ly(nebula.radius_ly), " in radius."),
                 inline(f"It formed via {nebula.formation_cause}."), _orbit(nebula)]
    return section("nebula", None, nebula.name, [facts("Nebula", rows), prose("summary", *sentences)])


def _quasar(quasar):
    horizon = build._radius_km(quasar.event_horizon_radius_km)
    rows = [
        fact("mass", "Black Hole Mass", build._mass_solar(quasar.black_hole_mass_solar)),
        fact("event_horizon", "Event Horizon Radius", horizon),
        fact("luminosity", "Luminosity", build._luminosity(quasar.luminosity_w)),
        fact("eddington", "Share of the Eddington Limit", qty("fraction", quasar.eddington_ratio, "frac")),
        fact("accretion", "Gas Swallowed", build._mass_solar(quasar.accretion_rate_solar_per_year), " a year"),
        fact("broad_line", "Broad-Line Region", build._held("ld", quasar.broad_line_region_light_days), " out"),
        fact("age", "Active For", _years(quasar.active_age_years)),
        fact("radio", "Radio", "Loud" if quasar.is_radio_loud else "Quiet"),
    ]
    if quasar.is_radio_loud:
        rows.append(fact("jets", "Jet Reach", _ly(quasar.jet_length_ly), " end to end"))
    sentences = [
        inline(f"{quasar.name} is the galaxy's active nucleus: a supermassive black hole swallowing gas, with an accretion "
               "disk that outshines every star in a galaxy like the Milky Way by a factor of ",
               qty("count", round(quasar.galaxy_luminosity_multiple)), "."),
        inline("Gas in its broad-line region orbits fast enough to smear its emission lines thousands of kilometres per "
               "second wide."),
        inline("Twin relativistic jets punch out of the galaxy entirely and feed its radio lobes." if quasar.is_radio_loud
               else "It is radio-quiet, with no large-scale jets."),
    ]
    return section("quasar", None, quasar.name, [facts("Quasar", rows), prose("summary", *sentences)])


def _rogue_surface_sentences(planet):
    """What its own heat does to a rogue planet's surface (`RoguePlanet.surface_sentences`, as quantities)."""
    regime = planet.surface_regime
    temperature = build._kelvin(planet.surface_temperature_k, companions=False)
    sentences = [inline("It is ", build._gy(planet.age_gy), " old. With no star to warm it, its only heat is its own: ",
                        qty("flux", planet.internal_heat_flux_w_m2), " leaking from its interior, so it glows at ",
                        build._kelvin(planet.effective_temperature_k, companions=False), ".")]
    if regime in ("gas-giant", "brown-dwarf"):
        sentences.append(inline("It has no solid surface, and at the 1 bar level the temperature is ", temperature, "."))
        return sentences
    if regime == "hydrogen-envelope":
        sentences.append(inline("A thick hydrogen envelope, opaque to heat at high pressure, blankets it, so the ground "
                                "beneath ", build._pascal(planet.surface_pressure_pa), " of gas sits at ", temperature, "."))
        if planet.has_liquid_water and not planet.ice_shell_thickness_km:
            sentences.append(inline("There it holds a liquid ocean about ", build._km(planet.ocean_depth_km), " deep."))
    elif regime == "bare-rock":
        sentences.append(inline("Its bare rock surface has frozen to ", temperature, "."))
    else:
        sentences.append(inline("Its surface has frozen to ", temperature, ", and whatever air it had lies on the ground "
                                "as frost."))
    if planet.ice_shell_thickness_km and planet.has_liquid_water:
        sentences.append(inline("Under an ice shell about ", build._km(planet.ice_shell_thickness_km), " thick, its own "
                                "heat keeps a liquid ocean about ", build._km(planet.ocean_depth_km), " deep."))
    elif planet.ice_shell_thickness_km:
        sentences.append(inline("Its water is frozen solid, an ice layer about ", build._km(planet.ice_shell_thickness_km),
                                " thick down to the rock."))
    return sentences


def _rogue_planet(planet):
    kind = planet.kind_label.lower()
    if getattr(planet, "planet_class", None):
        kind += f", class {planet.planet_class}"
    elif planet.spectral_type:
        kind += f", type {planet.spectral_type}"
    what = f"{planet.spectral_type} brown dwarf" if planet.spectral_type else "planet"
    rows = [
        fact("type", "Type", f"Rogue {kind}"),
        fact("radius", "Radius", build._radius_km(planet.radius_km)),
        fact("composition", "Composition", planet.composition),
        _orbit_fact(planet),
    ]
    sentences = [inline(f"{planet.name} is a free-floating {what} adrift in interstellar space, bound to no star."),
                 *_rogue_surface_sentences(planet)]
    if planet.has_moons:
        sentences.append(inline("A smaller companion body, likely captured after ejection, still orbits it."))
    sentences.append(_orbit(planet))
    return section("rogue-planet", None, planet.name, [facts("Rogue planet", rows), prose("summary", *sentences)])


def _interstellar_comet(comet):
    rows = [
        fact("type", "Type", "Interstellar Comet"),
        fact("nucleus", "Nucleus Diameter", build._km(comet.nucleus_diameter_km)),
        fact("speed", "Excess Speed", build._kms(comet.velocity_kms)),
        fact("activity", "Activity", "Active" if comet.is_active else "Dormant"),
        _orbit_fact(comet),
    ]
    activity = ("It shows an active coma and tail as its surface ices sublimate." if comet.is_active
                else "It shows no activity and passes through inert.")
    sentences = [inline(f"{comet.name} is a small icy body on a hyperbolic, unbound trajectory through interstellar space."),
                 inline("It travels at ", build._kms(comet.velocity_kms), " relative to any star it passes."),
                 inline(activity), _orbit(comet)]
    return section("interstellar-comet", None, comet.name,
                   [facts("Interstellar comet", rows), prose("summary", *sentences),
                    model.composition("Composition", summary=comet.get_composition_summary())])


def _asteroid_field(field):
    family = tuning.ASTEROID_FIELD_COMPOSITIONS[field.composition_family]["description"]
    rows = [
        fact("class", "Class", field.field_class),
        fact("density", "Density", field.density),
        fact("family", "Family", f"{field.composition_family.capitalize()} ({family})"),
        fact("radius", "Radius", _ly(field.radius_ly)),
        build.meteorite_group_fact(field),
        _orbit_fact(field),
    ]
    parts = [{"name": component, "level": level} for component, level in field.composition]
    sentences = [inline(f"{field.name} is a field of asteroid debris drifting in open interstellar space, about ",
                        _ly(field.radius_ly), " in radius."), _orbit(field)]
    return section("asteroid-field", None, field.name,
                   [facts("Asteroid field", rows), prose("summary", *sentences),
                    model.composition("Composition", parts, format_composition_summary(field.composition))])


def _supernova_remnant(remnant):
    class_name = tuning.NEBULA_CLASSES[remnant.remnant_class]["name"]
    rows = [
        fact("class", "Class", f"{remnant.remnant_class} ({class_name.lower()})"),
        fact("morphology", "Morphology", remnant.morphology),
        fact("progenitor", "Progenitor", remnant.progenitor_type),
        fact("age", "Age", _years(remnant.age_years)),
        fact("radius", "Radius", _ly(remnant.radius_ly)),
        *_gas_facts(remnant),
        _orbit_fact(remnant),
    ]
    sentences = [inline(f"{remnant.name} is {_article(remnant.morphology)} {remnant.morphology} supernova remnant, the "
                        f"expanding wreckage of a {remnant.progenitor_type} supernova ", _years(remnant.age_years),
                        " ago. It now reaches ", _ly(remnant.radius_ly), " from its centre."), _orbit(remnant)]
    children = []
    centre = remnant.compact_remnant
    if centre is not None:
        sentences.append(inline(f"At its centre lies the collapsed core of the progenitor star, {centre.name}."))
        children.append(_compact(centre))
    elif getattr(remnant, "compact_remnant_name", None):    # a stored remnant: the centre has a page of its own
        sentences.append(inline(f"At its centre lies the collapsed core of the progenitor star, {remnant.compact_remnant_name}."))
    elif remnant.progenitor_type == "core-collapse":
        sentences.append(inline("Whatever compact remnant the collapse left behind, if any, has not been found in the "
                                "debris."))
    return section("supernova-remnant", None, remnant.name, [facts("Supernova remnant", rows), prose("summary", *sentences)],
                   children)


def _compact(remnant):
    if isinstance(remnant, BlackHole):
        return section("black-hole", None, remnant.name, build._black_hole_blocks(remnant))
    return section("neutron-star", None, remnant.name, build._neutron_star_blocks(remnant))


_BUILDERS = (
    (Nebula, _nebula), (Quasar, _quasar), (RoguePlanet, _rogue_planet), (InterstellarComet, _interstellar_comet),
    (AsteroidField, _asteroid_field), (SupernovaRemnant, _supernova_remnant), (BlackHole, _compact),
    (NeutronStar, _compact),
)


def describe_phenomenon(phenomenon):
    """
    The Document of a generated phenomenon: one section named for it, with a facts table and its sentences. Raises
    `ValueError` for any other kind of object.
    """
    for cls, builder in _BUILDERS:
        if isinstance(phenomenon, cls):
            top = builder(phenomenon)
            return model.document("phenomenon", None, phenomenon.name, [top])
    raise ValueError(f"no description for a {type(phenomenon).__name__}")
