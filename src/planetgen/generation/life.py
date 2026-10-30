# planetgen/generation/life.py

"""
Planet Life Classes and Evolutionary Data
==========================================

This module determines which kind of life a planet or moon supports, which of the classes that fit it lives there, how fast its
biosphere would evolve, how far it has got, and (for habitable, non-moon worlds) its evolutionary timeline narrative.

Life data is intentionally NOT computed during `Planet` construction (see
`physics.planets.generate_planet_properties`) -- every planet and moon is
generated in a "no life data yet" state. Instead, `StarSystem.__init__`
calls `apply_life_data` on every planet and moon in the system, after all
of them have been generated and after `Star.adjust_age_for_planets` has
finalized the star's age, so evolutionary timelines are computed against
the star's final age rather than a provisional pre-adjustment one.

The class switch (GEN.889, docs/design/habitability/life-classes.md): `physics.life_classes` scores the body for each life class from its
stored conditions; the class that lives on it is the supported one with the most net primary production (ties within a factor of 2 are
drawn, weighted by it); the schedule of the five hard steps is drawn at the pace its temperature allows and with the first step its class
needs (`physics.life_timeline`); a class that splits water builds the oxygen its photosynthesis has made by the star's age
(`physics.biogenic_oxygen`); and the complexity of the class on the final air (`physics.life_complexity`) caps the milestones and gives
the animal grade. All the draws of a body are made from a stream of its own, seeded by one draw of the run's (GEN.872), so a change in
the life model no longer re-rolls the names and flavour text of every later body.
"""

from planetgen import tuning
from planetgen.generation import evolution
from planetgen.physics import aerosols, atmosphere, biogenic_oxygen, body_colour, habitability_world, life_air, life_classes, life_complexity, life_supply, life_timeline, planet_classes, radiation, surface_liquid_world
from planetgen.util import draw
from planetgen.util import log

def apply_life_data(planet, keep_stored=False):
    """
    Determines and assigns the life of a planet or moon: the class that lives on it (`life_class`, `life_class_variant`), the classes it
    supports (`life_class_scores`, `life_class_best`), its schedule of the stages of life (`physics.life_timeline`, a draw of five hard steps at
    the pace its temperature allows, the stored `evolution.SCHEDULE_FIELDS`), the oxygen its life has built into the air by now
    (`physics.biogenic_oxygen`, GEN.855), its complexity and animal grade (`physics.life_complexity`) and the evolutionary timeline.

    Called once per planet/moon by `StarSystem.__init__`, after all planets
    and moons in the system have been generated and the star's age has been
    finalized by `Star.adjust_age_for_planets`.

    Every draw is made on a stream of the body's own; the run's stream moves by one draw however many are made (GEN.872), and by none for an edit (`keep_stored`).

    Args:
        planet (Planet): The planet or moon to apply life data to. Its class, schedule, pace, spectra, `life_stage` and the two ages are set in
                         place, and the air takes the oxygen the life has made.
        keep_stored (bool): Keep the class and the schedule the body already has when the conditions still support the class
                            (`validation.reapply_life` after an edit), so a world is not given other life by a change it does not
                            care about (DB.100). Only a class that no longer fits is drawn again, with a schedule of its own.
    """
    if keep_stored:      # an edit makes no draw of the run's: what it must draw again comes from a stream of the body's own
        stream = draw.Stream(f"life-edit:{planet.is_moon}:{planet.mass!r}:{planet.radius!r}:{planet.distance!r}:{planet.life_class_id}")
    else:
        stream = draw.Stream(f"life:{draw.random()!r}")
    with draw.bound(stream):
        _apply(planet, keep_stored)


def _apply(planet, keep_stored):
    _set_description(planet)
    scores = life_classes.score_classes(planet)
    chosen, kept = _choose(planet, scores, keep_stored)
    if chosen is None:
        _lifeless(planet, scores, "no class is supported")
        return
    row = life_classes.describe(chosen.id)
    if kept:
        log.debug(f"Life class: kept {row.id}, still supported")
    else:
        log.choice("Life class", row.id, f"the supported class with the most net primary production ({chosen.npp:.3g} of Earth's); "
                   f"{sum(1 for s in scores.values() if s.score > 0)} classes are supported")
    _set_class(planet, row)
    # The schedule is drawn once; the life sentence, the timeline and the stored columns all read it.
    if not (kept and all(getattr(planet, field, None) is not None for field in evolution.SCHEDULE_FIELDS[:2]) and planet.evolutionary_pace):
        _draw_schedule(planet, row)
    if not _life_has_begun(planet):
        _lifeless(planet, scores, f"class {row.id} is supported but its first step comes after the star's age")
        return
    _build_air(planet, row)
    final = life_classes.score_classes(planet)             # the conditions after the oxygen the life has made
    if not (final.get(row.id) and final[row.id].score > 0):
        # The life has used up what it lived on (the donor or hydrogen it drew down, the air it changed): the class that lives here is still one
        # the world supports, at the score it was chosen at, so the stored classes always hold the one that lives (TEST.486).
        final = {**final, row.id: chosen}
    _set_best(planet, final)
    result = final[row.id]
    anoxic = _anoxic_animals(planet, row, keep_stored)
    complexity = life_complexity.assess(row, result.npp, atmosphere.partial_kpa(planet, "o2"), planet.surface_temperature, anoxic)
    _set_evolutionary_state(planet, evolution.get_evolutionary_state(
        planet.star, tuple(getattr(planet, field) for field in evolution.SCHEDULE_FIELDS), max_stage=complexity.top_stage))
    planet.life_complexity = complexity.index
    animals = planet.life_stage in ("multicellularity", "technological_civilization")
    planet.life_animal_grade = complexity.grade if animals else 0
    log.choice("Life complexity", f"{complexity.index:.3g} x Earth's net primary production",
               f"oxidant yield {complexity.oxidant_yield:.2f} of 1, temperature rate factor {complexity.rate_factor:.2f}, "
               f"animal grade {complexity.grade} of 3, highest stage {complexity.top_stage}")
    habitability_world.rescore(planet)  # GEN.207: Suitability reads the life that was just drawn


def _lifeless(planet, scores, why):
    """A world with no life: no class, no schedule, no stage, no complexity, and the score of a lifeless world redone once everything is final
    (GEN.217). `scores` are still stored as the classes the conditions support."""
    _set_best(planet, scores)
    _set_class(planet, None)
    _set_schedule(planet, None, None)
    _set_evolutionary_state(planet, None)
    planet.life_animal_grade = planet.life_complexity = None
    habitability_world.rescore(planet)
    log.debug(f"Life class: none, {why}")


def _life_has_begun(planet):
    """Whether the first step of the body's schedule (abiogenesis) is at or before its star's age. A world that has not reached it has no life
    yet, whatever its conditions support, so it is given no class (GEN.1001); a forced civilization is always reached."""
    first = evolution.get_evolutionary_state(planet.star, tuple(getattr(planet, field) for field in evolution.SCHEDULE_FIELDS),
                                             max_stage=evolution.MILESTONE_KEYS[0])
    return first["stage"] is not None


def _choose(planet, scores, keep_stored):
    """`(ClassScore or None, kept)`: the class that lives on the body. With `keep_stored` the class it has is kept while it is still supported; else
    the dominant one (`life_classes.dominant`). None whatever the scores when the star is younger than `tuning.LIFE_MIN_STAR_AGE_GY`."""
    age = getattr(getattr(planet, "star", None), "age", None)
    if age is not None and age < tuning.LIFE_MIN_STAR_AGE_GY:
        return None, False      # a crust and oceans take about 0.1 Gyr: a star younger than that has no life anywhere (GEN.997)
    if keep_stored:
        stored = life_classes.class_id(getattr(planet, "life_class", None), getattr(planet, "life_class_variant", None))
        if stored in scores and scores[stored].score > 0:
            return scores[stored], True
    return life_classes.dominant(scores, draw), False


def _anoxic_animals(planet, row, keep_stored):
    """Whether a world of an anaerobic class with a donor supply gets animals that have lost aerobic respiration (route R4, GEN.884): a draw of
    `tuning.LIFE_COMPLEXITY["anoxic_animal_chance"]` on the worlds with the donors; an edit keeps the answer the world has."""
    if keep_stored and getattr(planet, "life_animal_grade", None) is not None:
        return bool(planet.life_animal_grade)
    if row.oxidant == "o2" or life_supply.supplies(habitability_world.values_from_body(planet))["donor_supply"] < tuning.LIFE_COMPLEXITY["anoxic_animal_donor_min"]:
        return False
    return draw.random() < tuning.LIFE_COMPLEXITY["anoxic_animal_chance"]


def _set_best(planet, scores):
    """Stores the classes the body supports: the best one (the largest likelihood), its score and the score of every class at or above the line."""
    top = life_classes.best(scores)
    setattr(planet, life_classes.SCORES_FIELD, {class_id: result.score for class_id, result in scores.items() if result.score > 0} or None)
    planet.life_class_best = top.id if top else None
    planet.life_class_best_score = (top.score if top else 0) if scores else None


def _set_description(planet):
    """The class's own description. It claims no life and no planet, so a moon reads it as it stands; what the index
    allows is said by the life sentences."""
    planet.description = tuning.PLANET_CLASSES[planet.planet_class]["description"]


def _draw_schedule(planet, row):
    """
    Draws the schedule of `planet` of life class `row` (GEN.750): five waiting times at the pace its surface temperature allows
    (`life_timeline.rate_factor`), the first step scaled by the class's (`life.first_step_scale`: anoxygenic photosynthesis comes before the
    oxygenic kind), and the wait for complex cells starting once the class has its oxidant: at once for one that breathes something else, and for
    one that breathes oxygen when the air already holds it (an abiotic oxygen air) or a delay after the first photosynthesis when the class makes
    it and the world can make enough (`biogenic_oxygen.profile_of`). A forced technological civilization draws a timeline that finishes in the
    star's age.
    """
    factor = life_timeline.rate_factor(getattr(planet, "surface_temperature", None))
    needs_oxygen = row.oxidant == "o2"
    at_start = (not needs_oxygen) or atmosphere.partial_kpa(planet, "o2") * 1000.0 >= biogenic_oxygen.EUKARYOTE_MINIMUM_PA
    delay = biogenic_oxygen.profile_of(planet).delay_gyr if needs_oxygen and row.releases_oxygen else None
    waits = None
    if planet.star.system_config.INTELLIGENT_LIFE is True:
        waits = life_timeline.forced_waits(planet.star.age, factor, delay, at_start)
    if waits is None:
        waits = life_timeline.unit_waits()
    waits = (waits[0], waits[1] * row.life["first_step_scale"], *waits[2:])
    ages = life_timeline.schedule(waits, factor, delay, at_start)
    log.choice("Evolutionary schedule", ", ".join("never" if age is None else f"{age:.3g} Gyr" for age in ages),
               f"hard steps at temperature rate factor {factor:.3g}, oxygen "
               f"{'not needed' if not needs_oxygen else 'from the start' if at_start else 'after ' + str(delay) + ' Gyr of photosynthesis'}")
    _set_schedule(planet, ages, life_timeline.pace_of(factor))


def _build_air(planet, row):
    """
    Gives the air what the life of `planet` has put into it by its star's age: the oxygen of a class that splits water (GEN.855) and the
    rest of the class's ledger, the methane and the hydrogen drawdown (`physics.life_air`, GEN.857, GEN.858, GEN.859). When either changes the air,
    everything that reads it is made again once (`_rebuild`).
    """
    if planet.body_type == "g":
        return
    rose = _build_oxygen(planet, row)
    changed = life_air.apply_gases(planet, row, rose)
    if rose or changed:
        _rebuild(planet)
    if changed:
        log.choice("Life gases", f"{atmosphere.partial_kpa(planet, 'ch4'):.3g} kPa of methane",
                   f"class {row.id}: methane and hydrogen of the ledger at the star's age")


def _build_oxygen(planet, row):
    """
    Gives the air the oxygen the life of `planet` has made by its star's age (GEN.855), when its class splits water, its photosynthesis
    is old enough for the delay (`biogenic_oxygen.Profile.delay_gyr`) and the air holds less than the world's balance of source and
    sinks: the air takes that. An air that already holds the oxygen is left as it is. Returns whether the air took oxygen (the rise
    has happened: it burned the reduced gases, the methane among them).
    """
    if not row.releases_oxygen:
        return False
    photosynthesis_age = getattr(planet, evolution.SCHEDULE_FIELDS[1], None)
    if photosynthesis_age is None:
        return False
    profile = biogenic_oxygen.profile_of(planet)
    target = biogenic_oxygen.oxygen_mol(profile, planet.star.age - photosynthesis_age)
    if not biogenic_oxygen.apply(planet, target):
        return False
    log.choice("Biogenic oxygen", f"{atmosphere.partial_kpa(planet, 'o2'):.2f} kPa", f"oxygenic photosynthesis {planet.star.age - photosynthesis_age:.2f} Gyr old, "
               f"delay {profile.delay_gyr:.2f} Gyr")
    return True


def _rebuild(planet):
    """What reads the air, made again after life has changed it: the cloud, haze and dust (on the body's own stream), the ozone and doses (on a
    stream of the body's own, so the run's draws do not move), the colour of the atmosphere, the surface liquid (a different temperature and pressure
    leave a different liquid stable, GEN.1035), the habitability score and the class (derived from the air as it is now, GEN.980)."""
    from planetgen.physics import planets as physics_planets          # the heat flow the score reads; planets imports this chain
    aerosols.generate_aerosols(planet)
    surface_liquid_world.generate(planet)
    with draw.bound(f"life-air:{planet.is_moon}:{planet.mass!r}:{planet.radius!r}:{planet.distance!r}"):
        radiation.generate_dose(planet)
    body_colour.generate(planet)
    habitability_world.generate(planet, physics_planets.internal_heat_flux_w_m2(planet))
    physics_planets.apply_class(planet)               # GEN.980: the class is an output of the finished conditions, and the air has changed them


def _set_class(planet, row):
    """Stores the life class that lives on a body and the colours its life shows (`None`: no class)."""
    planet.life_class = row.number if row else None
    planet.life_class_variant = (row.variant or None) if row else None
    planet.reflection_spectrum_visible = (row.life["reflection_spectrum_visible"] or None) if row else None
    planet.reflection_spectrum_non_visible = (row.life["reflection_spectrum_non_visible"] or None) if row else None


def _set_schedule(planet, ages, pace):
    """Stores a schedule (one age for each `evolution.SCHEDULE_FIELDS`, `None` for a stage never reached) and its pace word (`None`:
    no schedule at all)."""
    planet.evolutionary_pace = pace
    for field, age in zip(evolution.SCHEDULE_FIELDS, ages if ages is not None else (None,) * len(evolution.SCHEDULE_FIELDS)):
        setattr(planet, field, age)


def _set_evolutionary_state(planet, state):
    """Stores a `get_evolutionary_state` on the body (`None`: no timeline; the schedule stays what `apply_life_data`
    drew for the class)."""
    planet.life_stage = state["stage"] if state else None
    planet.life_milestone_age_gyr = state["milestone_age_gyr"] if state else None
    planet.life_system_age_gyr = state["system_age_gyr"] if state else None


def _flavor_lists():
    return {"habitable": tuning.HABITABLE_FLAVOR, "planet": tuning.PLANET_FLAVOR, "orbital": tuning.ORBITAL_FLAVOR}


def flavor_key(category, text):
    """The `category:index` key of the flavour line `text` in a `tuning` flavour list (`habitable`, `planet` or
    `orbital`)."""
    return f"{category}:{_flavor_lists()[category].index(text)}"


def flavor_line(key):
    """The flavour line a `flavor_key` names, or `None` for no key."""
    if key is None:
        return None
    category, index = key.split(":")
    return _flavor_lists()[category][int(index)]


def decide_flavor_text(planet):
    """
    Rolls for, and if selected assigns, this planet's (or moon's) flavor
    text, updating the shared `system_config` flavor bookkeeping in the same
    step. Must be called after `apply_life_data(planet)` for the same
    planet/moon, since the habitable/multicellular-life check below reads
    `planet.life_stage`.

    Called once per planet/moon by `StarSystem.__init__`, in the same pass
    as `apply_life_data`, rather than at render time -- `to_paragraph_list()`
    can be called more than once (e.g. printed and then saved), and rolling
    here instead of there means repeated renders read a decision already
    made rather than re-rolling and double-counting `system_flavor_count`.

    Args:
        planet (Planet): The planet or moon to decide flavor text for. Its
                         `flavor_key`/`flavor_text_count` attributes are set
                         in place; `planet.system_config.system_flavor_count`/
                         `recent_flavor_texts` are updated when a flavor is
                         selected.
    """
    if not (draw.random() < tuning.FLAVOR_CHANCE_PLANET
            and planet.system_config.system_flavor_count < tuning.MAX_FLAVOR_TOTAL):
        log.debug("Planet flavor text: none (roll failed FLAVOR_CHANCE_PLANET or MAX_FLAVOR_TOTAL reached)")
        return

    selected_flavor = None
    system_config = planet.system_config

    # Filter out recently used flavor texts
    available_habitable_flavor = [f for f in tuning.HABITABLE_FLAVOR if f not in system_config.recent_flavor_texts]
    available_planet_flavor = [f for f in tuning.PLANET_FLAVOR if f not in system_config.recent_flavor_texts]
    available_orbital_flavor = [f for f in tuning.ORBITAL_FLAVOR if f not in system_config.recent_flavor_texts]

    # If all flavors of a category have been recently used, reset to the full list
    if not available_habitable_flavor:
        available_habitable_flavor = tuning.HABITABLE_FLAVOR
    if not available_planet_flavor:
        available_planet_flavor = tuning.PLANET_FLAVOR
    if not available_orbital_flavor:
        available_orbital_flavor = tuning.ORBITAL_FLAVOR

    # Check for habitable and multicellular/technological life
    is_habitable = planet.life_class is not None
    has_multicellular_life = planet.life_stage in ("multicellularity", "technological_civilization")

    if is_habitable and has_multicellular_life:
        selected_flavor = draw.choice(available_habitable_flavor)
        flavor_category = "habitable (multicellular/technological life present)"
        flavor_list = "habitable"
    elif planet.body_type == "t" and planet_classes.classify(planet) != "A":
        selected_flavor = draw.choice(available_planet_flavor)
        flavor_category = "planet (terrestrial body, not class A)"
        flavor_list = "planet"
    elif planet.body_type == "g" or planet_classes.classify(planet) == "A":
        selected_flavor = draw.choice(available_orbital_flavor)
        flavor_category = "orbital (gas giant or class A body)"
        flavor_list = "orbital"
    else:
        flavor_category = None

    if selected_flavor:
        log.choice("Planet flavor text", selected_flavor, f"category: {flavor_category}")
        planet.flavor_key = flavor_key(flavor_list, selected_flavor)
        system_config.system_flavor_count += 1
        planet.flavor_text_count += 1

        # Add the selected flavor text to the recent list
        system_config.recent_flavor_texts.append(selected_flavor)
        # Keep the recent_flavor_texts list to a maximum size
        if len(system_config.recent_flavor_texts) > tuning.MAX_RECENT_FLAVOR_TEXTS:
            system_config.recent_flavor_texts.pop(0) # Remove the oldest entry
