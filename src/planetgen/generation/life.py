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
of them have been generated, so evolutionary timelines are computed against
the star's age, which is fixed when the star is made (GEN.1546).

The class switch (GEN.889, docs/design/habitability/life-classes.md): `physics.life_classes` scores the body for each life class from its
stored conditions; the class that lives on it is the supported one with the most net primary production (ties within a factor of 2 are
drawn, weighted by it); the schedule of the five hard steps is drawn at the pace its temperature allows and with the first step its class
needs (`physics.life_timeline`); a class that splits water builds the oxygen its photosynthesis has made by the star's age
(`physics.biogenic_oxygen`); and the complexity of the class on the final air (`physics.life_complexity`) caps the milestones and gives
the animal grade. All the draws of a body are made from a stream of its own, seeded by one draw of the run's (GEN.872), so a change in
the life model no longer re-rolls the names and flavor text of every later body.

The adaptation traits of the life (GEN.876, `physics/adaptation_traits.py`) are drawn first, from the stress the body has and on a stream of their own, because
the class scores read them: a hard world's life is tougher, so it supports classes an easy world's would not. A world that ends up with no life stores none.
"""

import copy
from dataclasses import dataclass

from planetgen import tuning
from planetgen.generation import evolution
from planetgen.physics import adaptation_draw, adaptation_record, adaptation_traits, aerosols, air_record, atmosphere, biogenic_oxygen, body_colour, habitability_world, life_air, life_classes, life_complexity, flavor_tree, life_supply, life_timeline, ocean_oxygen, planet_type, radiation, surface_liquid_world, trace_air
from planetgen.physics import carbide_deposits, constants, hydrothermal_deposits, sedimentary_deposits
from planetgen.util import draw
from planetgen.util import log

CARBON_GATE = "co2_partial_pa"
"""str: The name of the gate of a life class that holds the carbon dioxide its carbon fixers need (`life_classes` rows, GEN.1581); a dominant class that fails it on the air it made is dropped (GEN.1736)."""

def apply_life_data(planet, keep_stored=False):
    """
    Determines and assigns the life of a planet or moon: the class that lives on it (`life_class`, `life_class_variant`), the classes it
    supports (`life_class_scores`, `life_class_best`), its schedule of the stages of life (`physics.life_timeline`, a draw of five hard steps at
    the pace its temperature allows, the stored `evolution.SCHEDULE_FIELDS`), the oxygen its life has built into the air by now
    (`physics.biogenic_oxygen`, GEN.855), its complexity and animal grade (`physics.life_complexity`) and the evolutionary timeline.

    Called once per planet/moon by `StarSystem.__init__`, after all planets
    and moons in the system have been generated and the star's age is
    fixed from the star's own state (GEN.1546).

    Every draw is made on a stream of the body's own; the run's stream moves by one draw however many are made (GEN.872), and by none for an edit (`keep_stored`).

    Args:
        planet (Planet): The planet or moon to apply life data to. Its class, schedule, pace, spectra, `life_stage` and the two ages are set in
                         place, and the air takes the oxygen the life has made.
        keep_stored (bool): Keep the class and the schedule the body already has when the conditions still support the class
                            (`validation.reapply_life` after an edit), so a world is not given other life by a change it does not
                            care about (DB.100). Only a class that no longer fits is drawn again, with a schedule of its own.
    """
    if keep_stored:      # an edit makes no draw of the run's: what it must draw again comes from a stream of the body's own
        key = f"life-edit:{planet.is_moon}:{planet.mass!r}:{planet.radius!r}:{planet.distance!r}:{planet.life_class_id}"
    else:
        key = f"life:{draw.random()!r}"
    with draw.bound(draw.Stream(key)):
        before = _trial_copy(planet)
        barred = set()
        while (frozen := _apply(planet, keep_stored, key, barred)) is not None:
            # The life of the dominant class took the liquid it lives in (GEN.1298): the world is made again without that class, from the body as it was
            barred.add(frozen)
            vars(planet).clear()
            vars(planet).update(vars(_trial_copy(before)))
    _remake_ground_if_thawed(planet, before)
    carbide_deposits.generate(planet)            # GEN.1423: moissanite and cohenite on a carbon-rich body, on streams of their own
    hydrothermal_deposits.generate(planet)       # GEN.1381: the hydrothermal and uranium deposits, which read the air the life pass left, on streams of their own
    sedimentary_deposits.generate(planet)        # GEN.1378: the deposits that read the life timeline, on streams of their own (no draw of the run's)
    from planetgen.physics import planets       # planets imports this chain
    planets.apply_class(planet)                  # GEN.1745 (Boss, 11:25 AM EST 2026-11-15): no class until the evolution is finished; the class is read once, here, from the finished body
    if not (keep_stored and getattr(planet, "air_history", None)):       # an edit keeps the history the body has (its air solve's inputs are not stored, so it cannot be worked again)
        planet.air_history = _history_of(planet)
    planet.life_adaptations = _adaptations_of(planet)        # GEN.1808: what the life had to become to live where it does, read from the finished body


def _adaptations_of(planet):
    """The record of the adaptations the life of `planet` needed (`adaptation_record.record_of`), made from the finished body: the classes that live on it, the age its life has reached, none when no class has a record."""
    living = getattr(planet, life_classes.LIVING_FIELD, None) or []
    classes = [entry["life_class"] for entry in living]
    if not any(class_id in adaptation_record.CLASSES for class_id in classes):
        return None
    star_age = getattr(getattr(planet, "star", None), "age", None)
    began = [age for age in (one.started(planet) for one in (getattr(getattr(planet, "build", None), "air_living", None) or [])) if age is not None]
    age = None if star_age is None or not began else max(star_age - min(began), 0.0)
    return adaptation_record.record_of(planet, classes, age)


def _history_of(planet):
    """The packed record of the air history of `planet` (`air_record.pack`), made from the finished body, or `None` when it gets none (`air_record.eligible`). The life it reads is what `_build_air` kept on the
    build record: the classes, the air as they found it, the air the oxygen's rise left and the age of the rise. It draws only on the streams of the body's own (the events and the extinctions, `air_events`), so the run's draws do not move."""
    build = planet.build
    living = build.air_living or []
    begin = min((age for age in (one.started(planet) for one in living) if age is not None), default=None)
    record = air_record.record_of(planet, living, build.air_found, build.air_risen, build.air_rise_gyr, begin)
    return None if record is None else air_record.pack(record)


def _remake_ground_if_thawed(planet, before):
    """
    GEN.1385: the air the life made can warm the body past the freezing point of the ice its ground was drawn with (a frozen hydrosphere at 240 K that the
    life's gases take to 408 K), and a ground of ice over a surface that hot is no ground. The water and the ground are made again at the temperature
    the life left, on the body's own stream (`planets.refresh_at_final_age`, the same step the life pass runs before it), only when the ground held ice
    that no longer stands: a body whose temperature the life did not move past that point keeps what it had, draw for draw.
    """
    from planetgen.physics import planets, surface_composition
    if planet.body_type in ("g", "a") or (planet.surface_temperature or 0.0) <= (before.surface_temperature or 0.0):
        return
    if surface_composition.mix_of(planet).get(surface_composition.ICE) and planet.surface_temperature >= constants.CELSIUS_ZERO_K:
        planets.refresh_at_final_age(planet)        # GEN.1590: the cover the thaw left; the class is read once, from the finished body (`apply_life_data`)


@dataclass
class _Cand:
    """A life class on a body while its life is settled (GEN.1209): its class, its score before the life, its share of the food production, and the
    five waiting times of its schedule (`waits`, drawn on a stream of the class's own) or its ages (`ages`, kept from the stored body)."""

    row: object
    score: object
    share: float = 1.0
    waits: tuple | None = None
    ages: tuple | None = None


def _apply(planet, keep_stored, key, barred):
    """Settles the life of `planet`; returns the id of the dominant class when its own air took its solvent away (the caller makes the world again
    without it; `barred` holds the classes it has barred), else `None`."""
    planet.life_traits = adaptation_traits.text_of(adaptation_draw.draw(habitability_world.values_from_body(planet), draw.Stream(f"{key}:traits")))   # GEN.876
    scores = {class_id: result for class_id, result in life_classes.score_classes(planet).items() if class_id not in barred}
    chosen, kept = _choose(planet, scores, keep_stored)
    if chosen is None:
        _lifeless(planet, scores, "no class is supported")
        return None
    row = life_classes.describe(chosen.id)
    if kept:
        if log.debug_log_active():
            log.debug(f"Life class: kept {row.id}, still supported.")
    else:
        log.choice("Life class", row.id, f"the supported class with the most net primary production ({chosen.npp:.3g} of Earth's); "
                   f"{sum(1 for s in scores.values() if s.score > 0)} classes are supported")
    _set_class(planet, row)
    headline = _Cand(row, chosen)
    # The schedule is drawn once; the life sentence, the timeline and the stored columns all read it.
    if kept and all(getattr(planet, field, None) is not None for field in evolution.SCHEDULE_FIELDS[:2]) and planet.evolutionary_pace:
        headline.ages = tuple(getattr(planet, field) for field in evolution.SCHEDULE_FIELDS)
    else:
        headline.waits = _draw_waits(planet, row)
    cands = [headline, *_minor_candidates(planet, scores, chosen, keep_stored, key)]
    cands = _redate(planet, cands)
    if not _begun(planet, headline):
        # The class with the most production has not got going by the star's age (GEN.1001), but another that the conditions support has (`_redate` kept
        # only those): the first of them, by production, is the one that lives here and heads the list.
        if len(cands) == 1:
            _lifeless(planet, scores, f"class {row.id} is supported but its first step comes after the star's age")
            return None
        cands = _redate(planet, cands[1:])
        headline = cands[0]
        row, chosen = headline.row, headline.score
        _set_class(planet, row)
        log.choice("Life class", row.id, "the first class that has begun by the star's age; the one with the most production has not")
    cands = _settle(planet, cands, key)
    final = life_classes.score_classes(planet)             # the conditions after the oxygen the life has made
    if not (final.get(row.id) and final[row.id].score > 0):
        if row.id in final and "solvent" in final[row.id].failed:
            # Its own air (the greenhouse gas it drew down) took the liquid it lives in: a class cannot live without its solvent (GEN.1298)
            if log.debug_log_active():
                log.debug(f"Life class: {row.id} dropped, its own air left no {row.solvent['medium']} for it.")
            return row.id
        if row.id in final and CARBON_GATE in final[row.id].failed:
            # Its own air (the oxygen it raised, and the carbon sink asked again at the temperature the oxygen left, `biogenic_oxygen.rebalance_co2`) holds less carbon dioxide
            # than the class needs to fix carbon: a class passed the gate before its air, and a class cannot live on carbon its own air no longer holds (GEN.1736)
            if log.debug_log_active():
                log.debug(f"Life class: {row.id} dropped, its own air left it too little carbon dioxide.")
            return row.id
        # The life has used up what it lived on (the donor or hydrogen it drew down, the air it changed): the class that lives here is still one
        # the world supports, at the score it was chosen at, so the stored classes always hold the one that lives (TEST.486).
        final = {**final, row.id: chosen}
    _set_best(planet, final)
    result = final[row.id]
    anoxic = _anoxic_animals(planet, row, keep_stored)
    ocean_o2 = _ocean_o2_kpa(planet)
    complexity = life_complexity.assess(row, result.npp, atmosphere.partial_kpa(planet, "o2"), _clock_k(planet, row), anoxic, ocean_o2)
    if ocean_o2 is not None:
        log.choice("Sealed ocean oxygen", f"{ocean_o2:.3g} kPa", "the lid's oxidant, the water's radiolysis and the faults' radicals over the reductants of the floor (GEN.1824)")
    _set_evolutionary_state(planet, evolution.get_evolutionary_state(
        planet.star, tuple(getattr(planet, field) for field in evolution.SCHEDULE_FIELDS), max_stage=complexity.top_stage))
    planet.life_complexity = complexity.index
    planet.life_animal_grade = _animal_grade(complexity, planet.life_stage)
    log.choice("Life complexity", f"{complexity.index:.3g} x Earth's net primary production",
               f"oxidant yield {complexity.oxidant_yield:.2f} of 1, temperature rate factor {complexity.rate_factor:.2f}, "
               f"animal grade {complexity.grade} of 3, highest stage {complexity.top_stage}")
    setattr(planet, life_classes.LIVING_FIELD, _living_entries(planet, cands, final, key))
    _trace_gases(planet)                # GEN.524: the trace gases of table m read the life stage just stored
    habitability_world.rescore(planet)  # GEN.207: Suitability reads the life that was just drawn
    kinds = habitability_world.living_kinds(habitability_world.values_from_body(planet))
    for entry in getattr(planet, life_classes.LIVING_FIELD) or ():
        entry["kinds"] = list(kinds.get(entry["life_class"], (0, 0, 0, 0)))


def _minor_candidates(planet, scores, chosen, keep_stored, key):
    """
    The other classes that may live beside the dominant one (GEN.1209): any class the conditions support, whatever its solvent, one for each class
    number (the variants are one class: the best-scoring). Each has a schedule of its own, drawn on a stream keyed by the class, so the dominant class's
    draws are the same whether others live there or not. An edit (`keep_stored`) keeps the minor classes the body has while they are still supported and
    adds none.
    """
    taken = life_classes.describe(chosen.id).number
    found = {}
    if keep_stored:
        stored = {entry["life_class"]: entry for entry in (getattr(planet, life_classes.LIVING_FIELD, None) or [])}
        for class_id, entry in stored.items():
            row = life_classes.BY_ID.get(class_id)
            if row is not None and row.number != taken and scores.get(class_id) is not None and scores[class_id].score > 0 and entry.get("ages"):
                found[row.number] = _Cand(row, scores[class_id], ages=tuple(entry["ages"]))
        return sorted(found.values(), key=lambda c: (-c.score.npp, c.row.id))
    for class_id, result in scores.items():
        row = life_classes.BY_ID[class_id]
        if result.score <= 0 or result.npp <= 0.0 or row.number == taken:
            continue
        if row.number not in found or (result.npp, class_id) > (found[row.number].score.npp, found[row.number].row.id):
            found[row.number] = _Cand(row, result)
    out = sorted(found.values(), key=lambda c: (-c.score.npp, c.row.id))
    for cand in out:
        cand.waits = _draw_waits(planet, cand.row, draw.Stream(f"{key}:class:{cand.row.id}"))
    return out


def _share_out(cands):
    """Each class's share of the body's food production: its net primary production over the sum (one class has all of it)."""
    total = sum(c.score.npp for c in cands)
    for cand in cands:
        cand.share = cand.score.npp / total if total > 0.0 else 1.0 / len(cands)


def _date(planet, cands):
    """The ages of the five steps of each class that has only its waiting times: at the pace the body's temperature allows, the wait for complex cells
    starting once the class has its oxidant: at once for one that breathes something else or when the air already holds oxygen, else when the oxygen
    of all the classes together has risen (`life_air.oxygen_rise_gyr`; never if it does not). A class's own clock runs at where it lives (`_clock_k`: the surface, or the cloud layer's
    temperature for a class of the acid cloud, GEN.1307)."""
    if all(c.ages is not None for c in cands):
        return
    factors = {c.row.id: life_timeline.rate_factor(_clock_k(planet, c.row)) for c in cands}
    living = [life_air.Living(c.row, c.share, (c.ages or life_timeline.schedule(c.waits, factors[c.row.id], 0.0, True))[1]) for c in cands]
    rise = life_air.oxygen_rise_gyr(planet, living)
    for cand, one in zip(cands, living):
        if cand.ages is not None:
            continue
        at_start = _oxidant_at_start(planet, cand.row)
        delay = None if rise is None else max(0.0, rise - one.producers_gyr)
        cand.ages = life_timeline.schedule(cand.waits, factors[cand.row.id], delay, at_start)
        log.choice("Evolutionary schedule", ", ".join("never" if age is None else f"{age:.3g} Gyr" for age in cand.ages),
                   f"class {cand.row.id}, hard steps at temperature rate factor {factors[cand.row.id]:.3g}, oxygen "
                   f"{'not needed' if cand.row.oxidant != 'o2' else 'from the start' if at_start else 'after the rise at ' + str(rise) + ' Gyr'}")


def _oxidant_at_start(planet, row):
    """Whether a class of `row` has its oxidant from the start: it breathes something other than oxygen, or the air already holds enough for
    complex cells (`biogenic_oxygen.EUKARYOTE_MINIMUM_PA`)."""
    return row.oxidant != "o2" or atmosphere.partial_kpa(planet, "o2") * 1000.0 >= biogenic_oxygen.EUKARYOTE_MINIMUM_PA


def _begun(planet, cand):
    """Whether the first step of the class's schedule (abiogenesis) is at or before its star's age. A class that has not reached it has no life yet,
    whatever the conditions support (GEN.1001); a forced civilization is always reached."""
    first = evolution.get_evolutionary_state(planet.star, cand.ages, max_stage=evolution.MILESTONE_KEYS[0])
    return first["stage"] is not None


def _living_of(cands):
    """The classes that live on the body as the air reads them (GEN.1213): `life_air.Living` for each, with its share of the food production and the
    age its producers began (the second step of its schedule)."""
    return [life_air.Living(c.row, c.share, c.ages[1]) for c in cands]


def _settle(planet, cands, key):
    """
    Settles which classes live on the body and gives it the air they make together (GEN.1209, GEN.1213). A minor class takes hold only if its first step is
    inside the star's age (`_begun`) and the air that the classes make does not poison it: the conditions score it above zero on the final air. A class
    the air poisons is dropped, the one with the smallest share first, and the air made again from the rest, on a trial copy of the body, until the classes
    that are left all live in the air they make. The dominant class is never dropped here (it has its own rule in `_apply`). One class makes its air on
    the body directly, as before. Returns the classes that live, with their ages and shares as they stand.
    """
    while len(cands) > 1:
        trial = _trial_copy(planet)
        _air_on(trial, cands, key)
        final = life_classes.score_classes(trial)
        poisoned = [c for c in cands[1:] if not (final.get(c.row.id) and final[c.row.id].score > 0)]
        if not poisoned:
            _air_on(planet, cands, key)
            return cands
        smallest = min(poisoned, key=lambda c: (c.share, c.row.id))
        if log.debug_log_active():
            log.debug(f"Life class: {smallest.row.id} dropped, the air of the other classes poisons it.")
        cands = _redate(planet, [c for c in cands if c is not smallest])
    _store_schedule(planet, cands)
    _build_air(planet, cands)
    return cands


def _air_on(body, cands, key):
    """Stores the dominant class's schedule on `body` and builds the air its classes make, on the body's own air stream."""
    _store_schedule(body, cands)
    with draw.bound(draw.Stream(f"{key}:air")):
        _build_air(body, cands)


def _store_schedule(planet, cands):
    """Stores the dominant class's schedule and its pace on the body (the columns the life sentence, the timeline and the air read). One kept from the
    stored body (an edit) is stored already."""
    if cands[0].waits is not None:
        _set_schedule(planet, cands[0].ages, life_timeline.pace_of(life_timeline.rate_factor(_clock_k(planet, cands[0].row))))


def _redate(planet, cands):
    """The classes with their shares and ages worked out for the classes that are there: the shares of the food production, then the ages (the oxygen rise
    moves with who lives), then a minor class whose first step comes after the star's age is dropped, and again until none is. The dominant class stays."""
    while True:
        _share_out(cands)
        for cand in cands:
            if cand.waits is not None:
                cand.ages = None
        _date(planet, cands)
        kept = [cands[0], *[c for c in cands[1:] if _begun(planet, c)]]
        if len(kept) == len(cands):
            return cands
        cands = kept

def _trial_copy(planet):
    """A copy of `planet` to try an air on: its attributes copied, with the containers (the air, the spectra, the layers) copied too, so building the air on
    it leaves the body as it was. The star and the moons are shared (nothing in the air pass changes them). Its build record is copied the same way."""
    clone = copy.copy(planet)       # the clone has a record of its own (`Planet.__copy__`)
    _copy_containers(planet, clone, skip=("star", "moons"))
    _copy_containers(planet.build, clone.build)
    return clone


def _copy_containers(source, clone, skip=()):
    """Gives `clone` its own copy of every dict, list and set `source` holds, so a pass on the clone leaves the source's as they were."""
    for name, value in vars(source).items():
        if name not in skip and isinstance(value, (dict, list, set)):
            # GEN.1734: a list of bodies (the moons the planet was born with, `build.born_moons`) is shared like `moons`: each body holds the siblings of its host, so a deep copy of
            # one walks every body of the system in one chain (a recursion of 160 deep systems) and copies bodies the air pass never changes
            setattr(clone, name, copy.copy(value) if _holds_bodies(value) else copy.deepcopy(value))


def _holds_bodies(container):
    """Whether a list, set or dict has a planet, moon or other body in it (an object with a `body_type`)."""
    items = container.values() if isinstance(container, dict) else container
    return any(getattr(item, "body_type", None) is not None for item in items)


def _living_entries(planet, cands, final, key):
    """The stored list of the classes that live on the body (`life_living`, the rows of `planet_life`): the dominant class first, then the rest by share.
    Each has its share, its schedule, the stage its schedule reached by the star's age (capped by its complexity on the final air), its animal grade
    and the complexity index."""
    entries = []
    for position, cand in enumerate(cands):
        row = cand.row
        if position == 0:
            stage, milestone, system_age = planet.life_stage, planet.life_milestone_age_gyr, planet.life_system_age_gyr
            grade, complexity_index = planet.life_animal_grade, planet.life_complexity
        else:
            result = final.get(row.id) or cand.score
            anoxic = _anoxic_animals(planet, row, False, draw.Stream(f"{key}:class:{row.id}:anoxic"))
            complexity = life_complexity.assess(row, result.npp, atmosphere.partial_kpa(planet, "o2"), _clock_k(planet, row), anoxic, _ocean_o2_kpa(planet))
            state = evolution.get_evolutionary_state(planet.star, cand.ages, max_stage=complexity.top_stage)
            stage, milestone, system_age = state["stage"], state["milestone_age_gyr"], state["system_age_gyr"]
            grade = _animal_grade(complexity, stage)
            complexity_index = complexity.index
        entries.append({"life_class": row.id, "share": cand.share, "stage": stage, "milestone_age_gyr": milestone, "system_age_gyr": system_age,
                        "animal_grade": grade, "complexity": complexity_index, "ages": list(cand.ages)})
    return entries


def _clock_k(planet, row):
    """The temperature the clock of class `row` runs at on `planet`: the surface's, the water's of an ice-covered ocean for a class that lives in it (GEN.1822), or the acid cloud layer's for a class that lives there (`life_classes.clock_temperature_k`, GEN.1307)."""
    if getattr(row, "habitat", None) is None and not (life_classes.lives_in_habitat_water(row) and getattr(planet, "hydrosphere", None) == "ice-covered ocean"):
        return getattr(planet, "surface_temperature", None)
    return life_classes.clock_temperature_k(row, life_classes.values_from_body(planet))


def _animal_grade(complexity, stage):
    """The animal grade of life at `stage`: the complexity's once the life is multicellular, else 0."""
    return complexity.grade if stage in ("multicellularity", "technological_civilization") else 0


def _lifeless(planet, scores, why):
    """A world with no life: no class, no schedule, no stage, no complexity, and the score of a lifeless world redone once everything is final
    (GEN.217). `scores` are still stored as the classes the conditions support."""
    _set_best(planet, scores)
    _set_class(planet, None)
    setattr(planet, life_classes.LIVING_FIELD, None)
    _set_schedule(planet, None, None)
    _set_evolutionary_state(planet, None)
    planet.life_animal_grade = planet.life_complexity = planet.life_traits = None
    _trace_gases(planet)                # GEN.524: a lifeless world loses the trace gases a drawn life had given it
    habitability_world.rescore(planet)
    log.debug(f"Life class: none, {why}")


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


def _ocean_o2_kpa(planet):
    """The oxygen, kPa, of the ocean of a body under an ice lid (`ocean_oxygen`, GEN.1824), the one its animals breathe; `None` for any other body."""
    return ocean_oxygen.ocean_o2_kpa(habitability_world.values_from_body(planet))


def _anoxic_animals(planet, row, keep_stored, stream=None):
    """Whether a world of an anaerobic class with a donor supply gets animals that have lost aerobic respiration (route R4, GEN.884): a draw of
    `tuning.LIFE_COMPLEXITY["anoxic_animal_chance"]` on the worlds with the donors; an edit keeps the answer the world has."""
    if keep_stored and getattr(planet, "life_animal_grade", None) is not None:
        return bool(planet.life_animal_grade)
    if row.oxidant == "o2" or life_supply.supplies(habitability_world.values_from_body(planet))["donor_supply"] < tuning.LIFE_COMPLEXITY["anoxic_animal_donor_min"]:
        return False
    return (draw if stream is None else stream).random() < tuning.LIFE_COMPLEXITY["anoxic_animal_chance"]


def _set_best(planet, scores):
    """Stores the classes the body supports: the best one (the largest likelihood), its score and the score of every class at or above the line."""
    top = life_classes.best(scores)
    setattr(planet, life_classes.SCORES_FIELD, {class_id: result.score for class_id, result in scores.items() if result.score > 0} or None)
    planet.life_class_best = top.id if top else None
    planet.life_class_best_score = (top.score if top else 0) if scores else None


def _draw_waits(planet, row, stream=None):
    """
    Draws the five waiting times of the schedule of `planet` of life class `row` (GEN.750): one exponential wait for each step, the first scaled by the
    class's (`life.first_step_scale`: anoxygenic photosynthesis comes before the oxygenic kind). `stream` is the class's own (a minor class); the
    dominant class draws on the body's. A request for a technological civilization raises the odds of waits that finish in the star's age (GEN.1130, `life_timeline.raised_waits`), for the
    dominant class, with the delay the class's own oxygen takes.
    """
    factor = life_timeline.rate_factor(_clock_k(planet, row))
    waits = None
    if stream is None and planet.star.system_config.INTELLIGENT_LIFE is True:
        at_start = _oxidant_at_start(planet, row)
        delay = biogenic_oxygen.profile_of(planet, producers=[(1.0, biogenic_oxygen.band_of(row))]).delay_gyr if row.oxidant == "o2" and row.releases_oxygen else None
        waits = life_timeline.raised_waits(planet.star.age, factor, delay, at_start)
    if waits is None:
        waits = life_timeline.unit_waits(stream)
    return (waits[0], waits[1] * row.life["first_step_scale"], *waits[2:])


def _build_air(planet, cands):
    """
    Gives the air what the life of `planet` has put into it by its star's age, summed over the classes that live there (GEN.1213): the oxygen of the
    classes that split water (GEN.855) and the rest of their ledgers, the methane and the hydrogen drawdown (`physics.life_air`, GEN.857, GEN.858,
    GEN.859). When either changes the air, everything that reads it is made again once (`_rebuild`). The oxygen is asked first, on the air the life found; when it is not there and the methane
    that the life then made has changed the air, it is asked again at the temperature that methane left (the oxygen rise of a world its own life thaws, GEN.1745 slice 3).
    """
    if not planet_type.has_ground(planet):      # a giant, or a planet under a deep envelope (GEN.1334): no ground for life to change the air of
        return
    living = _living_of(cands)
    was = (planet.surface_temperature, planet.atmospheric_pressure)
    build = planet.build
    build.air_living, build.air_found, build.air_risen, build.air_rise_gyr = living, biogenic_oxygen.air_moles(planet), None, None       # GEN.1745 slice 4b: what the history of the air reads
    rose = _build_oxygen(planet, living)
    if rose:
        build.air_risen = biogenic_oxygen.air_moles(planet)
    changed = life_air.apply_gases(planet, living, rose)
    if changed and not rose and _build_oxygen(planet, living):
        build.air_risen = biogenic_oxygen.air_moles(planet)
        # GEN.1745 slice 3: the methane the life made warmed the world, and the oxygen law reads the temperature of its sea (the producers' temperature factor); the air that was asked for oxygen when it was
        # cold is asked again at the temperature the life left it, as the history reads the law at each age. The rise burns the methane, so the gases of the rows that go on are put in again.
        rose = True
        life_air.apply_gases(planet, living, True)
    if rose or changed:
        _rebuild(planet, remake_water=(planet.surface_temperature, planet.atmospheric_pressure) != was)
    if changed:
        log.choice("Life gases", f"{atmosphere.partial_kpa(planet, 'ch4'):.3g} kPa of methane",
                   f"class {cands[0].row.id}: methane and hydrogen of the ledger at the star's age")


def _build_oxygen(planet, living):
    """
    Gives the air the oxygen the life of `planet` has made by its star's age (GEN.855, GEN.1213), when a class that splits water has
    producers old enough for the rise (`life_air.oxygen_state`: the delay the summed source of the oxygenic classes takes to fill the store, which holds
    the methane the earlier classes made) and the air holds less than the world's balance of source and sinks: the air takes that. An air that already
    holds the oxygen is left as it is. Returns whether the air took oxygen (the rise has happened: it burned the reduced gases, the methane among them).
    """
    rise, target = life_air.oxygen_state(planet, living)
    if rise is None or not biogenic_oxygen.apply(planet, target):
        return False
    planet.build.air_rise_gyr = rise
    log.choice("Biogenic oxygen", f"{atmosphere.partial_kpa(planet, 'o2'):.2f} kPa", f"oxygenic photosynthesis from the first producers; the oxygen rose at {rise:.2f} Gyr")
    if biogenic_oxygen.rebalance_co2(planet):      # GEN.1309 A: the carbon sink is asked again at the temperature the burnt methane left
        log.choice("Carbon balance", f"{atmosphere.partial_kpa(planet, 'co2'):.3g} kPa of carbon dioxide", f"the sink solved again at {planet.surface_temperature:.0f} K after the oxygen rise")
    return True


def _air_stream_key(planet):
    """The stream the doses are made again on after life changed the air: the body's own, so the run's draws do not move."""
    return f"life-air:{planet.is_moon}:{planet.mass!r}:{planet.radius!r}:{planet.distance!r}"


def _trace_gases(planet):
    """Gives the air its trace gases (`trace_air.apply`, GEN.524) and, when they changed it, makes the doses and the ultraviolet index of the air again: they
    read the air, and the trace gases shift it by parts in a thousand. The same stream as `_rebuild` takes, so the draws do not move."""
    before = atmosphere.gas_partials_pa(planet)
    trace_air.apply(planet)
    if planet_type.has_ground(planet) and getattr(planet, "uv_surface_index", None) is not None and atmosphere.gas_partials_pa(planet) != before:
        with draw.bound(_air_stream_key(planet)):
            radiation.generate_dose(planet)


def _rebuild(planet, remake_water=False):
    """What reads the air, made again after life has changed it: the cloud, haze and dust (on the body's own stream), the ozone and doses (on a
    stream of the body's own, so the run's draws do not move), the colour of the atmosphere, the surface liquid (a different temperature and pressure
    leave a different liquid stable, GEN.1035) and the habitability score. The class is not made here: it is read once, from the finished body (`apply_life_data`, GEN.1745). When the life moved the
    temperature or the pressure (`remake_water`), the water is made again first at the new ones (`planets.refresh_at_final_age`, on the body's own stream):
    the state of the hydrosphere, the ice, the ocean cover and the ground were drawn at the air the life found, and a liquid ocean stored at 432 K under 63 kPa
    is not one (GEN.1590)."""
    from planetgen.physics import planets as physics_planets          # the heat flow the score reads; planets imports this chain
    aerosols.generate_aerosols(planet)
    surface_liquid_world.generate(planet)
    with draw.bound(_air_stream_key(planet)):
        radiation.generate_dose(planet)
    body_colour.generate(planet)
    if remake_water:
        physics_planets.refresh_at_final_age(planet)       # the hydrosphere, the cover and the ground, and the habitability score that reads them
    else:
        habitability_world.generate(planet, physics_planets.internal_heat_flux_w_m2(planet))


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


def flavor_line(key):
    """The flavor line a stored `flavor_key` (a line id of `physics/data/flavor_tree.json`) names, or `None` for no key."""
    return flavor_tree.text_of(key)


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
    system_config = planet.system_config
    with draw.bound(draw.Stream(flavor_tree.seed_of("roll", planet))):       # the body's own stream: the roll does not depend on the order the bodies were made in
        rolled = draw.random()
    if not (rolled < tuning.FLAVOR_CHANCE_PLANET and system_config.system_flavor_count < tuning.MAX_FLAVOR_TOTAL):
        log.debug("Planet flavor text: none (roll failed FLAVOR_CHANCE_PLANET or MAX_FLAVOR_TOTAL reached).")
        return

    # The tree (GEN.1430) gathers every line whose gates pass on the body's stored conditions (a ground to stand on or none, molten all over, multicellular or technological life),
    # never the class letter (GEN.1399: the class is decided last), takes out the ids already used in the system and draws one by weight.
    option = flavor_tree.draw_line("planet", planet, system_config.used_flavor_ids)
    if option is None:
        log.debug("Planet flavor text: none (the tree offers no line this body has not already shown in its system).")
        return
    log.choice("Planet flavor text", option.text, f"line: {option.id}")
    planet.flavor_key = option.id
    system_config.system_flavor_count += 1
    planet.flavor_text_count += 1
    system_config.used_flavor_ids.append(option.id)


def decide_flavor_text_alone(body, scope="planet"):
    """
    GEN.1433: the flavor line of a body that belongs to no system's list of lines, a rogue planet or one of its moons (`starless_host.build_rogue_body`), an asteroid belt (`scope="belt"`), or one of the other kinds of object (a star, comet, black hole, neutron star, nebula, remnant, quasar, interstellar comet, asteroid field: the scope is the kind). The same roll
    (`tuning.FLAVOR_CHANCE_PLANET`) and the same tree walk as `decide_flavor_text`, on the body's own stream, with no system bookkeeping to share: there is no
    system total to fill and no ids already shown, so what one rogue shows never changes what another does. Sets `flavor_key` and `flavor_text_count` in place.
    """
    if getattr(body, "flavor_key", None) is not None:      # decided already (stored, or drawn when it was made): the same body never draws twice
        return
    with draw.bound(draw.Stream(flavor_tree.seed_of(f"roll:{scope}", body))):
        rolled = draw.random()
    if rolled >= tuning.FLAVOR_CHANCE_PLANET:
        return
    option = flavor_tree.draw_line(scope, body)
    if option is not None:
        body.flavor_key = option.id
        if hasattr(body, "flavor_text_count"):
            body.flavor_text_count += 1       # a planet and a moon count their lines; a belt keeps only the key
