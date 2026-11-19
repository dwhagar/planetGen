# planetgen/physics/air_history.py

"""
The clocks of a world's air history (GEN.1745, slice 1; design: reports/air-history-design-gen1739-2026-11-15.md, plan: reports/gen1745-build-plan-2026-11-15.md)
====================================================================================================================================================

Today a world's air is solved once, at its star's present age (`air_final.make`), and the laws of the air that depend on time (the star's light, the interior's heat, the water lost to the runaway
greenhouse, the reductants that hold the oxygen back) are each read once, at that age. The history of a world's air is those laws read at a series of ages instead of once. This module is the first
part of it: the **step schedule** (`schedule`) and the **clocks** every step reads (`Step`, `run`): the star's light, the heat the interior gives off, the reductant sink of the oxygen law over its
source, and the water lost so far. They are the existing laws, evaluated at the age of the step: `runaway_history.light_at` and `loss_kg` (with their trace), `interior_heat.heat_at`,
`biogenic_oxygen.reductant_history_factor` and `profile_of`. Nothing is drawn from the ambient stream and nothing is stored; nothing in the generator reads this module yet (the later slices of GEN.1745 feed the final
solve from it), so the output of a run is the one it was.

**The steps.** Logarithmic from `FIRST_STEP_GYR` (10 kyr) to `LOG_UNTIL_GYR` (1 Gyr), `LOG_STEPS_PER_DECADE` to the decade, then every `MAX_STEP_GYR` (the step the oxygen law's own integral takes,
`biogenic_oxygen.RISE_STEP_GYR`: the heat flow it reads changes by under 10 percent over a step after the first Gyr), cut at every age the caller names (the ages of the life stages and of the oxygen
rise) and ending at the star's age. These are numerical settings, not physical values.

**The temperature path (slice 1b).** `run(body, cuts, air_input)` with `air_input = (setting, start_k, entry_state)` of the body's final air solve (`air_final.make`) gives each step the surface
temperature of the body's air at that age. `aged_setting` is the solver input the body had at that age (the star's light and the interior's heat of the age, the rest the final one); the steps where the light or the
heat has moved by `ANCHOR_LIGHT_STEP` or `ANCHOR_HEAT_STEP` since the last one are solved in full (`_anchor_indices`, at most `MAX_ANCHORS`: a full solve is 60 to 600 ms), from the draws the final air began with
(`surface_air.Air.entry_state`), and the steps between are read on the line between their two anchors in the starlight absorbed. Against a full solve at every step on generated worlds the temperature is within
about 1 K on average and 5 K at most over the range where it matters for the air (100 Myr on, under 400 K), and up to 20 K in the first 100 Myr (a steam or magma ocean air, which the later slices give their own law) and across an
ice-albedo branch (`reports/gen1745-slice1b-measure-2026-11-16.md`).

**The oxygen trajectory (slice 2a).** `oxygen_trajectory(body, living, ages)` gives the oxygen, mol, the air holds at each age as the classes of `living` raise it: the air as the life found it until the rise, then the plateau
or the balance of the source against the reductant sink of that age (the laws of `biogenic_oxygen` and `life_air.oxygen_state`, read at each age and not once at the star's age). Its last value is the oxygen `biogenic_oxygen.apply` puts into the air, to the last
digit (`reports/gen1745-slice2a-check-2026-11-17.md`: 9 worlds on which the oxygen rose in 120 generated systems, every one equal). Nothing reads it yet; slice 2c swaps it in for `life._build_oxygen`.

**The ledger gases (slice 2b).** `gas_trajectory(body, living, ages, found, risen, rise_gyr)` gives the air, mol of each gas, at each age as the ledger rows of `living` make it: the methane the classes that have begun hold together
(and the carbon dioxide it is made from), the hydrogen, hydrogen sulfide and the like each draws down, counted at each age (`life_air.ledger_air`, `life_air._producing` with an age) and not once at the star's age. Its last value is the air
`life_air.apply_gases` stores. Nothing reads it yet; slice 2c swaps it in for `life_air.apply_gases`.

**The oxygen law's label** of a step (`oxygen_law`): `"none"` while the reductant sink `R(t)` is above the source `S` (no free oxygen can stand), `"rise"` once it is under `S` (the store the first
oxygen goes to can fill, `biogenic_oxygen.rise_age_gyr`), `"balance"` once it is under `S / SECOND_RISE_RATIO` (the second rise, GEN.1707); `None` for a body whose life makes no oxygen.
"""

import dataclasses
import math
from dataclasses import dataclass

from planetgen.generation.build_record import record_of
from planetgen.physics import biogenic_oxygen, constants, interior_heat, life_air, runaway_history, surface_air
from planetgen.util import draw

FIRST_STEP_GYR = 1.0e-5
"""float: The age, Gyr, of the first step of a history (10 kyr: the start of the air as the body left its disk; a numerical setting). (a numerical setting of the history's steps, defined by the integration: GEN.1739)"""
LOG_STEPS_PER_DECADE = 4
"""int: The steps to a decade of age before `LOG_UNTIL_GYR` (a numerical setting). (a numerical setting of the history's steps, defined by the integration: GEN.1739)"""
LOG_UNTIL_GYR = 1.0
"""float: The age, Gyr, up to which the steps are logarithmic; after it they are `MAX_STEP_GYR` apart. (a numerical setting of the history's steps, defined by the integration: GEN.1739)"""
MAX_STEP_GYR = biogenic_oxygen.RISE_STEP_GYR
"""float: The longest step, Gyr, after `LOG_UNTIL_GYR`: the step of the oxygen law's own integral (GEN.1463)."""
SAME_AGE_GYR = 1.0e-9
"""float: Two ages closer than this, Gyr, are one step. (a numerical setting of the history's steps, defined by the integration: GEN.1739)"""
ANCHOR_LIGHT_STEP = 0.03
"""float: A step is solved in full (an anchor) once the starlight its air absorbs has changed by this share since the last anchor; the temperature moves by a few kelvin over it. (a numerical setting of the history's steps, defined by the integration: GEN.1739)"""
ANCHOR_HEAT_STEP = 0.25
"""float: ... or the interior's heat flow has changed by this share since the last anchor. (a numerical setting of the history's steps, defined by the integration: GEN.1739)"""
MAX_ANCHORS = 12
"""int: The most full air solves of one history, the first and the last included (a full solve is 60 to 600 ms; the cost guard of GEN.1745). (a numerical setting of the history's steps, defined by the integration: GEN.1739)"""

NO_OXYGEN, RISE, BALANCE = "none", "rise", "balance"
OXYGEN_LAWS = (NO_OXYGEN, RISE, BALANCE)
"""tuple: The labels of `Step.oxygen_law`."""


@dataclass(frozen=True)
class Step:
    """
    One age of a world's history and what the clocks read there.

    Attributes:
        age_gyr (float): The age of the step, Gyr; the last step is the star's age.
        luminosity_ratio (float or None): The star's luminosity at that age over its luminosity now (`runaway_history.light_at`); `None` for a star whose past is not followed.
        heat_w_m2 (float or None): The heat the body's interior gives off at that age (`interior_heat.heat_at`, the body's own flux), W/m2; `None` for a body with no thermal history stored.
        heat_ratio (float or None): That heat over the heat it gives off at the star's age.
        reductant_ratio (float): `R(t) / R(now)` of the oxygen law (`biogenic_oxygen.reductant_history_factor`); 1 for a body with no thermal history.
        water_lost_kg (float): The water the body had lost to the runaway greenhouse by that age (`runaway_history.loss_kg`'s trace); 0 for a body that lost none.
        oxygen_law (str or None): The label of the oxygen law at that age (`OXYGEN_LAWS`); `None` when the body's oxygenic source is nothing.
        temperature_k (float or None): The surface temperature of the body's air at that age, K (`temperatures`); `None` when the run was given no air to solve.
        anchor (bool): Whether `temperature_k` is a full air solve at that age (else it is read between the two anchors round it).
    """

    age_gyr: float
    luminosity_ratio: float | None
    heat_w_m2: float | None
    heat_ratio: float | None
    reductant_ratio: float
    water_lost_kg: float
    oxygen_law: str | None
    temperature_k: float | None = None
    anchor: bool = False


def schedule(age_gyr, cuts=()):
    """
    The ages, Gyr, of the steps of a history that ends at `age_gyr`: logarithmic from `FIRST_STEP_GYR` to `LOG_UNTIL_GYR`, then `MAX_STEP_GYR` apart, with every age of `cuts` that falls inside
    the history added, and `age_gyr` last. Ascending, with no two the same age. Empty for a star with no age.
    """
    if not age_gyr or age_gyr <= 0.0:
        return []
    ages = set()
    step_ratio = 10.0 ** (1.0 / LOG_STEPS_PER_DECADE)
    age = FIRST_STEP_GYR
    while age < min(LOG_UNTIL_GYR, age_gyr):
        ages.add(age)
        age *= step_ratio
    age = LOG_UNTIL_GYR
    while age < age_gyr:
        ages.add(age)
        age += MAX_STEP_GYR
    ages.update(cut for cut in cuts if cut is not None and 0.0 < cut < age_gyr)
    ages.add(age_gyr)
    kept = []
    for age in sorted(ages):
        if not kept or age - kept[-1] > SAME_AGE_GYR:
            kept.append(age)
        elif age == age_gyr:
            kept[-1] = age_gyr            # the star's age is always the last step, exactly
    return kept


def _oxygen_law(profile, reductant_ratio):
    """The label of the oxygen law at a step: `R(t) = profile.reductant_mol_yr * reductant_ratio` against the source (`None` when there is none)."""
    if profile.source_mol_yr <= 0.0:
        return None
    reductant = profile.reductant_mol_yr * reductant_ratio
    if reductant < profile.source_mol_yr / biogenic_oxygen.SECOND_RISE_RATIO:
        return BALANCE
    return RISE if reductant < profile.source_mol_yr else NO_OXYGEN


def _water_lost_at(trace, age_gyr):
    """The water lost by `age_gyr` from the `(age_yr, kg)` trace of `runaway_history.loss_kg` (the loss at the end of the last step that ended by then; 0 before the first)."""
    lost = 0.0
    for age_yr, kilograms in trace:
        if age_yr > age_gyr * 1.0e9 * (1.0 + 1.0e-12):
            break
        lost = kilograms
    return lost


def aged_setting(setting, body, age_gyr):
    """
    The solver input (`surface_air.Setting`) of `body` as it was at `age_gyr`: the same body, the starlight of the star's luminosity at that age (`runaway_history.light_at`; the albedo law and the
    clouds act on it as before) and the heat of the interior at that age (`interior_heat.heat_at`: the heat released, the flow, and the heat released by the sealing time or the age if younger). Everything else
    (the water, the delivered volatiles, the haze and far-ultraviolet draws) is the final one, which slice 2 of GEN.1745 replaces with the history's own stores. `None` for a body or star the clocks do not cover.
    """
    star = getattr(body, "star", None)
    light = runaway_history.light_at(star, age_gyr)
    heat = interior_heat.heat_at(body, age_gyr)
    distance_au = record_of(body).star_distance_au or body.distance
    if light is None or heat is None or not distance_au or not setting.starlight_w_m2:
        return None
    flux = light[0] / (4.0 * math.pi * (distance_au * constants.AU_M) ** 2)
    scale = flux / setting.starlight_w_m2
    return dataclasses.replace(setting, age_gy=age_gyr, released_j=heat.released_j, heat_flux_w_m2=heat.flux_w_m2,
                               released_at_seal_j=interior_heat.released_by(body, min(age_gyr, 0.5)),
                               starlight_w_m2=setting.starlight_w_m2 * scale, absorbed_flux_w_m2=setting.absorbed_flux_w_m2 * scale)


def _anchor_indices(steps, scale_of, heat_of):
    """The indices of the steps solved in full: the first step with a light and a heat, the last, and each step where the light or the heat has moved by `ANCHOR_*_STEP` since the last one (at most `MAX_ANCHORS`: the ones left out are the ones
    nearest their neighbours)."""
    chosen = [0]
    for index in range(1, len(steps) - 1):
        last = chosen[-1]
        light_moved = abs(scale_of[index] / scale_of[last] - 1.0) if scale_of[last] else 1.0
        heat_moved = abs(heat_of[index] / heat_of[last] - 1.0) if heat_of[last] else 1.0
        if light_moved > ANCHOR_LIGHT_STEP or heat_moved > ANCHOR_HEAT_STEP:
            chosen.append(index)
    chosen.append(len(steps) - 1)
    while len(chosen) > MAX_ANCHORS:
        gaps = [(max(abs(math.log(max(scale_of[chosen[k + 1]], 1e-30) / max(scale_of[chosen[k - 1]], 1e-30))),
                     abs(math.log(max(heat_of[chosen[k + 1]], 1e-30) / max(heat_of[chosen[k - 1]], 1e-30)))), k) for k in range(1, len(chosen) - 1)]
        turning = {k for _, k in gaps if (scale_of[chosen[k]] - scale_of[chosen[k - 1]]) * (scale_of[chosen[k + 1]] - scale_of[chosen[k]]) < 0.0}
        kept_apart = [(gap, k) for gap, k in gaps if k not in turning] or gaps         # an anchor at a turning point of the light (the star's pre-main-sequence bump) is the last to go: the steps round it are read between it and its neighbours
        del chosen[min(kept_apart)[1]]
    return chosen


def _solved_temperature(setting, start_k, entry_state):
    """The surface temperature of `setting` solved in full from the draws the final air began with (`surface_air.Air.entry_state`), or `None` when the solve gives no air."""
    stream = draw.Stream(0)
    stream.setstate(entry_state)
    with draw.bound(stream):
        air = surface_air.solve(setting, start_k)
    return None if air is None else air.temperature_k


def temperatures(body, steps, air_input):
    """
    `steps` with the surface temperature of the body's air at each age: a full `surface_air.solve` at the anchors (`_anchor_indices`) and, between two anchors, the temperature read on the line between them in the
    starlight the body takes in (the one quantity of the aged setting that moves the temperature a step moves; the age when two anchors have the same light). `air_input` is `(setting, start_k, entry_state)` of the final
    solve of the body (`air_final.make`). Steps with no clocks (no star light or heat), and every step when an anchor gives no air, keep `temperature_k` `None`. Each solve reads its draws from `entry_state`, so every
    anchor is the air the body has now, solved under the light and heat of its age.
    """
    setting, start_k, entry_state = air_input
    aged = [aged_setting(setting, body, step.age_gyr) for step in steps]
    if not steps or any(item is None for item in aged):
        return list(steps)
    scale_of = [item.absorbed_flux_w_m2 for item in aged]
    heat_of = [item.heat_flux_w_m2 for item in aged]
    anchors = _anchor_indices(steps, scale_of, heat_of)
    solved = {}
    for index in anchors:
        temperature = _solved_temperature(aged[index], start_k, entry_state)
        if temperature is None:
            return list(steps)
        solved[index] = temperature
    result = []
    for index, step in enumerate(steps):
        if index in solved:
            result.append(dataclasses.replace(step, temperature_k=solved[index], anchor=True))
            continue
        below = max(k for k in anchors if k < index)
        above = min(k for k in anchors if k > index)
        span = scale_of[above] - scale_of[below]
        share = (scale_of[index] - scale_of[below]) / span if span else (math.log(step.age_gyr / steps[below].age_gyr) / math.log(steps[above].age_gyr / steps[below].age_gyr))
        share = min(max(share, 0.0), 1.0)
        result.append(dataclasses.replace(step, temperature_k=solved[below] + share * (solved[above] - solved[below])))
    return result


def run(body, cuts=(), air_input=None):
    """
    The steps of the history of a rocky `body` (`schedule` of its star's age, with the ages in `cuts` added): at each, the star's light, the interior's heat, the oxygen law's reductant ratio and
    label and the water lost. With `air_input` (`(setting, start_k, entry_state)` of the body's final air solve, `temperatures`) each step also has the surface temperature of its air. Reads the stored body and draws
    nothing from the ambient stream. Empty for a body with no star or no star age.
    """
    star = getattr(body, "star", None)
    star_age = getattr(star, "age", None)
    if star is None or not star_age:
        return []
    heat_now = interior_heat.heat_at(body, star_age)
    flux_now = None if heat_now is None else heat_now.own_flux_w_m2
    trace = []
    distance_au = record_of(body).star_distance_au or body.distance
    mass_kg, radius_km = getattr(body, "mass", None), getattr(body, "radius", None)
    if distance_au and mass_kg and radius_km:
        runaway_history.loss_kg(star, distance_au, mass_kg, radius_km, trace=trace)
    profile = biogenic_oxygen.profile_of(body)
    steps = []
    for age in schedule(star_age, cuts):
        light = runaway_history.light_at(star, age)
        heat = interior_heat.heat_at(body, age)
        own = None if heat is None else heat.own_flux_w_m2
        ratio = biogenic_oxygen.reductant_history_factor(body, age)
        steps.append(Step(age_gyr=age,
                          luminosity_ratio=None if light is None else light[0] / star.luminosity,
                          heat_w_m2=own,
                          heat_ratio=None if own is None or not flux_now else own / flux_now,
                          reductant_ratio=ratio,
                          water_lost_kg=_water_lost_at(trace, age),
                          oxygen_law=_oxygen_law(profile, ratio)))
    return temperatures(body, steps, air_input) if air_input is not None else steps


def oxygen_cuts(body, living):
    """The ages, Gyr, a history of the oxygen of `body` is cut at: the age the oxygen rises (`life_air.oxygen_state`) and the ages the producers of each class of `living` began; `[]` for a body whose life makes no oxygen."""
    rise, _ = life_air.oxygen_state(body, living)
    if rise is None:
        return []
    return [rise, *[age for age in (one.started(body) for one in living) if age is not None]]


def oxygen_trajectory(body, living, ages):
    """
    The oxygen, mol, that the air of `body` holds at each of `ages` (ascending, Gyr) as the classes of `living` raise it (GEN.1745 slice 2a): `[(age_gyr, mol), ...]`. It is the oxygen law of
    `biogenic_oxygen` read at each age instead of once at the star's age: the air as the life found it (`have`, the oxygen of an abiotic air, none for most) until the rise (`life_air.oxygen_state`);
    from the rise, the plateau (`biogenic_oxygen.PLATEAU_PA`, never more than the balance) while the reductant sink `R(age)` (`profile.reductant_mol_yr` times `reductant_history_factor`) is over the source `S` over
    `SECOND_RISE_RATIO`; and the balance of `S - R(age)` against oxidative weathering (`biogenic_oxygen.balance_pa`) once it is under, held to the fire limit. At the star's age the reductant ratio is 1, so the last
    value is what `biogenic_oxygen.stage_target_mol` gives and `biogenic_oxygen.apply` puts into the air (the burnt air's own oxygen if that is more).
    Must be called on the air as the life found it, before `apply` has burnt its reduced gases and added the oxygen. Reuse: `life_air.oxygen_state`, `_oxygenic` and `standing_methane_mol`, `biogenic_oxygen.profile_of`,
    `burn`, `balance_pa`, `reductant_history_factor`, `_oxygen_mol_at` and `_fire_limit_mol`.
    """
    moles = biogenic_oxygen.air_moles(body)
    have = moles.get("O2", 0.0)
    rise, _ = life_air.oxygen_state(body, living)
    if rise is None:
        return [(age, have) for age in ages]
    pairs = life_air._oxygenic(body, living)
    extra = 2.0 * life_air.standing_methane_mol(body, living, rise)
    profile = biogenic_oxygen.profile_of(body, extra_store_mol=extra, producers=[(share, band) for _, share, band in pairs])
    burnt, _ = biogenic_oxygen.burn(moles)
    held = burnt.get("O2", 0.0)
    fire = biogenic_oxygen._fire_limit_mol(body, burnt)
    plateau = biogenic_oxygen._oxygen_mol_at(body, burnt, biogenic_oxygen.PLATEAU_PA)
    area = biogenic_oxygen.area_m2(body)
    trajectory = []
    for age in ages:
        if age < rise:
            trajectory.append((age, have))
            continue
        reductant = profile.reductant_mol_yr * biogenic_oxygen.reductant_history_factor(body, age)
        surplus = profile.source_mol_yr - reductant
        steady = min(biogenic_oxygen._oxygen_mol_at(body, burnt, biogenic_oxygen.balance_pa(surplus, area)), fire) if surplus > 0.0 else 0.0
        target = 0.0 if profile.steady_mol <= 0.0 else steady if reductant <= profile.source_mol_yr / biogenic_oxygen.SECOND_RISE_RATIO else min(steady, plateau)
        trajectory.append((age, max(held, min(target, fire))))
    return trajectory


def gas_trajectory(body, living, ages, found, risen=None, rise_gyr=None):
    """
    The air, mol of each gas, that `body` holds at each of `ages` (ascending, Gyr) as the classes of `living` make it (GEN.1745 slice 2b): `[(age_gyr, moles), ...]`. It is the ledger gases of
    `life_air` (the methane the classes that have begun hold together, the carbon dioxide it is made from, and the hydrogen, hydrogen sulfide and the like each draws down) counted at each age
    instead of once at the star's age: a class counts from the age its producers began (`life_air._producing`), and a methane row that ends at the oxygen's rise stops there (`rise_gyr`).
    `found` is the air as the life found it, before the oxygen was made; `risen` is the air the oxygen rise leaves (the oxygen put in and the reduced gases burnt, `biogenic_oxygen.apply`), read
    from the rise on (`None` or `rise_gyr` `None`: the oxygen has not risen, so the air is `found` at every age). The last value is the air `life_air.apply_gases` stores (`life_air.ledger_air`,
    to the last digit; `reports/gen1745-slice2b-check-2026-11-18.md`). The temperature of the air at an age is the final air's, as in slice 2a: the history's own temperature path
    (`temperatures`) reaches the gases in slice 3. The air only changes where a class begins or the oxygen rises, so each distinct state is worked once.
    Reuse: `life_air.ledger_air`, `_producing`, `_with_ledger_gases`, `summed_methane_flux_mol_yr`.
    """
    if getattr(body, "body_type", "t") == "g" or not found:       # a giant, or no air: the ledger gives it nothing (`life_air.apply_gases`)
        return [(age, found) for age in ages]
    states = {}
    trajectory = []
    for age in ages:
        rose = rise_gyr is not None and risen is not None and age >= rise_gyr
        key = (rose, tuple((life_air._producing(body, one, rose, age), life_air._producing(body, one, False, age)) for one in living))
        if key not in states:
            states[key] = life_air.ledger_air(body, living, rose, risen if rose else found, age)[0]
        trajectory.append((age, states[key]))
    return trajectory
