# planetgen/physics/air_branches.py

"""
The air-state oddities: the nightside collapse of a locked planet's air (W19) and the other branch of a bistable climate (W20) (GEN.1051, lever L7 of docs/design/systems/weird-worlds.md)
=========================================================================================================================================================================================

Both rows act on the air the solver has just made (`planets.calculate_atmospheric_conditions` keeps it for them as `body._air_context`) and may replace it:

* **W19, the nightside collapse** (`physics/nightside_collapse.py`, `docs/design/composition/nightside-collapse-model.md`). A rocky planet that the star's tide has locked (the estimate the clouds
  use, `cloud_albedo.tidally_locked`; the spin is drawn after the air) and whose air is thin loses it onto its dark side. The main gas of the solved air goes through `collapse` at the planet's
  gravity and starlight; when it keeps under half of its gas (`Collapse.collapsed`) the planet has no air (the day side is bare and the surface temperature is the starlight's alone), as
  the airless branch of the solver. The pressure the gas keeps is a tiny share of the air (`Collapse.residual_pressure_pa`) and is not stored; the cap on the night side is not stored yet.
* **W20, the other branch.** A world with surface water has the ice-albedo feedback (GEN.972, `climate_branches`), and in a band of starlight it holds two stable climates: a warm one and a
  snowball. The solver starts from the airless temperature and finds one of them. The row is the history draw between them: when it rolls and the air is cool enough to be near the window, the air is solved
  again from the opposite starts in turn (`COLD_STARTS_K` for a warm answer, `WARM_STARTS_K` for a cold one: the solver's answer depends on its path, so one start can miss the other branch, GEN.1715), and the first answer that is a distinct
  one is kept (`MIN_BRANCH_GAP_K` away, and a different ice cover); the search stops when `AGREEING_STARTS_TO_STOP` starts have ended on the first answer. It is the same world, so the second solve repeats every draw of the first (`Air.entry_state`: the inventory, the efficiencies,
  the sink cycle and the late accretion) and only the start temperature differs (GEN.1646). Nothing is drawn from the run's stream.

Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.75 (rows NEB7 and NEB8).
"""

import dataclasses

from planetgen import tuning
from planetgen.physics import climate_branches, constants, nightside_collapse, oddities, saturation, surface_air
from planetgen.util import draw
from planetgen.physics import unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")

CONTEXT_ATTRIBUTE = "_air_context"
CHOICE_ATTRIBUTE = "_air_choice"

COLD_STARTS_K = (200.0, 150.0, 220.0, 260.0)
"""tuple: The temperatures the second solve starts from, in turn, to find a snowball (the first under the `climate_branches.ICE_COVERED_BELOW_K` ramp); it stops at the first start that ends on a distinct answer
(GEN.1715: the air solve ends on a branch that depends on the path as much as on the start, so one cold start of 200 K found the snowball of 1 world in 200 keys where 4 starts find it in about 10). Source tier T6, tweak level L0."""
WARM_START_K = 320.0
"""float: The temperature the solve starts from to find a warm world (above the `climate_branches.ICE_FREE_ABOVE_K` ramp), and the first start of `WARM_STARTS_K`. Source tier T6, tweak level L0."""
WARM_STARTS_K = (WARM_START_K, 280.0, 340.0)
"""tuple: The same for a second solve that looks for the warm world, tried in turn as `COLD_STARTS_K` are. Source tier T6, tweak level L0."""
MIN_BRANCH_GAP_K = 20.0
"""float: How far apart two answers must be to be two branches (the ramp is 65 K wide). Source tier T6, tweak level L1."""
MIN_ICE_GAP = 0.1
"""float: How different the Bond albedo of the two answers must be (the ice-free and ice-covered values differ by 0.32). Source tier T6, tweak level L1."""
AGREEING_STARTS_TO_STOP = 2
"""int: How many starts may end on the first answer (less than `MIN_BRANCH_GAP_K` and `MIN_ICE_GAP` from it) before the search stops, the world having one branch (the second solve is the cost of the oddity, 12 to 20 seconds
a world with a dense air, and of 132 worlds that rolled W20 129 had every start end on the first answer). A start that ends elsewhere (neither the same nor a distinct answer) does not count, so the search goes on
through them. Source tier T6, tweak level L1."""
NEAR_WINDOW_MAX_K = 330.0
"""float: A world whose first answer is hotter than this has no snowball branch to find (the ice line cannot be at the equator at that flux): the second solve is not made. Source tier T6, tweak level L1."""


def formation_start_k(no_air_temperature_k, *, steam, ice_law):
    """
    The temperature the first solve of a body starts from, the branch the body is stored on (GEN.1716: the answer is chosen by this rule on the physics and not by a number a caller passes).
    A body holding water past its runaway limit (`steam`, `air_final.past_runaway_with_water`) is on the steam branch: it starts at `tuning.RUNAWAY_START_TEMPERATURE_K`. A body with the ice law
    (`ice_law`) cooled from its formation, so it enters from the warm side (GEN.1696): `WARM_START_K`. Any other body starts from its airless temperature.
    """
    if steam:
        return max(no_air_temperature_k, tuning.RUNAWAY_START_TEMPERATURE_K)      # GEN.1238: the steam branch, not the cool one
    if ice_law:
        return max(no_air_temperature_k, WARM_START_K)
    return no_air_temperature_k


@dataclasses.dataclass
class Context:
    """What the air step keeps of itself for the rows: the setting it solved, the air it got, whether it stands past the runaway limit and the starlight."""

    setting: object
    air: object
    steam: bool
    starlight_w_m2: float
    distance_au: float
    albedo: float


def clear(body):
    """Forgets the context and the choice of a body (the step that made them is about to run again)."""
    setattr(body, CONTEXT_ATTRIBUTE, None)
    setattr(body, CHOICE_ATTRIBUTE, None)


def main_gas(partials_pa):
    """The formula of the gas with the largest partial pressure that has a saturation curve, or `None`."""
    for formula, pressure in sorted(partials_pa.items(), key=lambda item: -item[1]):
        if saturation.has(formula) and pressure > 0.0:
            return formula
    return None


@oddities.register("collapse_air")
def collapse_air(body, row):
    """
    The effect of oddity W19: the air of a locked planet freezes onto its dark side. Declines (False) when the planet is not locked, has no air, has no main gas the model knows or keeps
    half of its gas. Otherwise the air is lost: `_air_choice` is the string `"lost"` and the caller makes the airless state.
    """
    context = getattr(body, CONTEXT_ATTRIBUTE, None)
    if context is None or context.air is None or not context.air.partials_pa or getattr(body, "is_moon", False):
        return False
    law = surface_air.cloud_law_of(body, context.starlight_w_m2, context.distance_au)
    if not law.locked:
        return False
    formula = main_gas(context.air.partials_pa)
    if formula is None:
        return False
    gravity = constants.G * body.mass / (body.radius * 1000.0) ** 2
    found = nightside_collapse.collapse(formula, context.air.pressure_pa, gravity, context.starlight_w_m2, _kg_to_mearth(body.mass), context.albedo)
    if not found.collapsed:
        return False
    setattr(body, CHOICE_ATTRIBUTE, "lost")
    return None


@oddities.register("other_branch")
def other_branch(body, row):
    """
    The effect of oddity W20: a bistable world sits on its other branch. See the module text; declines (False) for a world with no ice law, a hot one, one the second start leaves on the
    same branch, one past its runaway limit, or one whose second solve finds no air.
    """
    context = getattr(body, CONTEXT_ATTRIBUTE, None)
    if context is None or context.air is None or context.setting.albedo_law is None or getattr(body, "is_moon", False) or getattr(body, CHOICE_ATTRIBUTE, None) is not None:
        return False
    if context.steam:
        return False        # GEN.1716: past the runaway limit no ocean can stand; the cool air a cold start finds there is the greenhouse table having no limit (`air_final.past_runaway_with_water`), not a branch
    first = context.air
    if first.temperature_k > NEAR_WINDOW_MAX_K or first.albedo is None:
        return False
    middle = 0.5 * (climate_branches.ICE_COVERED_BELOW_K + climate_branches.ICE_FREE_ABOVE_K)
    starts = COLD_STARTS_K if first.temperature_k > middle else WARM_STARTS_K
    other = None
    agreeing = 0
    for start in starts:
        with draw.bound(oddities.stream(f"{row.id}:air", oddities.body_key(body))):
            if first.entry_state is not None:     # GEN.1646: the other branch is the same world, so every draw of the first solve is repeated and only the start temperature differs
                draw.setstate(first.entry_state)
            candidate = surface_air.solve(context.setting, start)
        if candidate is None or not candidate.partials_pa or candidate.albedo is None:
            continue
        if abs(candidate.temperature_k - first.temperature_k) >= MIN_BRANCH_GAP_K and abs(candidate.albedo - first.albedo) >= MIN_ICE_GAP:
            other = candidate
            break
        if abs(candidate.temperature_k - first.temperature_k) < MIN_BRANCH_GAP_K and abs(candidate.albedo - first.albedo) < MIN_ICE_GAP:
            agreeing += 1
            if agreeing >= AGREEING_STARTS_TO_STOP:
                break
    if other is None:
        return False
    setattr(body, CHOICE_ATTRIBUTE, other)
    return None
