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
  again from the opposite start (a snowball start for a warm answer, a warm one for a cold answer), and the other answer is kept if it is a distinct
  one (`MIN_BRANCH_GAP_K` away, and a different ice cover). It is the same world, so the second solve repeats every draw of the first (`Air.entry_state`: the inventory, the efficiencies,
  the sink cycle and the late accretion) and only the start temperature differs (GEN.1646). Nothing is drawn from the run's stream.

Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.75 (rows NEB7 and NEB8).
"""

import dataclasses

from planetgen.physics import climate_branches, constants, nightside_collapse, oddities, saturation, surface_air
from planetgen.util import draw

CONTEXT_ATTRIBUTE = "_air_context"
CHOICE_ATTRIBUTE = "_air_choice"

COLD_START_K = 200.0
"""float: The temperature the second solve starts from to find a snowball (below the `climate_branches.ICE_COVERED_BELOW_K` ramp). Source tier T6, tweak level L0."""
WARM_START_K = 320.0
"""float: The temperature the second solve starts from to find a warm world (above the `climate_branches.ICE_FREE_ABOVE_K` ramp). Source tier T6, tweak level L0."""
MIN_BRANCH_GAP_K = 20.0
"""float: How far apart two answers must be to be two branches (the ramp is 65 K wide). Source tier T6, tweak level L1."""
MIN_ICE_GAP = 0.1
"""float: How different the Bond albedo of the two answers must be (the ice-free and ice-covered values differ by 0.32). Source tier T6, tweak level L1."""
NEAR_WINDOW_MAX_K = 330.0
"""float: A world whose first answer is hotter than this has no snowball branch to find (the ice line cannot be at the equator at that flux): the second solve is not made. Source tier T6, tweak level L1."""


@dataclasses.dataclass
class Context:
    """What the air step keeps of itself for the rows: the setting it solved, the air it got and the starlight."""

    setting: object
    air: object
    start_k: float
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
    found = nightside_collapse.collapse(formula, context.air.pressure_pa, gravity, context.starlight_w_m2, body.mass / constants.EARTH_MASS_TO_KG, context.albedo)
    if not found.collapsed:
        return False
    setattr(body, CHOICE_ATTRIBUTE, "lost")
    return None


@oddities.register("other_branch")
def other_branch(body, row):
    """
    The effect of oddity W20: a bistable world sits on its other branch. See the module text; declines (False) for a world with no ice law, a hot one, one the second start leaves on the
    same branch, or one whose second solve finds no air.
    """
    context = getattr(body, CONTEXT_ATTRIBUTE, None)
    if context is None or context.air is None or context.setting.albedo_law is None or getattr(body, "is_moon", False) or getattr(body, CHOICE_ATTRIBUTE, None) is not None:
        return False
    first = context.air
    if first.temperature_k > NEAR_WINDOW_MAX_K or first.albedo is None:
        return False
    middle = 0.5 * (climate_branches.ICE_COVERED_BELOW_K + climate_branches.ICE_FREE_ABOVE_K)
    start = COLD_START_K if first.temperature_k > middle else WARM_START_K
    with draw.bound(oddities.stream(f"{row.id}:air", oddities.body_key(body))):
        if first.entry_state is not None:     # GEN.1646: the other branch is the same world, so every draw of the first solve is repeated and only the start temperature differs
            draw.setstate(first.entry_state)
        other = surface_air.solve(context.setting, start)
    if other is None or not other.partials_pa or other.albedo is None:
        return False
    if abs(other.temperature_k - first.temperature_k) < MIN_BRANCH_GAP_K or abs(other.albedo - first.albedo) < MIN_ICE_GAP:
        return False
    setattr(body, CHOICE_ATTRIBUTE, other)
    return None
