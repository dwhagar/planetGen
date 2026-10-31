# planetgen/generation/helium_core_draw.py

"""
The helium-core white dwarf of a star drawn as a white dwarf (GEN.996, part 2)
===============================================================================

A white dwarf under 0.45 solar masses is the core of a giant that a companion stripped on the red giant branch before the helium flash
(`physics/helium_core_white_dwarf.py`). The tree gives a companion white dwarf that core where the pair's orbit says it overflowed
(`multistar._Builder.helium_core`). The white dwarf a star is drawn as (`Star`, `sample_living_star`) is made before its tree, so
it needs the same answer without one: this module draws, for a star of this birth mass and [Fe/H] that has become a white dwarf,
the pair it would have had (the chance of a companion of the `star-count` table, the companion's birth mass and the period of the
tree's own laws, `multistar.period_days`) and, when the pair's axis is where the lobe was reached with a core of 0.18 to 0.45 solar masses,
the white dwarf is that core. No rate is set: the share is what the generator's laws leave in the window (about 6 percent of the white
dwarfs of 1 to 1.5 solar mass progenitors; the 20 pc sample has 7.4 percent). The tree then gives such a white dwarf its companion
on the axis its core says (`multistar._bind_helium_core_companion`), so the star and its tree agree.
"""

import math

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.physics import constants, helium_core_white_dwarf, stability
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw

ZERO_AGE_RADIUS_EXPONENT = 0.8
"""float: The exponent of the zero-age main-sequence radius in solar masses (R ~ M^0.8 solar radii) the draw reads the closest a pair can be
from: its floor is the two radii touching. Source tier T6, tweak level L0 (the usual low-mass main-sequence relation)."""


def draw_core_mass_sol(initial_mass_sol, fe_h, rng=draw):
    """
    The helium core, solar masses, that a star of this birth mass and [Fe/H] left as a white dwarf, or None when it left an ordinary one: no
    companion drew, or the pair's axis is outside the window in which the giant's lobe was reached before the helium flash.
    """
    if initial_mass_sol > tuning.HELIUM_CORE_PROGENITOR_MAX_SOL:
        return None
    weights, inside_weight = multistar.close_trend(initial_mass_sol, fe_h)
    if weights is None:
        weights = multistar.star_count_weights(initial_mass_sol)
    if rng.random() >= 1.0 - weights[0] / sum(weights):
        return None
    companion = multistar.companion_mass(initial_mass_sol, rng, fe_h)
    total = initial_mass_sol + companion
    floor_au = ((initial_mass_sol ** ZERO_AGE_RADIUS_EXPONENT + companion ** ZERO_AGE_RADIUS_EXPONENT)
                * constants.SOLAR_RADIUS_M / constants.AU_M)
    period_d = multistar.period_days("star", total, floor_au, math.inf, initial_mass_sol >= multistar.MASSIVE_PRIMARY_MIN_MASS_SOL,
                                     inside_weight, rng.random())
    if period_d is None:
        return None
    axis = stability.axis_au(period_d / 365.25, total)
    return helium_core_white_dwarf.overflow_core_mass_sol(initial_mass_sol, companion, axis, 0.0, fe_h)


def with_helium_core(state, initial_mass_sol, age_gy, rng=draw):
    """The state of a drawn star, as `sample_living_star` gives it, with its white dwarf made a helium core when the draw says so
    (`draw_core_mass_sol`); a star that is not a white dwarf, or an ordinary white dwarf, comes back as it is."""
    if state["yerkes_class"] != "VII":
        return state
    core = draw_core_mass_sol(initial_mass_sol, state["metallicity_fe_h"], rng)
    if core is None:
        return state
    return {**state, **evolution.helium_core_state(initial_mass_sol, age_gy, state["metallicity_fe_h"], core, rng)}
