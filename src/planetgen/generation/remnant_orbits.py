# planetgen/generation/remnant_orbits.py

"""
The planets of a white dwarf and the disk of a planetary nebula's central pair (GEN.318)
=======================================================================================

* **Mass loss** (`widen_for_mass_loss`). The planets a white dwarf is given are drawn as the ones that formed around its
  progenitor; `StarSystem._apply_star_history` has already removed those inside the progenitor's engulfment radius
  (`physics.mass_loss.engulfment_radius_au`), and the ones left are moved out by what the star lost
  (`physics.mass_loss.expansion_factor`, 1.85 for a star of one solar mass that leaves 0.54). A belt's edges move with
  it. A planet is reconciled with its new distance as one moved by the spacing check is
  (`validation.reconcile_moved_planet`): its zone, climate and orbit follow, and its class is derived again from the new conditions.
* **A belt through the giant branch and the mass loss** (`engulf_belt`, `widen_for_mass_loss`, GEN.1072). A belt the progenitor's envelope reached is
  not lost whole: the part inside the engulfment radius is and the outer part keeps its share of the mass (`physics.belt_history.engulf`). The mass
  loss widens its edges as it does any orbit, and with giants left in the system the sudden change throws a share of the belt inward
  (`belt_history.mass_loss`); that scattered mass is what falls on the white dwarf (`pollution_rate_g_s`). Each event leaves its line on the belt
  (`history`, a tuple of sentences); the history and the scattered mass are worked out when the system is built and are not stored.
* **A planetary nebula's close pair** (`is_post_agb_pair`, `post_agb_disk`). A close pair whose white dwarf is still as
  hot as a planetary nebula's central star (`tuning.PLANETARY_NEBULA_CENTRAL_MIN_TEMPERATURE_K`) is a post-common-envelope
  pair inside its ejecta: it gets a circumbinary disk descriptor (`disk_kind` "post-agb", on the primary, as the pair's
  disk is the primary's) and no planets. The note gives no rate, so the disk is a flag on the pair, not a population.
"""

import math

from planetgen import tuning
from planetgen.generation import validation
from planetgen.physics import belt_history, mass_loss
from planetgen.physics import unit_registry as units
from planetgen.util import draw
from planetgen.util import log

_kg_to_msun = units.converter("kg", "Msun")

POST_AGB_DISK_KIND = "post-agb"
"""str: The `disk_kind` of a planetary nebula's circumbinary disk."""


def _members(star):
    return getattr(star, "stars", None) or [star]


def is_post_agb_pair(star):
    """Whether `star` is a close pair (a `BinaryStarProxy`) with a white dwarf as hot as a planetary nebula's central
    star."""
    members = getattr(star, "stars", None)
    if not members or len(members) < 2:
        return False
    return any(getattr(member, "yerkes_class", None) in ("VII", "D")
               and member.temperature >= tuning.PLANETARY_NEBULA_CENTRAL_MIN_TEMPERATURE_K for member in members)


def post_agb_disk(primary, inner_au):
    """Gives `primary` the disk descriptor of a planetary nebula's close pair: a mass log-uniform in
    `POST_AGB_DISK_MASS_RANGE_SOLAR`, the circumbinary stability limit `inner_au` for an inner edge and an outer edge
    log-uniform in the young stars' range (`DISK_OUTER_RADIUS_RANGE_AU`), with no bodies in it."""
    low, high = tuning.POST_AGB_DISK_MASS_RANGE_SOLAR
    near, far = tuning.DISK_OUTER_RADIUS_RANGE_AU
    primary.disk_kind = POST_AGB_DISK_KIND
    primary.disk_mass_solar = math.exp(draw.uniform(math.log(low), math.log(high)))
    primary.disk_inner_au = inner_au
    primary.disk_outer_au = max(inner_au, math.exp(draw.uniform(math.log(near), math.log(far))))
    log.choice("Planetary nebula pair", "circumbinary disk, no planets",
               f"a close pair with a white dwarf of {max(m.temperature for m in _members(primary)):.0f} K")


def _belt_state(belt):
    """The `belt_history.BeltState` of `belt` after the events so far: formed from its edges the first time, and the mass it had then kept on
    `formed_mass_kg`."""
    state = getattr(belt, "belt_state", None)
    if state is None:
        belt.formed_mass_kg = belt.mass_kg
        state = belt_history.BeltState.formed(belt.lower_limit, belt.upper_limit)
    return state


def _apply_belt_state(belt, state):
    """Writes `state` back to `belt`: its limits, its distance (kept inside them), its mass, the mass thrown inward and its history."""
    belt.belt_state = state
    belt.lower_limit, belt.upper_limit = state.lower_au, state.upper_au
    if not belt.lower_limit <= belt.distance <= belt.upper_limit:
        belt.distance = 0.5 * (belt.lower_limit + belt.upper_limit)
    belt.mass_kg = belt_history.mass_kg(state, belt.formed_mass_kg)
    belt.scattered_mass_kg = belt_history.scattered_mass_kg(state, belt.formed_mass_kg)
    belt.history = state.history


def engulf_belt(belt, radius_au):
    """
    The part of `belt` inside the envelope's reach `radius_au` is lost, in place, and the part outside keeps its share of the mass
    (`belt_history.engulf`). A belt entirely inside is left as it is for the caller to remove (its inner edge is inside the radius).
    """
    if radius_au <= belt.lower_limit or radius_au >= belt.upper_limit:
        return
    _apply_belt_state(belt, belt_history.engulf(_belt_state(belt), radius_au))


def pollution_rate_g_s(bodies):
    """The rate, grams per second, at which the rock the belts of `bodies` threw inward falls on the white dwarf (`belt_history.pollution_rate_g_s`)."""
    return belt_history.pollution_rate_g_s(sum(getattr(body, "scattered_mass_kg", 0.0) for body in bodies if body.body_type == 'a'))


def formation_mass_kg(star):
    """The mass, kg, `star` had when its planets formed: its members' (`mass_loss.member_formation_mass_kg`)."""
    return sum(mass_loss.member_formation_mass_kg(member) for member in _members(star))


def widen_for_mass_loss(planets, star, explicit_request=False):
    """
    Moves the planets and belts of `planets` out, in place, by the mass `star` lost: each orbit grows by `(M_formation +
    M_p) / (M_now + M_p)`. A star that lost nothing is left as it is, and so is a forced-type white dwarf (no model
    progenitor) when a habitable world or a `SLOTS` list is asked for (`explicit_request`), whose planets were placed
    where the request needs them.

    Args:
        planets (list): The `Planet` and `AsteroidBelt` instances around `star`, after the engulfed ones are removed.
        star (Star): The star or close pair's proxy they orbit.
        explicit_request (bool): Whether a habitable world or `SLOTS` is asked for.

    Returns:
        bool: Whether any orbit moved (the caller spaces them again: the Hill radii are over the lighter star).
    """
    if explicit_request and any(getattr(member, "yerkes_class", None) in ("VII", "D")
                                and member.initial_mass_sol is None for member in _members(star)):
        return False
    formed = formation_mass_kg(star)
    if formed <= star.mass * (1.0 + 1e-9) or not planets:
        return False
    giants_remain = any(body.body_type == 'g' for body in planets)
    for body in planets:
        if body.body_type == 'a':
            factor = mass_loss.expansion_factor(formed, star.mass)
            body.distance *= factor
            state = belt_history.mass_loss(_belt_state(body), _kg_to_msun(formed), _kg_to_msun(star.mass), giants_remain)
            _apply_belt_state(body, state)
            continue
        body.distance *= mass_loss.expansion_factor(formed, star.mass, body.mass)
        validation.reconcile_moved_planet(body)
    log.debug(f"White dwarf planets widened by {mass_loss.expansion_factor(formed, star.mass):.3g} for the mass the "
              f"star lost ({formed:.3g} to {star.mass:.3g} kg)")
    return True
