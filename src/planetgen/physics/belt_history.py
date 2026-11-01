# planetgen/physics/belt_history.py

"""
What happens to a belt over a system's life: giant-planet events, the star's giant branch and mass loss, a companion, and the dust that falls on a white dwarf (GEN.1064; docs/design/systems/belt-evolution.md sections 2.5 and 4 steps 3 and 6, belt-sources.md BR4, BR5, BR12 to BR14)
===================================================================================================================================================================================================================================================================

A belt is the surviving part of what formed, and each event that acts on it is a pure function from a `BeltState` to a new one that keeps what is left, what was scattered inward and a line of text the page can show. Nothing calls it yet
(the belt wiring, GEN.1069 and later, runs the events a system's history holds in order). Mass is a share of the belt as formed; the ledger (`belt_ledger.py`) and the grinding (`belt_mass.py`) turn it into kilograms.

1. **A giant-planet event** (`giant_event`). A Grand Tack-type migration through the belt leaves `GRAND_TACK_SURVIVAL` of it (2e-3 to 1e-2, [R3150, R4603], READ), a giant-planet instability alone
   `INSTABILITY_SURVIVAL` (0.1 to 0.5; 50 to 90 percent lost, [R4603]). The share is drawn log-uniform over the range (the widest range until a population synthesis narrows it) from a stream, so a system's history is reproducible.
2. **The giant branch** (`engulf`). A belt that overlaps the radius the star's envelope reaches (`mass_loss.engulfment_radius_au`) is not lost whole: the part inside is, and the part outside keeps its share of the mass
   (the surface density goes as `a^-SURFACE_DENSITY_SLOPE`, so the mass between two distances goes as `a^(2 - slope)`). Only a belt entirely inside is lost.
3. **Mass loss** (`mass_loss`). When the star loses mass slowly every orbit widens by `M_before / M_after` (adiabatic invariance; [R1031], READ: a planet beyond 5 au survives and widens), the belt with them.
   If giants remain, the sudden change moves their resonances through the belt and a share `MASS_LOSS_SCATTER_SHARE` times the fraction of mass lost is thrown inward; with no giants nothing destabilises it.
4. **A companion** (`limit`). A close companion truncates a belt outside the stability limit ([R1990], in `orbits.holman_wiegert_critical_semimajor_axis`) and a circumbinary belt starts beyond the inner one
   (`holman_wiegert_circumbinary_a_crit_au`); the part outside the limit is lost.
5. **White dwarf pollution** (`pollution_rate_g_s`, `is_polluted`). The scattered rock falls on the white dwarf over `DELIVERY_MYR`; the rate is the scattered mass over that time. [R2415] (READ) measure rates from a few 1e5 to a few 1e8 g/s;
   a white dwarf is polluted when the rate is at least `MIN_POLLUTION_RATE_G_S` (their lowest, which is also about what can be detected). So a belt of 1e18 to 1e21 kg scattered over 100 Myr spans the observed range.

The Solar System is a test target (a main-belt-mass belt scattered at a white dwarf makes a rate at the top of the observed range). The fractions of white dwarfs polluted (27 to 50 percent accreting now, 56 percent with metals, BR14) are checks on the population this gives, not inputs.

Limits. Each event is applied to a belt as one annulus; the shares of mass inside a distance assume a single power-law surface density; the scattered share of a mass drop is a stand-in (BHS3) with no source behind it; the rain of comets
and the planets that also scatter rock (the other half of white dwarf pollution) are not here; a belt that survives a Grand Tack with the survival drawn is not tied to the system's own migration history (the origin table, GEN.1054, decides whether the event happened).

Use:

    state = BeltState.formed(lower_au, upper_au)
    state = giant_event(state, GRAND_TACK, stream)           # or INSTABILITY
    state = engulf(state, engulfment_radius_au)
    state = mass_loss(state, mass_before_solar, mass_after_solar, giants_remain)
    state = limit(state, inner_au, outer_au, "a companion at 120 au")
    state.history                                              # the lines the page shows
"""

from dataclasses import dataclass, replace

from planetgen.util.checks import finite_domain
from planetgen.util.random import log_uniform

GRAND_TACK = "grand tack"
INSTABILITY = "instability"

GRAND_TACK_SURVIVAL = (2.0e-3, 1.0e-2)
"""tuple: The share of a belt left after a giant migrates through it ([R4603], and [R3150], READ: depleted by a factor of several hundred). [Stand-in BR4 of 6.95, row BHS1, tier T2, tweak L1.]"""

INSTABILITY_SURVIVAL = (0.1, 0.5)
"""tuple: The share of a belt left after a giant-planet instability alone ([R4603], READ: 50 to 90 percent lost). [Stand-in BR5 of 6.95, row BHS1, tier T2, tweak L1.]"""

SURFACE_DENSITY_SLOPE = 1.5
"""float: The belt's solids go as `a^-slope` (the minimum-mass nebula, as `scatter.DISC_SLOPE`), so the mass between two distances goes as `a^(2 - slope)`. [Stand-in BHS2, tier T6, tweak L1.]"""

MASS_LOSS_SCATTER_SHARE = 0.5
"""float: The share of a belt thrown inward by the sudden mass loss of the star, per unit fraction of the star's mass lost, if giants remain. [Stand-in BHS3, tier T9, tweak L3.]"""

DELIVERY_MYR = 100.0
"""float: How long the scattered rock takes to fall on the white dwarf ([R2415]'s cooling ages are 20 to 200 Myr; the middle in the log). [Stand-in BHS4, tier T6, tweak L1.]"""

MIN_POLLUTION_RATE_G_S = 1.0e5
"""float: The lowest accretion rate seen on a polluted white dwarf ([R2415], READ: a few 1e5 to a few 1e8 g/s). [Stand-in BHS5, tier T2, tweak L0.]"""

_SECONDS_PER_MYR = 3.15576e13


def _slope_exponent():
    return 2.0 - SURFACE_DENSITY_SLOPE


def _share_between(lower_au, upper_au, low_au, high_au):
    """The share of the mass of a belt from `lower_au` to `upper_au` that lies between `low_au` and `high_au` (a single power-law surface density)."""
    k = _slope_exponent()
    low, high = max(lower_au, low_au), min(upper_au, high_au)
    if high <= low:
        return 0.0
    return (high ** k - low ** k) / (upper_au ** k - lower_au ** k)


@dataclass(frozen=True)
class BeltState:
    """
    A belt after some of its history.

    Attributes:
        lower_au, upper_au (float): Its limits now (zero width once it is lost).
        mass_fraction (float): The share of the mass it formed with that it still has.
        scattered_fraction (float): The share of the formed mass thrown inward by its events so far (the pollution source).
        history (tuple): One sentence per event, in order.
    """

    lower_au: float
    upper_au: float
    mass_fraction: float
    scattered_fraction: float
    history: tuple

    @classmethod
    def formed(cls, lower_au, upper_au):
        if not 0.0 < lower_au < upper_au:
            raise ValueError("the belt needs 0 < lower < upper")
        return cls(lower_au, upper_au, 1.0, 0.0, ())

    @property
    def lost(self):
        return self.mass_fraction <= 0.0 or self.upper_au <= self.lower_au


def _lost(state, line):
    return replace(state, upper_au=state.lower_au, mass_fraction=0.0, history=state.history + (line,))


def _percent(share):
    return f"{100.0 * share:.0f}" if share >= 0.01 else f"{100.0 * share:.1f}"


def giant_event(state, kind, stream):
    """
    A giant planet migrates through the belt (`GRAND_TACK`) or the giants go unstable (`INSTABILITY`): the belt keeps a drawn share of its mass.

    Args:
        state (BeltState): The belt.
        kind (str): `GRAND_TACK` or `INSTABILITY`.
        stream: An object with `random()` (a `draw.Stream`); one draw is taken whatever the state, so the run's stream does not depend on the history.
    """
    if kind not in (GRAND_TACK, INSTABILITY):
        raise ValueError("the event is a grand tack or an instability")
    low, high = GRAND_TACK_SURVIVAL if kind == GRAND_TACK else INSTABILITY_SURVIVAL
    survival = log_uniform(low, high, rng=stream)
    if state.lost:
        return state
    what = "A giant planet migrated through the belt" if kind == GRAND_TACK else "The giant planets went unstable"
    line = f"{what}; {_percent(survival)} percent of the belt's mass was left."
    return replace(state, mass_fraction=state.mass_fraction * survival, history=state.history + (line,))


def engulf(state, radius_au):
    """
    The star's envelope reaches `radius_au` (`mass_loss.engulfment_radius_au`): the belt inside it is lost, the part outside keeps its share of the mass.
    """
    if not radius_au > 0.0:
        raise ValueError("the radius must be positive")
    if state.lost or radius_au <= state.lower_au:
        return state
    if radius_au >= state.upper_au:
        return _lost(state, f"The star swelled out to {radius_au:.2f} au and swallowed the whole belt.")
    kept = _share_between(state.lower_au, state.upper_au, radius_au, state.upper_au)
    line = f"The star swelled out to {radius_au:.2f} au; the part of the belt inside it was lost and the outer {_percent(kept)} percent of its mass survived."
    return replace(state, lower_au=radius_au, mass_fraction=state.mass_fraction * kept, history=state.history + (line,))


@finite_domain()
def mass_loss(state, mass_before_solar, mass_after_solar, giants_remain):
    """
    The star loses mass: every orbit, the belt's included, widens by `mass_before / mass_after`; if giants remain, a share of the belt is thrown inward.

    Args:
        state (BeltState): The belt.
        mass_before_solar, mass_after_solar (float): The star's mass before and after.
        giants_remain (bool): Whether the system still has giants to destabilise it.
    """
    if not (mass_before_solar >= mass_after_solar > 0.0):
        raise ValueError("a star loses mass: before at least after, after positive")
    if state.lost or mass_after_solar == mass_before_solar:
        return state
    factor = mass_before_solar / mass_after_solar
    lost_fraction = 1.0 - mass_after_solar / mass_before_solar
    thrown = min(1.0, MASS_LOSS_SCATTER_SHARE * lost_fraction) if giants_remain else 0.0
    line = f"The star lost {_percent(lost_fraction)} percent of its mass and the belt widened to {state.lower_au * factor:.1f} to {state.upper_au * factor:.1f} au"
    line += f"; {_percent(thrown)} percent of it was thrown inward by the giants' resonances." if thrown > 0.0 else "."
    return replace(state, lower_au=state.lower_au * factor, upper_au=state.upper_au * factor, mass_fraction=state.mass_fraction * (1.0 - thrown),
                   scattered_fraction=state.scattered_fraction + state.mass_fraction * thrown, history=state.history + (line,))


def limit(state, inner_au, outer_au, reason):
    """
    A companion limits the belt to `inner_au` to `outer_au` (use 0 and the Holman and Wiegert limit for a companion beyond the belt, the circumbinary limit and infinity for a close pair): the part outside is lost.

    Args:
        state (BeltState): The belt.
        inner_au, outer_au (float): The stable range.
        reason (str): Words for the history line ("a companion at 120 au").
    """
    if not 0.0 <= inner_au < outer_au:
        raise ValueError("the stable range needs 0 <= inner < outer")
    if state.lost or (inner_au <= state.lower_au and outer_au >= state.upper_au):
        return state
    if inner_au >= state.upper_au or outer_au <= state.lower_au:
        return _lost(state, f"The belt lay outside the stable range set by {reason} and was lost.")
    kept = _share_between(state.lower_au, state.upper_au, inner_au, outer_au)
    lower, upper = max(state.lower_au, inner_au), min(state.upper_au, outer_au)
    line = f"{reason[:1].upper()}{reason[1:]} limits the belt to {lower:.1f} to {upper:.1f} au; {_percent(kept)} percent of its mass was inside."
    return replace(state, lower_au=lower, upper_au=upper, mass_fraction=state.mass_fraction * kept, history=state.history + (line,))


@finite_domain()
def pollution_rate_g_s(scattered_mass_kg, delivery_myr=DELIVERY_MYR):
    """The accretion rate of the rock a belt scattered onto a white dwarf, grams per second: its mass over the time it takes to fall."""
    if not (scattered_mass_kg >= 0.0 and delivery_myr > 0.0):
        raise ValueError("the mass must not be negative and the delivery time must be positive")
    return scattered_mass_kg * 1000.0 / (delivery_myr * _SECONDS_PER_MYR)


def is_polluted(rate_g_s):
    """Whether a white dwarf accreting at this rate shows heavy elements (at least the lowest rate measured on one)."""
    return rate_g_s >= MIN_POLLUTION_RATE_G_S


def scattered_mass_kg(state, formed_mass_kg):
    """The mass thrown inward by the belt's events, kg, for a belt that formed with `formed_mass_kg`."""
    return state.scattered_fraction * formed_mass_kg


def mass_kg(state, formed_mass_kg):
    """The mass the belt has left of what it formed with, kg (before the collisional grinding, which `belt_mass` does)."""
    return state.mass_fraction * formed_mass_kg


def history_text(state):
    """The belt's history as one paragraph, or an empty string for a belt nothing has happened to."""
    return " ".join(state.history)
