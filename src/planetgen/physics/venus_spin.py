# planetgen/physics/venus_spin.py

"""
The Venus-type spin: the slow, non-synchronous (often retrograde) rotation of a thick-air planet (GEN.1038)
==========================================================================================================

A rocky planet with a dense air is not driven to synchronous rotation. The star heats the air unevenly; the air's thermal bulge is pulled by the star the other way
from the planet's solid tide, and the two torques balance at spin rates that are **not** the orbital rate. Add the friction between the core and the mantle, which
drives the tilt of the axis to 0 or 180 degrees, and a planet can end up turning slowly backwards (Venus: 243 days retrograde on a 225 day orbit). The model is
[R0896] (*Nature* 411:767; 2003, *Icarus* 163:1 and 24) as generalised to eccentric orbits and tilt by [R0962] (*Int. J.
Astrobiol.* 14:233, [R0962], read in full). This module is that closed-form result as a pure function; nothing calls it yet (GEN.1048 wires it into the
planet's spin and tilt, `weird-worlds.md` W2).

**The strength of the thermal tide.** The two tides are summed into one number, `ws`, the spin offset at which they cancel (their equation 57, circular orbit):

    ws / n  =  1.9255 * (L / L_Venus) * (R / R_Venus) * (a / a_Venus)^2.5 * (M / M_Venus)^-1.5 * (m_Venus / m)

for a planet of mass `m` and radius `R` at distance `a` from a star of mass `M` and luminosity `L` (their equations 55 and 68: the ground-heated thermal torque
over the gravitational torque, in orbital rates `n`). The coefficient is Venus itself (their Table 1: `ws/n` 1.9255, from the 243 day retrograde spin), so the
unknown tide parameters of the planet (the Love number, the two time lags, the scale height and the share of starlight that reaches the ground) are taken to be
Venus's for every thick-air rocky planet, as the paper does ("V-type planets"). A thin air has no such tide: the strength is cut in proportion to the surface
pressure below `FULL_TIDE_PRESSURE_BAR`.

**The equilibria.** With the viscous (linear in frequency) tide model of the paper, `d omega/dt = -K0 g(omega)`, and `g = 0` where (their equation 54, the signs
of `omega - n`, `2 omega - n`, `2 omega - 3 n` on the four ranges of the spin)

    omega = (1 + 6 e^2) n + ws [ (1 - 21 e^2 / 2) s1 - e^2 s2 + 9 e^2 s3 ]

which is one candidate per range, kept when it lies in its own range (equations 62 to 67):

    omega_2+ = n + ws + e^2 (6 n - 5 ws / 2)       above 3 n / 2
    omega_1+ = n + ws + e^2 (6 n - 41 ws / 2)      between n and 3 n / 2
    omega_1- = n - ws + e^2 (6 n + ws / 2)         between n / 2 and n
    omega_2- = n - ws + e^2 (6 n + 5 ws / 2)       below n / 2 (negative: the spin runs backwards)

No state is synchronous unless the thermal tide is nil and the orbit round (the sign jump of the torque at `omega = n` repels otherwise), and at most three coexist. Venus has `omega_2-` = -0.9250 n, 242.9 days retrograde
(Table 1 of the paper: -242.9). A retrograde state exists when `ws / n >= (1 + 6 e^2) / (1 - 5 e^2 / 2)` (about `1 + 17 e^2 / 2`).

**Which state a planet reaches.** The state depends on the tilt the planet is born with (their sections 4.3 and 4.4): the lowest tilts go to the fastest
state, the highest to the slowest. In the two runs of the paper that show the split, two states divide at 60 degrees (a thick-air planet at 0.6 AU M, circular: tilts
over 60 degrees reach the retrograde state), and three states divide at 40 and 100 degrees (0.42 AU M, `e` 0.1). The tilt the planet already has (`spin`, `impact_spin`) therefore decides
the end state with no further draw, and the share of planets that end retrograde is whatever share of them is born tilted past the split. Core-mantle friction then
leaves the axis at 0 or 180 degrees (Venus 177.4), which a caller reads as `final_obliquity_deg`.

**Whether it got there.** The e-folding time of the slow-rotation tide is `tau = m a^6 / (M^2 R^3)` times Venus's 2.3 Gyr over Venus's value of that (equation 48 with
Table 1): a planet whose `tau` is longer than the age has not reached an equilibrium and the module says so (`reached`).

**What it does not do.** Tilt is not integrated (the split rule stands for the paper's simulations). The planet is taken to have no large moon (a moon's tide
holds a planet's spin, as the Earth's), and spin-orbit resonances such as Mercury's 3:2 are `tidal_spin` and `spin_torques`. The paper's eccentricity limit is 0.3
(`PAPER_MAX_ECCENTRICITY`): beyond it the result is computed at 0.3 and flagged. Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.90.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")

VENUS_STRENGTH = 1.9255
"""float: `ws / n` of Venus ([R0962: Table 1], from the 243 day retrograde spin). Source tier T1, tweak level L0."""
VENUS_DAMPING_TIME_GY = 2.3
"""float: The time Venus's slow-rotation tide takes to reach its equilibrium, Gy (the same table). Source tier T1, tweak level L0."""
VENUS_MASS_KG = 4.8675e24
"""float: Venus's mass, kg. Source tier T1, tweak level L0."""
VENUS_RADIUS_M = 6.0518e6
"""float: Venus's mean radius, m. Source tier T1, tweak level L0."""
VENUS_DISTANCE_AU = 0.723332
"""float: Venus's semi-major axis, AU. Source tier T1, tweak level L0."""
FULL_TIDE_PRESSURE_BAR = 10.0
"""float: The surface pressure from which the thermal tide has its full Venus-like strength (the design note's "air above about 10 bar"); below it the strength falls in
proportion to the pressure. Source tier T6, tweak level L1."""
PAPER_MAX_ECCENTRICITY = 0.3
"""float: The paper's limit of validity ("moderate eccentricity, e < 0.3"). Source tier T1, tweak level L0."""
TWO_STATE_SPLIT_DEG = 60.0
"""float: With two end states, the birth tilt above which a planet reaches the slower one (the paper's run at 0.6 AU M). Source tier T3, tweak level L1."""
THREE_STATE_SPLITS_DEG = (40.0, 100.0)
"""tuple: With three end states, the birth tilts that divide the fastest, middle and slowest (the paper's run at 0.42 AU M, `e` 0.1). Source tier T3, tweak level L1."""

LABELS = ("2+", "1+", "1-", "2-")
"""tuple: The four states from the fastest to the slowest."""


class VenusSpinError(ValueError):
    """An input the model cannot use."""


@dataclass(frozen=True)
class Equilibrium:
    """
    One spin state the tides leave a planet in.

    Attributes:
        label (str): `"2+"`, `"1+"`, `"1-"` or `"2-"`.
        ratio (float): `omega / n`, signed: negative is retrograde.
    """

    label: str
    ratio: float

    @property
    def retrograde(self):
        """Whether the planet turns backwards."""
        return self.ratio < 0.0


@dataclass(frozen=True)
class VenusSpin:
    """
    The spin of a thick-air planet.

    Attributes:
        strength (float): `ws / n`, thermal over gravitational tide at the planet.
        damping_time_gy (float): The e-folding time of its slow-rotation tide, Gy.
        reached (bool): Whether the age is longer than that time.
        thin_air (bool): Whether the air is too thin for any thermal tide (strength zero).
        eccentricity_clamped (bool): Whether the eccentricity was above `PAPER_MAX_ECCENTRICITY` and the states were computed at that.
        equilibria (tuple): The `Equilibrium` states, fastest first.
        end_state (Equilibrium): The state the birth tilt leads to (the fastest when the tides have not had time: `reached` is False).
        final_obliquity_deg (float): 0 or 180 when `reached` (the equivalent of a retrograde spin is 180), else the birth tilt.
        period_days (float): The sidereal rotation period of the end state, days (positive), or `None` when `reached` is False.
        retrograde (bool): Whether the end state turns backwards.
    """

    strength: float
    damping_time_gy: float
    reached: bool
    thin_air: bool
    eccentricity_clamped: bool
    equilibria: tuple
    end_state: Equilibrium
    final_obliquity_deg: float
    period_days: float
    retrograde: bool


def _check(value, name):
    if not value > 0.0:
        raise VenusSpinError(f"{name} must be positive")


def strength_over_orbital_rate(mass_kg, radius_m, distance_au, star_mass_kg, star_luminosity_w, surface_pressure_pa):
    """`ws / n`: the thermal tide over the gravitational tide, scaled from Venus (the module docstring), with the pressure cut."""
    for value, name in ((mass_kg, "mass"), (radius_m, "radius"), (distance_au, "distance"), (star_mass_kg, "star mass"),
                        (star_luminosity_w, "star luminosity")):
        _check(value, name)
    if surface_pressure_pa <= 0.0:
        return 0.0
    cut = min(1.0, surface_pressure_pa / (FULL_TIDE_PRESSURE_BAR * 1.0e5))
    return (VENUS_STRENGTH * cut * (star_luminosity_w / constants.SOLAR_LUMINOSITY) * (radius_m / VENUS_RADIUS_M)
            * (distance_au / VENUS_DISTANCE_AU) ** 2.5 * (_kg_to_msun(star_mass_kg)) ** -1.5 * (VENUS_MASS_KG / mass_kg))


def damping_time_gy(mass_kg, radius_m, distance_au, star_mass_kg):
    """The e-folding time of the slow-rotation tide, Gy: `m a^6 / (M^2 R^3)` scaled to Venus's 2.3 Gy."""
    for value, name in ((mass_kg, "mass"), (radius_m, "radius"), (distance_au, "distance"), (star_mass_kg, "star mass")):
        _check(value, name)
    return (VENUS_DAMPING_TIME_GY * (mass_kg / VENUS_MASS_KG) * (distance_au / VENUS_DISTANCE_AU) ** 6
            * (constants.SOLAR_MASS_TO_KG / star_mass_kg) ** 2 * (VENUS_RADIUS_M / radius_m) ** 3)


def equilibria(strength, eccentricity):
    """
    The spin states, fastest first, for a thermal over gravitational tide `strength` (`ws / n`) on an orbit of `eccentricity` (at most `PAPER_MAX_ECCENTRICITY`): the
    candidates of equations 64 to 67 that lie in their own range of spin rate.
    """
    if strength < 0.0:
        raise VenusSpinError("the strength cannot be negative")
    if not 0.0 <= eccentricity <= PAPER_MAX_ECCENTRICITY:
        raise VenusSpinError(f"the eccentricity must be 0 to {PAPER_MAX_ECCENTRICITY}")
    e2, ws = eccentricity ** 2, strength
    # (label, ratio, lowest and highest spin of its range in orbital rates)
    candidates = (("2+", 1.0 + ws + e2 * (6.0 - 2.5 * ws), 1.5, math.inf),
                  ("1+", 1.0 + ws + e2 * (6.0 - 20.5 * ws), 1.0, 1.5),
                  ("1-", 1.0 - ws + e2 * (6.0 + 0.5 * ws), 0.5, 1.0),
                  ("2-", 1.0 - ws + e2 * (6.0 + 2.5 * ws), -math.inf, 0.5))
    return tuple(Equilibrium(label, ratio) for label, ratio, low, high in candidates if low <= ratio < high)


def retrograde_possible(strength, eccentricity):
    """Whether the slowest state runs backwards: `ws / n` at least `(1 + 6 e^2) / (1 - 5 e^2 / 2)`."""
    return strength >= (1.0 + 6.0 * eccentricity ** 2) / (1.0 - 2.5 * eccentricity ** 2)


def end_state(states, birth_obliquity_deg):
    """
    The state a planet born with a tilt of `birth_obliquity_deg` ends in: the only one, or the fastest, middle or slowest by the tilt (the module docstring).
    """
    if not states:
        raise VenusSpinError("no equilibrium state")
    if not 0.0 <= birth_obliquity_deg <= 180.0:
        raise VenusSpinError("the tilt must be 0 to 180 degrees")
    if len(states) == 1:
        return states[0]
    if len(states) == 2:
        return states[1] if birth_obliquity_deg >= TWO_STATE_SPLIT_DEG else states[0]
    low, high = THREE_STATE_SPLITS_DEG
    return states[0] if birth_obliquity_deg < low else states[1] if birth_obliquity_deg < high else states[2]


def venus_type_spin(mass_kg, radius_m, surface_pressure_pa, distance_au, eccentricity, star_mass_kg, star_luminosity_w, age_gy, birth_obliquity_deg):
    """
    The spin a rocky planet ends in under its star's gravitational and thermal tides and core-mantle friction (the module docstring).

    Args:
        mass_kg, radius_m (float): The planet.
        surface_pressure_pa (float): Its air.
        distance_au, eccentricity (float): Its orbit (semi-major axis; eccentricity above `PAPER_MAX_ECCENTRICITY` is computed at that and flagged).
        star_mass_kg, star_luminosity_w (float): The star.
        age_gy (float): The age the tides have had.
        birth_obliquity_deg (float): The tilt it was born with, 0 to 180.

    Returns:
        VenusSpin
    """
    clamped = eccentricity > PAPER_MAX_ECCENTRICITY
    e = min(max(eccentricity, 0.0), PAPER_MAX_ECCENTRICITY)
    strength = strength_over_orbital_rate(mass_kg, radius_m, distance_au, star_mass_kg, star_luminosity_w, surface_pressure_pa)
    tau = damping_time_gy(mass_kg, radius_m, distance_au, star_mass_kg)
    states = equilibria(strength, e)
    reached = age_gy >= tau
    chosen = end_state(states, birth_obliquity_deg)
    if not reached:
        return VenusSpin(strength, tau, False, surface_pressure_pa <= 0.0, clamped, states, chosen, birth_obliquity_deg, None, False)
    mean_motion = math.sqrt(constants.G * star_mass_kg / (distance_au * constants.AU_M) ** 3)
    period_days = 2.0 * math.pi / (abs(chosen.ratio) * mean_motion) / 86400.0
    return VenusSpin(strength, tau, True, surface_pressure_pa <= 0.0, clamped, states, chosen, 180.0 if chosen.retrograde else 0.0, period_days,
                     chosen.retrograde)
