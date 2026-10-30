# planetgen/physics/rogue_kick.py

"""
How fast an ejected planet leaves, and how fast it moves through the Galaxy (GEN.1112, RP6)
============================================================================================

The rogue planets had the speed of a star on a generic galactic orbit and no ejection kick (`rogue-planet-review.md` row 14). A planet thrown out of its system by a
scattering planet leaves the star with a speed at infinity set by how much orbital energy the encounters gave it beyond the little it needed, and it carries that speed,
added to its star's, into the Galaxy. This module gives that speed from the scatterer and the orbit it scatters on, the direction as a held input, and the velocity the rogue has
now relative to the local standard of rest, using the drift (`rogue_drift`) for the years between. Nothing calls it yet (GEN.1111 creates the rogue rows, and RP1 gives each
ejection its scatterer and orbit).

**The speed.** A body on an orbit at distance `a` from a star has the binding energy `GM / 2a`, which is `v_esc(a)^2 / 4` for the escape speed `v_esc = sqrt(2GM/a)` there. A
scattering encounter changes its energy by about `v_orb dv`, where `dv` is the most the scatterer can deflect it, the smaller of the scatterer's own escape speed `v_esc,p` and
the orbital speed `v_orb` (the Safronov number `Theta = v_esc,p^2 / 2 v_orb^2` above 1 means the scatterer ejects and below 1 means it only stirs). So the fractional energy
step of one encounter is `STEP = ENERGY_STEP_SCALE min(1, v_esc,p / v_orb)`. A random walk in energy reaches ejection when it first crosses zero, and the energy it has
beyond zero is then exponentially distributed with a mean of one step: `E_x = -STEP (v_esc^2 / 4) ln(1 - u)`, so

    v_inf = v_esc(a) sqrt(STEP (-ln(1 - u)) / 2)

for a uniform number `u`. With the step scaled to 0.3 (`ENERGY_STEP_SCALE`, tuned) Jupiter ejecting from 5 AU gives a median of 6 km/s and a tail to 15, a giant at 30 AU a median of 2.5 km/s; the design note's stand-in median
of 3 km/s for the whole population and Sayed, Profumo and Smyth's 0.7 to 10.7 km/s (abstract only, `rogue-planet-drift.md` section 2) bracket it. The speed means something only where the
scatterer ejects at all (`ejects`: Safronov number above `EJECTION_SAFRONOV_NUMBER`); an Earth at 1 AU only stirs.

**The direction** is a held input (Boss has not decided it; GEN.1113): `isotropic_direction` is the default, scattering having no preferred direction in the star's frame.

**The velocity through the Galaxy.** The rogue's velocity is its host star's plus the drift's relative velocity (`rogue_drift.drift`): the speed at infinity is the speed relative
to the star at ejection, and the Galaxy's tide turns it over the epicycle after. `rogue_velocity` adds the host's velocity relative to the local standard of rest (the host's
own, from the population's velocity dispersion, drawn by the caller) to that relative velocity.

Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.87.
"""

import math

from planetgen.physics import constants, rogue_drift
from planetgen.util.checks import finite_domain

ENERGY_STEP_SCALE = 0.3
"""float: The fractional orbital-energy step of one scattering encounter in units of `min(1, v_esc,p / v_orb)`: most encounters are distant and deflect much less than the most
the scatterer can. Tuned so that giants scattering between 5 and 30 AU give a median of 2.5 to 6 km/s, the design note's 3 km/s and inside Sayed et al.'s 0.7 to 10.7 km/s. Source tier T8, tweak level L3."""
EJECTION_SAFRONOV_NUMBER = 1.0
"""float: The Safronov number above which a scatterer ejects what it scatters and below which it only stirs (Safronov 1972; Ford and Rasio 2008 use the same cut). Source tier T6, tweak level L0."""
_KM = 1.0e3


class RogueKickError(ValueError):
    """An input the model cannot use."""


def escape_speed_km_s(host_mass_kg, distance_au):
    """The escape speed from a star at a distance, km/s."""
    if not host_mass_kg > 0.0 or not distance_au > 0.0:
        raise RogueKickError("the mass and the distance must be positive")
    return math.sqrt(2.0 * constants.G * host_mass_kg / (distance_au * constants.AU_M)) / _KM


def circular_speed_km_s(host_mass_kg, distance_au):
    """The circular orbital speed round a star at a distance, km/s."""
    return escape_speed_km_s(host_mass_kg, distance_au) / math.sqrt(2.0)


def planet_escape_speed_km_s(planet_mass_kg, planet_radius_m):
    """The escape speed from the surface of the scattering planet, km/s."""
    if not planet_mass_kg > 0.0 or not planet_radius_m > 0.0:
        raise RogueKickError("the planet's mass and radius must be positive")
    return math.sqrt(2.0 * constants.G * planet_mass_kg / planet_radius_m) / _KM


def safronov_number(planet_mass_kg, planet_radius_m, host_mass_kg, distance_au):
    """`Theta = v_esc,p^2 / (2 v_orb^2) = (M_p / M_*)(a / R_p)`: above 1 a scatterer ejects, below 1 it stirs."""
    return (planet_escape_speed_km_s(planet_mass_kg, planet_radius_m) / circular_speed_km_s(host_mass_kg, distance_au)) ** 2 / 2.0


def ejects(planet_mass_kg, planet_radius_m, host_mass_kg, distance_au):
    """Whether the scatterer ejects what it scatters (its Safronov number is above `EJECTION_SAFRONOV_NUMBER`)."""
    return safronov_number(planet_mass_kg, planet_radius_m, host_mass_kg, distance_au) > EJECTION_SAFRONOV_NUMBER


def energy_step(planet_mass_kg, planet_radius_m, host_mass_kg, distance_au):
    """The fractional orbital-energy step of one encounter, `ENERGY_STEP_SCALE min(1, v_esc,p / v_orb)`."""
    ratio = planet_escape_speed_km_s(planet_mass_kg, planet_radius_m) / circular_speed_km_s(host_mass_kg, distance_au)
    return ENERGY_STEP_SCALE * min(1.0, ratio)


@finite_domain()
def ejection_speed_km_s(u, host_mass_kg, distance_au, planet_mass_kg, planet_radius_m):
    """
    The speed at infinity relative to the star of a body ejected by a scatterer, km/s.

    Args:
        u (float): A uniform random number, 0 to 1 (excluding 1).
        host_mass_kg (float): The star's mass.
        distance_au (float): Where the scattering happens: the scatterer's distance from the star.
        planet_mass_kg, planet_radius_m (float): The scatterer's.
    """
    if not 0.0 <= u < 1.0:
        raise RogueKickError("u must be at least 0 and under 1")
    step = energy_step(planet_mass_kg, planet_radius_m, host_mass_kg, distance_au)
    return escape_speed_km_s(host_mass_kg, distance_au) * math.sqrt(step * -math.log1p(-u) / 2.0)


def isotropic_direction(u_cos, u_phi):
    """A direction uniform on the sphere `(outward, along the rotation, up)` from two uniform numbers, 0 to 1."""
    if not (0.0 <= u_cos <= 1.0 and 0.0 <= u_phi <= 1.0):
        raise RogueKickError("the random numbers must be from 0 to 1")
    z = 2.0 * u_cos - 1.0
    r = math.sqrt(max(0.0, 1.0 - z * z))
    phi = 2.0 * math.pi * u_phi
    return (r * math.cos(phi), r * math.sin(phi), z)


def kick(speed_km_s, direction):
    """The kick `(u, v, w)` of a speed in a direction, in `rogue_drift`'s frame."""
    return rogue_drift.kick_from_direction(speed_km_s, direction)


def rogue_velocity(host_velocity_km_s, kick_km_s, time_myr, galaxy=rogue_drift.SOLAR_NEIGHBOURHOOD):
    """
    The rogue's velocity relative to the local standard of rest, `time_myr` after the kick: the host star's velocity plus the velocity relative to the star that the drift gives.

    Args:
        host_velocity_km_s (sequence): The host star's `(U, V, W)` relative to the local standard of rest now, km/s: outward, along the rotation, up.
        kick_km_s (sequence): The kick at ejection.
        time_myr (float): The time since ejection.
    """
    if len(host_velocity_km_s) != 3 or any(not math.isfinite(c) for c in host_velocity_km_s):
        raise RogueKickError("a host velocity is three finite components")
    here = rogue_drift.drift(kick_km_s, time_myr, galaxy)
    return (host_velocity_km_s[0] + here.vx_km_s, host_velocity_km_s[1] + here.vy_km_s, host_velocity_km_s[2] + here.vz_km_s)


def relative_speed_km_s(velocity, reference=(0.0, 0.0, 0.0)):
    """The speed of a body with a velocity relative to an observer with the reference velocity (both relative to the same standard of rest), km/s."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(velocity, reference)))
