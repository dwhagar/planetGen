# planetgen/physics/obliquity_variation.py

"""
How far a planet's tilt wanders under the pull of the other bodies of its system (GEN.1454)
==========================================================================================

A planet's spin axis precesses about its orbit's normal because the star (and each moon) pulls on the planet's equatorial bulge, and the orbit's normal itself nods because the other
planets (and a companion star) pull on the orbit. The tilt of the axis to the orbit, the **obliquity** `eps`, stays put while the axis follows the orbit smoothly, and wanders when the
two motions are in step. This module works out how far, from the same quantities the secular exchange of GEN.1450 (`secular_chain`) already holds: the planets' masses and orbits, a
companion star, and the moons.

**The spin side.** The axis turns about the orbit normal at `alpha cos(eps)`, with the precession constant

    alpha = (3 / 2) (E / omega) sum_i G m_i / (a_i^3 (1 - e_i^2)^(3/2))

over the bodies `i` that pull on the bulge (the star, each moon, a companion), `omega` the spin rate and `E = (C - A) / C` the dynamical ellipticity of the planet. For a world that has
reached the shape its spin gives it, `E = J2 / (C / M R^2)` with `J2 = k_f q / 3`, `q = omega^2 R^3 / G M` and the fluid Love number from the Darwin-Radau relation,
`k_f = (3/2) / (1 + (5/2) (1 - (3/2) C / M R^2)^2)` (the Earth gives `E` = 0.00319 against the 0.00327 measured). A large moon raises `alpha` several-fold, which is what steadies the Earth.

**The orbit side.** The nodes of the planets' orbits regress with the frequencies `s_k` of the nodal eigenmodes of the Laplace-Lagrange matrix (`secular_chain.matrices`), the
planet `j` carrying the amplitude `I_jk` (radians) of each mode. A companion star in an inclined orbit adds one mode of its own: the planet's orbit precesses about the companion's
orbit at `s = -(3/4) n (M_c / M) (a / a_c)^3 (1 - e_c^2)^(-3/2) cos(i)` with the amplitude `sin(i)` of the mutual inclination `i` (the exchange with the siblings and the companion's
mode are taken as separate; they couple only through a mutual inclination, and a close companion that the generator aligned with the planets (`star_orbit_alignment`) leaves it
small).

**What follows.** With `X = cos(eps)`, the Hamiltonian of the axis near a mode's frequency is the pendulum of [R2571] with the resonance at
`alpha X = -s_k` and a half-width in `X` of

    W_k = 2 sqrt(2 I_k |s_k| sin(eps_k) / alpha),      cos(eps_k) = -s_k / alpha.

* **Away from every resonance** (`X` outside every `[X_k - W_k, X_k + W_k]`) the axis follows the orbit: the obliquity oscillates by `I_k |s_k| / |alpha X + s_k|` for each mode
  (the orbit's own nod when `alpha` is small against `|s_k|`, nothing when it is large), added in quadrature.
* **Inside one resonance** it librates over the zone `[X_k - W_k, X_k + W_k]`.
* **Where two zones overlap** (`W_j + W_k >= |X_j - X_k|`, Chirikov's criterion) the motion is chaotic over the union of the zones that chain together, and the obliquity over a long
  time covers it. This is the Moon-less Earth and Mars of [R5227]: a planet with few or no moons beside a heavy neighbour, whose `alpha` is small enough
  that several `alpha X = -s_k` fall inside its range of tilt.

`analyse` returns the `Variation` of a planet; `chaotic_tilt_deg` draws the tilt a chaotic planet has at a given time from its range. Nothing here is fitted: every number is a
quantity of the system or a constant of the theory.

**What it does not do.** The theory is the linear one: the amplitudes are small angles (the exchange holds the same limits, `secular_chain.LINEAR_INCLINATION_LIMIT_RAD`); a
mode's resonance is taken as isolated before the overlap test (the true boundary of a chaotic zone is rougher); the spin rate and the moons are taken as they are now (a planet that
has tidally despun, or whose big moon formed late, has had other values; `alpha` follows the spin down as `omega^1`, so the past of a slow spinner is more chaotic than its present
shows); the companion's mode is added without the coupling to the siblings' modes; and the dynamical ellipticity is the fluid one (a stiff planet carries more, Mars 0.0054 against 0.0041, which moves `alpha` and the resonances with it). The theory is also the lower
bound of the chaos: the Moon-less Earth of Laskar and others [R5227] wanders from 0 to 85 degrees through combination terms of the modes that the isolated pendulums here do not have, and comes
out as a 17 degree swing, and Mars comes out resonant (18 to 37 degrees at its measured `alpha`) where the full theory finds it chaotic from 0 to 60. The Earth with its Moon comes out at
22.6 to 24.3 degrees against the 22.1 to 24.5 measured. A body locked to its star or resonant with it takes the Cassini state instead (`tidal_spin.cassini_tilt_deg`), not this. The registered stand-ins are in `docs/design/research/gap-filling-plan.md`, section 6.271.
"""

import math
from dataclasses import dataclass

import numpy as np

from planetgen.physics import constants, secular_chain

ARCSEC_PER_RAD = 180.0 * 3600.0 / math.pi
"""float: Arcseconds in a radian, for reading the rates in the units the literature prints them in."""

MODE_FLOOR_RAD_YR = 1.0e-12
"""float: A nodal frequency below this is the invariable plane's (zero) mode, which does not move the orbit normal against the spin, and is dropped."""
MODE_AMPLITUDE_FLOOR_RAD = 1.0e-9
"""float: A mode with an amplitude below this (a plane that is already flat) forces nothing and is dropped."""


@dataclass(frozen=True)
class Torquer:
    """A body that pulls on the planet's equatorial bulge: its mass (kg), its semi-major axis (m) and its eccentricity about the planet (the star, a moon, a companion)."""

    mass_kg: float
    distance_m: float
    eccentricity: float = 0.0


@dataclass(frozen=True)
class Mode:
    """A nodal mode that drives the planet's orbit normal: the frequency `s` (radians per year, negative for regression) and the amplitude `I` of the orbit's inclination (radians)."""

    frequency_rad_yr: float
    amplitude_rad: float


@dataclass(frozen=True)
class Resonance:
    """One mode's resonance with the spin: the obliquity `eps_k` (degrees) where `alpha cos(eps) = -s`, and the half-width `W_k` of its zone in `cos(eps)`."""

    mode: Mode
    obliquity_deg: float
    cos_center: float
    cos_half_width: float


@dataclass(frozen=True)
class Variation:
    """
    How far a planet's obliquity wanders.

    Attributes:
        minimum_deg (float): The smallest obliquity over a long time, degrees.
        maximum_deg (float): The largest.
        chaotic (bool): Whether the tilt is in a chaotic zone (two or more overlapping resonances), so any value in the range is reached.
        resonant (bool): Whether the tilt lies inside a resonance's zone, chaotic or not.
        forced_deg (float): The smooth oscillation of the tilt away from the resonances (the quadrature sum of the modes').
        resonances (tuple): The `Resonance` of every mode that a tilt can reach (the resonance falls at `|cos| <= 1`).
        precession_rad_yr (float): The precession constant `alpha`.
    """

    minimum_deg: float
    maximum_deg: float
    chaotic: bool
    resonant: bool
    forced_deg: float
    resonances: tuple
    precession_rad_yr: float


def fluid_love_number(moment_factor):
    """The Darwin-Radau fluid Love number `k_f = (3/2) / (1 + (5/2) (1 - (3/2) C / M R^2)^2)` of a body of the moment factor `C / M R^2` (0.917 for the Earth's 0.3307)."""
    if not 0.0 < moment_factor <= 0.4:
        raise ValueError("the moment factor must be between 0 and 0.4")
    return 1.5 / (1.0 + 2.5 * (1.0 - 1.5 * moment_factor) ** 2)


def dynamical_ellipticity(mass_kg, radius_m, spin_rate_rad_s, moment_factor):
    """`E = (C - A) / C` of a planet whose figure is the fluid one its spin gives it: `J2 / (C / M R^2)` with `J2 = k_f q / 3`, `q = omega^2 R^3 / G M`."""
    if not (mass_kg > 0.0 and radius_m > 0.0 and spin_rate_rad_s >= 0.0):
        raise ValueError("mass and radius must be positive and the spin rate not negative")
    q = spin_rate_rad_s ** 2 * radius_m ** 3 / (constants.G * mass_kg)
    return fluid_love_number(moment_factor) * q / 3.0 / moment_factor


def precession_constant_rad_yr(ellipticity, spin_rate_rad_s, torquers):
    """`alpha = (3/2) (E / omega) sum G m / (a^3 (1 - e^2)^(3/2))`, radians per year (the rate the axis turns about the orbit normal is `alpha cos(eps)`)."""
    if not spin_rate_rad_s > 0.0:
        return 0.0
    total = 0.0
    for body in torquers:
        if not (body.mass_kg > 0.0 and body.distance_m > 0.0 and 0.0 <= body.eccentricity < 1.0):
            raise ValueError("a torquer needs a positive mass and distance and an eccentricity under 1")
        total += constants.G * body.mass_kg / (body.distance_m ** 3 * (1.0 - body.eccentricity ** 2) ** 1.5)
    return 1.5 * ellipticity / spin_rate_rad_s * total * constants.SECONDS_PER_YEAR


def all_sibling_modes(a_au, masses_kg, host_mass_kg, inclination_rad, node_rad):
    """
    The nodal modes that drive each planet of a star's planets (semi-major axes in AU, masses in kg, inclinations and nodes of the orbits in radians, to a common plane), one tuple of
    `Mode` for each planet: the frequencies of the Laplace-Lagrange nodal matrix and the amplitude each carries in that planet's inclination. The invariable plane's own mode (frequency
    zero) is left out. The matrix is made once for all the planets.
    """
    count = len(a_au)
    if count < 2:
        return tuple(() for _ in range(count))
    _a, B = secular_chain.matrices(a_au, masses_kg, host_mass_kg)
    frequencies, vectors, inverse = secular_chain._modes(B)
    start = np.sin(np.asarray(inclination_rad, float)) * np.exp(1j * np.asarray(node_rad, float))
    coefficients = inverse @ start
    out = []
    for index in range(count):
        modes = []
        for k, frequency in enumerate(frequencies):
            amplitude = float(abs(vectors[index, k] * coefficients[k]))
            if abs(frequency) > MODE_FLOOR_RAD_YR and amplitude > MODE_AMPLITUDE_FLOOR_RAD:
                modes.append(Mode(float(frequency), math.asin(min(amplitude, 1.0))))
        out.append(tuple(modes))
    return tuple(out)


def companion_mode(a_au, host_mass_kg, companion_mass_kg, companion_a_au, companion_eccentricity, mutual_inclination_rad):
    """
    The nodal mode a companion star drives in a planet at `a_au` round a host: the planet's orbit precesses about the companion's orbit at `-(3/4) n (M_c / M) (a / a_c)^3
    (1 - e_c^2)^(-3/2) cos(i)` with the amplitude `i` of the mutual inclination. `None` when the orbits are coplanar (nothing to drive) or the companion is too close to treat as a
    perturbation (inside five times the planet's distance).
    """
    if not (a_au > 0.0 and host_mass_kg > 0.0 and companion_mass_kg > 0.0 and companion_a_au > 5.0 * a_au and 0.0 <= companion_eccentricity < 1.0):
        return None
    amplitude = abs(math.sin(mutual_inclination_rad))
    if amplitude <= MODE_AMPLITUDE_FLOOR_RAD:
        return None
    n = math.sqrt(constants.G * host_mass_kg / (a_au * constants.AU_TO_M) ** 3) * constants.SECONDS_PER_YEAR
    rate = -0.75 * n * companion_mass_kg / host_mass_kg * (a_au / companion_a_au) ** 3 / (1.0 - companion_eccentricity ** 2) ** 1.5 * math.cos(mutual_inclination_rad)
    return Mode(rate, math.asin(min(amplitude, 1.0)))


def resonances_of(alpha_rad_yr, modes):
    """The `Resonance` of each mode whose condition `alpha cos(eps) = -s` can be met (`|s| <= alpha`); empty when `alpha` is not positive."""
    if not alpha_rad_yr > 0.0:
        return ()
    out = []
    for mode in modes:
        center = -mode.frequency_rad_yr / alpha_rad_yr
        if abs(center) >= 1.0:
            continue
        sine = math.sqrt(1.0 - center * center)
        half_width = 2.0 * math.sqrt(2.0 * mode.amplitude_rad * abs(mode.frequency_rad_yr) * sine / alpha_rad_yr)
        out.append(Resonance(mode, math.degrees(math.acos(center)), center, half_width))
    return tuple(sorted(out, key=lambda r: r.cos_center))


def _chains(resonances):
    """The resonances grouped into chains whose zones overlap (Chirikov: the sum of two half-widths reaches the distance of the centers), each chain a list in order of `cos`."""
    chains = []
    for resonance in resonances:
        if chains and resonance.cos_half_width + chains[-1][-1].cos_half_width >= resonance.cos_center - chains[-1][-1].cos_center:
            chains[-1].append(resonance)
        else:
            chains.append([resonance])
    return chains


def _clamped_cos(value):
    return max(-1.0, min(1.0, value))


def analyse(alpha_rad_yr, obliquity_deg, modes):
    """
    The `Variation` of a planet of the precession constant `alpha` and the present obliquity (degrees) under the nodal `modes` (`Mode` list).

    A planet with no modes (a lone planet round a lone star) has none: its range is its tilt. A planet whose `alpha` is not positive (no torque on the bulge, or no spin) has the tilt
    fixed in space, so the orbit's own nod is the whole of it.
    """
    if not 0.0 <= obliquity_deg <= 180.0:
        raise ValueError("the obliquity must be between 0 and 180 degrees")
    modes = tuple(modes)
    if not modes:
        return Variation(obliquity_deg, obliquity_deg, False, False, 0.0, (), max(alpha_rad_yr, 0.0))
    if not alpha_rad_yr > 0.0:
        swing = math.degrees(math.sqrt(sum(m.amplitude_rad ** 2 for m in modes)))
        return Variation(max(0.0, obliquity_deg - swing), min(180.0, obliquity_deg + swing), False, False, swing, (), max(alpha_rad_yr, 0.0))
    x = math.cos(math.radians(obliquity_deg))
    resonances = resonances_of(alpha_rad_yr, modes)
    for chain in _chains(resonances):
        low = _clamped_cos(min(r.cos_center - r.cos_half_width for r in chain))
        high = _clamped_cos(max(r.cos_center + r.cos_half_width for r in chain))
        if low <= x <= high:
            ceiling, floor = math.degrees(math.acos(low)), math.degrees(math.acos(high))      # a larger cos is a smaller tilt
            return Variation(floor, ceiling, len(chain) >= 2, True, 0.0, resonances, alpha_rad_yr)
    swing = math.degrees(math.sqrt(sum((m.amplitude_rad * abs(m.frequency_rad_yr) / abs(alpha_rad_yr * x + m.frequency_rad_yr)) ** 2 for m in modes)))
    return Variation(max(0.0, obliquity_deg - swing), min(180.0, obliquity_deg + swing), False, False, swing, resonances, alpha_rad_yr)


def chaotic_tilt_deg(variation, u):
    """
    The obliquity of a planet in a chaotic zone at some time, degrees, for the unit draw `u`: the time spent is even in `cos(eps)` over the zone (the axis wanders over the zone's phase
    space, and the area of the sphere between two tilts is proportional to the difference of the cosines), so the tilt is `acos` of a uniform in the cosines of the range's ends.
    """
    if not 0.0 <= u <= 1.0:
        raise ValueError("u must be between 0 and 1")
    nearest, farthest = math.cos(math.radians(variation.minimum_deg)), math.cos(math.radians(variation.maximum_deg))
    return math.degrees(math.acos(nearest + u * (farthest - nearest)))
