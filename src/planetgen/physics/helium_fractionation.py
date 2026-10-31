# planetgen/physics/helium_fractionation.py

"""
Helium enrichment of an envelope that loses hydrogen first (GEN.1045, W16 of docs/design/systems/weird-worlds.md, registry 6.77)
==================================================================================================================================

Audit of `envelope_escape` (GEN.822): it follows the *mass* of a hydrogen-helium envelope that the XUV removes (the hydro-based rate of Kubyshkina
et al. 2018), and its composition never changes: the envelope is solar at every step, so a planet that has lost 95 percent of its envelope
still has a solar mix and the helium-dominated air of Hu et al. 2015 and Malsky et al. 2023 cannot come out of it. This module adds the missing
step as a pure function of the loss history. It changes nothing in `envelope_escape` and no generated output.

Hydrogen drags helium out only while the escape flux is large enough; a smaller flux lets hydrogen leave and helium stay. With the homopause
(where molecular diffusion takes over from mixing) at `HOMOPAUSE_RADIUS_OVER_TRANSIT` times the transit radius and at a temperature T:

    b = 1.04e18 T^0.732 cm^-1 s^-1          binary diffusion of H in He      Hu 2015 eq. 9
    b' = b / (1 + x)                        with the ionised share x (row HF2)
    phi_DL = G M (m_He - m_H) b' / (R^2 k T)     the diffusion-limited number flux of hydrogen          Hu 2015 eq. 12
    phi_He / X_He = phi_H / X_H - phi_DL         the fluxes of the two gases                            Hu 2015 eq. 13

A total mass flux Phi = 4 pi R^2 (phi_H m_H + phi_He m_He) at or under `4 pi R^2 X_H m_H phi_DL` is all hydrogen (eq. 14); over it the two
leave in the proportions of eq. 15. The mixing ratios X are by number. `evolve` takes the envelope a step at a time at a given total loss, so the
hydrogen share falls and the threshold falls with it: a thin envelope that is still losing mass is helium-rich.

The total loss comes from outside: a history of `(age, envelope)` pairs, which `envelope_escape.envelope_left_earth` gives at a set of ages.

What the papers give and the module reproduces: a planet must lose at least half of its starting envelope to reach a helium mass fraction of 0.4
and be about 2 Gyr old (Malsky et al. 2023, Introduction and B.1); with a thin envelope (10^-3 of the planet's mass or less) hydrogen falls by
orders of magnitude in a few billion years (Hu et al. 2015); a 10 Earth-mass planet with 0.5 percent envelope has lost just over half of it
after 5 Gyr at 1e9 g/s (Malsky et al. 2023, B.4).

Use:

    history = [(age, envelope_escape.envelope_left_earth(core, f_env0, insolation, host_mass, luminosity, age)) for age in ages]
    envelope, y = helium_fraction_from_history(core, history, insolation)
    envelope_kg, y = evolve(envelope_kg, y_start, mass_kg, radius_m, rate_kg_s, duration_s)
"""

import math

from planetgen.physics import bulk_structure, constants, envelope_escape, unit_registry as units

_mearth_to_kg = units.converter("Mearth", "kg")

SOLAR_ENVELOPE_HELIUM = 0.245
"""float: The helium share of the mass of a hydrogen-helium envelope at the start, `Y / (X + Y)` of the solar mix (Malsky et al. 2023, A: X = 0.74,
Y = 0.24; row HF3)."""

HOMOPAUSE_TEMPERATURE_K = 1.0e4
"""float: The temperature at the homopause, K. Malsky et al. 2023 run 3,000 and 10,000 K and find helium-rich planets at both; Hu et al. 2015 use
10^4 K, the weaker fractionation (the diffusion-limited flux goes as T^-0.27), row HF2."""

IONISED_SHARE = 0.1
"""float: The hydrogen ionised share `x` at the homopause (Hu et al. 2015 Sec. 2.1 and Malsky and Rogers 2020: 0.1), row HF2."""

HOMOPAUSE_RADIUS_OVER_TRANSIT = 1.1
"""float: The homopause radius over the transit radius ("approximately 10 percent larger", Malsky et al. 2023 A.1), row HF2."""

DIFFUSION_A_CM = 1.04e18
"""float: `b = A T^0.732` in cm^-1 s^-1 (as used by Hu et al. 2015 and Malsky et al. 2023), row HF1."""

DIFFUSION_EXPONENT = 0.732
"""float: The temperature exponent of `b`, row HF1."""

HELIUM_ATOM_MASS_RATIO = 3.9726
"""float: The mass of a helium atom over a hydrogen atom."""

MIN_HELIUM_FOR_ENHANCED = 0.4
"""float: A planet is helium-enhanced from this helium mass fraction (Malsky et al. 2023, Introduction, row HF3)."""

STEP_SHARE = 0.02
"""float: The largest share of the envelope one step of `evolve` removes (numerical)."""

_M_H = constants.HYDROGEN_ATOM_MASS_KG
_M_HE = HELIUM_ATOM_MASS_RATIO * _M_H
_GY_S = 1.0e9 * constants.SECONDS_PER_YEAR
_R_EARTH_M = constants.EARTH_RADIUS_KM * 1.0e3


def binary_diffusion_m(temperature_k):
    """`b`, the binary diffusion coefficient of hydrogen in helium, m^-1 s^-1 (1 cm^-1 is 100 m^-1)."""
    return DIFFUSION_A_CM * temperature_k ** DIFFUSION_EXPONENT * 100.0


def effective_diffusion_m(temperature_k, ionised_share=IONISED_SHARE):
    """`b'`, m^-1 s^-1: the coefficient with ion-neutral drag, `b / (1 + x)`: at 10^4 K and x = 0.1 it is 8.0e20 cm^-1 s^-1, the value Hu et al.
    2015 give (row HF2)."""
    return binary_diffusion_m(temperature_k) / (1.0 + ionised_share)


def diffusion_limited_flux(mass_kg, radius_m, temperature_k=HOMOPAUSE_TEMPERATURE_K, ionised_share=IONISED_SHARE):
    """`phi_DL`, hydrogen atoms per m^2 per s (Hu et al. 2015 eq. 12): the most hydrogen flux that leaves helium behind. `radius_m` is the
    homopause radius."""
    return (constants.G * mass_kg * (_M_HE - _M_H) * effective_diffusion_m(temperature_k, ionised_share)
            / (radius_m ** 2 * constants.BOLTZMANN * temperature_k))


def number_fractions(helium_mass_fraction):
    """`(X_H, X_He)`, the mixing ratios by number of an envelope of hydrogen and helium of the helium mass fraction `y`."""
    hydrogen = (1.0 - helium_mass_fraction) / _M_H
    helium = helium_mass_fraction / _M_HE
    total = hydrogen + helium
    return hydrogen / total, helium / total


def critical_mass_flux_kg_s(mass_kg, radius_m, helium_mass_fraction, temperature_k=HOMOPAUSE_TEMPERATURE_K, ionised_share=IONISED_SHARE):
    """The total mass loss, kg/s, at and under which only hydrogen leaves (Hu et al. 2015 eq. 14): `4 pi R^2 X_H m_H phi_DL`."""
    x_h, _ = number_fractions(helium_mass_fraction)
    return 4.0 * math.pi * radius_m ** 2 * x_h * _M_H * diffusion_limited_flux(mass_kg, radius_m, temperature_k, ionised_share)


def partition(total_kg_s, mass_kg, radius_m, helium_mass_fraction, temperature_k=HOMOPAUSE_TEMPERATURE_K, ionised_share=IONISED_SHARE):
    """`(hydrogen kg/s, helium kg/s)` of a total mass loss `total_kg_s` (Hu et al. 2015 eqs. 13 to 15): all hydrogen at and under
    `critical_mass_flux_kg_s`, the two in the proportions of eq. 15 over it."""
    if total_kg_s <= 0.0:
        return 0.0, 0.0
    x_h, x_he = number_fractions(helium_mass_fraction)
    area = 4.0 * math.pi * radius_m ** 2
    phi_dl = diffusion_limited_flux(mass_kg, radius_m, temperature_k, ionised_share)
    helium_over_x = (total_kg_s / area - x_h * _M_H * phi_dl) / (x_h * _M_H + x_he * _M_HE)
    if helium_over_x <= 0.0:
        return total_kg_s, 0.0
    helium = area * helium_over_x * x_he * _M_HE
    return total_kg_s - helium, helium


def evolve(envelope_kg, helium_mass_fraction, mass_kg, radius_m, rate_kg_s, duration_s, temperature_k=HOMOPAUSE_TEMPERATURE_K,
           ionised_share=IONISED_SHARE):
    """
    The envelope and its helium mass fraction after `duration_s` of a total loss `rate_kg_s` (held for the interval). Steps of at most
    `STEP_SHARE` of the envelope; the envelope is gone (`(0.0, helium fraction)`) if the loss outlasts it. `radius_m` is the homopause radius.
    """
    hydrogen = envelope_kg * (1.0 - helium_mass_fraction)
    helium = envelope_kg * helium_mass_fraction
    elapsed = 0.0
    while elapsed < duration_s and hydrogen + helium > 0.0:
        total = hydrogen + helium
        y = helium / total
        h_rate, he_rate = partition(rate_kg_s, mass_kg, radius_m, y, temperature_k, ionised_share)
        step = duration_s - elapsed
        if rate_kg_s > 0.0:
            step = min(step, STEP_SHARE * total / rate_kg_s)
        hydrogen = max(hydrogen - h_rate * step, 0.0)
        helium = max(helium - he_rate * step, 0.0)
        elapsed += step
        if h_rate <= 0.0 and he_rate <= 0.0:
            break
    total = hydrogen + helium
    return total, (helium / total if total > 0.0 else helium_mass_fraction)


def helium_fraction_from_history(core_mass_earth, history, insolation, start=SOLAR_ENVELOPE_HELIUM, temperature_k=HOMOPAUSE_TEMPERATURE_K,
                                 ionised_share=IONISED_SHARE):
    """
    The envelope left, Earth masses, and its helium mass fraction at the end of a loss `history` of `(age_gy, envelope_earth)` pairs in rising age
    (the first is the start; `envelope_escape.envelope_left_earth` at a few ages gives one). The loss between two ages is held as a constant total
    rate and split by `partition`; the envelope follows the history and the split only sets the composition. Returns `(envelope_earth, y)`; a
    history that reaches zero is a bare core, `(0.0, start)`.
    """
    y = start
    age, envelope = history[0]
    for next_age, left in history[1:]:
        if left <= 0.0:
            return 0.0, start
        lost = envelope - left
        if lost > 0.0:
            planet = core_mass_earth + 0.5 * (envelope + left)
            radius = bulk_structure.envelope_radius_earth(planet, 0.5 * (envelope + left) / planet, insolation, next_age)
            duration = (next_age - age) * _GY_S
            _, y = evolve(_mearth_to_kg(envelope), y, _mearth_to_kg(planet),
                          radius * _R_EARTH_M * HOMOPAUSE_RADIUS_OVER_TRANSIT, _mearth_to_kg(lost) / duration, duration,
                          temperature_k, ionised_share)
        age, envelope = next_age, left
    return envelope, y


def is_helium_enhanced(helium_mass_fraction):
    """Whether the envelope is helium-enhanced (Y of 0.4 or more)."""
    return helium_mass_fraction >= MIN_HELIUM_FOR_ENHANCED
