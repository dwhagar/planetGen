# planetgen/physics/brown_dwarf.py

"""
A brown dwarf as a star (GEN.252): its mass, how bright and hot it is at an age, and the L, T or Y spectral type that
temperature makes.

A brown dwarf is too light to hold hydrogen fusion, so it only cools. Burrows et al. 2001 (Rev. Mod. Phys. 73:719, eqs. 1
and 2) give the old-age cooling as two power laws, `L = 4e-5 Lsun (t / 1 Gy)^-1.3 (M / 0.05 Msun)^2.64` and
`Teff = 1550 K (t / 1 Gy)^-0.32 (M / 0.05 Msun)^0.83` (solar metallicity). The first is already the rogue planets'
(`tuning.ROGUE_BROWN_DWARF_*`); the two together fix the radius, `R = sqrt(L / 4 pi sigma Teff^4)`, which is nearly
the same at every age (`t^-0.01`) and falls slowly with mass (`M^-0.34`): 0.088 solar radii at 0.05 solar masses.
A brown dwarf younger than a few hundred million years is larger and hotter than these old-age laws give, so the
temperature is capped at `tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K` (a late-M dwarf).

The spectral type is the mean relation of Kirkpatrick et al. 2021 (ApJS 253:7, Table 13, from 525 L, T and Y dwarfs
within 20 pc), inverted: L0 is 2,238 K, T0 1,255 K, Y0 461 K. Their relation has a scatter of 134 K among L dwarfs and 79 K
among T and Y dwarfs; that is measurement scatter and is not drawn, so a brown dwarf of a given temperature always has
the one type. Above L0 it is a late-M dwarf, named by the stellar scale (`stellar_evolution.spectral_letter_and_subclass`).

Masses follow the substellar mass function of the same census, `dN/dM ~ M^-0.6` (alpha 0.6 +- 0.1).

Left out: a young brown dwarf's larger radius and its gravity-lowered type, the L/T transition's cloud-clearing
plateau (about 18 K per subtype from L8.75 to T4.75, which the piecewise fit already flattens), metallicity and
dust opacity.
"""

import math

from planetgen import tuning
from planetgen.physics import constants, stellar_evolution as evolution
from planetgen.util import draw

REFERENCE_MASS_SOL = 0.05
"""float: The mass the Burrows laws are scaled to, solar masses."""


def mass_range_solar():
    """`(lowest, highest)` brown dwarf mass in solar masses: the deuterium limit to the hydrogen limit."""
    low, high = (jupiter * constants.JUPITER_MASS_TO_KG / constants.SOLAR_MASS_TO_KG
                 for jupiter in tuning.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER)
    return low, high


def sample_mass_solar(rng=draw):
    """A brown dwarf's mass in solar masses, from `dN/dM ~ M^-alpha` (`tuning.BROWN_DWARF_MASS_FUNCTION_SLOPE`) between
    `mass_range_solar`."""
    low, high = mass_range_solar()
    power = 1.0 - tuning.BROWN_DWARF_MASS_FUNCTION_SLOPE
    return (low ** power + rng.random() * (high ** power - low ** power)) ** (1.0 / power)


def luminosity_sol(mass_sol, age_gy):
    """Burrows et al. 2001's old-age luminosity, solar luminosities (the rogue brown dwarfs' law)."""
    return (tuning.ROGUE_BROWN_DWARF_LUMINOSITY_SOL * age_gy ** tuning.ROGUE_GIANT_COOLING_TIME_EXPONENT
            * (mass_sol / REFERENCE_MASS_SOL) ** tuning.ROGUE_BROWN_DWARF_MASS_EXPONENT)


def old_age_temperature_k(mass_sol, age_gy):
    """Burrows et al. 2001's old-age effective temperature, kelvin, before the young cap."""
    return (tuning.BROWN_DWARF_TEMPERATURE_AT_1_GY_K * age_gy ** tuning.BROWN_DWARF_TEMPERATURE_AGE_EXPONENT
            * (mass_sol / REFERENCE_MASS_SOL) ** tuning.BROWN_DWARF_TEMPERATURE_MASS_EXPONENT)


def effective_temperature_k(mass_sol, age_gy):
    """Effective temperature, kelvin: the old-age law, no hotter than `tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K`."""
    return min(old_age_temperature_k(mass_sol, age_gy), tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K)


def radius_sol(mass_sol):
    """Radius in solar radii, from the two Burrows laws at 1 Gy (`t^-0.01`, so age is left out): the luminosity
    is `4 pi R^2 sigma Teff^4`."""
    luminosity_w = luminosity_sol(mass_sol, 1.0) * constants.SOLAR_LUMINOSITY
    temperature = old_age_temperature_k(mass_sol, 1.0)
    radius_m = math.sqrt(luminosity_w / (constants.FOUR_PI * constants.STEFAN_BOLTZMANN_CONSTANT * temperature ** 4))
    return radius_m / constants.SOLAR_RADIUS_M


def _fit_temperature_k(coefficients, number):
    return sum(coefficient * number ** power for power, coefficient in enumerate(coefficients))


def type_temperature_k(number):
    """The mean effective temperature of the spectral subclass number `number` (0 = L0, 10 = T0, 20 = Y0), kelvin
    (`tuning.BROWN_DWARF_SPECTRAL_TYPE_FITS`)."""
    for _low, high, coefficients, _rms in tuning.BROWN_DWARF_SPECTRAL_TYPE_FITS:
        if number <= high:
            return _fit_temperature_k(coefficients, number)
    raise ValueError(f"spectral subclass number {number} is past the last fit")


def spectral_number(temperature_k):
    """The L, T, Y subclass number (0 = L0, 10 = T0, 20 = Y0) whose mean temperature is `temperature_k`
    (`tuning.BROWN_DWARF_SPECTRAL_TYPE_FITS` inverted), or None when the temperature is that of L0 or hotter. Colder
    than the last fit's end (22) is 22."""
    fits = tuning.BROWN_DWARF_SPECTRAL_TYPE_FITS
    if temperature_k >= _fit_temperature_k(fits[0][2], fits[0][0]):
        return None
    for low, high, coefficients, _rms in fits:
        if temperature_k >= _fit_temperature_k(coefficients, high) or high == fits[-1][1]:
            break
    c0, c1, c2 = coefficients
    if c2 == 0.0:
        number = (temperature_k - c0) / c1
    else:
        # The smaller root is the one inside the fit's range (the parabola turns over beyond it).
        number = (-c1 - math.sqrt(c1 * c1 - 4.0 * c2 * (c0 - temperature_k))) / (2.0 * c2)
    return min(max(number, low), high)


def spectral_type(temperature_k):
    """The spectral type of a brown dwarf of this effective temperature: `"M8"` (a late-M dwarf, above L0), `"L5"`,
    `"T2"` or `"Y0"`. The subclass is the whole number the mean relation reaches."""
    number = spectral_number(temperature_k)
    if number is None:
        letter, subclass = evolution.spectral_letter_and_subclass(temperature_k)
        return f"{letter}{subclass}"
    letter = "L" if number < 10.0 else "T" if number < 20.0 else "Y"
    return f"{letter}{int(number) % 10 if letter != 'Y' else int(number) - 20}"


def state(mass_sol, age_gy):
    """A brown dwarf of this mass and age: `luminosity_sol`, `temperature_k`, `radius_sol` (see `radius_sol`; the
    luminosity is the one the capped temperature and the radius make) and its `spectral_type`."""
    radius = radius_sol(mass_sol)
    temperature = effective_temperature_k(mass_sol, age_gy)
    luminosity = (radius ** 2) * (temperature / tuning.SUN_EFFECTIVE_TEMPERATURE_K) ** 4
    return {"luminosity_sol": luminosity, "temperature_k": temperature, "radius_sol": radius,
            "spectral_type": spectral_type(temperature)}
