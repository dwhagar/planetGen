# planetgen/physics/brown_dwarf.py

"""
A brown dwarf as a star (GEN.252, GEN.549): its mass, how bright and hot it is at an age, and the L, T or Y spectral type
that temperature makes.

A brown dwarf is too light to hold hydrogen fusion, so it only cools. Its temperature and radius at a mass and an age
come from the evolutionary models of Baraffe, Homeier, Allard and Chabrier 2015 (A&A 577:A42, `tuning.BROWN_DWARF_TRACKS`):
one cooling track per mass, interpolated in the logarithms of mass and age. The models stop where an object has cooled
below their range (about 1,500 K); a track that ends before 10 Gy continues at the old-age cooling of Burrows et al. 2001
(Rev. Mod. Phys. 73:719, eq. 2, `Teff` falling as `t^-0.32`) with its last radius, and a track that reaches 10 Gy holds its
last values: the objects at the hydrogen-burning limit settle on a plateau, which is why L dwarfs hotter than 1,500 K
exist in the numbers the 20 pc census counts. The temperature is capped at `tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K` (a
late-M dwarf). The luminosity is `4 pi R^2 sigma Teff^4` of the two.

The spectral type is the mean relation of Kirkpatrick et al. 2021 (ApJS 253:7, Table 13, from 525 L, T and Y dwarfs
within 20 pc), inverted: L0 is 2,238 K, T0 1,255 K, Y0 461 K. Their relation has a scatter of 134 K among L dwarfs and 79 K
among T and Y dwarfs; that is measurement scatter and is not drawn, so a brown dwarf of a given temperature always has
the one type. Above L0 it is a late-M dwarf, named by the stellar scale (`stellar_evolution.spectral_letter_and_subclass`).

Masses follow the substellar mass function of the same census, `dN/dM ~ M^-0.6` (alpha 0.6 +- 0.1).

Left out: the L/T transition's cloud-clearing plateau (about 18 K per subtype from L8.75 to T4.75, which the piecewise fit
already flattens), metallicity and dust opacity, and a gravity-lowered type for a young brown dwarf. Objects between 0.072
and 0.0764 solar masses are inside the 80 Jupiter mass limit of `ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER` but above the models'
own hydrogen-burning limit, so the models make them stable late-M and L dwarfs; they are kept, as the census's L dwarfs.
"""

import bisect
import math

from planetgen import tuning
from planetgen.physics import constants, stellar_evolution as evolution
from planetgen.util import draw


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


def _track_values(track, age_gy):
    """`(temperature, radius)` of one `BROWN_DWARF_TRACKS` row at `age_gy`: log-linear in age between its ages, the first
    values before them, and past the last age the radius held and the temperature held (a row that reaches 10 Gy) or
    continued at `BROWN_DWARF_TEMPERATURE_AGE_EXPONENT` (a row that ended sooner)."""
    ages = [row[0] for row in track]
    if age_gy <= ages[0]:
        return float(track[0][1]), float(track[0][2])
    if age_gy >= ages[-1]:
        last_age, last_temperature, last_radius = track[-1]
        if last_age >= tuning.BROWN_DWARF_TRACK_END_GY:
            return float(last_temperature), float(last_radius)
        return last_temperature * (age_gy / last_age) ** tuning.BROWN_DWARF_TEMPERATURE_AGE_EXPONENT, float(last_radius)
    upper = bisect.bisect_right(ages, age_gy)
    (age_low, temperature_low, radius_low), (age_high, temperature_high, radius_high) = track[upper - 1], track[upper]
    weight = math.log(age_gy / age_low) / math.log(age_high / age_low)
    return (math.exp((1.0 - weight) * math.log(temperature_low) + weight * math.log(temperature_high)),
            math.exp((1.0 - weight) * math.log(radius_low) + weight * math.log(radius_high)))


def _track(mass_sol, age_gy):
    """`(temperature, radius)` of the cooling track at `mass_sol` and `age_gy`: the two rows around the mass, each at the
    age, joined log-linearly in mass."""
    tracks = tuning.BROWN_DWARF_TRACKS
    masses = [row[0] for row in tracks]
    if not masses[0] <= mass_sol <= masses[-1]:
        raise ValueError(f"no cooling track for {mass_sol} solar masses (the models run {masses[0]} to {masses[-1]})")
    upper = min(max(bisect.bisect_left(masses, mass_sol), 1), len(masses) - 1)
    (mass_low, track_low), (mass_high, track_high) = tracks[upper - 1], tracks[upper]
    weight = math.log(mass_sol / mass_low) / math.log(mass_high / mass_low)
    temperature_low, radius_low = _track_values(track_low, age_gy)
    temperature_high, radius_high = _track_values(track_high, age_gy)
    return (math.exp((1.0 - weight) * math.log(temperature_low) + weight * math.log(temperature_high)),
            math.exp((1.0 - weight) * math.log(radius_low) + weight * math.log(radius_high)))


def effective_temperature_k(mass_sol, age_gy):
    """Effective temperature, kelvin: the cooling track's, no hotter than `tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K`."""
    return min(_track(mass_sol, age_gy)[0], tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K)


def radius_sol(mass_sol, age_gy):
    """Radius in solar radii at this mass and age: the cooling track's (a young brown dwarf is larger, an old one nearly
    Jupiter's size)."""
    return _track(mass_sol, age_gy)[1]


def luminosity_sol(mass_sol, age_gy):
    """Luminosity in solar luminosities: `R^2 (Teff / Teff_sun)^4` of the track's radius and the capped temperature."""
    temperature, radius = _track(mass_sol, age_gy)
    temperature = min(temperature, tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K)
    return radius ** 2 * (temperature / tuning.SUN_EFFECTIVE_TEMPERATURE_K) ** 4


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
    """A brown dwarf of this mass and age: `luminosity_sol`, `temperature_k` (capped), `radius_sol` and its
    `spectral_type`."""
    temperature, radius = _track(mass_sol, age_gy)
    temperature = min(temperature, tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K)
    return {"luminosity_sol": radius ** 2 * (temperature / tuning.SUN_EFFECTIVE_TEMPERATURE_K) ** 4,
            "temperature_k": temperature, "radius_sol": radius, "spectral_type": spectral_type(temperature)}
