# planetgen/physics/brown_dwarf.py

"""
A brown dwarf as a star (GEN.252, GEN.549): its mass, how bright and hot it is at an age, and the L, T or Y spectral type
that temperature makes.

A brown dwarf is too light to hold hydrogen fusion, so it only cools. Its temperature and radius at a mass and an age
come from the evolutionary models of [R0254] (A&A 577:A42, the 30 ages of `datasets/baraffe_isochrones.py`, GEN.948):
one cooling track per mass, interpolated in the logarithms of mass and age. The models stop where an object has cooled
below their range (about 1,500 K); a track that ends before 10 Gy continues at the old-age cooling of [R0606]
(Rev. Mod. Phys. 73:719, eq. 2, `Teff` falling as `t^-0.32`) with its last radius, and a track that reaches 10 Gy holds its
last values: the objects at the hydrogen-burning limit settle on a plateau, which is why L dwarfs hotter than 1,500 K
exist in the numbers the 20 pc census counts. The temperature is capped at the temperature of the
hydrogen-burning limit at the age (`substellar_ceiling.cap_k`, BHAC15: 2,950 K at the peak, 2,050 K at 1 Gy), the hottest a brown dwarf of that age can be. The luminosity is `4 pi R^2 sigma Teff^4` of the two.

The spectral type is the map of [R2374] (ApJS 253:7, Table 13, from 525 L, T and Y dwarfs within 20 pc,
`brown_dwarf_types`): L0 is 2,238 K, T0 1,255 K, Y0 461 K. The map has a scatter of 134 K among L dwarfs and 79 K among T
and Y dwarfs, and a dwarf's type is drawn with it (GEN.760) on a stream of its temperature, so a given dwarf always has the
one type but dwarfs of one temperature do not all share it. Above L0 it is a late-M dwarf, named by the stellar scale
(`stellar_evolution.spectral_letter_and_subclass`).

Masses follow the substellar mass function of the same census, `dN/dM ~ M^-0.6` (alpha 0.6 +- 0.1).

Left out: the radius of an old brown dwarf of the lightest masses (the table loses them while they are still young and large, so
the radius is held at that young value, 1.7 Jupiter radii at 13 Jupiter masses; the temperature continues), the L/T transition's cloud-clearing plateau (about 18 K per subtype from L8.75 to T4.75, which the piecewise fit
already flattens), metallicity and dust opacity, and a gravity-lowered type for a young brown dwarf. Objects between 0.072
and the top of the range (0.08 solar masses, where the stars' mass function begins: `mass_boundaries`, GEN.936) are above the
models' own hydrogen-burning limit, so the models make them stable late-M and L dwarfs; they are kept, as the census's L
dwarfs.
"""

from planetgen import tuning
from planetgen.datasets import baraffe_isochrones
from planetgen.physics import brown_dwarf_types, stellar_evolution as evolution, substellar_ceiling, unit_registry as units
from planetgen.util import draw

_kg_to_msun = units.converter("kg", "Msun")
_mjup_to_kg = units.converter("Mjup", "kg")


def mass_range_solar():
    """`(lowest, highest)` brown dwarf mass in solar masses: the deuterium limit to the hydrogen limit."""
    low, high = (_kg_to_msun(_mjup_to_kg(jupiter))
                 for jupiter in tuning.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER)
    return low, high


def sample_mass_solar(rng=draw):
    """A brown dwarf's mass in solar masses, from `dN/dM ~ M^-alpha` (`tuning.BROWN_DWARF_MASS_FUNCTION_SLOPE`) between
    `mass_range_solar`."""
    low, high = mass_range_solar()
    power = 1.0 - tuning.BROWN_DWARF_MASS_FUNCTION_SLOPE
    return (low ** power + rng.random() * (high ** power - low ** power)) ** (1.0 / power)


def _track(mass_sol, age_gy):
    """`(temperature, radius)` of the cooling track at `mass_sol` and `age_gy`: the Baraffe table, continued past the age it loses the
    mass (`baraffe_isochrones.cooling_track`)."""
    return baraffe_isochrones.cooling_track(mass_sol, age_gy)


def effective_temperature_k(mass_sol, age_gy):
    """Effective temperature, kelvin: the cooling track's, no hotter than `substellar_ceiling.cap_k` of the mass and age."""
    return min(_track(mass_sol, age_gy)[0], substellar_ceiling.cap_k(mass_sol, age_gy))


def radius_sol(mass_sol, age_gy):
    """Radius in solar radii at this mass and age: the cooling track's (a young brown dwarf is larger, an old one nearly
    Jupiter's size)."""
    return _track(mass_sol, age_gy)[1]


def luminosity_sol(mass_sol, age_gy):
    """Luminosity in solar luminosities: `R^2 (Teff / Teff_sun)^4` of the track's radius and the capped temperature."""
    temperature, radius = _track(mass_sol, age_gy)
    temperature = min(temperature, substellar_ceiling.cap_k(mass_sol, age_gy))
    return radius ** 2 * (temperature / tuning.SUN_EFFECTIVE_TEMPERATURE_K) ** 4


def spectral_type(temperature_k):
    """
    The spectral type of a brown dwarf of this effective temperature: `"M8"` (a late-M dwarf), `"L5"`, `"T2"` or `"Y0"`
    (GEN.760). The map of [R2374] (`brown_dwarf_types`) gives the type's mean temperature with a scatter
    of 134 K among L dwarfs and 79 K among T and Y dwarfs, and the type is drawn with it (`brown_dwarf_types.drawn_number`)
    on a stream of the temperature alone, so the same dwarf has the same type wherever it is read from its stored
    temperature. A temperature that reads as L0 or hotter after the offset is a late-M dwarf, named by the stellar scale.
    """
    number, read_at = brown_dwarf_types.drawn_number(
        temperature_k, draw.Stream(f"brown-dwarf-type:{float(temperature_k)!r}"))
    if number is None:
        letter, subclass = evolution.spectral_letter_and_subclass(read_at)
        return f"{letter}{subclass}"
    letter = "L" if number < 10.0 else "T" if number < 20.0 else "Y"
    return f"{letter}{int(number) % 10 if letter != 'Y' else int(number) - 20}"


def state(mass_sol, age_gy):
    """A brown dwarf of this mass and age: `luminosity_sol`, `temperature_k` (capped), `radius_sol` and its
    `spectral_type`."""
    temperature, radius = _track(mass_sol, age_gy)
    temperature = min(temperature, substellar_ceiling.cap_k(mass_sol, age_gy))
    return {"luminosity_sol": radius ** 2 * (temperature / tuning.SUN_EFFECTIVE_TEMPERATURE_K) ** 4,
            "temperature_k": temperature, "radius_sol": radius, "spectral_type": spectral_type(temperature)}
