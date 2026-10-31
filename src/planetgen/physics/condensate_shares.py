# planetgen/physics/condensate_shares.py

"""
What has condensed at a disk temperature, and the water share of the solids (GEN.950, part 1).

Lodders 2003 (ApJ 591, 1220) computes the equilibrium condensation of a gas of solar composition at 10^-4 bar. Table 10 gives the
temperatures at which the ices come out (water at 182 K, ammonia hydrate at 131 K, methane at 41 K) and Table 11 the mass that
each condensate makes, in percent of the nebula, for the Solar System abundances: rock 0.489 (silicates and oxides 0.303, metal and
troilite 0.186), water ice 0.571, ammonia ice 0.097 and methane ice 0.330. Rock is whole below the troilite line (704 K) and 22.75
percent of the oxygen is in it before the water ice condenses, so the rest of the oxygen is water.

So the water share of the solids is a step in temperature: none above 182 K, 0.571 / (0.489 + 0.571) = 0.539 between 182 K and 131 K,
0.494 with the ammonia, 0.384 with the methane. A host star of another composition changes the three ice masses: water in proportion
to the oxygen left after the rock has taken its fixed amount (`ROCK_OXYGEN_SHARE` of the Sun's, per atom of rock-forming metal),
ammonia to the nitrogen and methane to the carbon. The caller gives each as the star's ratio to the rock-forming metals over the
Sun's (`oxygen`, `carbon`, `nitrogen`; 1 is the Sun).

The temperature of a disk at a distance is the blackbody one, `disk_temperature_k`, the same law `formation.snow_line_au` inverts for
its fixed line.

Limits. The table is for a gas of solar C/O at 10^-4 bar with every carbon atom in methane and every nitrogen atom in ammonia; the
kinetic models of Lodders's Table 10 (graphite kept, methane and ammonia not formed) are not read here. The ice mass of a comet
(Marboeuf et al. 2014, ice over rock near 1) is half of this equilibrium (about 2); which one planetesimals follow is open. The snow line
(`formation.snow_line_au`), the ice boost of the surface density (`WATER_ICE_SOLID_BOOST`), the temperature below which a ring is cold
and the water share an icy body can reach (`bulk_draw`) read it (GEN.950); the star's oxygen over its rock-forming metals is not
passed on, so the ice masses are those of the Sun's composition.
"""

import math

from planetgen.physics import constants as physical_constants, unit_registry as units
from planetgen.util.checks import finite_domain

ICE_CONDENSATION_K = {"H2O": 182.0, "NH3": 131.0, "CH4": 41.0}
"""dict: The temperature (K) below which each ice condenses, complete equilibrium at 10^-4 bar, Solar System composition (Lodders 2003
Table 10: water ice 182, ammonia hydrate 131, pure methane ice 41; the methane clathrate at 78 K takes about 9 percent of the carbon
and is left out). [Stand-in CS1, tier T1, tweak L0.]"""

TABLE_11_MASS_PERCENT = {"rock": 0.489, "H2O": 0.571, "NH3": 0.097, "CH4": 0.330}
"""dict: The mass each condensate makes, in percent of the nebula, for the Solar System composition (Lodders 2003 Table 11; rock is
the silicates and oxides 0.303 plus metal and FeS 0.186). [Stand-in CS2, tier T1, tweak L0.]"""

ROCK_OXYGEN_SHARE = 0.2275
"""float: The share of the Sun's oxygen that is bound in rock before water ice condenses (Lodders 2003 Table 11 note: 22.75 percent);
the rest becomes water. [Stand-in CS3, tier T1, tweak L0.]"""

ROCK = "rock"
ICES = tuple(ICE_CONDENSATION_K)

_m_to_au = units.converter("m", "au")


def _check_ratios(oxygen, carbon, nitrogen):
    for name, value in (("oxygen", oxygen), ("carbon", carbon), ("nitrogen", nitrogen)):
        if not value >= 0.0:
            raise ValueError(f"the {name} ratio must be zero or more, not {value!r}")


@finite_domain()
def condensed_masses(temperature_k, oxygen=1.0, carbon=1.0, nitrogen=1.0):
    """
    The mass of each condensate that has formed at `temperature_k`, in the percent of the nebula of Lodders's Table 11.

    Args:
        temperature_k (float): The disk temperature, in K (rock is whole at any temperature below the troilite line at 704 K, and
            this is not for a hotter gas).
        oxygen (float): The star's O over rock-forming metal, over the Sun's. Below `ROCK_OXYGEN_SHARE` all of it is in rock and
            there is no water.
        carbon (float): The star's C over rock-forming metal, over the Sun's.
        nitrogen (float): The star's N over rock-forming metal, over the Sun's.

    Returns:
        dict: `"rock"` and each ice that has condensed (`"H2O"`, `"NH3"`, `"CH4"`) to its mass, in percent.
    """
    if not temperature_k > 0.0:
        raise ValueError("the temperature must be positive")
    _check_ratios(oxygen, carbon, nitrogen)
    table = TABLE_11_MASS_PERCENT
    scale = {
        "H2O": max(oxygen - ROCK_OXYGEN_SHARE, 0.0) / (1.0 - ROCK_OXYGEN_SHARE),
        "NH3": nitrogen,
        "CH4": carbon,
    }
    masses = {ROCK: table[ROCK]}
    for ice, line in ICE_CONDENSATION_K.items():
        if temperature_k < line:
            masses[ice] = table[ice] * scale[ice]
    return masses


def water_share_of_solids(temperature_k, oxygen=1.0, carbon=1.0, nitrogen=1.0):
    """Water ice over everything that has condensed (rock and the ices) at `temperature_k`; 0 above 182 K."""
    masses = condensed_masses(temperature_k, oxygen, carbon, nitrogen)
    return masses.get("H2O", 0.0) / sum(masses.values())


def ice_over_rock(temperature_k, oxygen=1.0, carbon=1.0, nitrogen=1.0):
    """The mass of all the ices that have condensed over the rock at `temperature_k` (1.17 for water alone, 2.04 for all three)."""
    masses = condensed_masses(temperature_k, oxygen, carbon, nitrogen)
    return (sum(masses.values()) - masses[ROCK]) / masses[ROCK]


def solid_boost(temperature_k, oxygen=1.0, carbon=1.0, nitrogen=1.0):
    """
    How many times the rock's surface density the solids make at `temperature_k`: `1 + ice_over_rock`. The jump of the minimum-mass
    nebula at the snow line (Hayashi 1981 puts it at 30 / 7 = 4.3 where Lodders's water alone gives 2.17).
    """
    return 1.0 + ice_over_rock(temperature_k, oxygen, carbon, nitrogen)


WATER_ICE_SOLID_BOOST = solid_boost(0.5 * (ICE_CONDENSATION_K["H2O"] + ICE_CONDENSATION_K["NH3"]))
"""float: How many times the rock's surface density the solids make just beyond the water line, between the water and ammonia lines (2.17 for
the Sun's composition; Hayashi 1981 has 30 / 7 = 4.3). The jump of the minimum-mass nebula at the snow line (`formation.mmsn_surface_density_gcm2`)."""


@finite_domain()
def disk_temperature_k(distance_au, luminosity_w):
    """The temperature a blackbody reaches in the light of a star of `luminosity_w` at `distance_au`: `(L / 16 pi sigma)^(1/4) / sqrt(d)`."""
    if not (distance_au > 0.0 and luminosity_w > 0.0):
        raise ValueError("the distance and the luminosity must be positive")
    distance_m = distance_au * physical_constants.AU_M
    return (luminosity_w / (16.0 * math.pi * physical_constants.STEFAN_BOLTZMANN_CONSTANT * distance_m ** 2)) ** 0.25


@finite_domain()
def ice_line_au(luminosity_w, ice="H2O"):
    """The distance at which a blackbody in the light of a star of `luminosity_w` is as cold as the condensation temperature of `ice`, in AU."""
    if ice not in ICE_CONDENSATION_K:
        raise ValueError(f"not an ice of the table: {ice!r}")
    distance_m = math.sqrt(luminosity_w / (16.0 * math.pi * physical_constants.STEFAN_BOLTZMANN_CONSTANT)) / ICE_CONDENSATION_K[ice] ** 2
    return _m_to_au(distance_m)
