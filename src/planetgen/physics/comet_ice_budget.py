# planetgen/physics/comet_ice_budget.py

"""
The ices a comet can hold, from where in the disc it formed (GEN.595, part 1; docs/design/composition/comet-composition-tuning.md section 5).

`comet_population` draws a comet's CO and CO2 from the Solar System's measured comets (a log-normal about the population's median) and, for an interstellar comet,
from the wider ranges its two objects force (CO/H2O to 2, CO2/H2O to 10). Those are the Solar System's shape with a per-system offset knob (GEN.520 D6). This module
derives the shape from the disc instead, from the snow lines and the carbon and oxygen budget of Oberg, Murray-Clay and Bergin 2011 (ApJL 743, L16, arXiv:1110.5567,
Table 1 and section 2, READ):

1. **Where the ices freeze.** The disc is a power law `T = T0 (r / 1 AU)^-q` with the average of a large sample of discs, T0 200 K and q 0.62 (Andrews and Williams
   2007, as Oberg et al. use it), and the star's luminosity scales T0 as `L^0.25`. The evaporation temperatures at disc densities are 135 K (H2O), 47 K (CO2) and 20 K
   (CO), which put the lines of a Sun-like star at 1.9, 10 and 41 AU (`ice_line_au`).
2. **What is there to freeze.** Table 1 gives the O and C in each carrier per 10^4 hydrogen: H2O 0.9, CO2 0.3 and CO 1.5 (Pontoppidan 2006, from the ice and gas of one
   disc, as Oberg et al. take them). So a body formed beyond the CO2 line has CO2/H2O 0.33 and beyond the CO line CO/H2O 1.67, in the star's own carbon and oxygen: the CO and
   CO2 scale with the star's carbon and the H2O with its oxygen after the rock has taken its fixed share (`condensate_shares.ROCK_OXYGEN_SHARE`).
3. **What a warm body keeps.** Above its line a molecule is not frozen but can be trapped in water ice: Oberg et al. cite Fayolle et al. 2011 that it is hard to trap more than
   5 to 10 percent of the CO (`TRAPPED_FRACTION`); the same fraction is taken for CO2, which sublimes warmer.
4. **Where bodies formed.** A solids surface density that goes as `1 / r` puts equal masses of planetesimals in equal intervals of `ln r`, and `zone_shares` gives the share
   of a formation range in each of the four zones on that rule (the rule is a choice: [Stand-in CI5, tier T6, tweak L1]).

Tests against measurement (never inputs): 2I/Borisov had CO/H2O of 130 to 155 percent (Bodewits et al. 2020, arXiv:2004.08972, READ), which is the CO zone's 167 percent for a
Sun-like carbon and oxygen to within 30 percent; the Solar System's median comet (CO2/H2O 10.9, CO/H2O 2.5 percent, `comet_population`) falls between the CO2 zone and the
warm zone, 0.5 dex under the CO2 zone's 33 percent in CO2 and 0.7 dex under its trapped CO (the model's absolute level for comets of the Solar System is the stand-in that stays); and
3I/ATLAS, CO2/H2O 7.6 (Cordiner et al. 2025, arXiv:2508.18209, READ), is 23 times the most a disc of this budget gives, an excess no zone explains here (the interstellar bracket of
`comet_population` stays a registered stand-in).

Limits. The temperature profile is the average of a sample of discs, one for all stars (Oberg et al. 2011 also use it); the lines do not move as the disc cools or the
snow line walks out; transport (drift and diffusion of icy grains) and processing by cosmic rays and radiation are not in it; the trapped fraction is a bound, not a measurement
of comets; the H2O line here (135 K) is the disc-density value and `condensate_shares` has Lodders's 182 K at 10^-4 bar, which the wiring leaves be (comets form beyond the star's own snow line, `formation_ratios`). Wired (train 7): a comet of a star's system forms at a distance drawn between the star's snow line and the edge of its disc (`formation_ratios`), and
`comet_population.draw_volatiles` takes its median CO and CO2 from the zone (the scatter stays the population's). The interstellar bracket stays a registered stand-in.
"""

import math

from planetgen.physics import condensate_shares
from planetgen.util.checks import finite_domain

ICE_LINES_K = {"H2O": 135.0, "CO2": 47.0, "CO": 20.0}
"""dict: The evaporation temperature of each ice at midplane densities (Oberg, Murray-Clay and Bergin 2011 Table 1, the adopted values; the ranges are 120-150, 42-52
and 18-22 K). [Stand-in CI1, tier T2, tweak L0.]"""

DISK_T0_K = 200.0
DISK_INDEX = 0.62
"""float: The temperature at 1 AU and the power-law index of the average protoplanetary disc round a Sun-like star (Andrews and Williams 2007, as Oberg et al. use them). [Stand-in CI1, tier T2, tweak L0.]"""

LUMINOSITY_EXPONENT = 0.25
"""float: The disc temperature goes as the star's luminosity to this power (radiative equilibrium: T goes as L^(1/4)). [Stand-in CI2, tier T5, tweak L0.]"""

SOLAR_LUMINOSITY_W = 3.828e26
"""float: The luminosity the disc profile is for (IAU nominal solar luminosity)."""

NEBULA_PER_1E4_H = {"H2O": 0.9, "CO2": 0.3, "CO": 1.5}
"""dict: The molecules of each carrier per 10^4 hydrogen atoms in the midplane of a typical disc (Oberg et al. 2011 Table 1: nO of H2O 0.9, of CO2 0.6 (two atoms: 0.3 molecules), of CO 1.5; from
Pontoppidan 2006). [Stand-in CI3, tier T2, tweak L0.]"""

TRAPPED_FRACTION = 0.075
"""float: The share of the CO (and CO2) above its snow line that is trapped in water ice and so in the body (Fayolle et al. 2011: at most 5 to 10 percent of the CO, as Oberg et al. cite it;
2ND-HAND). The CO2 value is the same by assumption. [Stand-in CI4, tier T2, tweak L1.]"""

ZONES = ("dry", "water", "carbon dioxide", "carbon monoxide")


@finite_domain()
def disk_temperature_k(distance_au, luminosity_w=SOLAR_LUMINOSITY_W):
    """The midplane temperature at a distance round a star of this luminosity: `T0 (L / Lsun)^0.25 r^-q`."""
    if not (distance_au > 0.0 and luminosity_w > 0.0):
        raise ValueError("the distance and the luminosity must be positive")
    return DISK_T0_K * (luminosity_w / SOLAR_LUMINOSITY_W) ** LUMINOSITY_EXPONENT * distance_au ** -DISK_INDEX


@finite_domain()
def ice_line_au(ice, luminosity_w=SOLAR_LUMINOSITY_W):
    """The distance, AU, inside which an ice is not frozen: where the disc is at its evaporation temperature (1.9, 10 and 41 AU for H2O, CO2 and CO round the Sun)."""
    if ice not in ICE_LINES_K:
        raise ValueError(f"no ice line for {ice!r}")
    if not luminosity_w > 0.0:
        raise ValueError("the luminosity must be positive")
    return (DISK_T0_K * (luminosity_w / SOLAR_LUMINOSITY_W) ** LUMINOSITY_EXPONENT / ICE_LINES_K[ice]) ** (1.0 / DISK_INDEX)


@finite_domain()
def nebular_ratios(carbon=1.0, oxygen=1.0):
    """
    The CO2 and the CO over the H2O when every one of them is frozen (the budget of a body formed beyond the CO line), as `(co2, co)`.

    Args:
        carbon (float): The star's C over rock-forming metal, over the Sun's.
        oxygen (float): The star's O over rock-forming metal, over the Sun's; below `ROCK_OXYGEN_SHARE` all of it is in rock and the water is zero.
    """
    if not (carbon >= 0.0 and oxygen >= 0.0):
        raise ValueError("the carbon and oxygen ratios must be zero or more")
    rock = condensate_shares.ROCK_OXYGEN_SHARE
    water = NEBULA_PER_1E4_H["H2O"] * max(oxygen - rock, 0.0) / (1.0 - rock)
    if water <= 0.0:
        raise ValueError("no oxygen is left for water")
    return NEBULA_PER_1E4_H["CO2"] * carbon / water, NEBULA_PER_1E4_H["CO"] * carbon / water


@finite_domain()
def zone_of(distance_au, luminosity_w=SOLAR_LUMINOSITY_W):
    """The zone a body formed in: `dry` inside the water line, `water` to the CO2 line, `carbon dioxide` to the CO line, `carbon monoxide` beyond it."""
    if not distance_au > 0.0:
        raise ValueError("the distance must be positive")
    if distance_au < ice_line_au("H2O", luminosity_w):
        return "dry"
    if distance_au < ice_line_au("CO2", luminosity_w):
        return "water"
    return "carbon dioxide" if distance_au < ice_line_au("CO", luminosity_w) else "carbon monoxide"


@finite_domain()
def ice_ratios(distance_au, luminosity_w=SOLAR_LUMINOSITY_W, carbon=1.0, oxygen=1.0, trapped=TRAPPED_FRACTION):
    """
    The CO2/H2O and CO/H2O of the ice of a body formed at a distance, as `(co2, co)`, or None for a body inside the water line (no water ice to hold them).

    Beyond the CO line both are frozen whole; between the CO2 and CO lines the CO2 is whole and `trapped` of the CO is held in the water ice; between the water and CO2 lines
    `trapped` of each is.
    """
    if not 0.0 <= trapped <= 1.0:
        raise ValueError("the trapped fraction must be in [0, 1]")
    zone = zone_of(distance_au, luminosity_w)
    if zone == "dry":
        return None
    co2, co = nebular_ratios(carbon, oxygen)
    if zone == "carbon monoxide":
        return co2, co
    if zone == "carbon dioxide":
        return co2, trapped * co
    return trapped * co2, trapped * co


@finite_domain()
def formation_ratios(luminosity_w, snow_line_au, outer_au, carbon, oxygen, u):
    """
    The CO2/H2O and CO/H2O of the ice of one comet, as `(co2, co)`: formed at a distance drawn log-uniformly (the equal masses in equal intervals of ln r of
    `zone_shares`, by the uniform `u` in [0, 1)) between the star's snow line, or the water line of the disc model past it, and `outer_au`, the edge of the disc.

    Args:
        luminosity_w (float): The star's luminosity.
        snow_line_au (float): The star's own snow line (`Star.snow_line_au`): comets form beyond it.
        outer_au (float): The outer edge of the disc.
        carbon, oxygen (float): As `nebular_ratios`.
        u (float): The uniform draw that places the comet.
    """
    if not 0.0 <= u < 1.0:
        raise ValueError("u must be in [0, 1)")
    inner = max(snow_line_au, ice_line_au("H2O", luminosity_w) * (1.0 + 1e-9))
    if not outer_au > inner:
        return ice_ratios(inner, luminosity_w, carbon, oxygen)
    return ice_ratios(inner * (outer_au / inner) ** u, luminosity_w, carbon, oxygen)


@finite_domain()
def zone_shares(inner_au, outer_au, luminosity_w=SOLAR_LUMINOSITY_W):
    """The share of the planetesimals formed between two distances in each zone, as a dict over `ZONES` adding to 1: equal masses form in equal intervals of ln r (a solids surface density as 1 / r)."""
    if not 0.0 < inner_au < outer_au:
        raise ValueError("the formation range needs 0 < inner < outer")
    edges = [0.0, *(ice_line_au(ice, luminosity_w) for ice in ("H2O", "CO2", "CO")), math.inf]
    total = math.log(outer_au / inner_au)
    shares = {}
    for zone, low, high in zip(ZONES, edges[:-1], edges[1:]):
        lo, hi = max(inner_au, low), min(outer_au, high)
        shares[zone] = math.log(hi / lo) / total if hi > lo else 0.0
    return shares
