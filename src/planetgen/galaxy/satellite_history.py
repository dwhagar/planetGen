# planetgen/galaxy/satellite_history.py

"""
The star formation history of a satellite galaxy, and the young and intermediate stars it holds (GEN.526;
docs/design/galaxy/intergalactic-and-halo-contents.md 5.5).

A satellite's row (`galaxy.components.Component`) holds its stellar mass, its half-mass radius and the mean age of its
stars. That is enough to say how much of it formed recently, with no share typed in:

1. **The history.** Stars formed from `tuning.SATELLITE_HISTORY_SPAN_GY` ago to now at a rate that rises or falls
   exponentially, `exp(-lambda t)` (`t` the time since the first stars), the one-parameter family whose mean stellar age
   can be set to the row's: a rising rate (negative `lambda`) is a satellite whose stars are younger than the middle of
   its history (the Large Cloud), a falling one an older one (the dwarf spheroidals), and a flat one has half the span as
   its mean age. The stellar mass sets the amplitude, so the present rate is the share of the mass formed in the
   history's last instant.
2. **The gas.** The Schmidt-Kennicutt law (`tuning.STAR_FORMATION_LAW_EXPONENT`, `STAR_FORMATION_LAW_COEFFICIENT`)
   inverts to the gas the present rate needs: the rate over the area of the half-mass radius gives a surface density, and
   that the gas per area. The gas disk is taken as long as the stellar one, as the Milky Way's disk is
   (`tuning.GAS_TO_STELLAR_SCALE_LENGTH_RATIO`).
3. **The threshold.** Below `tuning.SATELLITE_STAR_FORMATION_THRESHOLD_MSUN_PC2` a satellite forms no stars, so its young
   share is nothing; the stars its history put in the last 3 Gy are still there, and stay its intermediate share.

The young share is the history's last `tuning.STELLAR_POPULATION_AGE_RANGES_GY["satellite_young"]` window, the
intermediate share the window of `satellite_intermediate`, and the old share the rest.
"""

import dataclasses
import math

from planetgen import tuning

_YOUNG_WINDOW = tuning.STELLAR_POPULATION_AGE_RANGES_GY["satellite_young"]
_INTERMEDIATE_WINDOW = tuning.STELLAR_POPULATION_AGE_RANGES_GY["satellite_intermediate"]
_PC2_PER_KPC2 = 1.0e6
_LAMBDA_LIMIT_PER_GY = 12.0
"""float: The steepest rise or fall the family is searched over, per Gy: far steeper than any mean age asks for."""


@dataclasses.dataclass(frozen=True)
class SatelliteHistory:
    """
    A satellite's history as `history` works it out.

    Attributes:
        young_share, intermediate_share, old_share (float): The shares of its stars in the three populations
            (`satellite_young`, `satellite_intermediate`, `satellite`); they add up to 1.
        rate_per_gy (float): `lambda`, the exponent of the rate over the history, per Gy (negative: rising).
        present_rate_msun_yr (float): The present star formation rate the history gives, solar masses per year.
        gas_surface_density (float): The gas per area that rate needs by the Schmidt-Kennicutt law, solar masses per pc2.
        gas_mass_solar (float): That density over the area of the half-mass radius.
        forming (bool): Whether the gas is at or above the threshold, so the satellite forms stars now.
    """
    young_share: float
    intermediate_share: float
    old_share: float
    rate_per_gy: float
    present_rate_msun_yr: float
    gas_surface_density: float
    gas_mass_solar: float
    forming: bool


def _fraction_formed_before(rate_per_gy, span_gy, age_gy):
    """The share of the stars of an exponential history that are younger than `age_gy`: the integral of `exp(-lambda t)`
    over the last `age_gy` of the span, over the whole."""
    age = min(max(age_gy, 0.0), span_gy)
    if abs(rate_per_gy * span_gy) < 1.0e-9:
        return age / span_gy
    # The history is exp(-lambda t), t from 0 (the oldest stars) to the span; the last `age` is t from span - age.
    return (math.exp(-rate_per_gy * (span_gy - age)) - math.exp(-rate_per_gy * span_gy)) / (
        1.0 - math.exp(-rate_per_gy * span_gy))


def mean_age_gy(rate_per_gy, span_gy):
    """The mean age of the stars of an exponential history of this exponent, Gy."""
    x = rate_per_gy * span_gy
    if abs(x) < 1.0e-4:
        return 0.5 * span_gy + rate_per_gy * span_gy * span_gy / 12.0      # the series about a flat history
    mean_since_first_gy = 1.0 / rate_per_gy - span_gy / math.expm1(x)
    return span_gy - mean_since_first_gy


def rate_for_mean_age(mean_age, span_gy):
    """
    The exponent `lambda` (per Gy) whose history has stars of this mean age, found by bisection (the mean age falls as
    `lambda` falls, from the span at a steep fall to nothing at a steep rise). An age outside the range the search
    covers takes the steepest history.
    """
    low, high = -_LAMBDA_LIMIT_PER_GY, _LAMBDA_LIMIT_PER_GY
    if mean_age <= mean_age_gy(low, span_gy):
        return low
    if mean_age >= mean_age_gy(high, span_gy):
        return high
    for _ in range(80):
        middle = 0.5 * (low + high)
        if mean_age_gy(middle, span_gy) < mean_age:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def history(mass_solar, half_mass_radius_pc, mean_age):
    """
    The history of a satellite of this stellar mass, half-mass radius (pc) and mean stellar age (Gy): see the module
    docstring.

    Returns:
        SatelliteHistory
    """
    span = tuning.SATELLITE_HISTORY_SPAN_GY
    rate = rate_for_mean_age(mean_age, span)
    young = _fraction_formed_before(rate, span, _YOUNG_WINDOW[1])
    intermediate = _fraction_formed_before(rate, span, _INTERMEDIATE_WINDOW[1]) - young
    # The present rate per year: the share of the mass formed in the last year, which is the history's density at t = span.
    density_per_gy = (rate * math.exp(-rate * span) / (1.0 - math.exp(-rate * span))
                      if abs(rate * span) >= 1.0e-9 else 1.0 / span)
    present_rate = mass_solar * density_per_gy / 1.0e9
    area_pc2 = math.pi * float(half_mass_radius_pc) ** 2
    rate_per_area = present_rate / (area_pc2 / _PC2_PER_KPC2)
    surface = (rate_per_area / tuning.STAR_FORMATION_LAW_COEFFICIENT) ** (1.0 / tuning.STAR_FORMATION_LAW_EXPONENT)
    forming = surface >= tuning.SATELLITE_STAR_FORMATION_THRESHOLD_MSUN_PC2
    if not forming:
        young = 0.0
    return SatelliteHistory(
        young_share=young, intermediate_share=intermediate, old_share=1.0 - young - intermediate, rate_per_gy=rate,
        present_rate_msun_yr=present_rate, gas_surface_density=surface, gas_mass_solar=surface * area_pc2, forming=forming)
