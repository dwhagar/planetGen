# planetgen/physics/saturation_time.py

"""
How long a star stays saturated in X-rays, by mass and by initial rotation (GEN.493; docs/design/composition/host-activity-tuning.md section 4).

Row p6 of the atmosphere recipes draws a star's initial-rotation percentile and takes the saturation time from three solar-mass points of
Tu et al. 2015 (5.7, 23 and 226 Myr at the 10th, 50th and 90th percentile); the mass dependence was missing (every mass used the solar one, GEN.523).
[R2183] model the rotation and X-ray evolution of 0.1 to 1.25 solar mass stars with
the spin-down of a model and give, in their Fig. 10, the age at which a star falls below the saturation threshold against mass for the slow,
medium and fast rotators, the observed 5th, 50th and 95th percentiles of the rotation distribution at 150 Myr (READ). `ANCHORS_MYR` is that figure
read off: the pixel positions of the three curves on the published page, converted with the axes (a logarithmic axis up to 1000 Myr and a linear
one above), taken at thirteen masses. At one solar mass the curves give 5, 23 and 297 Myr, against Tu et al.'s 5.7, 23 and 226 Myr at the 10th,
50th and 90th percentiles: the same numbers where the two overlap. Lower mass stars stay saturated for gigayears (0.3 solar masses: 1.0, 1.5 and
1.9 Gyr) because they feel a smaller spin-down torque and must spin down further to cross a lower saturation threshold.

`saturation_time_myr(mass, percentile)` interpolates log-log in mass between the anchors and linearly in the logarithm of the time in the percentile
between 5, 50 and 95 (the same slope outside, to the 2nd and 98th percentile the draw reaches). The percentile of the observed distribution is uniform
by construction, so the row's draw is a uniform percentile and carries no weights of its own.

Limits. The anchors are figure read-offs (tier T3, tweak L0 as read; a few percent from the pixel scale, more on the steep 1.1 to 1.2 solar mass edge of
the fast track); the model is for solar metallicity and the percentiles of the 150 Myr distribution observed in clusters; fully convective stars at
0.1 solar masses carry the least certain spin-down. Wired (GEN.493): `activity.generate_activity` draws each star's percentile and takes its saturation age here, and row p6 of the atmosphere
recipes reads it through the recipe function `saturation_time_myr`.
"""

import math

from planetgen.util.checks import finite_domain

PERCENTILES = (5.0, 50.0, 95.0)
"""tuple: The percentiles of the 150 Myr rotation distribution the slow, medium and fast curves of [R2183] are."""

ANCHORS_MYR = {
    0.10: (3150.0, 4460.0, 4700.0),
    0.15: (1960.0, 2730.0, 3320.0),
    0.20: (1480.0, 2050.0, 2510.0),
    0.30: (1030.0, 1490.0, 1880.0),
    0.40: (412.0, 824.0, 1230.0),
    0.50: (120.0, 491.0, 869.0),
    0.60: (46.6, 300.0, 679.0),
    0.70: (29.2, 154.0, 570.0),
    0.80: (15.7, 38.4, 486.0),
    0.90: (7.2, 29.7, 387.0),
    1.00: (5.2, 22.8, 297.0),
    1.10: (4.7, 16.9, 130.0),
    1.20: (2.5, 13.1, 17.4),
}
"""dict: The age, Myr, at which a star of this mass (solar masses) falls below the X-ray saturation threshold, for the slow (5th percentile), medium
(50th) and fast (95th percentile) rotators, read off Fig. 10 of [R2183] (the fast value at 0.1 solar masses is read by eye where the
curve runs under the medium one). [Stand-in RS1, tier T3, tweak L0.]"""

_MASSES = sorted(ANCHORS_MYR)


def _log_time_at_percentile(row, percentile):
    """Log of the time for one row of three, linear in the percentile between the points and on the nearer segment's slope outside them."""
    logs = [math.log(value) for value in row]
    if percentile <= PERCENTILES[1]:
        low, high = 0, 1
    else:
        low, high = 1, 2
    slope = (logs[high] - logs[low]) / (PERCENTILES[high] - PERCENTILES[low])
    return logs[low] + slope * (percentile - PERCENTILES[low])


@finite_domain()
def saturation_time_myr(mass_solar, percentile=50.0):
    """
    The age, Myr, at which a star of this mass and rotation percentile (2 to 98 of the initial distribution) falls out of X-ray saturation.

    Masses outside 0.1 to 1.2 solar masses take the nearest anchor (the model does not reach them).
    """
    if not 0.0 <= percentile <= 100.0:
        raise ValueError(f"a percentile is 0 to 100, not {percentile!r}")
    mass = min(max(mass_solar, _MASSES[0]), _MASSES[-1])
    upper = next((m for m in _MASSES if m >= mass), _MASSES[-1])
    lower = max((m for m in _MASSES if m <= mass), default=_MASSES[0])
    at_low = _log_time_at_percentile(ANCHORS_MYR[lower], percentile)
    if upper == lower:
        return math.exp(at_low)
    at_high = _log_time_at_percentile(ANCHORS_MYR[upper], percentile)
    weight = math.log(mass / lower) / math.log(upper / lower)
    return math.exp(at_low + weight * (at_high - at_low))
