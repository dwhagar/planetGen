# planetgen/physics/chemically_peculiar.py

"""
Ap and Am stars: the chemically peculiar branch of the A-star draw (GEN.1027, after GEN.944).

A main-sequence star in the temperature range of the late-B to early-F stars is, with a share of its own, chemically peculiar:

- **Ap** (and Bp): a strong, ordered magnetic field (a dipole of 300 gauss to tens of kilogauss at the pole, Auriere et al. 2007) with abundance
  patches, which slows the star's rotation (periods from half a day to decades, most of days).
- **Am**: a metallic-line star with no detected magnetic field, a slow rotator (equatorial speed under about 100 km/s) that is usually in a close
  binary.

`kind_of` names the kind from the star's effective temperature and one uniform number, so a star keeps its kind from that number
(`stellar_abundances.branch_of` does the same for the alpha branch); `rotation_period_hours` and `polar_field_gauss` give the slow spin and the field
from further uniform numbers.

`of_star` draws the three numbers from a stream of the star's own (keyed by its mass, temperature and age, so the same star is the same kind and has the
same field however it is reached: a drawn star, a stored bright star's row, the sheet of a stored star). `spin.star_rotation_period_hours` takes the
period from it for a main-sequence dwarf, so the star and its row agree, and the sheet shows the kind and the field; the field is worked out from the
stored mass, temperature and age and is not a column.
"""

import math
from collections import namedtuple

from planetgen.util import draw

AP = "Ap"
AM = "Am"
KINDS = (AP, AM)

AP_TEMPERATURE_RANGE_K = (7500.0, 13000.0)
"""tuple: The temperatures, K, in which a main-sequence star may be Ap or Bp (late B to early F). [Stand-in PEC1, tier T3, tweak L0: the magnetic chemically
peculiar stars are found from B8 to F0 (Kochukhov and Bagnulo 2006; Auriere et al. 2007).]"""
AM_TEMPERATURE_RANGE_K = (7000.0, 9500.0)
"""tuple: The temperatures, K, in which a main-sequence star may be Am. [Stand-in PEC2, tier T3, tweak L0: Am stars are the A1 to F2 dwarfs, about 7,000 to
9,500 K (Abt 2009; Smith 1996).]"""
AP_SHARE = 0.10
"""float: The share of the stars in `AP_TEMPERATURE_RANGE_K` that are Ap or Bp. [Stand-in PEC3, tier T3, tweak L1: about 10 percent of A stars (Wolff 1968), 5 to
10 percent in the magnetic surveys (Auriere et al. 2007: 8 percent).]"""
AM_SHARE = 0.33
"""float: The share of the stars in `AM_TEMPERATURE_RANGE_K` that are Am. [Stand-in PEC4, tier T3, tweak L1: about a third of A stars between A4 and F1 (Smith
1996; Abt 2009 has 30 to 40 percent).]"""
AP_PERIOD_MEDIAN_DAYS = 4.5
"""float: The median rotation period of an Ap star, days. [Stand-in PEC5, tier T6, tweak L1: the periods run from half a day to decades with the bulk at days
(Mathys 2017); read second-hand.]"""
AP_PERIOD_SIGMA_LN = 1.3
"""float: The ln-width of the Ap rotation periods. [Stand-in PEC6, tier T6, tweak L1: spans the half-day to 30-day bulk within about two widths.]"""
AP_PERIOD_RANGE_DAYS = (0.5, 30.0)
"""tuple: The Ap periods kept, days. [Stand-in PEC7, tier T6, tweak L0: the long-period tail of decades is left out of the draw.]"""
AM_SPEED_RANGE_KMS = (10.0, 90.0)
"""tuple: The equatorial speeds, km/s, an Am star is drawn between, uniform. [Stand-in PEC8, tier T6, tweak L1: Am stars rotate slowly, v sin i under about 100
km/s with a mean near 50 (Abt and Morrell 1995).]"""
AP_POLAR_FIELD_GAUSS = (300.0, 3.0e4)
"""tuple: The polar field of an Ap star, gauss, drawn log-uniform between. [Stand-in PEC9, tier T3, tweak L0: dipole fields of 300 G to 34 kG (Auriere et al. 2007; Babcock's star
at 34 kG); the lower edge is the survey's detection floor.]"""

SOLAR_RADIUS_KM = 695700.0

DWARF_CLASS = "V"
"""str: The Yerkes class of the stars that may be peculiar: the main-sequence dwarfs."""

Peculiarity = namedtuple("Peculiarity", "kind period_hours field_gauss")
"""namedtuple: A peculiar star's kind (`AP` or `AM`), rotation period (hours) and polar field (gauss, None for an Am star)."""


def kind_of(temperature_k, u):
    """`AP`, `AM` or None for a main-sequence star of this effective temperature and uniform number `u` in [0, 1): the lowest `AP_SHARE` of the uniform numbers
    of a star in the Ap range are Ap, the next `AM_SHARE` of one in the Am range are Am."""
    ap_low, ap_high = AP_TEMPERATURE_RANGE_K
    am_low, am_high = AM_TEMPERATURE_RANGE_K
    if ap_low <= temperature_k <= ap_high and u < AP_SHARE:
        return AP
    if am_low <= temperature_k <= am_high and AP_SHARE <= u < AP_SHARE + AM_SHARE:
        return AM
    return None


def rotation_period_hours(kind, radius_rsun, u):
    """The rotation period, hours, of a peculiar star of `kind` and radius `radius_rsun` from the uniform number `u` (an Ap star's from its log-normal period
    law held to `AP_PERIOD_RANGE_DAYS`, an Am star's from its equatorial speed)."""
    if kind == AP:
        low, high = (math.log(days) for days in AP_PERIOD_RANGE_DAYS)
        centre, width = math.log(AP_PERIOD_MEDIAN_DAYS), AP_PERIOD_SIGMA_LN
        from statistics import NormalDist
        cut_low, cut_high = NormalDist().cdf((low - centre) / width), NormalDist().cdf((high - centre) / width)
        z = NormalDist().inv_cdf(cut_low + min(max(u, 1.0e-12), 1.0 - 1.0e-12) * (cut_high - cut_low))
        return 24.0 * math.exp(centre + width * z)
    if kind == AM:
        low, high = AM_SPEED_RANGE_KMS
        speed = low + u * (high - low)
        return 2.0 * math.pi * radius_rsun * SOLAR_RADIUS_KM / speed / 3600.0
    raise ValueError(f"{kind!r} is not a chemically peculiar kind ({', '.join(KINDS)})")


def polar_field_gauss(kind, u):
    """The polar magnetic field, gauss, of a peculiar star of `kind` from the uniform number `u`: an Ap star's log-uniform over `AP_POLAR_FIELD_GAUSS`, an Am
    star's None (no field is detected in them)."""
    if kind == AP:
        low, high = AP_POLAR_FIELD_GAUSS
        return low * (high / low) ** u
    if kind == AM:
        return None
    raise ValueError(f"{kind!r} is not a chemically peculiar kind ({', '.join(KINDS)})")


def of_star(yerkes_class, mass_kg, radius_km, temperature_k, age_gy):
    """
    The `Peculiarity` of a star, or None for one that is not chemically peculiar: only a main-sequence dwarf (`DWARF_CLASS`) in the temperature range of a
    kind may be. The kind, the period and the field come from three uniform numbers of a stream keyed by the star's mass, temperature and age alone, not
    from the draws of the run, so a star is the same kind with the same period and field every time it is built or read.
    """
    if yerkes_class != DWARF_CLASS or None in (mass_kg, radius_km, temperature_k, age_gy):
        return None
    stream = draw.Stream(f"peculiar:{mass_kg!r}:{temperature_k!r}:{age_gy!r}")
    u_kind, u_period, u_field = stream.random(), stream.random(), stream.random()
    kind = kind_of(temperature_k, u_kind)
    if kind is None:
        return None
    return Peculiarity(kind, rotation_period_hours(kind, radius_km / SOLAR_RADIUS_KM, u_period), polar_field_gauss(kind, u_field))
