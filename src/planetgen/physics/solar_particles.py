# planetgen/physics/solar_particles.py

"""
Solar energetic particle events: the event set, the yearly spectrum and the dose it brings (GEN.747)
=====================================================================================================

The stellar-particle dose of a body (`radiation.sep_dose_msv_yr`) used a free-space constant of 20 mSv/yr and an exponential in the
proton's rigidity (`R0` 0.2 GV). This module replaces both with the particles of the events that were measured.

**The events.** 47 Band-function fits of the event-integrated fluence of the solar proton events of 1997 to 2015 (`EVENTS`): 21 events that
reached the ground (the ground-level enhancements, GLE, with four of them in a second, storm-shock fit, `esp`) and 26 that did not (sub-GLE; as fitted and tabulated in Vainio et al. 2017, A&A 604:A8, Tables A.1 and A.2; typed second-hand into
`reports/secondhand-citations-2026-10-24-files/sep-band-parameters-vainio2017.csv`). Each is the integral rigidity spectrum of
`surface_dose.BandSpectrum`, `J(>R) = J0 R^-g1 exp(-R/R0)` up to `R1 = (g2 - g1) R0` and a steeper power law above, `J` in cm-2 for the
whole event at 1 AU and `R` in GV. The spectra are hard or soft, the fluences span five orders of magnitude, and a few events carry the dose.

**The rate.** The events are all the events of `RECORD_YEARS` years that the fits cover and that were large enough to fit. They are the
events of 100 pfu (>10 MeV peak flux) and above: NOAA SWPC's list of 255 events of 10 pfu or more in 42 years (1976 to 2017;
`NOAA_CUMULATIVE_COUNTS`) has 2.6 a year above 100 pfu and the set has 2.4 (`events_per_year_above`, tested), and the smaller ones add under
a few percent of the fluence. So the yearly mean spectrum of the Sun at 1 AU is the sum of the 47 spectra over the record (`yearly_fluence`),
with no model of the size distribution between them.

**The dose.** A proton of kinetic energy `T` crossing the body of a reference person deposits its stopping power in water, `S(T)`
(PSTAR, `proton_range` table `H2O liquid`), and the equivalent dose is `wR` (2 for protons, ICRP 103) times that. The person is taken at the
depth of the blood-forming organs, `BFO_DEPTH_G_CM2` of water (5 g/cm2, NASA's standard): a proton needs 78 MeV to reach it, so the soft
storm-shock protons that dominate the fluence below 30 MeV do not count (the skin depth of 1 g/cm2 would give 60 Sv a year). The sum over the
set is `FREE_MSV_YR`, the mean yearly dose in free space at 1 AU from the Sun with no field and no air. A column of air and a magnetic cutoff
remove the protons below their rigidity: `dose_share(rigidity_gv, geomagnetic_cutoff_gv)` is the share of the free dose left, with the
cutoff averaged over the latitudes (area weighted, `Rc cos^4(latitude)`, the dipole's Stormer law).

**What it does not do.** The star's flare rate and the distance scale it (`radiation.free_sep_msv_yr`, the activity module GEN.493); the dose is
the primary protons' only (no cascade, no helium or heavier ions, no energy lost in the column: a proton that crosses it is taken to deposit its
incident `S(T)`), and the 18 years are 1.6 solar cycles of one star. An event of 1e10 cm-2 above 30 MeV, bigger than any in the record (about 1 in 130 years, Usoskin and Kovaltsov 2012), adds about 8 percent of the mean fluence above 30 MeV (1e10 cm-2 in 130 years against 9e8 a year) and is left out. These are stand-ins and approximations stated in the gap-filling plan
(section 6.15).
"""

import bisect
import math

from planetgen.physics import proton_range
from planetgen.physics.surface_dose import BandSpectrum

EVENTS = (
    ('GLE', '55', '1997-Nov-06', 8.15e8, 0.28, 5.38, 0.116),
    ('GLE', '56', '1998-May-02', 8.98e6, 1.31, 6.51, 0.196),
    ('GLE', '57', '1998-May-06', 1.64e6, 1.92, 7.46, 0.202),
    ('GLE', '58', '1998-Aug-24', 2.1e5, 2.98, 5.27, 0.677),
    ('GLE', '58 esp', '1998-Aug-24', 1.63e4, 5.26, 7.74, 0.983),
    ('GLE', '59', '2000-Jul-14', 2.94e9, 0.51, 7.46, 0.123),
    ('GLE', '59 esp', '2000-Jul-14', 6.01e7, 3.24, 7.85, 0.226),
    ('GLE', '60', '2001-Apr-15', 5.22e7, 1.39, 5.69, 0.260),
    ('GLE', '61', '2001-Apr-18', 8.39e6, 1.85, 5.02, 0.237),
    ('GLE', '62', '2001-Nov-04', 2.14e9, 0.24, 6.67, 0.093),
    ('GLE', '62 esp', '2001-Nov-04', 4.78e8, 2.36, 11.2, 0.129),
    ('GLE', '63', '2001-Dec-26', 2.17e7, 1.81, 7.86, 0.180),
    ('GLE', '64', '2002-Aug-24', 5.06e6, 2.36, 6.70, 0.225),
    ('GLE', '65', '2003-Oct-28', 8.44e9, 0.01, 6.48, 0.089),
    ('GLE', '65 esp', '2003-Oct-28', 1.12e8, 2.81, 8.92, 0.171),
    ('GLE', '66', '2003-Oct-29', 7.62e7, 2.04, 6.86, 0.206),
    ('GLE', '67', '2003-Nov-02', 2.27e6, 3.50, 7.01, 0.321),
    ('GLE', '68', '2005-Jan-17', 3.51e7, 2.65, 8.29, 0.162),
    ('GLE', '69', '2005-Jan-20', 3.81e8, 0.72, 5.78, 0.204),
    ('GLE', '70', '2006-Dec-13', 1.33e8, 1.05, 5.80, 0.177),
    ('GLE', '71', '2012-May-17', 1.03e7, 1.36, 6.96, 0.219),
    ('subGLE', '', '1997-Nov-04', 2.97e6, 1.62, 8.89, 0.170),
    ('subGLE', '', '1998-Nov-14', 2.22e7, 1.39, 8.20, 0.128),
    ('subGLE', '', '1998-Sep-30', 3.64e5, 3.86, 10.9, 0.214),
    ('subGLE', '', '2000-Nov-09', 3.64e10, 0.19, 10.9, 0.082),
    ('subGLE', '', '2001-Apr-03', 1.35e9, 0.61, 7.98, 0.066),
    ('subGLE', '', '2001-Apr-12', 3.63e5, 2.42, 6.52, 0.331),
    ('subGLE', '', '2001-Aug-16', 5.47e8, 0.23, 9.12, 0.105),
    ('subGLE', '', '2001-Nov-22', 2.44e11, 0.04, 6.33, 0.038),
    ('subGLE', '', '2001-Oct-22', 4.56e5, 1.98, 4.47, 0.313),
    ('subGLE', '', '2001-Sep-24', 4.71e8, 1.93, 12.3, 0.110),
    ('subGLE', '', '2002-Apr-21', 6.13e8, 1.36, 10.7, 0.105),
    ('subGLE', '', '2004-Nov-01', 1.04e8, 0.28, 6.17, 0.068),
    ('subGLE', '', '2004-Nov-10', 2.41e8, 1.13, 5.99, 0.065),
    ('subGLE', '', '2005-Jun-16', 1.19e8, 0.04, 4.36, 0.090),
    ('subGLE', '', '2005-Sep-07', 5.94e6, 3.10, 5.55, 0.317),
    ('subGLE', '', '2006-Dec-06', 4.18e7, 2.32, 8.97, 0.170),
    ('subGLE', '', '2012-Jan-23', 1.96e5, 4.96, 9.65, 0.403),
    ('subGLE', '', '2012-Jan-27', 2.22e7, 1.94, 6.04, 0.161),
    ('subGLE', '', '2012-Mar-07', 1.53e9, 1.06, 7.73, 0.112),
    ('subGLE', '', '2012-Mar-13', 3.33e7, 1.30, 7.45, 0.087),
    ('subGLE', '', '2013-May-22', 5.07e5, 3.77, 8.25, 0.325),
    ('subGLE', '', '2014-Feb-25', 1.74e6, 3.46, 5.49, 0.620),
    ('subGLE', '', '2014-Jan-06', 6.35e6, 1.08, 8.34, 0.180),
    ('subGLE', '', '2014-Jan-07', 5.05e6, 3.06, 11.7, 0.154),
    ('subGLE', '', '2014-Sep-01', 1.87e6, 1.54, 4.57, 0.363),
    ('subGLE', '', '2015-Oct-29', 1.1e9, 0.51, 4.51, 0.048),
)
"""tuple: `(kind, event, date, J0 cm-2, gamma1, gamma2, R0 GV)` of the 47 fits: kind GLE or subGLE, the event number of the GLE ("58 esp" is the second
fit of the storm-shock particles of GLE 58) or "" for a sub-GLE, the date of the event, and the four Band parameters (Vainio et al. 2017,
Tables A.1 and A.2). The errors of the paper are not kept. Source tier T2 (the paper is read second-hand and the table typed from it)."""
RECORD_YEARS = 18.0
"""float: The years the set covers: the events run from 1997-11-04 to 2015-10-29 (17.98 years). Source tier T3, tweak level L0."""

NOAA_YEARS = 42.0
NOAA_CUMULATIVE_COUNTS = ((10.0, 255), (30.0, 177), (100.0, 109), (300.0, 72), (1000.0, 40), (3000.0, 18), (10000.0, 9))
"""tuple: `(peak flux above 10 MeV, pfu; events at or above it)` in the NOAA SWPC list of 1976 to 2017 (`NOAA_YEARS`), counted from the 255
events of 10 pfu or more (reports/best-guess-standins-2026-10-24.md section 17). Source tier T1 (the list is NOAA's; the counts are computed)."""
TAIL_ANCHOR_PFU = 1.2e5
"""float: The peak flux above 10 MeV, pfu, paired with an annual fluence above 30 MeV of 1e10 cm-2 (`TAIL_ANCHOR_RATE_PER_YEAR`); the pairing is
McCracken's for the Carrington event, recalled (reports/best-guess-standins-2026-10-24.md section 17). [Stand-in SE4, tier T3, tweak level L1.]"""
TAIL_ANCHOR_RATE_PER_YEAR = 0.0077
"""float: The events a year with an annual fluence above 30 MeV of 1e10 cm-2: four in the 600 years of annual 10Be data, 0.0077 (+0.0073, -0.0045,
90 percent) with Poisson errors (Usoskin and Kovaltsov 2012, ApJ 757, 92, section 4, READ). Nothing above 2e10 in 600 years (median 0.0012 a year,
one-sided 90 percent under 0.0027) and nothing above 5e10 in 11 400 years; the Carrington event itself is not supported by the isotope data.
[Stand-in SE4, tier T1, tweak level L0.]"""
EXTREME_SLOPE = 4.0
"""float: The cumulative power-law slope of the events past `TAIL_ANCHOR_PFU`: Usoskin and Kovaltsov 2012 fit the tail above 7e9 cm-2 by
p proportional to F30^-4 (an exponential fits as well; section 4, READ). [Stand-in SE4, tier T1, tweak level L1.]"""
TAIL_SLOPE = math.log(NOAA_CUMULATIVE_COUNTS[-1][1] / NOAA_YEARS / TAIL_ANCHOR_RATE_PER_YEAR) / math.log(
    TAIL_ANCHOR_PFU / NOAA_CUMULATIVE_COUNTS[-1][0])
"""float: The cumulative slope of the number of events from 10 000 pfu to `TAIL_ANCHOR_PFU`, about 1.34: the straight line in log-log between the
last NOAA node (9 events in 42 years) and the isotope-record rate at the anchor (Usoskin and Kovaltsov 2012). Replaces the 1.9 of GEN.747, which was
fitted to one Carrington-class event in 500 years. [Stand-in SE4, tier T3, tweak level L1.]"""

BFO_DEPTH_G_CM2 = 5.0
"""float: The depth of water at which the dose is taken, g/cm2: the blood-forming organs (NASA's standard 5 cm; recalled). Source tier T3."""
PROTON_WEIGHT = 2.0
"""float: ICRP 103 radiation weight of a proton (`surface_dose.WEIGHT_PROTON_PION`)."""
GY_PER_MEV_PER_G = 1.602176634e-10
"""float: One MeV per gram, in gray."""

GRID_POINTS = 600
GRID_MIN_MEV = 10.0
GRID_MAX_MEV = 20000.0


def _spectra():
    return [BandSpectrum(j0, gamma1, gamma2, r0) for _, _, _, j0, gamma1, gamma2, r0 in EVENTS]


SPECTRA = _spectra()
"""list: The `surface_dose.BandSpectrum` of each of `EVENTS`."""


def yearly_fluence(rigidity_gv):
    """The mean number of protons per cm2 and year above `rigidity_gv` at 1 AU: the 47 spectra summed over the record."""
    return sum(spectrum.flux(rigidity_gv) for spectrum in SPECTRA) / RECORD_YEARS


def _stopping_mev_cm2_g(energy_mev):
    return proton_range.stopping_power_mev_cm2_g(energy_mev, "H2O liquid")


def _build_dose_table():
    """The dose, mSv/yr at 1 AU, of the protons above each rigidity of a log grid (cumulative from the top), for protons that reach the organs."""
    grid = [GRID_MIN_MEV * (GRID_MAX_MEV / GRID_MIN_MEV) ** (i / (GRID_POINTS - 1)) for i in range(GRID_POINTS)]
    rigidities = [proton_range.rigidity_gv(t) for t in grid]
    flux = [yearly_fluence(r) for r in rigidities]
    organ_mev = proton_range.energy_mev(BFO_DEPTH_G_CM2, "H2O liquid")
    scale = PROTON_WEIGHT * GY_PER_MEV_PER_G * 1e3
    cumulative = [0.0] * GRID_POINTS
    for i in range(GRID_POINTS - 2, -1, -1):
        middle = math.sqrt(grid[i] * grid[i + 1])
        step = (flux[i] - flux[i + 1]) * scale * _stopping_mev_cm2_g(middle) if middle >= organ_mev else 0.0
        cumulative[i] = cumulative[i + 1] + step
    return rigidities, cumulative, proton_range.rigidity_gv(organ_mev)


_RIGIDITIES_GV, _DOSE_ABOVE, ORGAN_RIGIDITY_GV = _build_dose_table()
FREE_MSV_YR = _DOSE_ABOVE[0]
"""float: The mean yearly dose of the protons that reach the organs, mSv/yr, in free space at 1 AU from the Sun with no field and no air (about 265)."""


def dose_above_msv_yr(rigidity_gv):
    """The yearly dose, mSv/yr at 1 AU, of the protons above `rigidity_gv` (log-log between the grid's points; nil above its top)."""
    if rigidity_gv <= ORGAN_RIGIDITY_GV:
        return FREE_MSV_YR
    if rigidity_gv >= _RIGIDITIES_GV[-1]:
        return 0.0
    i = bisect.bisect_right(_RIGIDITIES_GV, rigidity_gv) - 1
    low, high = _DOSE_ABOVE[i], _DOSE_ABOVE[i + 1]
    if high <= 0.0:
        return low * (1.0 - (rigidity_gv - _RIGIDITIES_GV[i]) / (_RIGIDITIES_GV[i + 1] - _RIGIDITIES_GV[i]))
    x = math.log(rigidity_gv / _RIGIDITIES_GV[i]) / math.log(_RIGIDITIES_GV[i + 1] / _RIGIDITIES_GV[i])
    return low * (high / low) ** x


def dose_share(rigidity_gv, geomagnetic_cutoff_gv=0.0, steps=90):
    """
    The share of `FREE_MSV_YR` that gets past a column that stops protons below `rigidity_gv` and a dipole whose equatorial vertical cutoff is
    `geomagnetic_cutoff_gv`: the dose of the protons above `max(rigidity_gv, Rc cos^4(latitude))` averaged over the latitudes with the
    weight `cos(latitude)` (equal areas of a sphere), over the free dose. No column and no field give 1.
    """
    floor = max(rigidity_gv, ORGAN_RIGIDITY_GV)
    if geomagnetic_cutoff_gv <= floor:
        return min(dose_above_msv_yr(floor) / FREE_MSV_YR, 1.0)
    total = 0.0
    weights = 0.0
    for i in range(steps + 1):
        latitude = (math.pi / 2.0) * i / steps
        weight = math.cos(latitude) * (0.5 if i in (0, steps) else 1.0)
        total += weight * dose_above_msv_yr(max(floor, geomagnetic_cutoff_gv * math.cos(latitude) ** 4))
        weights += weight
    return min(total / weights / FREE_MSV_YR, 1.0)


def events_per_year_above(peak_flux_pfu):
    """The events a year of the NOAA list at or above a peak flux (>10 MeV, pfu): straight lines in log-log between `NOAA_CUMULATIVE_COUNTS`,
    `TAIL_SLOPE` from 10 000 pfu to `TAIL_ANCHOR_PFU` and `EXTREME_SLOPE` beyond; the count of the first node (255 in 42 years) held below
    10 pfu."""
    nodes = NOAA_CUMULATIVE_COUNTS
    if peak_flux_pfu <= nodes[0][0]:
        return nodes[0][1] / NOAA_YEARS
    if peak_flux_pfu >= TAIL_ANCHOR_PFU:
        return TAIL_ANCHOR_RATE_PER_YEAR * (peak_flux_pfu / TAIL_ANCHOR_PFU) ** -EXTREME_SLOPE
    if peak_flux_pfu >= nodes[-1][0]:
        return nodes[-1][1] / NOAA_YEARS * (peak_flux_pfu / nodes[-1][0]) ** -TAIL_SLOPE
    for (p0, n0), (p1, n1) in zip(nodes, nodes[1:]):
        if peak_flux_pfu <= p1:
            return math.exp(math.log(n0) + math.log(peak_flux_pfu / p0) / math.log(p1 / p0) * math.log(n1 / n0)) / NOAA_YEARS
    return nodes[-1][1] / NOAA_YEARS
