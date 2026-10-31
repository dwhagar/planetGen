# planetgen/physics/radiation.py

"""
Surface radiation dose, UV and galactic hazards (GEN.87).

The rules are docs/design/habitability/activity-magnetism-radiation-hydrosphere.md
section 4 (design defaults). A rocky body's yearly dose, mSv/yr, is

    surface = gcr * helio_mult + sep + ground

1. **Cosmic rays (`dose_gcr_msv_yr`, GEN.453).** `D = D0(X) f_B(M, X) f_sol` under a column of air `X = P0 / g`
   (g/cm2). `D0` is the no-field curve: the proton-primary air showers of Atri et al. 2013 (553, 263, 17.2, 3.02 and 0.65
   mSv/yr at X of 100, 200, 500, 700 and 1036), log-log between them; below 100 g/cm2 the primaries thinned by collisions plus
   the shower building up behind them (`thin_air_dose_msv_yr`, GEN.746): the Moon at X 0 (500 mSv/yr, Chang'E-4), Mars at X 22.6 (232
   mSv/yr and its measured fall with depth, Hassler 2014 and Guo 2018) and the shower's peak at the first run, 553 at 100; above 1036 the larger of the
   table's last slope and the muon tail (Crouch). `f_B` is the dipole's cut of that dose, the ratio of Atri's and Griessmeier's
   runs to the no-field run (it depends on the column), looked up at the field's equivalent moment `M = Rc0 / 14.34 GV`
   (a body of another radius has the cutoff of an Earth-sized body of that moment); `f_sol` is the solar-modulation
   spread (0.5 to 1.64 of the mean, Mars, a seeded draw). The cutoff rigidity is `Rc0 = 14.34 GV (M_dip f_dip / 7.768e22)
   (R_E / R_p)^2` (Stormer: moment over radius squared).
2. **Compressed heliosphere (`dose_helio_mult`).** `1 + 1.5 exp(-X / 100)
   c`, at most 2.5, `c` running from 0 in open space to 1 when a cloud
   presses the star's heliosphere to 1 AU or less (`heliosphere_compression`).
   Open space at generation; `store.refresh_containment` updates it when
   the system is found inside a nebula or remnant.
3. **Stellar particles (`dose_sep_msv_yr`, GEN.747).** A free-space `FREE_MSV_YR sqrt(N33 / N33_sun) / d_AU^2` (`FREE_MSV_YR` about 265 mSv/yr) times the
   square root of the flares' energy weight for the star's own frequency slope (`flare_energy_weight`, `flare_alpha`, GEN.805: a shallow distribution puts
   more of the energy in the big flares), the dose to the organs of the protons of the 47 measured solar events of 1997 to 2015 summed over their 18
   years (`solar_particles`), times the share of that dose that gets past the column's proton cutoff (`proton_range`: NIST PSTAR) and the dipole's
   `Rc0 cos^4(latitude)`, averaged over the latitudes (`solar_particles.dose_share`; 1 without air or field).
4. **Ground (`dose_ground_msv_yr`).** The gamma dose in air at 1 m over the surface rock and the radon of its uranium, from the
   body's own stored potassium, thorium and uranium (`radiogenic`, GEN.459), at its star's age: the same concentrations its interior
   heat reads, so there is no independent dose draw (`a_rad` is gone). The surface rock is the area-weighted rock of the body's stored
   ground (`surface_composition.land_rock_ppm`, GEN.441); the dose is the open-air, full-time effective dose
   of a reference person (radiation-model.md R12). Radon has no pressure gate (R16): its flux is drawn between 1 and
   `sqrt(P0 / P)` times the 1-bar value.
5. **UV (`uv_surface_index`).** The Setlow DNA-damage dose at the ground as a multiple of Earth's (Earth is 1 by definition, GEN.87), from
   the two-stream spectrum below the body's own air, ozone column (`ozone_column_du`, GEN.413: the Kozakis fit by oxygen and host star,
   `physics/ozone.py`), haze, dust and cloud, lit by the star's own spectrum with its stored UV-C share (`physics/body_radiation.py`, GEN.443); a
   body with no air reads the star's unfiltered light (GEN.445; docs/design/composition/radiation-model.md section 7). There is no proxy: no
   ratio to the Sun, no ozone exponent.
6. **Ozone loss (`ozone_loss_flag`, GEN.805, GEN.808).** True when there is ozone and either the system's lethal galactic events exceed one
   per 100 Myr or the planet's flare irradiation index (`flare_irradiation_index`, `N33 / d_AU^2`: Earth 0.003, Proxima b about 2,100) times
   the field's shielding `f_B` reaches `FLARE_OZONE_LIMIT` (300; `f_B` is 1 without a dipole, 0.3 for a magnetopause under 3 radii and 0.1
   above: the design note section 2.3). A world that loses its ozone keeps `OZONE_LOSS_REMAINING` of it (Segura et al. 2010: 94 percent
   removed by one large flare; Howard et al. 2018: 90 percent by repeated flaring): the stored `ozone_column_du` is the column that is left, the
   UV index, the habitability score and the colour of the sky and ground are worked out from it. A loss event that removes air (stripping,
   impact, escape) is not applied to the stored ozone column here: the final-age refresh (GEN.839) works the dose and the ozone again from the air that is left.
7. **Galactic hazard (`lethal_event_rate_per_gyr`, on the star).** A
   supernova within 8 pc, 1.5 per Gyr at 8 kpc from the centre, times
   `exp((8 kpc - R) / 1.45 kpc)` (250x at the centre), plus gamma-ray bursts
   (2 per Gyr at the Sun, falling outward) and nuclear-activity episodes (10
   per Gyr inside 2.9 kpc), GEN.189.
8. **Trapped particles (`dose_belt_msv_yr`, GEN.983).** A moon in the belts of its parent's field takes their dose at its orbit, from the half of the sky
   the body does not block, cut by its own air (`radiation_belts.moon_exposure_msv_yr`, `apply_parent_belts`). A planet's own belts lie above its air
   and give its ground 0; a giant has no surface.

A planet of a pulsar or a black hole is rated by GEN.89, not here.
"""

import math

from planetgen import tuning
from planetgen.physics import body_radiation, constants, magnetism, mixture, ozone, proton_range, radiation_belts, radiogenic, solar_particles, surface_composition
from planetgen.util import draw
from planetgen.util.random import log_uniform

BODY_FIELDS = (
    "surface_dose_msv_yr", "dose_gcr_msv_yr", "dose_sep_msv_yr", "dose_ground_msv_yr", "dose_helio_mult",
    "uv_surface_index", "ozone_loss_flag", "ozone_column_du", "dose_belt_msv_yr",
)
"""tuple: What `generate_dose` sets on a planet or moon, as stored in
`planets` and `moons`."""

GCR_MOON_MSV_YR = 500.0
"""float: The no-field dose with no air, mSv/yr (Chang'E-4 LND, 1369 uSv/day, a solar-minimum bound; Zhang 2020). Source tier T1."""
GCR_MARS_G_CM2 = 22.6
GCR_MARS_MSV_YR = 232.0
"""float: Mars at Gale (MSL RAD, 0.64 mSv/day, GCR only; Hassler 2014, Table 3; Guo 2015 gives 255 mSv/yr with a range of 128 to 420) under the
column of its mean pressure, 840 Pa over 3.71 m/s2 = 22.6 g/cm2 (Guo 2018; the 21 of the earlier stand-in was recalled). Source tier T1."""
GCR_MARS_LOG_SLOPE_PER_G_CM2 = -0.19 / (25.7 - 18.9)
"""float: d ln D / dX at Mars: the modelled dose equivalent falls by 126 uSv/day, 19 percent of its mean, between the column of 700 Pa (18.9 g/cm2)
and of 950 Pa (25.7 g/cm2) and "declines continuously and nonlinearly" with depth (Guo et al. 2018, arXiv:1712.06987, read for GEN.746; the paper's
dose-against-depth figures are not tabulated). Source tier T2."""
GCR_PEAK_G_CM2 = 100.0
"""float: Where the shower dose peaks, g/cm2: "about 100" (Guo 2018) and 122.7 +- 6.3 (the Pfotzer maximum, Sarkar 2017). The curve below the runs has
its peak at the first run of Atri 2013 (`GCR_SHOWER_POINTS`), so it joins them at its highest point. Source tier T2."""
GCR_SHOWER_POINTS = ((100.0, 553.0), (200.0, 263.0), (500.0, 17.2), (700.0, 3.02), (1036.0, 0.65))
"""tuple: `(X g/cm2, mSv/yr)`: total biological dose from muons, electrons, photons and neutrons in proton-primary CORSIKA
air showers with no magnetic field (Atri, Hariharan, Fryer et al. 2013, Table XV, via Griessmeier et al. 2016, Table 2).
Helium and heavier primaries are absent, so these are lower bounds."""
GCR_SEA_LEVEL_G_CM2 = GCR_SHOWER_POINTS[-1][0]

FIELD_COLUMNS_G_CM2 = (100.0, 200.0, 500.0, 700.0, 1036.0)
"""tuple: The columns the field's cut was run at (Atri 2013)."""
FIELD_FULL_MOMENTS = (0.0, 0.1, 0.15, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 6.0, 10.0)
"""tuple: The dipole moments, Earth's = 1, of the runs at 100 and 1036 g/cm2."""
FIELD_DOSE_100 = (553.0, 527.0, 510.0, 405.0, 257.0, 216.0, 172.0, 53.0, 15.0, 5.7, 2.3)
FIELD_DOSE_1036 = (0.65, 0.48, 0.46, 0.44, 0.42, 0.39, 0.34, 0.28, 0.23, 0.19, 0.10)
"""tuple: Dose, mSv/yr, at each of `FIELD_FULL_MOMENTS`."""
FIELD_MIDDLE_MOMENTS = (0.0, 0.15, 0.5, 1.0)
FIELD_MIDDLE_DOSE = {
    200.0: (263.0, 260.0, 141.0, 97.0), 500.0: (17.2, 18.2, 12.0, 8.7), 700.0: (3.02, 3.38, 2.53, 1.84),
}
"""dict: Dose, mSv/yr, of the runs between those columns, at `FIELD_MIDDLE_MOMENTS` (Atri 2013, Table XV). A field cannot
add primaries, so a ratio the runs put above 1 (a weak field at 500 and 700 g/cm2, within the runs' scatter) is held at 1."""
FIELD_MAX_MOMENT = FIELD_FULL_MOMENTS[-1]
"""float: The strongest field the runs reach (10 Earth moments); a stronger one takes this ratio, an upper bound."""

"""float: The solar-cycle spread of the dose as fractions of its mean (Mars: 0.35 to 1.15 mSv/day around 0.7; a total spread
of 3.3, Guo 2015). `tuning.GCR_SOLAR_MODULATION_RANGE` is the draw."""

MUON_LOG_A = (-11.22, -14.10)
MUON_RATE_PER_HG_CM2 = (0.00262, 0.001213)
MUON_FLOOR = 2.17e-13
MUON_MEV_PER_G = 2.0
MUON_J_PER_MEV = 1.602176634e-13
SECONDS_PER_YEAR = 3.15576e7
"""The Crouch fit of the vertical muon intensity (cm-2 s-1 sr-1) at depth h in hg/cm2 (m of water): `exp(-11.22 - 0.00262 h)
+ exp(-14.10 - 0.001213 h) + 2.17e-13` (docs/design/composition/radiation-model.md section 5), valid from about 1 km of water; the dose is
pi I times 2 MeV cm2/g. It is 9 times below a published value at Titan (10 850 g/cm2), so it is a floor under the table's tail."""

EARTH_DIPOLE_A_M2 = 7.768e22
HELIO_MAX_EXTRA = 1.5
HELIO_FADE_G_CM2 = 100.0

N33_SUN_PER_YR = 10 ** -2.5
"""float: The Sun's rate of 1e33 erg flares (one in 316 years), the unit of a star's `flare_n33_per_yr`: the dose of another star's particles
scales with the square root of its rate over this (`free_sep_msv_yr`, a design default pending GEN.493)."""

EARTH_AGE_GY = 4.5
"""float: The age a star without one is given (the Sun's, rounded)."""
SV_PER_GY = 0.7
"""float: Effective dose per absorbed dose in air for gamma rays from the ground, UNSCEAR (RECALLED, part E 5; radiation-model.md
section 6)."""
HOURS_PER_YEAR = 8766.0
GAMMA_AIR_ATTENUATION_CM2_PER_G = 0.0637
"""float: Gamma dose under a dense gas falls as exp(-this * density g/cm3 * 100 cm) (radiation-model.md section 6). The Earth's
coefficients are valid for a density under about 0.03 g/cm3 (F9); beyond it the factor is an extrapolation."""
GAMMA_AIR_HEIGHT_CM = 100.0

RADON_EMANATION = 0.2
"""float: The share of radon that leaves the grains (radiation-model.md section 6, part 07 2.10)."""
RADON_BULK_DENSITY_KG_M3 = 1600.0
RADON_DIFFUSION_M2_S = 2e-6
"""float: The effective diffusion coefficient of radon in soil at one bar; its dependence on pressure is not reached, so the flux is
drawn between 1 and `sqrt(P0/P)` times the 1-bar value (R16, Rule applied: most permissive)."""
RADON_HALF_LIFE_S = 3.824 * 86400.0
RADON_MEAN_LIFE_S = RADON_HALF_LIFE_S / math.log(2)
"""float: 5.52 days."""
RADON_BOX_HEIGHT_KM = 1.0
RADON_REFERENCE_SCALE_HEIGHT_KM = 8.5
"""float: The box a radon flux is spread through is 1 km tall for the Earth's scale height and scales with it (design inference,
R16)."""
RADON_SV_PER_BQ_H_M3 = 9e-9
RADON_EQUILIBRIUM_OUTDOORS = 0.6
"""float: The dose of 9 nSv per Bq h/m3 of radon at 0.6 of equilibrium outdoors (UNSCEAR; part 07 2.10)."""

FLARE_OZONE_LIMIT = 300.0
"""float: The flare irradiation index times `f_B` above which flares keep a world's ozone down (the design note section 2.3, a design default):
Proxima b without a field is flagged and Earth is not."""
SHIELD_NO_FIELD, SHIELD_WEAK_FIELD, SHIELD_STRONG_FIELD = 1.0, 0.3, 0.1
"""float: `f_B`, the share of the flares' particles that reach the ozone layer: no dipole, a magnetopause under `WEAK_MAGNETOPAUSE_RP`, above it."""
WEAK_MAGNETOPAUSE_RP = 3.0
"""float: The standoff, in radii of the body, under which a field is too squeezed to shield well."""
OZONE_LOSS_REMAINING = 0.1
"""float: The share of its ozone a world keeps when flares hold the flag (Segura et al. 2010: 94 percent removed; Howard et al. 2018: 90 percent)."""
FLARE_ENERGY_SPAN = 1.0e3
"""float: The span of flare energies, over the 1e33 erg of `N33`, the energy weight is worked out over (1e33 to 1e36 erg: the largest superflares seen
on solar-type stars are 1e35 to 1e36). A design default."""

LETHAL_EVENT_LIMIT_PER_GYR = 10.0
"""float: One lethal event per 100 Myr: the rate above which galactic hazards
flag ozone loss (Boss's default, GEN.87)."""


def _log_linear(x, xs, ys):
    """`y` at `x` along the straight lines of ln y against x through `(xs, ys)` (`xs` rising), held at the ends."""
    if x <= xs[0]:
        return ys[0]
    for (x0, y0), (x1, y1) in zip(zip(xs, ys), zip(xs[1:], ys[1:])):
        if x <= x1:
            return math.exp(math.log(y0) + (x - x0) / (x1 - x0) * (math.log(y1) - math.log(y0)))
    return ys[-1]


def _muon_dose_msv_yr(column_g_cm2):
    """The muon dose, mSv/yr, behind `column_g_cm2` by the Crouch fit (a floor under the table's tail)."""
    depth = column_g_cm2 / 100.0
    intensity = (math.exp(MUON_LOG_A[0] - MUON_RATE_PER_HG_CM2[0] * depth)
                 + math.exp(MUON_LOG_A[1] - MUON_RATE_PER_HG_CM2[1] * depth) + MUON_FLOOR)
    return math.pi * intensity * MUON_MEV_PER_G * MUON_J_PER_MEV / 1e-3 * SECONDS_PER_YEAR * 1e3


def _fit_thin_air_curve():
    """
    The three numbers of the thin-air dose curve, `(primary_length, build_up_exponent, build_up_length)` (rule 1): the dose below the first
    shower run is the primaries thinned by nuclear collisions plus the shower that builds up behind them,

        D(X) = D_moon exp(-X / L_p) + K (X / X_1)^a exp(-(X - X_1) / L_s),       K = D_1 - D_moon exp(-X_1 / L_p),

    and the three numbers are the ones that make it pass through Mars (value `GCR_MARS_MSV_YR` and log slope
    `GCR_MARS_LOG_SLOPE_PER_G_CM2`, at `GCR_MARS_G_CM2`) and peak at the first run `(X_1, D_1)` of Atri 2013. The Moon is the value at X = 0
    (the shower is nil there). For a trial `L_p` the shower's `a` and `1 / L_s` follow from the Mars value and the peak by two linear
    equations, and the Mars slope is then a function of `L_p` alone, found by bisection.
    """
    x1, d1 = GCR_PEAK_G_CM2, GCR_SHOWER_POINTS[0][1]
    xm, dm, slope_m = GCR_MARS_G_CM2, GCR_MARS_MSV_YR, GCR_MARS_LOG_SLOPE_PER_G_CM2 * GCR_MARS_MSV_YR

    def solve(length):
        primary = lambda x: GCR_MOON_MSV_YR * math.exp(-x / length)
        k = d1 - primary(x1)
        # peak: a / x1 - 1 / ls = primary(x1) / (length k);  Mars value: a ln(xm / x1) - (xm - x1) / ls = ln((dm - primary(xm)) / k)
        a11, a12, b1 = 1.0 / x1, -1.0, primary(x1) / (length * k)
        a21, a22, b2 = math.log(xm / x1), -(xm - x1), math.log((dm - primary(xm)) / k)
        det = a11 * a22 - a12 * a21
        exponent = (b1 * a22 - a12 * b2) / det
        inverse_length = (a11 * b2 - a21 * b1) / det
        shower = (dm - primary(xm))
        residual = -primary(xm) / length + shower * (exponent / xm - inverse_length) - slope_m
        return exponent, 1.0 / inverse_length, residual

    low, high = 10.0, 0.999 * xm / math.log(GCR_MOON_MSV_YR / dm)  # above this the primaries alone are more than Mars's dose
    for _ in range(80):
        middle = 0.5 * (low + high)
        if (solve(low)[2] > 0) == (solve(middle)[2] > 0):
            low = middle
        else:
            high = middle
    exponent, shower_length, _ = solve(0.5 * (low + high))
    return 0.5 * (low + high), exponent, shower_length


assert GCR_SHOWER_POINTS[0][0] == GCR_PEAK_G_CM2, "the thin-air curve joins the first shower run at the peak"
GCR_PRIMARY_LENGTH_G_CM2, GCR_BUILD_UP_EXPONENT, GCR_SHOWER_LENGTH_G_CM2 = _fit_thin_air_curve()
"""float: The thin-air curve's fitted numbers (28.3 g/cm2, 5.94 and 17.1 g/cm2). The primaries' length of 28 g/cm2 is not a free number: it
falls out of the Mars value and slope, and is the size of the nuclear collision length of the carbon-to-oxygen nuclei that carry most of the
dose equivalent (about 25 g/cm2 in air; the protons' is 90). Source tier T5, tweak level L0 for the length, T6 and L2 for the shape of the
build-up (the form is the Gaisser-Hillas shape of a shower, fitted by its three conditions)."""


def thin_air_dose_msv_yr(column_g_cm2):
    """Rule 1's `D0` below the first shower run (`_fit_thin_air_curve`): 500 mSv/yr (the Moon) at no air, down through Mars's 232 at 21 g/cm2
    to a minimum near 35 g/cm2 and up to the runs' 553 at 100."""
    x1, d1 = GCR_PEAK_G_CM2, GCR_SHOWER_POINTS[0][1]
    x = min(max(column_g_cm2, 0.0), x1)
    primary = GCR_MOON_MSV_YR * math.exp(-x / GCR_PRIMARY_LENGTH_G_CM2)
    if x <= 0.0:
        return primary
    k = d1 - GCR_MOON_MSV_YR * math.exp(-x1 / GCR_PRIMARY_LENGTH_G_CM2)
    return primary + k * (x / x1) ** GCR_BUILD_UP_EXPONENT * math.exp(-(x - x1) / GCR_SHOWER_LENGTH_G_CM2)


def no_field_dose_msv_yr(column_g_cm2):
    """Rule 1's `D0`: the dose, mSv/yr, under `column_g_cm2` with no magnetic field and the mean solar modulation."""
    x = max(column_g_cm2, 0.0)
    points = GCR_SHOWER_POINTS
    if x < points[0][0]:
        return thin_air_dose_msv_yr(x)
    columns = [math.log(column) for column, _ in points]
    doses = [dose for _, dose in points]
    if x <= points[-1][0]:
        return _log_linear(math.log(x), columns, doses)
    slope = math.log(doses[-1] / doses[-2]) / (columns[-1] - columns[-2])
    return max(doses[-1] * (x / points[-1][0]) ** slope, _muon_dose_msv_yr(x))


def _field_ratio_at_column(moment, column_index):
    """The dose with a dipole of `moment` Earth moments over the dose with none at one of `FIELD_COLUMNS_G_CM2`."""
    ends = (FIELD_DOSE_100, FIELD_DOSE_1036)
    if column_index in (0, len(FIELD_COLUMNS_G_CM2) - 1):
        doses = ends[0 if column_index == 0 else 1]
        return _log_linear(moment, FIELD_FULL_MOMENTS, [dose / doses[0] for dose in doses])
    x = FIELD_COLUMNS_G_CM2[column_index]
    doses = FIELD_MIDDLE_DOSE[x]
    if moment <= FIELD_MIDDLE_MOMENTS[-1]:
        ratio = _log_linear(moment, FIELD_MIDDLE_MOMENTS, [dose / doses[0] for dose in doses])
    else:
        # Beyond the runs at this column the strong-field shape is that of the end columns, set at this column.
        def end_ratio(m):
            low = _log_linear(m, FIELD_FULL_MOMENTS, [d / FIELD_DOSE_100[0] for d in FIELD_DOSE_100])
            high = _log_linear(m, FIELD_FULL_MOMENTS, [d / FIELD_DOSE_1036[0] for d in FIELD_DOSE_1036])
            f = math.log(x / FIELD_COLUMNS_G_CM2[0]) / math.log(FIELD_COLUMNS_G_CM2[-1] / FIELD_COLUMNS_G_CM2[0])
            return math.exp(math.log(low) + f * (math.log(high) - math.log(low)))
        ratio = doses[-1] / doses[0] * end_ratio(moment) / end_ratio(FIELD_MIDDLE_MOMENTS[-1])
    return min(ratio, 1.0)


def magnetic_factor(cutoff_gv, column_g_cm2):
    """Rule 1's `f_B`: the dose with a dipole of equatorial cutoff `cutoff_gv` over the dose with none, under
    `column_g_cm2`. The runs reach 100 to 1036 g/cm2 (a thinner column takes the ratio at 100, a thicker one at 1036)
    and 10 Earth moments (a stronger field takes that ratio)."""
    if cutoff_gv <= 0:
        return 1.0
    moment = min(cutoff_gv / EARTH_CUTOFF_GV, FIELD_MAX_MOMENT)
    columns = FIELD_COLUMNS_G_CM2
    x = min(max(column_g_cm2, columns[0]), columns[-1])
    ratios = [min(_field_ratio_at_column(moment, i), 1.0) for i in range(len(columns))]
    return _log_linear(math.log(x), [math.log(c) for c in columns], ratios)


def gcr_dose_msv_yr(column_g_cm2, cutoff_gv=0.0, modulation=1.0):
    """The cosmic-ray dose under `column_g_cm2`, with a dipole of equatorial cutoff `cutoff_gv` (0 for none) and the
    solar modulation `modulation` as a fraction of the mean (`tuning.GCR_SOLAR_MODULATION_RANGE`)."""
    return no_field_dose_msv_yr(column_g_cm2) * magnetic_factor(cutoff_gv, column_g_cm2) * modulation


def cutoff_rigidity_gv(dipole_moment_a_m2, radius_km):
    """`Rc0`, the vertical Stormer cutoff at the magnetic equator, GV, of a dipole moment (A m2, the dipole part) on a
    body of `radius_km`: `c mu0 M / (4 pi r^2) / 4`, moment over radius squared."""
    if not dipole_moment_a_m2 or radius_km <= 0:
        return 0.0
    return STORMER_GV_PER_UNIT * dipole_moment_a_m2 / (radius_km * 1e3) ** 2


STORMER_GV_PER_UNIT = constants.SPEED_OF_LIGHT_M_S * 1e-7 / 4.0 / 1e9
"""float: GV per (A m2 / m2) of the Stormer vertical cutoff, c (mu0 / 4 pi) / 4."""
EARTH_CUTOFF_GV = cutoff_rigidity_gv(EARTH_DIPOLE_A_M2, constants.EARTH_RADIUS_KM)
"""float: Earth's equatorial vertical cutoff, 14.34 GV (the earlier constant 14.9 was 4 percent high)."""


def heliosphere_compression(open_radius_au, compressed_radius_au):
    """0 in open space to 1 when the heliosphere is pressed to 1 AU or
    less (log scale)."""
    if not compressed_radius_au or compressed_radius_au >= open_radius_au or open_radius_au <= 1.0:
        return 0.0
    return min(1.0, math.log(open_radius_au / compressed_radius_au) / math.log(open_radius_au))


def helio_multiplier(column_g_cm2, compression):
    """Rule 2's cosmic-ray multiplier (1 to 2.5)."""
    return 1.0 + HELIO_MAX_EXTRA * math.exp(-column_g_cm2 / HELIO_FADE_G_CM2) * min(max(compression, 0.0), 1.0)


def proton_cutoff_gv(column_g_cm2, material="air"):
    """The rigidity, GV, of the proton whose range is `column_g_cm2` of `material` (`proton_range`: NIST PSTAR)."""
    return proton_range.rigidity_gv(proton_range.energy_mev(column_g_cm2, material))


def sep_dose_msv_yr(free_msv_yr, column_g_cm2, cutoff_gv):
    """Rule 3's stellar-particle dose under `column_g_cm2` for a free-space dose and a dipole's equatorial `cutoff_gv`: the share of the protons
    of the measured events (`solar_particles.dose_share`) that gets past the column's proton cutoff and the field's, latitude-averaged."""
    return free_msv_yr * solar_particles.dose_share(proton_cutoff_gv(column_g_cm2), cutoff_gv)


def flare_energy_weight(flare_alpha):
    """
    The energy a flare distribution of slope `flare_alpha` (`dN/dE ~ E^-alpha`, `N33` flares above 1e33 erg) puts out over 1e33 to 1e36 erg
    (`FLARE_ENERGY_SPAN`), over what a slope of 2 does: 1 for a star with no slope, above 1 for a shallower one and below for a steeper
    (about 2.2 and 0.5 at slopes of 1.8 and 2.2).
    """
    if not flare_alpha:
        return 1.0
    exponent = 2.0 - flare_alpha
    if abs(exponent) < 1e-9:
        return 1.0
    return (FLARE_ENERGY_SPAN ** exponent - 1.0) / (exponent * math.log(FLARE_ENERGY_SPAN))


def free_sep_msv_yr(flare_n33_per_yr, distance_au, flare_alpha=None):
    """Rule 3's free-space dose at `distance_au`: the Sun's mean yearly dose at 1 AU from its measured events (`solar_particles.FREE_MSV_YR`) scaled by
    the square root of the star's rate of 1e33 erg flares over the Sun's, with the flares' energy weight for the star's slope (the dose goes as the
    square root of the energy they put out), and the inverse square of the distance."""
    if not flare_n33_per_yr or not distance_au:
        return 0.0
    return solar_particles.FREE_MSV_YR * math.sqrt(flare_n33_per_yr / N33_SUN_PER_YR * flare_energy_weight(flare_alpha)) / distance_au ** 2


def field_shielding(magnetopause_rp, dipole_class):
    """`f_B`: the share of the flares' particles that reach the ozone layer through the field: 1 for no dipole, `SHIELD_WEAK_FIELD` for a magnetopause
    under `WEAK_MAGNETOPAUSE_RP` radii, `SHIELD_STRONG_FIELD` above."""
    if dipole_class in (None, "none") or not magnetopause_rp:
        return SHIELD_NO_FIELD
    return SHIELD_WEAK_FIELD if magnetopause_rp < WEAK_MAGNETOPAUSE_RP else SHIELD_STRONG_FIELD


def ozone_lost(body):
    """Whether flares or galactic events keep a world's ozone down (rule 6): needs a stored ozone column, and the system's lethal events above the limit or
    the stored flare irradiation index times `field_shielding` at `FLARE_OZONE_LIMIT`."""
    if not (getattr(body, "ozone_column_du", None) or 0.0) > 0.0:
        return False
    hazard = getattr(body.star, "lethal_event_rate_per_gyr", None)
    if hazard is not None and hazard > LETHAL_EVENT_LIMIT_PER_GYR:
        return True
    index = getattr(body, "flare_irradiation_index", None) or 0.0
    return index * field_shielding(getattr(body, "magnetopause_rp", None), getattr(body, "dipole_class", None)) >= FLARE_OZONE_LIMIT


def gas_density_g_cm3(pressure_pa, molar_mass_kg_mol, temperature_k):
    """The ideal-gas density of the air at the ground, g/cm3 (0 for none)."""
    if not pressure_pa or not molar_mass_kg_mol or not temperature_k:
        return 0.0
    return pressure_pa * molar_mass_kg_mol / (constants.R * temperature_k) / 1000.0


def gamma_dose_msv_yr(rock, age_gy, gas_g_cm3=0.0):
    """
    The effective dose a reference person takes in the open air all year from the gamma rays of `rock` (a `radiogenic.Content` of
    the surface rock) at `age_gy`, mSv/yr: the dose rate in air at 1 m (`radiogenic.gamma_ngy_per_hour`) over 8766 h at
    `SV_PER_GY`, cut by exp(-0.0637 rho h) under a dense gas of `gas_g_cm3` (radiation-model.md section 6).
    """
    rate_ngy_h = radiogenic.gamma_ngy_per_hour(rock, age_gy)
    attenuation = math.exp(-GAMMA_AIR_ATTENUATION_CM2_PER_G * gas_g_cm3 * GAMMA_AIR_HEIGHT_CM)
    return rate_ngy_h * 1e-9 * HOURS_PER_YEAR * SV_PER_GY * 1e3 * attenuation


def radon_flux_bq_m2_s(rock, age_gy, diffusion_factor=1.0):
    """The radon leaving the ground, Bq/(m2 s): J = A(226Ra) f rho_b sqrt(lambda D_e), the radium in equilibrium with the uranium-238.
    `diffusion_factor` is the factor on D_e (1 at one bar; up to P0/P under thin air, R16)."""
    activity = radiogenic.activity_bq_per_kg(rock, age_gy)["U-238"]
    return (activity * RADON_EMANATION * RADON_BULK_DENSITY_KG_M3
            * math.sqrt(math.log(2) / RADON_HALF_LIFE_S * RADON_DIFFUSION_M2_S * diffusion_factor))


def radon_dose_msv_yr(rock, age_gy, scale_height_km, diffusion_factor=1.0):
    """The effective dose of the radon outdoors all year, mSv/yr: the flux spread through a box 1 km tall for the Earth's scale
    height (`RADON_BOX_HEIGHT_KM`, scaled), a mean life of 5.52 days, 9 nSv per Bq h/m3 at 0.6 of equilibrium (R16)."""
    if not scale_height_km or scale_height_km <= 0:
        return 0.0
    height_m = RADON_BOX_HEIGHT_KM * scale_height_km / RADON_REFERENCE_SCALE_HEIGHT_KM * 1000.0
    concentration = radon_flux_bq_m2_s(rock, age_gy, diffusion_factor) * RADON_MEAN_LIFE_S / height_m
    return concentration * HOURS_PER_YEAR * RADON_EQUILIBRIUM_OUTDOORS * RADON_SV_PER_BQ_H_M3 * 1e3


def ground_dose_msv_yr(content, age_gy, rock_ppm, gas_g_cm3=0.0, scale_height_km=None, diffusion_factor=1.0):
    """
    Rule 4: the gamma and the radon dose of a body's surface rock, mSv/yr: the bulk silicate `content` raised to the surface rock
    `rock_ppm`, `(K, Th, U)` ppm today (`radiogenic.surface_content`). Pass no `scale_height_km` for a body with no air (no radon).
    """
    rock = radiogenic.surface_content(content, rock_ppm)
    return (gamma_dose_msv_yr(rock, age_gy, gas_g_cm3)
            + radon_dose_msv_yr(rock, age_gy, scale_height_km, diffusion_factor))


def radon_diffusion_factor(pressure_pa, unit):
    """The factor on the radon diffusion coefficient, log-uniform from 1 to P0/P by `unit` in [0, 1) (R16); 1 at one bar or more."""
    if not pressure_pa or pressure_pa >= constants.STANDARD_ATMOSPHERE_PA:
        return 1.0
    return (constants.STANDARD_ATMOSPHERE_PA / pressure_pa) ** unit


def _clear(body):
    for field in BODY_FIELDS:
        setattr(body, field, None)


def _ground_dose(body, pressure_pa):
    """The body's ground dose from its stored content (`radiogenic.content_of`) and the rock of its stored ground
    (`surface_composition.land_rock_ppm`, GEN.441), mSv/yr; a body with no content (built without it) has no ground dose. Draws the
    radon diffusion factor."""
    diffusion_unit = draw.random()
    content = radiogenic.content_of(body)
    if content is None:
        return 0.0
    has_air = pressure_pa > 0.0 and (getattr(body, "land_fraction", None) or 0.0) > 0.0
    gas = gas_density_g_cm3(pressure_pa, getattr(body, "atm_molar_density", None), getattr(body, "surface_temperature", None))
    age = getattr(body.star, "age", None) or EARTH_AGE_GY
    return ground_dose_msv_yr(content, age, surface_composition.land_rock_ppm(body), gas,
                              getattr(body, "scale_height", None) if has_air else None,
                              radon_diffusion_factor(pressure_pa, diffusion_unit))


def generate_dose(body):
    """Sets `BODY_FIELDS` on a planet or moon (call after its atmosphere,
    water and magnetic field and once `magnetism.update_exposure` has given
    its distance). A gas giant has no surface (every field `None`). Draws
    the solar modulation, the surface rock and the radon diffusion of the ground (`_ground_dose`) and the ozone column once each."""
    _clear(body)
    if body.body_type == "g":
        return
    star = body.star
    distance = getattr(body, "_star_distance_au", None)
    pressure_pa = getattr(body, "atmospheric_pressure", None) or 0.0
    column = column_g_cm2(pressure_pa, getattr(body, "gravity", None))
    moment = getattr(body, "magnetic_moment_a_m2", None) or 0.0
    cutoff = cutoff_rigidity_gv(moment, body.radius)
    modulation = log_uniform(*tuning.GCR_SOLAR_MODULATION_RANGE)
    gcr = gcr_dose_msv_yr(column, cutoff, modulation)
    free = free_sep_msv_yr(getattr(star, "flare_n33_per_yr", None), distance, getattr(star, "flare_alpha", None))
    sep = sep_dose_msv_yr(free, column, cutoff)
    ground = _ground_dose(body, pressure_pa)
    p_o2_pa = mixture.partial_pressures_pa(body).get("O2", 0.0)
    body.ozone_column_du = ozone.draw_column_du(getattr(star, "temperature", None), p_o2_pa, getattr(body, "gravity", None))
    body.dose_gcr_msv_yr = gcr
    body.dose_sep_msv_yr = sep
    body.dose_ground_msv_yr = ground
    body.dose_helio_mult = 1.0
    body.dose_belt_msv_yr = belt_dose_msv_yr(body, column)
    body.surface_dose_msv_yr = gcr + sep + ground + body.dose_belt_msv_yr
    body.ozone_loss_flag = ozone_lost(body)
    if body.ozone_loss_flag:     # GEN.808: the ozone the flares leave is what filters the ultraviolet and tints the sky
        body.ozone_column_du *= OZONE_LOSS_REMAINING
    body.uv_surface_index = body_radiation.dna_dose_index(body)


def belt_dose_msv_yr(body, column):
    """The dose, mSv/yr, a moon takes from its parent's belts through an air of `column` g/cm2: the exposure `apply_parent_belts` stored, cut by the air
    (GEN.983). 0 for a planet, and for a moon until it has the exposure."""
    return (getattr(body, "_belt_exposure_msv_yr", None) or 0.0) * radiation_belts.air_cut(column)


def apply_parent_belts(moon, parent):
    """
    GEN.983: stores on a moon the dose it takes, bare, from the belts of its parent's field (`radiation_belts.moon_exposure_msv_yr`, no draws), and
    makes its stored belt dose and surface dose from it. Call once the parent's field and the moon's orbit and dose are set; `generate_dose` keeps the
    exposure and cuts it by the air as it then is. A moon with no stored dose (a body saved without one) is left as it is.
    """
    moon._belt_exposure_msv_yr = radiation_belts.moon_exposure_msv_yr(parent, moon.distance)
    old = getattr(moon, "dose_belt_msv_yr", None)
    if old is None:
        return
    new = belt_dose_msv_yr(moon, column_g_cm2(getattr(moon, "atmospheric_pressure", None) or 0.0, getattr(moon, "gravity", None)))
    moon.surface_dose_msv_yr += new - old
    moon.dose_belt_msv_yr = new


def update_ground_dose(body):
    """
    GEN.839: the ground dose and the surface dose that holds it, made again from the ground the body now stores (the volcanism gate follows
    the heat at the star's final age); the rest of the dose is kept. Draws the radon diffusion factor again. A body with no stored dose (a
    giant, or one saved before it) is left as it is.
    """
    old = getattr(body, "dose_ground_msv_yr", None)
    if body.body_type == "g" or old is None:
        return
    ground = _ground_dose(body, getattr(body, "atmospheric_pressure", None) or 0.0)
    body.surface_dose_msv_yr += ground - old
    body.dose_ground_msv_yr = ground


def column_g_cm2(pressure_pa, gravity_g):
    """X = P0 / g, g/cm2, from a stored pressure (Pa) and gravity (Earth
    g)."""
    gravity = (gravity_g or 0.0) * constants.EARTH_GRAVITY
    return pressure_pa / gravity * 0.1 if gravity > 0 and pressure_pa else 0.0


def total_dose_msv_yr(gcr, helio_mult, sep, ground, belt=0.0):
    """`surface_dose_msv_yr` from its parts."""
    return gcr * helio_mult + sep + ground + (belt or 0.0)
