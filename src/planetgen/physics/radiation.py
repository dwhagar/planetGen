# planetgen/physics/radiation.py

"""
Surface radiation dose, UV and galactic hazards (GEN.87).

The rules are docs/design/habitability/activity-magnetism-radiation-hydrosphere.md
section 4 (design defaults). A rocky body's yearly dose, mSv/yr, is

    surface = gcr * helio_mult + sep + ground

1. **Cosmic rays (`dose_gcr_msv_yr`, GEN.453).** `D = D0(X) f_B(M, X) f_sol` under a column of air `X = P0 / g`
   (g/cm2). `D0` is the no-field curve: the proton-primary air showers of Atri et al. 2013 (553, 263, 17.2, 3.02 and 0.65
   mSv/yr at X of 100, 200, 500, 700 and 1036), log-log between them; below 100 g/cm2 a stand-in through the Moon (X 0, 500
   mSv/yr, Chang'E-4) and Mars (X 21, 232 mSv/yr, Hassler 2014) until a shower run replaces it; above 1036 the larger of the
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
3. **Stellar particles (`dose_sep_msv_yr`).** A free-space `20 mSv/yr
   sqrt(N33 / N33_sun) / d_AU^2`, times `exp(-p_cut(X) / 0.2 GV)` (the
   rigidity of the proton whose range is X, from NIST PSTAR, `proton_range`) and the
   polar-cap share `1 - sqrt(1 - sqrt(R_eff / Rc0))` with `R_eff = max(0.5
   p_cut, 0.2 GV)` (1 without a field).
4. **Ground (`dose_ground_msv_yr`).** The gamma dose in air at 1 m over the surface rock and the radon of its uranium, from the
   body's own stored potassium, thorium and uranium (`radiogenic`, GEN.459), at its star's age: the same concentrations its interior
   heat reads, so there is no independent dose draw (`a_rad` is gone). The surface rock is a basaltic crust with a drawn felsic share
   (`radiogenic.SURFACE_RESERVOIRS`, a stand-in for the ground type of GEN.441); the dose is the open-air, full-time effective dose
   of a reference person (radiation-model.md R12). Radon has no pressure gate (R16): its flux is drawn between 1 and
   `sqrt(P0 / P)` times the 1-bar value.
5. **UV (`uv_surface_index`).** The Setlow DNA-damage dose at the ground as a multiple of Earth's (Earth is 1 by definition, GEN.87), from
   the two-stream spectrum below the body's own air, ozone column (`ozone_column_du`, GEN.413: the Kozakis fit by oxygen and host star,
   `physics/ozone.py`), haze, dust and cloud, lit by the star's own spectrum with its stored UV-C share (`physics/body_radiation.py`, GEN.443); a
   body with no air reads the star's unfiltered light (GEN.445; docs/design/composition/radiation-model.md section 7). There is no proxy: no
   ratio to the Sun, no ozone exponent.
6. **Ozone loss (`ozone_loss_flag`).** True when there is ozone and either
   the system's lethal galactic events exceed one per 100 Myr or particle
   events deliver 100 mSv/yr at the ozone layer (X = 10 g/cm2).
7. **Galactic hazard (`lethal_event_rate_per_gyr`, on the star).** A
   supernova within 8 pc, 1.5 per Gyr at 8 kpc from the centre, times
   `exp((8 kpc - R) / 1.45 kpc)` (250x at the centre), plus gamma-ray bursts
   (2 per Gyr at the Sun, falling outward) and nuclear-activity episodes (10
   per Gyr inside 2.9 kpc), GEN.189.

A planet of a pulsar or a black hole is rated by GEN.89, not here.
"""

import math

from planetgen import tuning
from planetgen.physics import body_radiation, constants, magnetism, ozone, proton_range, radiogenic
from planetgen.util import draw
from planetgen.util.random import log_uniform

BODY_FIELDS = (
    "surface_dose_msv_yr", "dose_gcr_msv_yr", "dose_sep_msv_yr", "dose_ground_msv_yr", "dose_helio_mult",
    "uv_surface_index", "ozone_loss_flag", "ozone_column_du",
)
"""tuple: What `generate_dose` sets on a planet or moon, as stored in
`planets` and `moons`."""

GCR_ANCHORS = ((0.0, 500.0), (21.0, 232.0))
"""tuple: `(X g/cm2, mSv/yr)`, the no-field dose below 100 g/cm2, a stand-in until a shower run (AD5) replaces it: the Moon
(Chang'E-4 LND, 1369 uSv/day, a solar-minimum bound; Zhang 2020) and Mars at Gale (MSL RAD, 0.64 mSv/day, GCR only, column
about 21 g/cm2; Hassler 2014, Table 3). The shape between them and from Mars to 100 g/cm2 is unknown (the Pfotzer
maximum is near 100 to 120 g/cm2); Guo 2015 gives 255 mSv/yr at Mars with a range of 128 to 420."""
GCR_SHOWER_POINTS = ((100.0, 553.0), (200.0, 263.0), (500.0, 17.2), (700.0, 3.02), (1036.0, 0.65))
"""tuple: `(X g/cm2, mSv/yr)`: total biological dose from muons, electrons, photons and neutrons in proton-primary CORSIKA
air showers with no magnetic field (Atri, Hariharan, Fryer et al. 2013, Table XV, via Griessmeier et al. 2016, Table 2).
Helium and heavier primaries are absent, so these are lower bounds."""
GCR_SEA_LEVEL_G_CM2 = GCR_SHOWER_POINTS[-1][0]

FIELD_COLUMNS_G_CM2 = (100.0, 200.0, 500.0, 700.0, 1036.0)
"""tuple: The columns the field's cut was run at (Atri 2013)."""
FIELD_FULL_MOMENTS = (0.0, 0.1, 0.15, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 6.0, 10.0)
"""tuple: The dipole moments, Earth's = 1, of the runs at 100 and 1036 g/cm2 (Griessmeier 2016, Table 2)."""
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
pi I times 2 MeV cm2/g. It is 9 times below Sittler 2019 at Titan (10 850 g/cm2), so it is a floor under the table's tail."""

EARTH_DIPOLE_A_M2 = 7.768e22
HELIO_MAX_EXTRA = 1.5
HELIO_FADE_G_CM2 = 100.0

SEP_FREE_MSV_YR = 20.0
SEP_R0_GV = 0.2
SEP_MIN_RIGIDITY_GV = 0.2
N33_SUN_PER_YR = 10 ** -2.5

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

OZONE_LAYER_G_CM2 = 10.0
OZONE_SEP_MSV_YR = 100.0

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


def no_field_dose_msv_yr(column_g_cm2):
    """Rule 1's `D0`: the dose, mSv/yr, under `column_g_cm2` with no magnetic field and the mean solar modulation."""
    x = max(column_g_cm2, 0.0)
    points = GCR_SHOWER_POINTS
    if x < points[0][0]:
        nodes = GCR_ANCHORS + points[:1]
        columns = [math.log1p(column) for column, _ in nodes]
        return _log_linear(math.log1p(x), columns, [dose for _, dose in nodes])
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
    """Rule 3's stellar-particle dose under `column_g_cm2` for a free-space
    dose and a dipole's `cutoff_gv`."""
    p_cut = proton_cutoff_gv(column_g_cm2)
    r_eff = max(0.5 * p_cut, SEP_MIN_RIGIDITY_GV)
    polar = 1.0 if cutoff_gv <= r_eff else 1.0 - math.sqrt(1.0 - math.sqrt(r_eff / cutoff_gv))
    return free_msv_yr * math.exp(-p_cut / SEP_R0_GV) * polar


def free_sep_msv_yr(flare_n33_per_yr, distance_au):
    """Rule 3's free-space dose at `distance_au`."""
    if not flare_n33_per_yr or not distance_au:
        return 0.0
    return SEP_FREE_MSV_YR * math.sqrt(flare_n33_per_yr / N33_SUN_PER_YR) / distance_au ** 2


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


def ground_dose_msv_yr(content, age_gy, felsic_fraction, gas_g_cm3=0.0, scale_height_km=None, diffusion_factor=1.0):
    """
    Rule 4: the gamma and the radon dose of a body's surface rock, mSv/yr: the bulk silicate `content` raised to the surface rock
    (`radiogenic.surface_content`, `felsic_fraction`). Pass no `scale_height_km` for a body with no air (no radon).
    """
    rock = radiogenic.surface_content(content, felsic_fraction)
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
    """The body's ground dose from its stored content (`radiogenic.content_of`), mSv/yr; a body with none (built without it) has no
    ground dose. Draws, in this order: whether the lid is mobile, the felsic share, the radon diffusion factor."""
    mobile = draw.random()
    felsic_unit = draw.random()
    diffusion_unit = draw.random()
    content = radiogenic.content_of(body)
    if content is None:
        return 0.0
    mass_earth = body.mass / constants.EARTH_MASS_TO_KG
    felsic = (tuning.RADIOGENIC_FELSIC_FRACTION_MAX * felsic_unit
              if mobile < magnetism.mobile_lid_chance(body, mass_earth) else 0.0)
    has_air = pressure_pa > 0.0 and (getattr(body, "land_fraction", None) or 0.0) > 0.0
    gas = gas_density_g_cm3(pressure_pa, getattr(body, "atm_molar_density", None), getattr(body, "surface_temperature", None))
    age = getattr(body.star, "age", None) or EARTH_AGE_GY
    return ground_dose_msv_yr(content, age, felsic, gas, getattr(body, "scale_height", None) if has_air else None,
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
    moment = (getattr(body, "magnetic_moment_a_m2", None) or 0.0) * magnetism.body_dipole_share(body)
    cutoff = cutoff_rigidity_gv(moment, body.radius)
    modulation = log_uniform(*tuning.GCR_SOLAR_MODULATION_RANGE)
    gcr = gcr_dose_msv_yr(column, cutoff, modulation)
    free = free_sep_msv_yr(getattr(star, "flare_n33_per_yr", None), distance)
    sep = sep_dose_msv_yr(free, column, cutoff)
    ground = _ground_dose(body, pressure_pa)
    p_o2_pa = (getattr(body, "p_o2_kpa", None) or 0.0) * 1000.0
    body.ozone_column_du = ozone.draw_column_du(getattr(star, "temperature", None), p_o2_pa, getattr(body, "gravity", None))
    body.dose_gcr_msv_yr = gcr
    body.dose_sep_msv_yr = sep
    body.dose_ground_msv_yr = ground
    body.dose_helio_mult = 1.0
    body.surface_dose_msv_yr = gcr + sep + ground
    body.uv_surface_index = body_radiation.dna_dose_index(body)
    hazard = getattr(star, "lethal_event_rate_per_gyr", None)
    ozone_sep = sep_dose_msv_yr(free, OZONE_LAYER_G_CM2, cutoff)
    body.ozone_loss_flag = bool(body.ozone_column_du > 0.0 and (
        (hazard is not None and hazard > LETHAL_EVENT_LIMIT_PER_GYR) or ozone_sep >= OZONE_SEP_MSV_YR))


def column_g_cm2(pressure_pa, gravity_g):
    """X = P0 / g, g/cm2, from a stored pressure (Pa) and gravity (Earth
    g)."""
    gravity = (gravity_g or 0.0) * constants.EARTH_GRAVITY
    return pressure_pa / gravity * 0.1 if gravity > 0 and pressure_pa else 0.0


def total_dose_msv_yr(gcr, helio_mult, sep, ground):
    """`surface_dose_msv_yr` from its parts."""
    return gcr * helio_mult + sep + ground
