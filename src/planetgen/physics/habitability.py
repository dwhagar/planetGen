# planetgen/physics/habitability.py

"""
The planetary habitability index (GEN.84): the reference maths of
docs/design/habitability-index.md section 7, every constant with its source.

Two layers, as Boss chose (2026-10-07 17:11Z): the PHI-4 display -- four
domains (Pressure, Temperature, Chemistry, Radiation), each a score in
[0, 1] and a colour tier (Blue, Green, Yellow, Red) for a human visitor --
and the three numbers behind it from the Xenobiology document: PHI_bio
(microbial life), PHI_cpx (complex, animal life) and Phi_tech (human
operability with equipment). `equipment` names the kit a human needs.

Everything here is a pure function of a `World` (the inputs GEN.85 to
GEN.88 will store for every planet and moon); GEN.89 computes and stores
the scores. Units: pressures kPa, temperatures deg C, doses Sv/yr (Gy/yr
for microbes, taken as equal), column mass g/cm^2, energy flux W/m^2,
phosphate umol/L, inventories in Earth units.
"""

import math
from dataclasses import dataclass, field

# --- Physiology (human visitor) ----------------------------------------------------------

ARMSTRONG_KPA = 6.3
"""float: Water boils at body temperature below this (Armstrong, 1939), and
the saturated water vapour in the lungs takes this much of any breath."""

INSPIRED_O2_MIN_KPA = 8.0
"""float: Least inspired O2 for consciousness (SaO2 above 75%; West, 2012)."""

INSPIRED_O2_CHRONIC_KPA = 50.0
"""float: Most inspired O2 for months without lung damage (Lorrain Smith
effect, 40 to 53 kPa, Lambertsen 1971; 50 as in the Xenobiology doc's
f(pO2))."""

INSPIRED_O2_ACUTE_KPA = 160.0
"""float: Acute CNS oxygen toxicity (1.6 bar; Lambertsen, 1971)."""

MASK_MIN_KPA = ARMSTRONG_KPA + INSPIRED_O2_MIN_KPA
"""float: 14.3 kPa, the least ambient pressure at which a pure-oxygen mask
still gives 8.0 kPa inspired O2 (the inspired-gas equation with F = 1);
below it a pressure suit has to add the rest."""

WET_BULB_CRIT_C = 31.0
"""float: Uncompensable wet-bulb temperature for people (Vecellio et al., 2022)."""

# --- Tier bands ---------------------------------------------------------------------------
#
# Each band is (red_low, yellow_low, green_low, blue_low, blue_high, green_high,
# yellow_high, red_high): Blue between blue_low and blue_high, Green out to the
# green edges, Yellow out to the yellow edges, Red beyond, and the score falls
# to 0 at the red edges. None is an open side. "log" bands are compared in
# log10 (pressures, doses, partial pressures spanning decades).

TIERS = ("Blue", "Green", "Yellow", "Red")

PRESSURE_BAND = (0.63, ARMSTRONG_KPA, MASK_MIN_KPA, 50.0, 250.0, 400.0, 10_000.0, 100_000.0)
"""Total pressure, kPa (log). Blue 50 to 250 (no ebullism, no extra work of
breathing: the PHI-4 doc and Atmospheric Toxicity agree); Green down to the
14.3 kPa pure-O2 mask limit and up to 400 kPa, where nitrogen narcosis
starts (320 to 400 kPa, Bennett & Rostain 2003); Yellow down to the
Armstrong limit (counterpressure suit) and up to 10 MPa (the PHI-4 doc's
hyperbaric edge, past the 7 MPa COMEX Hydra X record); Red beyond, falling to
0 a decade further out."""

TEMPERATURE_BAND = (-120.0, -50.0, -20.0, 0.0, 31.0, 45.0, 122.0, 200.0)
"""Dry-bulb temperature, deg C (linear). Blue 0 to 31, Green -20 to 45
(ordinary clothing or cooling garments), Yellow -50 (eutectic brines) to
122 (the limit of carbon-based life, Takai et al. 2008; the PHI-4 doc's Red
edge), Red beyond. A wet bulb at or above 31 deg C makes Blue or Green
Yellow (see `temperature_score`)."""

DOSE_BAND = (None, None, None, None, 0.05, 0.1, 10.0, 1000.0)
"""Surface dose, Sv/yr (log). Blue below 50 mSv/yr (occupational limit),
Green to 100 mSv/yr, Yellow to 10 Sv/yr (shelter underground; the PHI-4
doc's Red edge and the Eigen threshold the Xenobiology doc puts at 10
Gy/yr), Red beyond, 0 at 1000 Sv/yr."""

# Gas limits, kPa partial pressure: (chronic, acute), from the Atmospheric
# Toxicity table (Lambertsen 1971; Schwieterman et al. 2019).
GAS_LIMITS_KPA = {
    "CO2": (0.93, 5.0),
    "CO": (0.005, 0.1),
    "H2S": (0.001, 0.05),
    "SO2": (0.0005, 0.01),
}

FILTER_FACTOR = 10.0
"""float: A half-face filtering respirator's assigned protection factor
(NIOSH APF 10): a gas a filter can take out (CO with Hopcalite, H2S, SO2)
is Yellow up to ten times its acute limit, Red above."""

CO2_GREEN_KPA = 2.0
"""float: CO2 is Green up to 2.0 kPa, where progressive hypercapnic acidosis
starts (1.0 to 2.0 kPa, Atmospheric Toxicity; the Xenobiology doc's
pCO2_tox), Yellow to its 5.0 kPa acute limit (closed-loop amine scrubbers;
it can't be filtered), Red above."""

WATER_ACTIVITY_CRIT = 0.605
"""float: Lowest water activity for cell division (Stevenson et al., 2015)."""

WATER_ACTIVITY_GREEN = 0.75
"""float: Halite saturation (aw 0.75): water below it is a brine."""

CHAOTROPICITY_CRIT_KJ_KG = 73.8
"""float: Chaotropicity that dismantles hydration shells (Ball & Hallsworth, 2015)."""

# --- PHI_bio, PHI_cpx -------------------------------------------------------------------

K_W = 50.0
"""float: Steepness of f(aw) round 0.605: f is 0.12 at 0.565 and 0.88 at
0.645, the spread of measured growth limits (Stevenson et al., 2015). The
doc leaves k_w open; planetGen's choice."""

SIGMA_CHI_KJ_KG = 10.0
"""float: Chaotropicity fall-off past 73.8 kJ/kg (the Xenobiology doc)."""

SOLVENT_FACTOR = {"water": 1.0, "hydrocarbon": 0.2, "sulfuric": 0.1, None: 0.0}
"""dict: L_solv's factor by solvent. Water 1; liquid methane/ethane 0.2 and
concentrated sulfuric acid 0.1 (possible but unproven solvents, Bains et
al. 2024, Petkowski et al. 2020: planetGen's choice); no stable liquid 0,
a hard gate."""

X_REF_EARTH = 0.25
"""float: X_j,ref for carbon, nitrogen and hydrogen, in Earth inventories:
a quarter of Earth's gives tanh(4) = 0.999, a tenth of Earth's 0.38.
Planetgen's choice (the doc leaves it open)."""

K_P_UMOL_L = 0.1
"""float: Phosphate half-saturation, umol/L: the order of microbial and
phytoplankton phosphate uptake (0.01 to 0.1 umol/L). The doc's 1.0 umol/L
(Toner & Catling 2020) would score Earth's 2.3 umol/L ocean only 0.70."""

ENERGY_FLUX_MIN_W_M2 = 1.0e-6
ENERGY_FLUX_REF_W_M2 = 10.0
"""float: L_ener = log10(Phi / min) / log10(ref / min), clamped to [0, 1]:
0 below a microwatt per square metre, 1 from 10 W/m^2 (Earth's surface has
about 80 W/m^2 of photosynthetic light; a radiolytic subsurface 1e-4 to
1e-3). Log-scaled because the sources span ten decades; the doc's
tanh(Phi / Phi_ref) has no Phi_ref that scores both."""

DOSE_THRESH_GY = 10.0
SIGMA_DOSE_GY = 100.0
"""float: L_rad's threshold and fall-off, Gy/yr (Daly, 2012; the doc)."""

CO_TOX_KPA = 0.01
CO2_TOX_KPA = 2.0
"""float: T_metazoa's CO and CO2 scales, kPa (the Xenobiology doc)."""

METAZOA_O2_LOW_KPA = 8.0
METAZOA_O2_HIGH_KPA = 50.0
METAZOA_O2_SIGMA_LOW_KPA = 2.0
METAZOA_O2_SIGMA_HIGH_KPA = 20.0
"""float: f(pO2) is 1 between 8 and 50 kPa of ambient O2 and falls off as
a half Gaussian outside, sigma 2 kPa below and 20 kPa above (planetGen's
choice; the doc names only the band)."""

# --- Phi_tech ----------------------------------------------------------------------------

SIGMA_P_KPA = 2500.0
"""float: M_press fall-off above 250 kPa: 0.50 at 2 MPa, 0.02 at 10 MPa
(the PHI-4 Yellow/Red edge). Planetgen's choice."""

K_T_PER_K = 0.5
"""float: M_therm's wet-bulb steepness (Vecellio et al., 2022; the doc)."""

COLD_LIMIT_C = -50.0
K_COLD_PER_K = 0.2
COLD_FLOOR = 0.3
"""float: M_therm's cold side, which the doc leaves out: a logistic at -50
deg C (the PHI-4 Yellow edge), steepness 0.2 per K, floored at 0.3 --
heating a sealed habitat is cheap next to cooling one, so no cold makes a
base impossible (planetGen's choice)."""

VACUUM_PRESS_FACTOR = 0.2
"""float: M_press below the Armstrong limit: every base is a sealed
pressure vessel, as on the Moon (planetGen's choice; the doc's 0.2)."""

WATER_SOURCE_FACTOR = {"liquid": 1.0, "ice": 0.8, "hydrated": 0.6, "vapour": 0.4, None: 0.1}
"""dict: M_isru by the most accessible water: liquid at the surface, ice,
hydrated minerals (perchlorate brines counted here), atmospheric vapour
only, or none (imported). Replaces the doc's undefined tanh((H2O + alpha
ClO4) / E_specific)."""

SHIELD_EXPONENT = 0.2
"""float: M_rad = (0.05 / D)^0.2 above 50 mSv/yr: the shielding a base must
add grows with log(D); Mars's 0.25 Sv/yr scores 0.73 and Europa's surface
(about 2000 Sv/yr) 0.12. Replaces the doc's undefined lambda."""

# --- Surface dose from the atmosphere ---------------------------------------------------

GCR_SOFT_SV_YR = 0.45
GCR_SOFT_LENGTH_G_CM2 = 22.0
GCR_HARD_SV_YR = 0.05
GCR_HARD_LENGTH_G_CM2 = 200.0
GROUND_SV_YR = 0.002
"""float: Galactic cosmic-ray surface dose, D = 0.45 exp(-X/22) + 0.05
exp(-X/200) Sv/yr plus 2 mSv/yr from the ground, X the column mass P0/g.
Fitted to the Moon's surface, 0.50 Sv/yr (Zhang et al., 2020), Mars's,
about 0.25 Sv/yr at 16 g/cm^2 (Hassler et al., 2014), 50 to 100 mSv/yr at
0.05 bar (the radiation doc) and Earth's 2.4 mSv/yr at 1033 g/cm^2 (Atri,
2017). Stellar particle events, magnetic fields and a nebula's compressed
heliosphere add to it in GEN.87."""


def column_mass_g_cm2(pressure_kpa, gravity_ms2):
    """X = P0 / g, g/cm^2."""
    if gravity_ms2 <= 0:
        return 0.0
    return pressure_kpa * 1000.0 / gravity_ms2 / 10.0


def gcr_surface_dose_sv_yr(column_g_cm2):
    """The cosmic-ray dose under `column_g_cm2` of air, plus the ground's."""
    return (GCR_SOFT_SV_YR * math.exp(-column_g_cm2 / GCR_SOFT_LENGTH_G_CM2)
            + GCR_HARD_SV_YR * math.exp(-column_g_cm2 / GCR_HARD_LENGTH_G_CM2) + GROUND_SV_YR)


# --- A world -----------------------------------------------------------------------------

@dataclass
class World:
    """What the index reads. `gases` maps species to mole fraction; `solvent`
    is the stable surface or subsurface liquid ("water", "hydrocarbon",
    "sulfuric" or None); `habitat_dose_sv_yr` is the dose where that liquid
    is (the surface dose for surface water, about the ground's for an ocean
    under ice); `water_source` is what a base can draw water from."""
    pressure_kpa: float
    temperature_c: float
    gases: dict
    gravity_ms2: float = 9.81
    relative_humidity: float = 0.5
    solvent: str = "water"
    water_activity: float = 1.0
    chaotropicity_kj_kg: float = 0.0
    ph: float = None
    inventories: dict = field(default_factory=lambda: {"C": 1.0, "N": 1.0, "H": 1.0})
    phosphate_umol_l: float = 2.3
    energy_flux_w_m2: float = 80.0
    surface_dose_sv_yr: float = None
    habitat_dose_sv_yr: float = None
    water_source: str = "liquid"

    def partial_kpa(self, gas):
        return self.gases.get(gas, 0.0) * self.pressure_kpa

    @property
    def dose_sv_yr(self):
        if self.surface_dose_sv_yr is not None:
            return self.surface_dose_sv_yr
        return gcr_surface_dose_sv_yr(column_mass_g_cm2(self.pressure_kpa, self.gravity_ms2))


# --- PHI-4 display ---------------------------------------------------------------------------

_EDGE_SCORES = (0.0, 0.4, 0.75, 1.0)
"""The score at the outer red, yellow, green and blue edges: Blue scores 1,
Green 0.75 to 1, Yellow 0.4 to 0.75, Red below 0.4 down to 0 at the outer
red edge, linear in between."""


def band_score(x, band, log=False):
    """The PHI-4 score of `x` in `band` (see the band constants): 1 inside
    Blue, falling linearly (in log10 for a log band) through 0.75 at the
    green edge and 0.4 at the yellow edge to 0 at the outer red edge."""
    t = (lambda v: math.log10(v) if v and v > 0 else -math.inf) if log else (lambda v: v)
    low = [band[0], band[1], band[2], band[3]]
    high = [band[7], band[6], band[5], band[4]]
    value = t(x)
    for edges, sign in ((low, 1.0), (high, -1.0)):
        points = [(t(edge) * sign, score) for edge, score in zip(edges, _EDGE_SCORES) if edge is not None]
        if not points:
            continue
        v = value * sign
        if v >= points[-1][0]:
            continue
        if v <= points[0][0]:
            return points[0][1]
        for (x0, s0), (x1, s1) in zip(points, points[1:]):
            if x0 <= v <= x1:
                return s0 + (s1 - s0) * (v - x0) / (x1 - x0)
    return 1.0


def tier(score):
    """The colour of a domain score."""
    if score >= 1.0:
        return "Blue"
    if score >= 0.75:
        return "Green"
    if score >= 0.4:
        return "Yellow"
    return "Red"


def wet_bulb_c(temperature_c, relative_humidity):
    """Wet-bulb temperature, deg C (Stull, 2011; fitted for -20 to 50 deg C
    and 5 to 99% humidity, clamped to those)."""
    t = temperature_c
    rh = min(max(relative_humidity, 0.05), 0.99) * 100.0
    return (t * math.atan(0.151977 * math.sqrt(rh + 8.313659)) + math.atan(t + rh) - math.atan(rh - 1.676331)
            + 0.00391838 * rh ** 1.5 * math.atan(0.023101 * rh) - 4.686035)


def pressure_score(world):
    return band_score(world.pressure_kpa, PRESSURE_BAND, log=True)


def temperature_score(world):
    score = band_score(world.temperature_c, TEMPERATURE_BAND)
    if -20.0 <= world.temperature_c <= 50.0 and wet_bulb_c(world.temperature_c, world.relative_humidity) >= WET_BULB_CRIT_C:
        score = min(score, 0.75 - 1e-9)  # an uncompensable wet bulb needs cooling garments: Yellow
    return score


def inspired_o2_kpa(world):
    """Inspired O2 (the inspired-gas equation): F_O2 (P - 6.3 kPa)."""
    return world.gases.get("O2", 0.0) * max(0.0, world.pressure_kpa - ARMSTRONG_KPA)


def _gas_scores(world):
    scores = {}
    o2 = inspired_o2_kpa(world)
    # Low O2: a mask fixes it wherever a mask works (Green); where it doesn't,
    # the pressure domain already says so.
    o2_low = 1.0 if o2 >= INSPIRED_O2_MIN_KPA else 0.75
    o2_high = band_score(o2, (None, None, None, None, INSPIRED_O2_CHRONIC_KPA, INSPIRED_O2_CHRONIC_KPA,
                              INSPIRED_O2_ACUTE_KPA, INSPIRED_O2_ACUTE_KPA * 10), log=True)
    scores["O2"] = min(o2_low, o2_high)
    chronic, acute = GAS_LIMITS_KPA["CO2"]
    scores["CO2"] = band_score(world.partial_kpa("CO2"),
                               (None, None, None, None, chronic, CO2_GREEN_KPA, acute, acute * 10), log=True)
    for gas in ("CO", "H2S", "SO2"):
        chronic, acute = GAS_LIMITS_KPA[gas]
        scores[gas] = band_score(world.partial_kpa(gas), (None, None, None, None, chronic, chronic,
                                                          acute * FILTER_FACTOR, acute * FILTER_FACTOR * 10),
                                 log=True)
    return scores


def _water_score(world):
    if world.solvent != "water" or world.ph is None:
        return 1.0
    ph = band_score(world.ph, (0.0, 1.0, 5.0, 6.0, 8.5, 9.5, 11.5, 14.0))
    aw = band_score(world.water_activity, (0.5, WATER_ACTIVITY_CRIT, WATER_ACTIVITY_GREEN, 0.90, 1.0, 1.0, 1.0, 1.0))
    chi = band_score(world.chaotropicity_kj_kg, (None, None, None, None, CHAOTROPICITY_CRIT_KJ_KG,
                                                 CHAOTROPICITY_CRIT_KJ_KG, CHAOTROPICITY_CRIT_KJ_KG,
                                                 CHAOTROPICITY_CRIT_KJ_KG * 2))
    return min(ph, aw, chi)


def chemistry_score(world):
    return min(min(_gas_scores(world).values()), _water_score(world))


def radiation_score(world):
    return band_score(world.dose_sv_yr, DOSE_BAND, log=True)


def phi4(world):
    """`{domain: (score, tier)}` and `"PHI-4"`: the geometric mean of the four
    (equal weights; any domain at 0 makes it 0, the doc's gates)."""
    domains = {"Pressure": pressure_score(world), "Temperature": temperature_score(world),
               "Chemistry": chemistry_score(world), "Radiation": radiation_score(world)}
    result = {name: (round(score, 4), tier(score)) for name, score in domains.items()}
    product = math.prod(domains.values())
    result["PHI-4"] = round(product ** 0.25, 4) if product > 0 else 0.0
    return result


def equipment(world):
    """The kit a human visitor needs: "ideal", "breathing mask", "mask
    with scrubber", "sealed suit" (a pressure suit, or a sealed suit
    against a Red temperature or chemistry) or "full life support with
    radiation hardening", from the worst of the four domains and what made
    it so."""
    scores = {"Pressure": pressure_score(world), "Temperature": temperature_score(world),
              "Chemistry": chemistry_score(world), "Radiation": radiation_score(world)}
    gases = _gas_scores(world)
    if scores["Radiation"] < 0.4:
        return "full life support with radiation hardening"
    if scores["Pressure"] < 0.75 or min(scores.values()) < 0.4:
        return "sealed suit"
    if min(gases[g] for g in ("CO2", "CO", "H2S", "SO2")) < 1.0:
        return "mask with scrubber"
    if gases["O2"] < 1.0 or scores["Pressure"] < 1.0:
        return "breathing mask"
    return "ideal"


# --- PHI_bio, PHI_cpx, Phi_tech ------------------------------------------------------------

def l_solv(world):
    factor = SOLVENT_FACTOR.get(world.solvent, 0.0)
    if factor == 0.0:
        return 0.0
    f_aw = 1.0 / (1.0 + math.exp(-K_W * (world.water_activity - WATER_ACTIVITY_CRIT))) \
        if world.solvent == "water" else 1.0
    excess = max(0.0, world.chaotropicity_kj_kg - CHAOTROPICITY_CRIT_KJ_KG)
    return factor * f_aw * math.exp(-0.5 * (excess / SIGMA_CHI_KJ_KG) ** 2)


def l_chem(world):
    cnh = math.prod(math.tanh(world.inventories.get(j, 0.0) / X_REF_EARTH) for j in ("C", "N", "H"))
    return cnh ** (1.0 / 3.0) * world.phosphate_umol_l / (world.phosphate_umol_l + K_P_UMOL_L)


def l_ener(world):
    if world.energy_flux_w_m2 <= ENERGY_FLUX_MIN_W_M2:
        return 0.0
    span = math.log10(ENERGY_FLUX_REF_W_M2 / ENERGY_FLUX_MIN_W_M2)
    return min(1.0, math.log10(world.energy_flux_w_m2 / ENERGY_FLUX_MIN_W_M2) / span)


def l_rad(world):
    dose = world.habitat_dose_sv_yr if world.habitat_dose_sv_yr is not None else world.dose_sv_yr
    excess = max(0.0, dose - DOSE_THRESH_GY)
    return math.exp(-0.5 * (excess / SIGMA_DOSE_GY) ** 2)


def phi_bio(world):
    """The product of the four likelihoods (not the doc's geometric mean,
    which forgives a missing solvent: each is a requirement, and all must
    hold)."""
    return l_solv(world) * l_chem(world) * l_ener(world) * l_rad(world)


def t_metazoa(world):
    o2 = world.partial_kpa("O2")
    if o2 < METAZOA_O2_LOW_KPA:
        f_o2 = math.exp(-0.5 * ((METAZOA_O2_LOW_KPA - o2) / METAZOA_O2_SIGMA_LOW_KPA) ** 2)
    elif o2 > METAZOA_O2_HIGH_KPA:
        f_o2 = math.exp(-0.5 * ((o2 - METAZOA_O2_HIGH_KPA) / METAZOA_O2_SIGMA_HIGH_KPA) ** 2)
    else:
        f_o2 = 1.0
    return (math.exp(-(world.partial_kpa("CO") / CO_TOX_KPA) ** 2)
            * math.exp(-(world.partial_kpa("CO2") / CO2_TOX_KPA) ** 2) * f_o2)


def phi_cpx(world):
    return phi_bio(world) * t_metazoa(world)


def m_press(world):
    p = world.pressure_kpa
    if p < ARMSTRONG_KPA:
        return VACUUM_PRESS_FACTOR
    if p <= 250.0:
        return 1.0
    return math.exp(-(p - 250.0) / SIGMA_P_KPA)


def _logistic(x):
    """1 / (1 + exp(-x)) without overflowing at a very hot or cold extreme."""
    if x >= 0.0:
        return 1.0 / (1.0 + math.exp(-x))
    e = math.exp(x)
    return e / (1.0 + e)


def m_therm(world):
    hot = 1.0
    if world.temperature_c >= -20.0:
        twb = wet_bulb_c(min(world.temperature_c, 50.0), world.relative_humidity)
        if world.temperature_c > 50.0:
            twb += world.temperature_c - 50.0  # past Stull's fit, the wet bulb tracks the dry
        hot = _logistic(-K_T_PER_K * (twb - WET_BULB_CRIT_C))
    cold = _logistic(K_COLD_PER_K * (world.temperature_c - COLD_LIMIT_C))
    return hot * (COLD_FLOOR + (1.0 - COLD_FLOOR) * cold)


def m_isru(world):
    return WATER_SOURCE_FACTOR.get(world.water_source, WATER_SOURCE_FACTOR[None])


def m_rad(world):
    dose = world.dose_sv_yr
    return 1.0 if dose <= 0.05 else (0.05 / dose) ** SHIELD_EXPONENT


def phi_tech(world):
    """The geometric mean of the four mitigation factors (the doc's form:
    equipment can make up for one hard factor in part)."""
    return (m_press(world) * m_therm(world) * m_isru(world) * m_rad(world)) ** 0.25


def scores(world):
    """Every number the index gives a world, rounded for display."""
    out = {key: round(fn(world), 3) for key, fn in (
        ("L_solv", l_solv), ("L_chem", l_chem), ("L_ener", l_ener), ("L_rad", l_rad),
        ("PHI_bio", phi_bio), ("PHI_cpx", phi_cpx), ("Phi_tech", phi_tech))}
    out.update(phi4(world))
    out["equipment"] = equipment(world)
    return out
