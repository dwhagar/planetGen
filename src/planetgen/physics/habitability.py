# planetgen/physics/habitability.py

"""
The planetary habitability index (GEN.84): the reference maths of
docs/design/habitability/habitability-index.md section 7, every constant with its source.

Two layers, as Boss chose (2026-10-07 17:11Z): the PHI-6 display -- since
GEN.207 six domains (Pressure, Temperature, Air, Water, Radiation and
Suitability, the long-term one), each a score in [0, 1] and a colour tier
(Blue, Green, Yellow, Red) for a human visitor, and an overall score that
mixes the colour gates with their geometric mean (`phi4`) -- and the two
numbers behind the life scores from the Xenobiology document: PHI_bio
(microbial life) and PHI_cpx (complex, animal life), which `life_scores`
turns into the Microbial, Flora, Fauna and Sentient scores (GEN.649; the
document's Phi_tech, human operability, is gone). The kit a human needs is worked out in `physics/equipment.py`.

Everything here is a pure function of a `World` (the inputs GEN.85 to
GEN.88 will store for every planet and moon); GEN.89 computes and stores
the scores. Units: pressures kPa, temperatures deg C, doses Sv/yr (Gy/yr
for microbes, taken as equal), column mass g/cm^2, energy flux W/m^2,
phosphate umol/L, inventories in Earth units.
"""

import math
from collections import namedtuple
from dataclasses import dataclass, field

from planetgen.physics import constants, radiation
from planetgen.physics.protection_limits import limit

# --- Physiology (human visitor) ----------------------------------------------------------

ARMSTRONG_KPA = limit("armstrong_kpa")
"""float: Water boils at body temperature below this (Armstrong, 1939), and
the saturated water vapour in the lungs takes this much of any breath."""

INSPIRED_O2_MIN_KPA = limit("inspired_o2_min_kpa")
"""float: Least inspired O2 for consciousness (SaO2 above 75%; West, 2012)."""

INSPIRED_O2_CHRONIC_KPA = limit("inspired_o2_chronic_kpa")
"""float: Most inspired O2 for months without lung damage (Lorrain Smith
effect, 40 to 53 kPa, Lambertsen 1971; 50 as in the Xenobiology doc's
f(pO2))."""

INSPIRED_O2_ACUTE_KPA = limit("inspired_o2_acute_kpa")
"""float: Acute CNS oxygen toxicity (1.6 bar; Lambertsen, 1971)."""

MASK_MIN_KPA = ARMSTRONG_KPA + INSPIRED_O2_MIN_KPA
"""float: 14.3 kPa, the least ambient pressure at which a pure-oxygen mask
still gives 8.0 kPa inspired O2 (the inspired-gas equation with F = 1);
below it a pressure suit has to add the rest."""

WET_BULB_CRIT_C = limit("wet_bulb_c")
"""float: Uncompensable wet-bulb temperature for people (Vecellio et al., 2022)."""

# --- Tier bands ---------------------------------------------------------------------------
#
# Each band is (red_low, yellow_low, green_low, blue_low, blue_high, green_high,
# yellow_high, red_high): Blue between blue_low and blue_high, Green out to the
# green edges, Yellow out to the yellow edges, Red beyond, and the score falls
# to 0 at the red edges. None is an open side. "log" bands are compared in
# log10 (pressures, doses, partial pressures spanning decades).

TIERS = ("Blue", "Green", "Yellow", "Red")

POTENTIALLY_HABITABLE_TIERS = ("Blue", "Green")
"""tuple: The overall PHI-6 colours of a habitable world (Boss, 2026-10-23 6:02 PM EDT: "Habitable should now be based on
the PHI-6 score which must be Blue or Green overall"; the name predates it, GEN.220). `is_habitable` is the one place the
word is defined."""

POTENTIALLY_HABITABLE_MIN_SCORE = 0.75
"""float: The lowest overall PHI-6 score that is still Green (`tier`), so the lowest that counts as potentially habitable
(Boss, 2026-10-23: "habitable" is an overall Blue or Green); the stored `phi4` of a planet or moon is compared with it."""

PRESSURE_BAND = (limit("pressure_low_outer_kpa"), ARMSTRONG_KPA, limit("mask_min_kpa"), limit("pressure_no_gear_kpa"),
                 limit("pressure_blue_max_kpa"), limit("pressure_narcosis_kpa"), limit("pressure_hyperbaric_max_kpa"),
                 limit("pressure_outer_kpa"))
"""Total pressure, kPa (log). Blue 50 to 250 (no ebullism, no extra work of
breathing: the PHI-4 doc and Atmospheric Toxicity agree); Green down to the
14.3 kPa pure-O2 mask limit and up to 400 kPa, where nitrogen narcosis
starts (320 to 400 kPa, Bennett & Rostain 2003); Yellow down to the
Armstrong limit (counterpressure suit) and up to 10 MPa (the PHI-4 doc's
hyperbaric edge, past the 7 MPa COMEX Hydra X record); Red beyond, falling to
0 a decade further out."""

TEMPERATURE_BAND = (limit("cold_outer_c"), limit("cold_suit_c"), limit("cold_clothing_c"), limit("comfort_low_c"),
                    limit("comfort_high_c"), limit("heat_garment_c"), limit("heat_max_c"), limit("hot_outer_c"))
"""Dry-bulb temperature, deg C (linear). Blue 0 to 31, Green -20 to 45
(ordinary clothing or cooling garments), Yellow -50 (eutectic brines) to
122 (the limit of carbon-based life, Takai et al. 2008; the PHI-4 doc's Red
edge), Red beyond. A wet bulb at or above 31 deg C makes Blue or Green
Yellow (see `temperature_score`)."""

DOSE_BAND = (None, None, None, None, limit("dose_blue_sv_yr"), limit("dose_green_sv_yr"), limit("dose_yellow_sv_yr"),
             limit("dose_outer_sv_yr"))
"""Surface dose, Sv/yr (log), for a person. Blue to 20 mSv/yr (the ICRP
occupational limit averaged over five years), Green to 100 mSv/yr (past the
50 mSv single-year limit, `SINGLE_YEAR_DOSE_SV`, a visit is still allowed),
Yellow to 1 Sv/yr (shelter or shielding needed; 1 Sv in a year is the
threshold of radiation sickness), Red beyond, 0 at 1000 Sv/yr. Fills the gap
the PHI-4 doc leaves between 100 mSv/yr and 10 Gy/yr (the Eigen threshold of
`L_rad` is for microbes and stays there): activity-magnetism-radiation-
hydrosphere.md 4.6, Boss's defaults of 2026-10-10 (Blue 20 mSv/yr, 50 as the
single-year limit)."""

UV_BAND = (None, None, None, None, limit("uv_cover_index"), limit("uv_suit_index"), limit("uv_lethal_index"), limit("uv_outer_index"))
"""The DNA-weighted ultraviolet dose at the ground as a multiple of Earth's (log; `World.uv_index`, GEN.445). Blue below the index where
sun-cover clothing starts (3), Green to where a full-cover suit starts (30), Yellow to the unfiltered light of an ozone-free world (300:
shelter from the sun), Red beyond, 0 a decade further out. The score of the Radiation domain is the lower of this and the dose's."""

SINGLE_YEAR_DOSE_SV = limit("dose_single_year_sv")
"""float: ICRP's limit for any one year, 50 mSv."""

Bands = namedtuple("Bands", "armstrong o2_min o2_chronic o2_acute wet_bulb pressure temperature dose o2_high uv")
"""namedtuple: The limits the five core domains are scored against: the Armstrong pressure, the least inspired
O2, the chronic and acute inspired O2 edges, the uncompensable wet bulb, and the pressure, temperature and
dose bands with `o2_high`, the band of too much O2. `HUMAN` is the visitor's; `bands` builds a species'
(POP.12)."""


def bands(overrides=None):
    """
    The `Bands` of a set of limits: every `protection_limits.LIMIT_ROWS` value in use, except those `overrides`
    (`{limit name: value}`) replaces. Nothing else about a species changes the score.
    """
    def value(name):
        return limit(name) if overrides is None or name not in overrides else overrides[name]

    armstrong, chronic, acute = value("armstrong_kpa"), value("inspired_o2_chronic_kpa"), value("inspired_o2_acute_kpa")
    return Bands(
        armstrong=armstrong, o2_min=value("inspired_o2_min_kpa"), o2_chronic=chronic, o2_acute=acute,
        wet_bulb=value("wet_bulb_c"),
        pressure=(value("pressure_low_outer_kpa"), armstrong, value("mask_min_kpa"), value("pressure_no_gear_kpa"),
                  value("pressure_blue_max_kpa"), value("pressure_narcosis_kpa"), value("pressure_hyperbaric_max_kpa"),
                  value("pressure_outer_kpa")),
        temperature=(value("cold_outer_c"), value("cold_suit_c"), value("cold_clothing_c"), value("comfort_low_c"),
                     value("comfort_high_c"), value("heat_garment_c"), value("heat_max_c"), value("hot_outer_c")),
        dose=(None, None, None, None, value("dose_blue_sv_yr"), value("dose_green_sv_yr"), value("dose_yellow_sv_yr"),
              value("dose_outer_sv_yr")),
        o2_high=(None, None, None, None, chronic, chronic, acute, acute * 10),
        uv=(None, None, None, None, value("uv_cover_index"), value("uv_suit_index"), value("uv_lethal_index"), value("uv_outer_index")))


HUMAN = bands()
"""Bands: The human visitor's limits."""

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

LIEBIG_MIN_WEIGHT = 0.65
"""float: The share of PHI_bio's exponent that goes to the weakest likelihood
(Liebig's law of the minimum); the rest goes to the geometric mean. 0.65
minimises the error against the Master doc's own worked table, 0.046 RMS
(weights 0.55 to 0.75 give 0.046 to 0.050; habitability-index.md 7.1)."""

DOSE_THRESH_GY = 10.0
SIGMA_DOSE_GY = 100.0
"""float: L_rad's threshold and fall-off, Gy/yr (Daly, 2012; the doc)."""

SV_PER_GY = 3.06
"""float: Sv per Gy of the cosmic-ray and particle dose a surface takes:
Mars's surface is 0.236 Sv/yr in 0.077 Gy/yr (MSL/RAD, Hassler et al. 2014).
`L_rad` reads absorbed dose in Gy/yr and a person's bands read effective dose
in Sv/yr (`World.dose_gy_yr`); a gamma-only dose (weighting 1) is 3 times
higher in Gy than this says, which no world near the 10 Gy/yr line notices."""

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


# --- A world -----------------------------------------------------------------------------

@dataclass
class World:
    """What the index reads. `gases` maps species to mole fraction; `solvent`
    is the stable surface or subsurface liquid ("water", "hydrocarbon",
    "sulfuric" or None); `habitat_dose_sv_yr` is the dose where that liquid
    is (the surface dose for surface water, about the ground's for an ocean
    under ice); `water_source` is what a base can draw water from; `groundwater` is the liquid held underground by a world with none
    at its surface (`"none"`, `"reachable"` within `tuning.AQUIFER_REACH_KM` of the surface, or `"deep"`, GEN.654). The three dose parts (the cosmic rays
    with the heliosphere's multiplier, the stellar particles and the ground, Sv/yr; `None` when unknown, then the whole
    dose counts as cosmic rays), the dipole's `cutoff_gv`, the DNA-weighted `uv_index` and the airborne `dust_mg_m3` (`None` when unknown) are what the survival
    equipment (`physics/equipment.py`) reads beyond the score's own inputs."""
    pressure_kpa: float
    temperature_c: float
    gases: dict
    flux_variation: float = 0.0
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
    groundwater: str = "none"
    life_stage: str = "none"
    life_on_surface: bool = True
    dose_gcr_sv_yr: float = None
    dose_sep_sv_yr: float = None
    dose_ground_sv_yr: float = None
    dose_belt_sv_yr: float = None
    cutoff_gv: float = 0.0
    uv_index: float = 0.0
    dust_mg_m3: float = None

    def partial_kpa(self, gas):
        return self.gases.get(gas, 0.0) * self.pressure_kpa

    @property
    def dose_sv_yr(self):
        if self.surface_dose_sv_yr is not None:
            return self.surface_dose_sv_yr
        return self.bare_cosmic_ray_dose_sv_yr

    @property
    def bare_cosmic_ray_dose_sv_yr(self):
        """The dose of a world whose record carries none (a hand-built one, a body saved before GEN.87): the cosmic rays of the tables of
        GEN.453 (`radiation.gcr_dose_msv_yr`) under its air at the mean solar modulation, with no field, no stellar particles and no ground. The
        one curve for it; the model of the dose is `radiation.py`."""
        column = radiation.column_g_cm2(self.pressure_kpa * 1000.0, self.gravity_ms2 / constants.EARTH_GRAVITY)
        return radiation.gcr_dose_msv_yr(column) / 1000.0

    @property
    def habitat_dose_gy_yr(self):
        """The absorbed dose where the solvent is, Gy/yr (`SV_PER_GY`): the
        quantity `L_rad` and the microbes' thresholds are written in."""
        sv = self.habitat_dose_sv_yr if self.habitat_dose_sv_yr is not None else self.dose_sv_yr
        return sv / SV_PER_GY


# --- PHI-6 display ---------------------------------------------------------------------------

_EDGE_SCORES = (0.0, 0.4, 0.75, 1.0)
"""The score at the outer red, yellow, green and blue edges: Blue scores 1,
Green 0.75 to 1, Yellow 0.4 to 0.75, Red below 0.4 down to 0 at the outer
red edge, linear in between."""


def band_score(x, band, log=False):
    """The PHI-6 score of `x` in `band` (see the band constants): 1 inside
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


def is_habitable(overall_score):
    """Whether a body with this overall PHI-6 score (`phi4`, `None` for a body with no score) is habitable: its overall
    colour is one of `POTENTIALLY_HABITABLE_TIERS`. The API, the lists, the search, the map and the counts all call this."""
    return overall_score is not None and tier(overall_score) in POTENTIALLY_HABITABLE_TIERS


FALLBACK_TIERS = ("Yellow",)
"""tuple: The overall colours a request for a habitable world accepts when it must fall short of `POTENTIALLY_HABITABLE_TIERS`, provided
the Suitability colour (the sustainability factor) is at least `FALLBACK_SUITABILITY_TIER`."""

FALLBACK_SUITABILITY_TIERS = ("Blue", "Green", "Yellow")
"""tuple: The Suitability colours a Yellow fallback may have (not Red)."""


def habitable_request_result(overall_score, suitability_score):
    """
    How a world answers a request for a habitable world (Boss, 2026-10-23): `"aim"` when it is habitable (overall Blue
    or Green), `"fallback"` when the request must settle for less and the world is overall Yellow with Suitability at
    least Yellow ("technically habitable"), else `None`. Either score may be `None` (an unscored body): `None`.
    """
    if overall_score is None:
        return None
    if is_habitable(overall_score):
        return "aim"
    if (tier(overall_score) in FALLBACK_TIERS and suitability_score is not None
            and tier(suitability_score) in FALLBACK_SUITABILITY_TIERS):
        return "fallback"
    return None


def wet_bulb_c(temperature_c, relative_humidity):
    """Wet-bulb temperature, deg C (Stull, 2011; fitted for -20 to 50 deg C
    and 5 to 99% humidity, clamped to those)."""
    t = temperature_c
    rh = min(max(relative_humidity, 0.05), 0.99) * 100.0
    return (t * math.atan(0.151977 * math.sqrt(rh + 8.313659)) + math.atan(t + rh) - math.atan(rh - 1.676331)
            + 0.00391838 * rh ** 1.5 * math.atan(0.023101 * rh) - 4.686035)


def pressure_score(world, limits=HUMAN):
    return band_score(world.pressure_kpa, limits.pressure, log=True)


def temperature_extremes_c(world):
    """
    `(coldest, hottest)` surface temperature over the cycle of a world whose light swings (GEN.148). The
    equilibrium temperature goes as the fourth root of the insolation, so a flux that varies by the
    fraction `flux_variation` peak to peak moves T by `(1 +- f/2)^(1/4)`; a lone star (0) gives the one value.
    """
    f = min(max(world.flux_variation or 0.0, 0.0), 1.0)
    if f == 0.0:
        return world.temperature_c, world.temperature_c
    kelvin = world.temperature_c + 273.15
    return kelvin * (1.0 - f / 2.0) ** 0.25 - 273.15, kelvin * (1.0 + f / 2.0) ** 0.25 - 273.15


def temperature_score(world, limits=HUMAN):
    """Pressure's partner: the worse of the cycle's cold and hot extremes (they are one value for a lone star)."""
    cold, hot = temperature_extremes_c(world)
    score = min(band_score(cold, limits.temperature), band_score(hot, limits.temperature))
    if -20.0 <= hot <= 50.0 and wet_bulb_c(hot, world.relative_humidity) >= limits.wet_bulb:
        score = min(score, 0.75 - 1e-9)  # an uncompensable wet bulb needs cooling garments: Yellow
    return score


def inspired_o2_kpa(world, limits=HUMAN):
    """Inspired O2 (the inspired-gas equation): F_O2 (P - 6.3 kPa)."""
    return world.gases.get("O2", 0.0) * max(0.0, world.pressure_kpa - limits.armstrong)


def gas_inputs(world):
    """`[(gas, partial kPa, band)]` for the gases the filter rules read (CO2, CO, H2S, SO2): the value
    and the log band `band_score` scores it in. Shared with the describer (`protection`)."""
    chronic, acute = GAS_LIMITS_KPA["CO2"]
    inputs = [("CO2", world.partial_kpa("CO2"), (None, None, None, None, chronic, CO2_GREEN_KPA, acute, acute * 10))]
    for gas in ("CO", "H2S", "SO2"):
        chronic, acute = GAS_LIMITS_KPA[gas]
        inputs.append((gas, world.partial_kpa(gas), (None, None, None, None, chronic, chronic,
                                                     acute * FILTER_FACTOR, acute * FILTER_FACTOR * 10)))
    return inputs


O2_HIGH_BAND = (None, None, None, None, INSPIRED_O2_CHRONIC_KPA, INSPIRED_O2_CHRONIC_KPA,
                INSPIRED_O2_ACUTE_KPA, INSPIRED_O2_ACUTE_KPA * 10)
"""tuple: The log band of inspired O2 against too much of it (the low side is a mask, see `_gas_scores`)."""


def _gas_scores(world, limits=HUMAN):
    scores = {}
    o2 = inspired_o2_kpa(world, limits)
    # Low O2: a mask fixes it wherever a mask works (Green); where it doesn't,
    # the pressure domain already says so.
    o2_low = 1.0 if o2 >= limits.o2_min else 0.75
    scores["O2"] = min(o2_low, band_score(o2, limits.o2_high, log=True))
    for gas, kpa, band in gas_inputs(world):
        scores[gas] = band_score(kpa, band, log=True)
    return scores


def water_inputs(world):
    """`[(label, value, band)]` for the liquid water's acidity, activity and chaotropicity, `[]` when there
    is no water to drink. Shared with the describer (`protection`)."""
    if world.solvent != "water" or world.ph is None:
        return []
    return [("pH", world.ph, (0.0, 1.0, 5.0, 6.0, 8.5, 9.5, 11.5, 14.0)),
            ("water activity", world.water_activity,
             (0.5, WATER_ACTIVITY_CRIT, WATER_ACTIVITY_GREEN, 0.90, 1.0, 1.0, 1.0, 1.0)),
            ("chaotropicity", world.chaotropicity_kj_kg,
             (None, None, None, None, CHAOTROPICITY_CRIT_KJ_KG, CHAOTROPICITY_CRIT_KJ_KG,
              CHAOTROPICITY_CRIT_KJ_KG, CHAOTROPICITY_CRIT_KJ_KG * 2))]


AQUIFER_WATER_SCORE = 0.6
"""float: The Water domain of a world whose only liquid water is an aquifer within reach: Yellow, as the ocean under ice is."""


def water_score(world):
    """Water Chemistry: is the water drinkable. A world with no liquid water has none, which is Red
    (Boss, 2026-10-14); a liquid water body scores its worst of pH, water activity and chaotropicity."""
    if world.groundwater != "none":  # liquid water only underground: Yellow within a drill's reach, Red beyond it (GEN.654)
        return AQUIFER_WATER_SCORE if world.groundwater == "reachable" else 0.0
    if world.solvent != "water":
        return 0.0
    inputs = water_inputs(world)
    if not inputs:
        return 1.0
    return min(band_score(value, band) for _label, value, band in inputs)


NO_ATMOSPHERE_KPA = 0.05
"""float: Under this surface pressure a world has no air to speak of (it prints as 0.0 kPa): the Air domain
is Red with a score of 0, not Green for lack of toxins (GEN.215, Boss 2026-10-15)."""


def has_atmosphere(world):
    return world.pressure_kpa >= NO_ATMOSPHERE_KPA


def air_score(world, limits=HUMAN):
    """Air Chemistry: is the air breathable, the worst of oxygen and the filtered gases. A world with no
    atmosphere has no air to breathe: 0."""
    if not has_atmosphere(world):
        return 0.0
    return min(_gas_scores(world, limits).values())


def uv_score(world, limits=HUMAN):
    """The ground ultraviolet's score, or `None` for a world with no ultraviolet index."""
    if not world.uv_index or world.uv_index <= 0.0:
        return None
    return band_score(world.uv_index, limits.uv, log=True)


def radiation_score(world, limits=HUMAN):
    """Radiation: the lower of the dose's score and the ultraviolet's (a world with no ultraviolet index scores on its dose alone)."""
    dose = band_score(world.dose_sv_yr, limits.dose, log=True)
    uv = uv_score(world, limits)
    return dose if uv is None else min(dose, uv)


# --- Suitability (GEN.207) --------------------------------------------------------------------
# Long-term survival: enough water, food and shelter for a person to live on the
# surface (Exoplanetary Habitability.md, Table 2). The number is the lowest of three
# rungs, then held under the colour ceiling the other five domains set.

LIFE_STAGES = ("none", "microbial", "macroscopic", "biosphere")
"""tuple: The native life a person could eat from, least to most: nothing (an abiotic
surface), microbes only (chemosynthetic mats, plankton), multicellular algae, mosses and
the first land plants, a full biosphere of forests and animals."""

FOOD_SCORE = {"none": 0.2, "microbial": 0.5, "macroscopic": 0.85, "biosphere": 1.0}
"""dict: The food rung by `World.life_stage`: Table 2's Red (completely abiotic), Yellow
(chemosynthetic mats), Green (bryophytes, invertebrates; whose chemistry may not suit a
person) and Blue (forest canopies, complex food webs). planetGen's numbers, one per tier."""

WATER_ACCESS_SCORE = {"liquid": 0.9, "ice": 0.6, "aquifer": 0.6, "hydrated": 0.25, "vapour": 0.2, None: 0.0}
"""dict: The water rung by the most accessible water (`World.water_source`): Table 2's Green
(rivers, lakes and rain), Yellow (glacial ice, subglacial liquid, an aquifer within reach of a drill),
Red (perchlorate brine, hydrated minerals, vapour). Fresh liquid water that is also Blue Water Chemistry is Blue (1)."""

YELLOW_FLOOR = 0.4
"""float: The lowest score of Yellow: where a domain that cannot make the overall Red is held in the mean."""

CEILING_SCORE = {"Blue": 1.0, "Green": 0.9999, "Yellow": 0.7499, "Red": 0.3999}
"""dict: The highest score of each colour, the clamp a gate puts on a score."""


def water_access_score(world):
    """The water rung: how reachable the water is. Only a world with liquid water has any
    (the Water domain is Red without it)."""
    score = WATER_ACCESS_SCORE.get(world.water_source, 0.0)
    if world.water_source == "liquid" and water_score(world) >= 1.0:
        return 1.0
    return score


def food_score(world):
    """The food rung from the native life; microbes that live only in an ocean under ice (or in a deep aquifer) are out of
    reach, so they count as none. Flora and fauna are food wherever the world scores them: a world with flora, fauna or a
    civilization has food (Boss, 2026-10-24)."""
    stage = world.life_stage
    if not world.life_on_surface and stage == "microbial":
        stage = "none"
    return FOOD_SCORE[stage]


def shelter_score(world, limits=HUMAN):
    """The shelter rung: shelter from radiation and the open sky, which is the surface dose's
    colour (Table 2: on Red worlds only deep lava tubes or buried habitats, on Blue ones open ground)."""
    return radiation_score(world, limits)


def suitability_ceiling(core):
    """The best colour Suitability can have, from the five core domains' scores (`core`, a dict
    by name): Red if Radiation, Temperature or Pressure is Red or there is no liquid water (Water
    Red); Yellow if any of the five is Yellow, or Air is Red; Green if any is Green; else Blue."""
    colours = {name: tier(score) for name, score in core.items()}
    if any(colours[name] == "Red" for name in ("Radiation", "Temperature", "Pressure", "Water")):
        return "Red"
    if any(colour in ("Yellow", "Red") for colour in colours.values()):
        return "Yellow"
    if any(colour == "Green" for colour in colours.values()):
        return "Green"
    return "Blue"


def suitability_score(world, core=None, limits=HUMAN):
    """Suitability: the lowest of the water, food and shelter rungs, held under `suitability_ceiling`."""
    core = core if core is not None else core_scores(world, limits)
    lowest = min(water_access_score(world), food_score(world), shelter_score(world, limits))
    return min(lowest, CEILING_SCORE[suitability_ceiling(core)])


CORE_DOMAINS = ("Pressure", "Temperature", "Air", "Water", "Radiation")
"""tuple: The five domains of whether a person survives the surface, in the order the page lists them."""

DOMAINS = (*CORE_DOMAINS, "Suitability")
"""tuple: The six domains, in the order the page lists them."""


def core_scores(world, limits=HUMAN):
    """The five core domains' scores by name."""
    return {"Pressure": pressure_score(world, limits), "Temperature": temperature_score(world, limits),
            "Air": air_score(world, limits), "Water": water_score(world), "Radiation": radiation_score(world, limits)}


def overall_score(scores):
    """
    The overall score from the six domains' scores (`scores`, by name): their geometric mean
    with the colour gates applied first (Boss, 2026-10-11 and 2026-10-14):

    * Radiation and Temperature cap the overall colour at their own (Yellow caps at Yellow,
      Red at Red).
    * Pressure: Red makes the overall Red and Yellow caps it at Yellow, unless a full suit is
      needed anyway (Temperature or Radiation Red, or Air Red), when it changes nothing.
    * Air or Water at Yellow or Red caps the overall at Yellow (never lower on their own).
    * Suitability Red drops the overall to Yellow, never to Red.
    * Blue needs all six Blue.

    The ones that cannot make the overall Red (Air, Water, Suitability, and Pressure when a full
    suit is needed anyway) enter the mean no lower than the Yellow edge (`YELLOW_FLOOR`), so a
    zero among them does not drag the mean to 0.
    """
    colours = {name: tier(score) for name, score in scores.items()}
    suit_anyway = any(colours[name] == "Red" for name in ("Temperature", "Radiation", "Air"))
    soft = {"Air", "Water", "Suitability"} | ({"Pressure"} if suit_anyway else set())
    terms = [max(score, YELLOW_FLOOR) if name in soft else score for name, score in scores.items()]
    product = math.prod(terms)
    mean = product ** (1.0 / len(terms)) if product > 0 else 0.0
    caps = [colours["Radiation"], colours["Temperature"]]
    if not suit_anyway and colours["Pressure"] in ("Red", "Yellow"):
        caps.append(colours["Pressure"])
    for name in ("Air", "Water"):
        if colours[name] in ("Yellow", "Red"):
            caps.append("Yellow")
    if colours["Suitability"] == "Red":
        caps.append("Yellow")
    if any(colour != "Blue" for colour in colours.values()):
        caps.append("Green")
    return min(mean, CEILING_SCORE[max(caps, key=TIERS.index)])


def phi4(world, limits=HUMAN):
    """`{domain: (score, tier)}` for the six domains and `"PHI-4"`, the overall PHI-6 score
    (`overall_score`; the key and the stored name `phi4` are the ones it had with four domains). `limits` are the `Bands` of
    whoever is scored: the human visitor unless a species' (POP.12)."""
    core = core_scores(world, limits)
    domains = {**core, "Suitability": suitability_score(world, core, limits)}
    result = {name: (round(score, 4), tier(score)) for name, score in domains.items()}
    overall = overall_score(domains)
    result["PHI-4"] = round(overall, 4)
    return result


# --- PHI_bio, PHI_cpx -------------------------------------------------------------

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
    dose_gy = world.habitat_dose_gy_yr
    excess = max(0.0, dose_gy - DOSE_THRESH_GY)
    return math.exp(-0.5 * (excess / SIGMA_DOSE_GY) ** 2)


def phi_bio(world):
    """The Liebig blend of the four likelihoods (GEN.148; habitability-index.md
    7.1): the weakest term to the power `LIEBIG_MIN_WEIGHT` times the
    geometric mean of all four to the rest. Any term at 0 gives 0, so a
    missing solvent is still a hard gate (the doc's plain geometric mean
    forgave it), while one weak term no longer drags a world to the strict
    product's score."""
    terms = (l_solv(world), l_chem(world), l_ener(world), l_rad(world))
    lowest = min(terms)
    if lowest <= 0.0:
        return 0.0
    return lowest ** LIEBIG_MIN_WEIGHT * math.prod(terms) ** ((1.0 - LIEBIG_MIN_WEIGHT) / len(terms))


def toxin_factor(world):
    """The share of animal life that carbon monoxide and carbon dioxide leave (the toxic terms of `t_metazoa`): 1 for clean air."""
    return (math.exp(-(world.partial_kpa("CO") / CO_TOX_KPA) ** 2) * math.exp(-(world.partial_kpa("CO2") / CO2_TOX_KPA) ** 2))


def t_metazoa(world):
    """The human-comfort window of animal life: the oxygen band of a person's air and the toxins (`toxin_factor`). The life scores of
    a world read its life class's own oxygen floors instead (`life_complexity`, GEN.877)."""
    o2 = world.partial_kpa("O2")
    if o2 < METAZOA_O2_LOW_KPA:
        f_o2 = math.exp(-0.5 * ((METAZOA_O2_LOW_KPA - o2) / METAZOA_O2_SIGMA_LOW_KPA) ** 2)
    elif o2 > METAZOA_O2_HIGH_KPA:
        f_o2 = math.exp(-0.5 * ((o2 - METAZOA_O2_HIGH_KPA) / METAZOA_O2_SIGMA_HIGH_KPA) ** 2)
    else:
        f_o2 = 1.0
    return toxin_factor(world) * f_o2


def phi_cpx(world):
    return phi_bio(world) * t_metazoa(world)


def scores(world):
    """Every number the index gives a world, rounded for display."""
    out = {key: round(fn(world), 3) for key, fn in (
        ("L_solv", l_solv), ("L_chem", l_chem), ("L_ener", l_ener), ("L_rad", l_rad),
        ("PHI_bio", phi_bio), ("PHI_cpx", phi_cpx))}
    out.update(phi4(world))
    return out
