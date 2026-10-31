# planetgen/physics/habitability_world.py

"""
The habitability score of every planet and moon (GEN.89): turns a body's
stored air, water, temperature, radiation and light into the `World` that
`habitability.scores` reads (GEN.84) and stores the result.

The mapping, every rule a planetGen default (docs/design/habitability/habitability-index.md
section 8):

* **Air.** Total pressure, temperature and gravity as stored; each species
  as a mole fraction of the total (the balance gas counts toward the total
  only); relative humidity is the water vapour against its saturation
  pressure at that temperature (Magnus formula, clamped to 5 to 99%).
* **Solvent.** Water when the body has a surface ocean (a Hycean one too) or
  an ocean under its ice; liquid methane when the air holds saturated CH4
  below its critical 190 K (Titan); otherwise none, which zeroes `L_solv`.
  Sulfuric acid is never generated yet. The ocean's pH and water activity
  are its stored ones; a chloride brine is chaotropic, 200 kJ/kg per unit
  of water activity below 1 (70 to 140 kJ/kg for 0.65 to 0.30).
* **Inventories.** Carbon, nitrogen and hydrogen in Earth units: hydrogen is
  the water (`water_mass_fraction` times the mass, against Earth's
  1.4e21 kg of ocean); carbon and nitrogen are at least that (volatiles
  arrive with water), or the air's own carbon (CO2, CO, CH4) against
  7.5e19 kg or its nitrogen against 4.0e18 kg, when larger.
* **Phosphate.** The ocean's `phosphorus` class: 1000 umol/L "high" (soda
  lakes), 2.3 "limited" (Earth), 0.01 "starved" (ice-sealed); a methane sea
  takes the "limited" value (the term has no non-aqueous model).
* **Energy.** A surface liquid gets 6% of the starlight (Earth's 1361 W/m2
  gives its 80 W/m2 of usable light) or the chemical power of the body's
  heat flow, 1% of it (Earth's 0.09 W/m2 gives 9e-4), whichever is more; an
  ocean under ice gets the chemical power only, and the heat flow it counts is its own plus the tide of its orbit (`tidal_flux_w_m2`, GEN.737).
* **Dose.** The stored `surface_dose_msv_yr`, in Sv; the habitat of an ocean
  under ice sees only the ground's.
* **Water source** (for the base): "liquid" at a surface ocean, "ice" for
  an icy surface or shell, "vapour" for a steam world, "hydrated" for a dry
  body with water in its rock, else vapour in the air or none.
* **Life** (GEN.207). The food rung of Suitability reads the most advanced native
  life as `native_life` (`native_life`): none, microbial, macroscopic (the first
  flora and fauna) or a full biosphere, from the body's life chemistry and the
  milestone its evolutionary timeline reached. Life that lives only in an ocean
  under ice is out of reach of a person (`life_on_surface` false). Life is drawn after
  the score, so `rescore` runs again once it is known.
* **Compact hosts.** A pulsar, neutron star or black hole's planet is rated
  as if the surface took 1000 Sv/yr (`LETHAL_DOSE_SV_YR`): Red radiation,
  the heaviest survival equipment and a note naming the host. No X-ray binaries
  are generated yet; the rule is by the host's class.

A gas giant has no surface and no score (every field `None`), except the two energy inputs the sheet shows for every body (GEN.784):
`energy_flux_w_m2` is the energy a giant radiates through its 1 bar level, the starlight it absorbs (its Bond albedo, `giant_draw.BOND_ALBEDO`)
plus the heat of its own interior (`giant_draw.own_flux_w_m2`), which is `sigma T_eff^4`; `flux_variation` is the swing of its starlight that
`insolation` gives every planet.
"""

import math

from planetgen import tuning
from planetgen.physics import aerosols, aquifers, atmosphere, constants, equipment, giant_draw, habitability, life_scores, magnetism, outside_hours, protection_limits, radiation, regional_temperature
from planetgen.physics import unit_registry as units

_pa_to_kpa = units.converter("Pa", "kPa")
_kpa_to_pa = units.converter("kPa", "Pa")
_km_to_m = units.converter("km", "m")
_au_to_m = units.converter("au", "m")

TIER_INDEX = {"Blue": 0, "Green": 1, "Yellow": 2, "Red": 3}
TIER_NAMES = tuple(TIER_INDEX)
"""tuple: The tier colours by `tier_<domain>` number."""
"""dict: A domain's colour as stored (`tier_<domain>`): 0 is best."""

DOMAINS = tuple(name.lower() for name in habitability.DOMAINS)
"""tuple: The six domains as stored (`phi4_<domain>`, `tier_<domain>`): pressure, temperature, air, water,
radiation, suitability."""
CORE_DOMAINS = tuple(name.lower() for name in habitability.CORE_DOMAINS)
"""tuple: The five that say whether a person survives the surface (Suitability is the long-term sixth)."""

BODY_FIELDS = (
    "phi4", *(f"phi4_{domain}" for domain in DOMAINS), *(f"tier_{domain}" for domain in DOMAINS),
    *life_scores.FIELDS, "l_solv", "l_chem", "l_ener", "l_rad",
    *equipment.FIELDS, *regional_temperature.FIELDS, *outside_hours.FIELDS, "energy_flux_w_m2", "flux_variation",
)
"""tuple: What `generate` sets on a planet or moon, as stored in `planets` and
`moons` (`energy_flux_w_m2` and `flux_variation` are the inputs the score keeps: the
dose refresh re-scores a body from its columns alone)."""

SCORE_FIELDS = tuple(name for name in BODY_FIELDS if name not in ("energy_flux_w_m2", "flux_variation"))
"""tuple: The stored fields `compute` returns (all but the stored energy flux and flux variation)."""

LIKELIHOOD_FIELDS = ("phi_bio", "phi_cpx")
"""tuple: The microbial and complex-life likelihoods `compute` also returns and `generate` and `rescore` keep on the body
in memory: the life draw reads them (`life.apply_life_data`) and the life scores are made from them. They are not stored
(GEN.649), so a loaded body has none; the stored `life_*` scores say what they gave."""

LETHAL_DOSE_SV_YR = 1000.0
"""float: The dose a compact object's planet is rated at, Sv/yr: the outer
Red edge of the radiation domain, where its score is 0."""

LETHAL_HABITAT_DOSE_GY_YR = 10_000.0
"""float: The absorbed dose where a compact object's planet keeps its life, Gy/yr:
10 kGy/yr, the D10 of Deinococcus radiodurans (5 to 10 kGy; design note
activity-magnetism-radiation-hydrosphere.md 4.6), the hardest microbe known,
so no `L_rad` is left."""

COMPACT_HOSTS = {"NS": "neutron star", "BH": "black hole"}
"""dict: The `yerkes_class` of a compact host and what it is called."""

DEFAULT_HARSHNESS = 0.5
"""float: The harshness of a world with no score (a gas giant, a body saved
before the scores existed): the middle of the 0 to 1 range, so it moves
nothing (population-and-politics.md, tech-level trait hooks)."""

LIGHT_SHARE = 0.06
CHEMICAL_SHARE = 0.01
"""float: The share of starlight, and of the body's heat flow, that life can use."""

PHOSPHATE_UMOL_L = {"high": 1000.0, "limited": 2.3, "starved": 0.01, None: 0.0}
"""dict: An ocean's `phosphorus` class as phosphate, umol/L. The water of an aquifer, which has no ocean class, is `limited`
(rock weathering supplies it; DESIGN stand-in, GEN.654)."""
NON_AQUEOUS_PHOSPHATE_UMOL_L = 2.3

CHLORIDE_BRINE_CHAOTROPICITY = 200.0
"""float: kJ/kg per unit of water activity below 1 in a chloride brine."""

METHANE_SEA_MAX_K = 190.0
METHANE_SATURATION = 0.98
"""float: A methane sea needs CH4 at this share of its vapour pressure."""

HYDRATED_MIN_WATER_FRACTION = 1e-5

EARTH_OCEAN_KG = 1.4e21
EARTH_CARBON_KG = 7.5e19
EARTH_NITROGEN_KG = 4.0e18

MOLAR_KG = {
    "o2": 0.031998, "co2": 0.04401, "co": 0.02801, "n2": 0.028013, "ar": 0.039948, "h2": 0.002016,
    "h2o": 0.018015, "ch4": 0.01604, "h2s": 0.03408, "so2": 0.06407,
}
"""dict: The stored species' molar masses, kg/mol."""
BALANCE_MOLAR_KG = 0.028
CARBON_KG_MOL = 0.012011
NITROGEN_KG_MOL = 0.014007

_GAS_NAMES = {"o2": "O2", "co2": "CO2", "co": "CO", "h2s": "H2S", "so2": "SO2"}
"""dict: The stored partials the world reads by their lower-case key; every other gas of `protection_limits.TOXIC` (NH3, HF, HCl, Cl2, NO2, HCN, O3) is read by its formula (GEN.657)."""


def values_from_body(body):
    """A planet or moon (in memory) as the column-named dict `compute` reads."""
    values = {
        "planet_class": getattr(body, "planet_class", None),
        "body_type": getattr(body, "body_type", None),
        "atmospheric_pressure_pa": getattr(body, "atmospheric_pressure", None),
        "surface_temperature_k": getattr(body, "surface_temperature", None),
        "gravity_g": getattr(body, "gravity", None),
        "radius_km": getattr(body, "radius", None),
        "mass_kg": getattr(body, "mass", None),
        "hydrosphere": getattr(body, "hydrosphere", None),
        "ocean_class": getattr(body, "ocean_class", None),
        "ocean_ph": getattr(body, "ocean_ph", None),
        "water_activity": getattr(body, "water_activity", None),
        "phosphorus": getattr(body, "phosphorus", None),
        "water_mass_fraction": getattr(body, "water_mass_fraction", None),
        "ice_mass_fraction": getattr(body, "ice_mass_fraction", None),
        "rock_mg_si_ratio": getattr(body, "rock_mg_si_ratio", None),
        "aquifer_depth_km": getattr(body, "aquifer_depth_km", None),
        "aquifer_thickness_km": getattr(body, "aquifer_thickness_km", None),
        "surface_dose_msv_yr": getattr(body, "surface_dose_msv_yr", None),
        "dose_ground_msv_yr": getattr(body, "dose_ground_msv_yr", None),
        "dose_gcr_msv_yr": getattr(body, "dose_gcr_msv_yr", None),
        "dose_sep_msv_yr": getattr(body, "dose_sep_msv_yr", None),
        "dose_belt_msv_yr": getattr(body, "dose_belt_msv_yr", None),
        "dose_helio_mult": getattr(body, "dose_helio_mult", None),
        "uv_surface_index": getattr(body, "uv_surface_index", None),
        "tau_dust": getattr(body, "tau_dust", None),
        "scale_height_km": getattr(body, "scale_height", None),
        "magnetic_moment_a_m2": getattr(body, "magnetic_moment_a_m2", None),
        "rotation_period_hours": getattr(body, "rotation_period_hours", None),
        "daylight_hours": getattr(body, "daylight_hours", None),
        "axial_tilt_deg": getattr(body, "axial_tilt_deg", None),
        "period_years": getattr(body, "period", None),
        "is_moon": bool(getattr(body, "is_moon", False)),
        "energy_flux_w_m2": getattr(body, "energy_flux_w_m2", None),
        "native_life": native_life(getattr(body, "life_class", None), getattr(body, "life_stage", None),
                                   getattr(body, "life_system_age_gyr", None) is not None),
        "flux_variation": getattr(body, "flux_variation", None),
        "life_class": getattr(body, "life_class", None),
        "life_class_variant": getattr(body, "life_class_variant", None),
        "life_animal_grade": getattr(body, "life_animal_grade", None),
        "life_stage": getattr(body, "life_stage", None),
        "evolutionary_pace": getattr(body, "evolutionary_pace", None),
        "life_milestone_age_gyr": getattr(body, "life_milestone_age_gyr", None),
        "life_system_age_gyr": getattr(body, "life_system_age_gyr", None),
        "internal_heat_flux_w_m2": getattr(body, "internal_heat_flux_w_m2", None),
        "tidal_flux_w_m2": getattr(body, "tidal_flux_w_m2", None),
        "mantle_redox": getattr(body, "mantle_redox", None),
        "star_teff_k": getattr(getattr(body, "star", None), "temperature", None),
    }
    values[atmosphere.GAS_FIELD] = atmosphere.gas_partials_pa(body)
    return values


_STAGE_OF_MILESTONE = {
    "technological_civilization": "biosphere", "multicellularity": "macroscopic", "complex_cells": "microbial",
    "photosynthesis": "microbial", "abiogenesis": "microbial",
}
"""dict: `{evolutionary milestone: life stage}`."""


def native_life(life_class, milestone, has_timeline):
    """
    The `habitability.LIFE_STAGES` stage of a body's native life, from its life
    class (`None` for none), the latest evolutionary milestone its timeline
    reached (`life_stage` on the body, `None` for none) and whether it has a
    timeline (`life_system_age_gyr` is set). Every body with a class has a
    timeline (a moon or a cold world too); one without (a hand-built body) is microbial; a
    timeline that has reached no milestone yet is none.
    """
    if not life_class:
        return "none"
    if not has_timeline:
        return "microbial"
    return _STAGE_OF_MILESTONE.get(milestone, "none")


def saturation_kpa(temperature_c):
    """Saturated water vapour pressure over liquid water, kPa (Magnus
    formula, 6.1094 hPa at 0 deg C; temperatures clamped to -45 to 60)."""
    t = min(max(temperature_c, -45.0), 60.0)
    return 0.61094 * math.exp(17.625 * t / (t + 243.04))


def relative_humidity(values):
    """The air's water vapour against saturation at the surface temperature,
    clamped to 5 to 99% (the range Stull's wet-bulb formula is fitted to)."""
    temperature_c = values["surface_temperature_k"] - constants.CELSIUS_ZERO_K
    share = atmosphere.partial_kpa(values, "h2o") / saturation_kpa(temperature_c)
    return min(max(share, 0.05), 0.99)


def groundwater(values):
    """`"none"`, `"reachable"` or `"deep"`: the liquid water a world with none at its surface holds underground (GEN.654,
    `physics/aquifers.py`), reachable when the top of its water-bearing layer is within `tuning.AQUIFER_REACH_KM`."""
    if not (values.get("aquifer_thickness_km") or 0.0) or values.get("hydrosphere") not in aquifers.SURFACE_STATES:
        return "none"
    return "reachable" if (values.get("aquifer_depth_km") or 0.0) <= tuning.AQUIFER_REACH_KM else "deep"


def solvent(values):
    """"water", "hydrocarbon" or `None`: the stable liquid life would use."""
    if values.get("hydrosphere") in ("surface ocean", "brine seep", "hycean", "ice-covered ocean") or groundwater(values) != "none":
        return "water"
    temperature = values["surface_temperature_k"]
    methane = atmosphere.partial_kpa(values, "ch4")
    if (methane > 0.0 and temperature <= METHANE_SEA_MAX_K
            and methane >= METHANE_SATURATION * atmosphere.vapour_pressure_kpa("ch4", temperature)):
        return "hydrocarbon"
    return None


def water_source(values):
    """What a base can draw water from (`habitability.WATER_SOURCE_FACTOR`)."""
    state = values.get("hydrosphere")
    if state in ("surface ocean", "brine seep", "hycean"):
        return "liquid"
    if state in ("ice", "ice-covered ocean"):
        return "ice"
    if groundwater(values) == "reachable":
        return "aquifer"
    if state == "vapour":
        return "vapour"
    if (values.get("water_mass_fraction") or 0.0) >= HYDRATED_MIN_WATER_FRACTION:
        return "hydrated"
    return "vapour" if atmosphere.partial_kpa(values, "h2o") > 0.0 else None


def air_inventories_kg(values):
    """`(carbon_kg, nitrogen_kg)` in the air, from the species' partial
    pressures; the balance gas is taken as 28 g/mol and carries neither."""
    pressure_kpa = _pa_to_kpa(values.get("atmospheric_pressure_pa") or 0.0)
    gravity = (values.get("gravity_g") or 0.0) * constants.EARTH_GRAVITY
    radius_m = _km_to_m(values.get("radius_km") or 0.0)
    if pressure_kpa <= 0.0 or gravity <= 0.0 or radius_m <= 0.0:
        return 0.0, 0.0
    partials = atmosphere.partial_pressures_kpa(values)
    balance = max(0.0, pressure_kpa - sum(partials.values()))
    mean_molar = (sum(partials[gas] * MOLAR_KG[gas] for gas in partials) + balance * BALANCE_MOLAR_KG) / (
        sum(partials.values()) + balance)
    moles = _kpa_to_pa(pressure_kpa) * 4.0 * math.pi * radius_m ** 2 / (gravity * mean_molar)
    share = {gas: partials[gas] / pressure_kpa for gas in partials}
    carbon = moles * (share["co2"] + share["co"] + share["ch4"]) * CARBON_KG_MOL
    nitrogen = moles * 2.0 * share["n2"] * NITROGEN_KG_MOL
    return carbon, nitrogen


def inventories(values):
    """Carbon, nitrogen and hydrogen in Earth units."""
    hydrogen = (values.get("water_mass_fraction") or 0.0) * (values.get("mass_kg") or 0.0) / EARTH_OCEAN_KG
    carbon, nitrogen = air_inventories_kg(values)
    return {"H": hydrogen, "C": max(hydrogen, carbon / EARTH_CARBON_KG),
            "N": max(hydrogen, nitrogen / EARTH_NITROGEN_KG)}


def giant_energy_flux_w_m2(body):
    """The energy a giant radiates through its 1 bar level, W/m2 (GEN.784): `sigma T_eff^4` of the starlight it absorbs and the heat of its own
    interior (`giant_draw.effective_temperature_k`). A giant that is a moon has no starlight of its own here (`giant_draw.starlight_w_m2`) and
    radiates its own heat."""
    starlight = giant_draw.starlight_w_m2(body) or 0.0
    return constants.STEFAN_BOLTZMANN_CONSTANT * giant_draw.effective_temperature_k(starlight, giant_draw.own_flux_w_m2(body)) ** 4


def energy_flux_w_m2(values, luminosity_w, distance_au, internal_flux_w_m2):
    """The light or chemical power available to life, W/m2 (see the module
    docstring). `luminosity_w` and `distance_au` describe the star (either
    `None`: no light)."""
    chemical = CHEMICAL_SHARE * (internal_flux_w_m2 or 0.0)
    if values.get("hydrosphere") == "ice-covered ocean" or groundwater(values) != "none" or not luminosity_w or not distance_au:
        return chemical
    light = LIGHT_SHARE * luminosity_w / (4.0 * math.pi * _au_to_m(distance_au) ** 2)
    return max(light, chemical)


def world_from(values, lethal_host=None):
    """The `habitability.World` of a rocky body's `values` (see
    `values_from_body`); `lethal_host` floors the dose it takes."""
    pressure_kpa = (values.get("atmospheric_pressure_pa") or 0.0) / 1000.0
    medium = solvent(values)
    under_ice = values.get("hydrosphere") == "ice-covered ocean"
    underground = groundwater(values)
    surface_dose = (values.get("surface_dose_msv_yr") or 0.0) / 1000.0
    habitat_dose = (values.get("dose_ground_msv_yr") or 0.0) / 1000.0 if under_ice or underground != "none" else surface_dose
    if lethal_host:
        surface_dose = max(surface_dose, LETHAL_DOSE_SV_YR)
        habitat_dose = max(habitat_dose, LETHAL_HABITAT_DOSE_GY_YR * habitability.SV_PER_GY)
    parts = _dose_parts_sv(values) if not lethal_host else (None, None, None)
    aqueous = medium == "water"
    activity = values.get("water_activity") if aqueous and values.get("water_activity") is not None else 1.0
    brine = values.get("ocean_class") == "chloride brine"
    gases = {}
    if pressure_kpa > 0.0:
        for gas, name in _GAS_NAMES.items():
            gases[name] = atmosphere.partial_kpa(values, gas) / pressure_kpa
        air = atmosphere.gas_partials_pa(values)
        for formula in protection_limits.TOXIC:
            gases.setdefault(formula, air.get(formula, 0.0) / (1000.0 * pressure_kpa))
    return habitability.World(
        pressure_kpa=pressure_kpa,
        temperature_c=values["surface_temperature_k"] - constants.CELSIUS_ZERO_K,
        gases=gases,
        gravity_ms2=(values.get("gravity_g") or 0.0) * constants.EARTH_GRAVITY,
        flux_variation=values.get("flux_variation") or 0.0,
        relative_humidity=relative_humidity(values),
        solvent=medium,
        water_activity=activity,
        chaotropicity_kj_kg=CHLORIDE_BRINE_CHAOTROPICITY * (1.0 - activity) if aqueous and brine else 0.0,
        ph=values.get("ocean_ph") if aqueous else None,
        inventories=inventories(values),
        phosphate_umol_l=(PHOSPHATE_UMOL_L.get(values.get("phosphorus") or ("limited" if underground != "none" else None), 0.0)
                          if aqueous else NON_AQUEOUS_PHOSPHATE_UMOL_L if medium else 0.0),
        energy_flux_w_m2=values.get("energy_flux_w_m2") or 0.0,
        surface_dose_sv_yr=surface_dose,
        habitat_dose_sv_yr=habitat_dose,
        water_source=water_source(values),
        groundwater=underground,
        life_stage=values.get("native_life") or "none",
        life_on_surface=values.get("hydrosphere") != "ice-covered ocean" and underground != "deep",
        dose_gcr_sv_yr=parts[0], dose_sep_sv_yr=parts[1], dose_ground_sv_yr=parts[2],
        dose_belt_sv_yr=(values.get("dose_belt_msv_yr") or 0.0) / 1000.0 if parts[0] is not None else None,
        cutoff_gv=cutoff_gv(values), uv_index=values.get("uv_surface_index") or 0.0,
        dust_mg_m3=aerosols.dust_loading_mg_m3(values.get("tau_dust"), values.get("scale_height_km")),
    )


def _dose_parts_sv(values):
    """`(cosmic rays with the heliosphere's multiplier, stellar particles, ground)` in Sv/yr from the stored parts,
    or three `None` when a body was saved without them."""
    gcr, sep, ground = (values.get(name) for name in ("dose_gcr_msv_yr", "dose_sep_msv_yr", "dose_ground_msv_yr"))
    if gcr is None or sep is None or ground is None:
        return None, None, None
    multiplier = values.get("dose_helio_mult") or 1.0
    return gcr * multiplier / 1000.0, sep / 1000.0, ground / 1000.0


def cutoff_gv(values):
    """The equatorial cutoff rigidity of the body's dipole, GV (0 without a field): the stored dipole moment over the radius squared
    (`radiation.cutoff_rigidity_gv`)."""
    return radiation.cutoff_rigidity_gv(values.get("magnetic_moment_a_m2") or 0.0, values.get("radius_km") or 0.0)


def compute(values, compact_host=None):
    """
    The `SCORE_FIELDS` of a rocky body from its column-named `values`
    (`values_from_body`, or a `planets`/`moons` row plus the stored
    `energy_flux_w_m2`). `compact_host` is a `COMPACT_HOSTS` name (or any
    text) for a pulsar's, neutron star's or black hole's planet. A body
    without a temperature, or a gas giant, scores nothing (every field
    `None`).
    """
    if values.get("body_type") == "g" or values.get("surface_temperature_k") is None:
        return {name: None for name in (*SCORE_FIELDS, *LIKELIHOOD_FIELDS)}
    world = world_from(values, lethal_host=compact_host)
    # The life scores come first: the Suitability food rung is made from them (Flora and Fauna, GEN.649)
    stage = life_scores.stage_of(world.life_stage, values.get("life_stage"), values.get("evolutionary_pace") is not None)
    life = class_life_scores(values, world, stage)
    world.life_stage = life_scores.food_stage(life[0], life[1], life[2])
    result = habitability.scores(world)
    out = {"phi4": result["PHI-4"], "phi_bio": result["PHI_bio"], "phi_cpx": result["PHI_cpx"],
           "l_solv": result["L_solv"], "l_chem": result["L_chem"], "l_ener": result["L_ener"], "l_rad": result["L_rad"]}
    out.update(zip(life_scores.FIELDS, life))
    for domain in DOMAINS:
        score, colour = result[domain.capitalize()]
        out[f"phi4_{domain}"] = score
        out[f"tier_{domain}"] = TIER_INDEX[colour]
    out.update(equipment.compute(world))
    out.update(regional_temperature.compute(values))
    out.update(outside_hours.fields(world, out, {**values, **out}))
    return out


def class_life_scores(values, world, stage):
    """
    `(microbial, flora, fauna, sentient)` of the life class a body stores (`life_class`, `life_class_variant`), at the `stage` its timeline
    reached: the class's likelihood on the body is the microbial one, the complex-life one that times the factor its complexity index gives
    (`life_complexity`) and the toxins of the air, and the animal grade is the stored one (`life_animal_grade`, which includes route R4) or the
    one the index allows. No class is no life.
    """
    from planetgen.physics import life_classes, life_complexity      # life_classes reads this module
    class_id = life_classes.class_id(values.get("life_class"), values.get("life_class_variant"))
    if class_id is None or stage is None:
        return 0, 0, 0, 0
    result = life_classes.score_classes(values).get(class_id)
    if result is None:
        return 0, 0, 0, 0
    row = life_classes.BY_ID[class_id]
    complexity = life_complexity.assess(row, result.npp, atmosphere.partial_kpa(values, "o2"), values["surface_temperature_k"])
    grade = values.get("life_animal_grade")
    complex_likelihood = result.likelihood * life_complexity.complex_factor(complexity.index) * habitability.toxin_factor(world)
    return life_scores.scores(round(result.likelihood, 3), round(complex_likelihood, 3), stage, values.get("life_milestone_age_gyr"),
                              values.get("life_system_age_gyr"), light_producer=row.light, grade=complexity.grade if grade is None else grade)


def compact_host_name(star):
    """The `COMPACT_HOSTS` name of a star object, or `None`. A pulsar is
    called that."""
    kind = COMPACT_HOSTS.get(getattr(star, "yerkes_class", None))
    if kind == "neutron star" and getattr(star, "pulsar_type", "non-pulsing") != "non-pulsing":
        return "pulsar"
    return kind


def harshness(values):
    """
    How hard a world is on the people living on it, 0 (every domain Blue) to
    1 (a domain at 0), for the tech draw and facility ratings: one minus the
    worst of the five core domain scores (`phi4_<domain>`, not Suitability, which
    reads the life that the tech draw comes after), so it follows the colour tier
    of the worst domain. `DEFAULT_HARSHNESS` for a body with no score. `values`
    is a body or a `planets`/`moons` row, read by column name.
    """
    get = values.get if hasattr(values, "get") else (lambda name: getattr(values, name, None))
    scores = [get(f"phi4_{domain}") for domain in CORE_DOMAINS]
    if any(score is None for score in scores):
        return DEFAULT_HARSHNESS
    return round(1.0 - min(scores), 4)


def rescore(body):
    """Recomputes a scored body's `SCORE_FIELDS` once its life is known (`native_life`
    is part of the score's inputs and life is drawn after it). A gas giant, or a body
    that was not scored, is left alone."""
    if body.body_type == "g" or getattr(body, "surface_temperature", None) is None:
        return
    if getattr(body, "phi4", None) is None:
        return
    values = values_from_body(body)
    star = getattr(body, "star", None)
    for name, value in compute(values, compact_host_name(star)).items():
        setattr(body, name, value)


def _clear(body):
    for name in (*BODY_FIELDS, *LIKELIHOOD_FIELDS):
        setattr(body, name, None)


def generate(body, internal_flux_w_m2):
    """Sets `BODY_FIELDS` on a planet or moon. Call last, once its air,
    water and dose are set and `magnetism.update_exposure` has given its
    distance from its star; `internal_flux_w_m2` is its own heat flow, and the tide of its orbit (`tidal_flux_w_m2`, GEN.737) is added to it.
    The flux swing (`flux_variation`) is an input the body keeps, as the dose refresh keeps it: a body scored again after life changed its air
    keeps the swing the insolation gave it, and is scored with it."""
    swing = getattr(body, "flux_variation", None)
    _clear(body)
    if body.body_type == "g":
        body.flux_variation = swing
        body.energy_flux_w_m2 = giant_energy_flux_w_m2(body)
        return
    body.flux_variation = swing
    values = values_from_body(body)
    if values["surface_temperature_k"] is None:
        return
    star = getattr(body, "star", None)
    body.energy_flux_w_m2 = energy_flux_w_m2(
        values, getattr(star, "luminosity", None), getattr(body, "_star_distance_au", None),
        internal_flux_w_m2 + (getattr(body, "tidal_flux_w_m2", None) or 0.0))
    values["energy_flux_w_m2"] = body.energy_flux_w_m2
    for name, value in compute(values, compact_host_name(star)).items():
        setattr(body, name, value)
