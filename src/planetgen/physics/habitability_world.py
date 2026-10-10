# planetgen/physics/habitability_world.py

"""
The habitability score of every planet and moon (GEN.89): turns a body's
stored air, water, temperature, radiation and light into the `World` that
`habitability.scores` reads (GEN.84) and stores the result.

The mapping, every rule a planetGen default (docs/design/habitability-index.md
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
  ocean under ice gets the chemical power only. Tides are not counted.
* **Dose.** The stored `surface_dose_msv_yr`, in Sv; the habitat of an ocean
  under ice sees only the ground's.
* **Water source** (for the base): "liquid" at a surface ocean, "ice" for
  an icy surface or shell, "vapour" for a steam world, "hydrated" for a dry
  body with water in its rock, else vapour in the air or none.
* **Compact hosts.** A pulsar, neutron star or black hole's planet is rated
  as if the surface took 1000 Sv/yr (`LETHAL_DOSE_SV_YR`): Red radiation,
  the lowest equipment tier and a note naming the host. No X-ray binaries
  are generated yet; the rule is by the host's class.

A gas giant has no surface and no score (every field `None`).
"""

import math

from planetgen.physics import atmosphere, constants, habitability

TIER_INDEX = {"Blue": 0, "Green": 1, "Yellow": 2, "Red": 3}
"""dict: A domain's colour as stored (`tier_<domain>`): 0 is best."""

EQUIPMENT_NAMES = (
    "shirtsleeve", "breathing mask", "mask with scrubber", "sealed suit",
    "full life support with radiation hardening",
)
"""tuple: `equipment_tier` 0 to 4 and the kit it names (`habitability.equipment`)."""

EQUIPMENT_LABELS = (
    "Shirtsleeve", "Mask", "Mask and scrubber", "Sealed suit", "Full life support",
)
"""tuple: The same tiers as the page and the search show them."""

DOMAINS = ("pressure", "temperature", "chemistry", "radiation")

BODY_FIELDS = (
    "phi4", "phi4_pressure", "phi4_temperature", "phi4_chemistry", "phi4_radiation",
    "tier_pressure", "tier_temperature", "tier_chemistry", "tier_radiation",
    "phi_bio", "phi_cpx", "phi_tech", "l_solv", "l_chem", "l_ener", "l_rad",
    "equipment_tier", "hab_note", "energy_flux_w_m2",
)
"""tuple: What `generate` sets on a planet or moon, as stored in `planets` and
`moons` (`energy_flux_w_m2` is the one input the score keeps: the dose
refresh re-scores a body from its columns alone)."""

SCORE_FIELDS = tuple(name for name in BODY_FIELDS if name != "energy_flux_w_m2")
"""tuple: The fields `compute` returns (all but the stored energy flux)."""

LETHAL_DOSE_SV_YR = 1000.0
"""float: The dose a compact object's planet is rated at, Sv/yr: the outer
Red edge of the radiation domain, where its score is 0."""

COMPACT_HOSTS = {"NS": "neutron star", "BH": "black hole"}
"""dict: The `yerkes_class` of a compact host and what it is called."""

LIGHT_SHARE = 0.06
CHEMICAL_SHARE = 0.01
"""float: The share of starlight, and of the body's heat flow, that life can use."""

PHOSPHATE_UMOL_L = {"high": 1000.0, "limited": 2.3, "starved": 0.01, None: 0.0}
"""dict: An ocean's `phosphorus` class as phosphate, umol/L."""
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
        "surface_dose_msv_yr": getattr(body, "surface_dose_msv_yr", None),
        "dose_ground_msv_yr": getattr(body, "dose_ground_msv_yr", None),
        "energy_flux_w_m2": getattr(body, "energy_flux_w_m2", None),
    }
    for field in atmosphere.PARTIAL_PRESSURE_FIELDS:
        values[field] = getattr(body, field, None)
    return values


def saturation_kpa(temperature_c):
    """Saturated water vapour pressure over liquid water, kPa (Magnus
    formula, 6.1094 hPa at 0 deg C; temperatures clamped to -45 to 60)."""
    t = min(max(temperature_c, -45.0), 60.0)
    return 0.61094 * math.exp(17.625 * t / (t + 243.04))


def relative_humidity(values):
    """The air's water vapour against saturation at the surface temperature,
    clamped to 5 to 99% (the range Stull's wet-bulb formula is fitted to)."""
    temperature_c = values["surface_temperature_k"] - constants.CELSIUS_ZERO_K
    share = (values.get("p_h2o_kpa") or 0.0) / saturation_kpa(temperature_c)
    return min(max(share, 0.05), 0.99)


def solvent(values):
    """"water", "hydrocarbon" or `None`: the stable liquid life would use."""
    if values.get("hydrosphere") in ("surface ocean", "hycean", "ice-covered ocean"):
        return "water"
    temperature = values["surface_temperature_k"]
    methane = values.get("p_ch4_kpa") or 0.0
    if (methane > 0.0 and temperature <= METHANE_SEA_MAX_K
            and methane >= METHANE_SATURATION * atmosphere.vapour_pressure_kpa("ch4", temperature)):
        return "hydrocarbon"
    return None


def water_source(values):
    """What a base can draw water from (`habitability.WATER_SOURCE_FACTOR`)."""
    state = values.get("hydrosphere")
    if state in ("surface ocean", "hycean"):
        return "liquid"
    if state in ("ice", "ice-covered ocean"):
        return "ice"
    if state == "vapour":
        return "vapour"
    if (values.get("water_mass_fraction") or 0.0) >= HYDRATED_MIN_WATER_FRACTION:
        return "hydrated"
    return "vapour" if (values.get("p_h2o_kpa") or 0.0) > 0.0 else None


def air_inventories_kg(values):
    """`(carbon_kg, nitrogen_kg)` in the air, from the species' partial
    pressures; the balance gas is taken as 28 g/mol and carries neither."""
    pressure_kpa = (values.get("atmospheric_pressure_pa") or 0.0) / 1000.0
    gravity = (values.get("gravity_g") or 0.0) * constants.EARTH_GRAVITY
    radius_m = (values.get("radius_km") or 0.0) * 1000.0
    if pressure_kpa <= 0.0 or gravity <= 0.0 or radius_m <= 0.0:
        return 0.0, 0.0
    partials = {gas: values.get(f"p_{gas}_kpa") or 0.0 for gas in atmosphere.SPECIES}
    balance = max(0.0, pressure_kpa - sum(partials.values()))
    mean_molar = (sum(partials[gas] * MOLAR_KG[gas] for gas in partials) + balance * BALANCE_MOLAR_KG) / (
        sum(partials.values()) + balance)
    moles = pressure_kpa * 1000.0 * 4.0 * math.pi * radius_m ** 2 / (gravity * mean_molar)
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


def energy_flux_w_m2(values, luminosity_w, distance_au, internal_flux_w_m2):
    """The light or chemical power available to life, W/m2 (see the module
    docstring). `luminosity_w` and `distance_au` describe the star (either
    `None`: no light)."""
    chemical = CHEMICAL_SHARE * (internal_flux_w_m2 or 0.0)
    if values.get("hydrosphere") == "ice-covered ocean" or not luminosity_w or not distance_au:
        return chemical
    light = LIGHT_SHARE * luminosity_w / (4.0 * math.pi * (distance_au * constants.AU_TO_M) ** 2)
    return max(light, chemical)


def world_from(values, lethal_host=None):
    """The `habitability.World` of a rocky body's `values` (see
    `values_from_body`); `lethal_host` floors the dose it takes."""
    pressure_kpa = (values.get("atmospheric_pressure_pa") or 0.0) / 1000.0
    medium = solvent(values)
    under_ice = values.get("hydrosphere") == "ice-covered ocean"
    surface_dose = (values.get("surface_dose_msv_yr") or 0.0) / 1000.0
    habitat_dose = (values.get("dose_ground_msv_yr") or 0.0) / 1000.0 if under_ice else surface_dose
    if lethal_host:
        surface_dose = max(surface_dose, LETHAL_DOSE_SV_YR)
        habitat_dose = max(habitat_dose, LETHAL_DOSE_SV_YR)
    aqueous = medium == "water"
    activity = values.get("water_activity") if aqueous and values.get("water_activity") is not None else 1.0
    brine = values.get("ocean_class") == "chloride brine"
    gases = {}
    if pressure_kpa > 0.0:
        for gas, name in _GAS_NAMES.items():
            gases[name] = (values.get(f"p_{gas}_kpa") or 0.0) / pressure_kpa
    return habitability.World(
        pressure_kpa=pressure_kpa,
        temperature_c=values["surface_temperature_k"] - constants.CELSIUS_ZERO_K,
        gases=gases,
        gravity_ms2=(values.get("gravity_g") or 0.0) * constants.EARTH_GRAVITY,
        relative_humidity=relative_humidity(values),
        solvent=medium,
        water_activity=activity,
        chaotropicity_kj_kg=CHLORIDE_BRINE_CHAOTROPICITY * (1.0 - activity) if aqueous and brine else 0.0,
        ph=values.get("ocean_ph") if aqueous else None,
        inventories=inventories(values),
        phosphate_umol_l=(PHOSPHATE_UMOL_L.get(values.get("phosphorus"), 0.0) if aqueous
                          else NON_AQUEOUS_PHOSPHATE_UMOL_L if medium else 0.0),
        energy_flux_w_m2=values.get("energy_flux_w_m2") or 0.0,
        surface_dose_sv_yr=surface_dose,
        habitat_dose_sv_yr=habitat_dose,
        water_source=water_source(values),
    )


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
        return {name: None for name in SCORE_FIELDS}
    world = world_from(values, lethal_host=compact_host)
    result = habitability.scores(world)
    out = {"phi4": result["PHI-4"], "phi_bio": result["PHI_bio"], "phi_cpx": result["PHI_cpx"],
           "phi_tech": result["Phi_tech"], "l_solv": result["L_solv"], "l_chem": result["L_chem"],
           "l_ener": result["L_ener"], "l_rad": result["L_rad"]}
    limits = []
    for domain in DOMAINS:
        score, colour = result[domain.capitalize()]
        out[f"phi4_{domain}"] = score
        out[f"tier_{domain}"] = TIER_INDEX[colour]
        if colour != "Blue":
            limits.append(f"{domain.capitalize()} {colour}")
    out["equipment_tier"] = EQUIPMENT_NAMES.index(result["equipment"])
    if compact_host:
        out["hab_note"] = (f"Planet of a {compact_host}: lethal radiation, rated as {LETHAL_DOSE_SV_YR:,.0f} "
                           f"Sv/yr at the surface")
        out["equipment_tier"] = len(EQUIPMENT_NAMES) - 1
    else:
        out["hab_note"] = ", ".join(limits) if limits else "All four domains Blue"
    return out


def compact_host_name(star):
    """The `COMPACT_HOSTS` name of a star object, or `None`. A pulsar is
    called that."""
    kind = COMPACT_HOSTS.get(getattr(star, "yerkes_class", None))
    if kind == "neutron star" and getattr(star, "pulsar_type", "non-pulsing") != "non-pulsing":
        return "pulsar"
    return kind


def _clear(body):
    for name in BODY_FIELDS:
        setattr(body, name, None)


def generate(body, internal_flux_w_m2):
    """Sets `BODY_FIELDS` on a planet or moon. Call last, once its air,
    water and dose are set and `magnetism.update_exposure` has given its
    distance from its star; `internal_flux_w_m2` is its own heat flow."""
    _clear(body)
    if body.body_type == "g":
        return
    values = values_from_body(body)
    if values["surface_temperature_k"] is None:
        return
    star = getattr(body, "star", None)
    body.energy_flux_w_m2 = energy_flux_w_m2(
        values, getattr(star, "luminosity", None), getattr(body, "_star_distance_au", None), internal_flux_w_m2)
    values["energy_flux_w_m2"] = body.energy_flux_w_m2
    for name, value in compute(values, compact_host_name(star)).items():
        setattr(body, name, value)
