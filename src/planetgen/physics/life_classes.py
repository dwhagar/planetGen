# planetgen/physics/life_classes.py

"""
The life classes as a function of a finished body (GEN.886, docs/design/habitability/life-classes.md)
====================================================================================================

A life class is a way of making a living: what the life drinks (its solvent), eats (its energy), breathes and gives off. The
class number is the metabolism, a letter after it the pigment or habitat variant. `score_classes` reads the conditions a
finished planet or moon stores (its temperature and pressure, its air, its water, the light its star gives in each pigment's
band, its heat flow and dose) and gives every class a likelihood from 0 to 1 and a score from 0 to `tuning.LIFE_SCORE_MAX`:
"this body supports Class N life". It reads nothing the generator draws. The class switch (GEN.889, `generation/life.py`) draws
the class that lives on a world from these scores: the supported class with the largest net primary production (`ClassScore.npp`),
and stores the scores that reach the line (`life_class_scores`).

The rules are the rows of `physics/data/life_classes.json`, one per class, and every number in them carries its source tier
(T1 to T9) and tweak level (L0 to L4) beside it:

* `solvent`: the liquid the class needs. Water and a hydrocarbon sea are the index's own (`habitability_world.solvent`, scored by
  `habitability.l_solv`, so Class 1a on an Earth-like world is today's `PHI_bio`); liquid or supercritical carbon dioxide is read
  from the temperature and pressure; sulfuric acid and ammonia water are not modelled yet (GEN.888, GEN.887), so those classes
  score 0 and the sweep says which quantity they wait for;
* `energy`: starlight in the pigment's band (`band_fraction` of the star's photons over the Sun's share of 400 to 700 nm) or the
  chemical share of the body's own heat flow, on the index's log curve (`ENERGY_FLUX_MIN_W_M2` to `ENERGY_FLUX_REF_W_M2`);
* `gates`: conditions that must hold (`min`, `lt`, `equals`, `not_equals`, `in`, `not_in`) on named quantities;
* `windows`: a band a quantity should lie in; outside it the score falls as a half Gaussian of the rule's sigma;
* `supply`: the fuel and oxidant supplies the class needs (one entry or a list), each a share of Earth's from `physics/life_supply.py`
  (GEN.887), scored against the entry's `ref`; `placeholder` stands in only for a body whose values carry no such quantity;
* `life`: how the class lives: its `metabolism`, the `oxidant` it breathes and that oxidant's `oxidant_yield` (aerobic respiration is 1;
  an oxygen breather's yield follows the oxygen of the air, `physics/life_complexity.py`), the `first_step_scale` of its first
  light or chemistry (`physics/life_timeline.py`) and the colours its life shows.

The likelihood of a class is the index's Liebig blend (`habitability.LIEBIG_MIN_WEIGHT`) of the solvent, chemistry, energy and
dose terms, times the product of its windows and its supply. A likelihood under `tuning.HABITABILITY_LIFE_SCORE_MIN` is none.

Use:

    scores = score_classes(planet)          # {"1a": ClassScore, ..., "9": ClassScore}
    best(scores)                            # the class with the largest likelihood, or None
"""

import dataclasses
import functools
import json
import math
import os
from collections.abc import Mapping

from planetgen import tuning
from planetgen.physics import atmosphere, constants, elements, habitability, habitability_world, life_supply, strange_solvents, surface_liquid_world

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "life_classes.json")

TIERS = tuple(f"T{n}" for n in range(1, 10))
TWEAKS = tuple(f"L{n}" for n in range(5))
OPERATORS = ("min", "lt", "equals", "not_equals", "in", "not_in")
SOLVENT_MEDIA = ("water", "hydrocarbon", "co2", "sulfuric", "ammonia_water")
ENERGY_SOURCES = ("light", "chemical")

PLANCK_J_S = 6.62607015e-34
BOLTZMANN_J_K = 1.380649e-23
LIGHT_M_S = 299792458.0
SPECTRUM_NM = (100.0, 20000.0)
"""tuple: The wavelengths, nm, over which a star's photons are counted (a blackbody beyond them holds a negligible share)."""
SPECTRUM_STEPS = 600

CO2_TRIPLE_K, CO2_TRIPLE_PA = 216.59, 518500.0
CO2_CRITICAL_K, CO2_CRITICAL_PA = 304.13, 7377300.0
"""float: The triple and critical points of carbon dioxide (read, NIST)."""
OPTIONAL_QUANTITIES = ("acid_medium", "ammonia_water_mass_fraction", "silicon_rich_crust")
"""tuple: Quantities worked out from a body's values when it carries none (`derived_quantities`: GEN.888, GEN.894); a body that carries one in its
values is scored on that, and a rule that reads one the body lacks fails and names it."""
OXIDANTS = ("o2", "nitrate", "perchlorate", "sulfate", "co2", "acid", "none")
"""tuple: What a class breathes (`life.oxidant`); `"o2"` is the oxidant whose yield follows the air's oxygen."""
METABOLISMS = ("oxygenic", "anoxygenic", "hydrogenotroph", "chemolithotroph", "hydrogenator", "acidophile")
"""tuple: How a class makes its living (`life.metabolism`): only an `"oxygenic"` one splits water and so makes free oxygen."""
CO2_DENSE_SHARE = 0.5
"""float: The share of the air that carbon dioxide must be for its phase to be the surface fluid's (the same as the Class 6 rule)."""


BODY_FIELDS = ("life_class", "life_class_variant", "life_class_best", "life_class_best_score", "life_animal_grade", "life_complexity")
"""tuple: What the life pass stores on a body besides the schedule, the stage and the four scores (columns of `planets` and `moons`): the class that
lives on it (the number and the variant letter), the best-supported class and its score, the animal grade (0 to 3) and the complexity index."""

SCORES_FIELD = "life_class_scores"
"""str: The body attribute holding `{class id: score}` for every class that reaches the line (the child table `life_class_scores`)."""


class LifeClassError(ValueError):
    """A rule file that is refused, or a body that cannot be scored."""


class ClassScore:
    """One class's result for one body: its likelihood, its 0 to 9 score, the terms it is made of and what it lacked."""

    __slots__ = ("id", "likelihood", "score", "terms", "failed", "missing", "placeholders", "band_ratio", "npp")

    def __init__(self, class_id, likelihood, terms, failed, missing, placeholders, band_ratio=None, npp=0.0):
        self.id = class_id
        self.band_ratio = band_ratio
        self.npp = npp
        self.likelihood = likelihood
        self.score = score_of(likelihood)
        self.terms = terms
        self.failed = failed
        self.missing = missing
        self.placeholders = placeholders

    def __repr__(self):
        return f"ClassScore({self.id}, score={self.score}, likelihood={self.likelihood:.4f})"


def score_of(likelihood):
    """0 below the line `tuning.HABITABILITY_LIFE_SCORE_MIN`, else 1..`tuning.LIFE_SCORE_MAX` (the index's step, `life_scores`)."""
    if likelihood is None or likelihood < tuning.HABITABILITY_LIFE_SCORE_MIN:
        return 0
    return max(1, min(tuning.LIFE_SCORE_MAX, math.ceil(tuning.LIFE_SCORE_MAX * likelihood)))


# --- The star's light in a band --------------------------------------------------------------

def _photon_weight(wavelength_nm, temperature_k):
    """The photons per second per unit wavelength a blackbody emits, up to a constant: lambda^-4 / (exp(hc / lambda k T) - 1)."""
    wavelength_m = wavelength_nm * 1e-9
    x = PLANCK_J_S * LIGHT_M_S / (wavelength_m * BOLTZMANN_J_K * temperature_k)
    if x > 700.0:
        return 0.0
    return wavelength_m ** -4 / math.expm1(x)


def _integral(low_nm, high_nm, temperature_k, steps):
    width = (high_nm - low_nm) / steps
    return sum(_photon_weight(low_nm + (i + 0.5) * width, temperature_k) for i in range(steps)) * width


@functools.lru_cache(maxsize=4096)
def _band_fraction(temperature_k, low_nm, high_nm):
    total = _integral(SPECTRUM_NM[0], SPECTRUM_NM[1], temperature_k, SPECTRUM_STEPS * 4)
    if total <= 0.0:
        return 0.0
    steps = max(24, int((high_nm - low_nm) / 5.0))
    return min(1.0, _integral(low_nm, high_nm, temperature_k, steps) / total)


def band_fraction(temperature_k, band_nm):
    """The share of a blackbody star's photons that fall in the band `(low, high)` nm: the star's light a pigment can use. 0 for a
    body with no star temperature (or under half a kelvin, a black hole's, which rounds to none)."""
    rounded = round(float(temperature_k or 0.0), 0)
    if rounded <= 0.0:
        return 0.0
    return _band_fraction(rounded, float(band_nm[0]), float(band_nm[1]))


# --- The phase of carbon dioxide -------------------------------------------------------------

def co2_saturation_pa(temperature_k):
    """The vapour pressure of carbon dioxide, Pa, between its triple and critical points: ln p linear in 1/T through the two points
    (within a few percent of the measured curve, which is all a phase test needs)."""
    slope = math.log(CO2_CRITICAL_PA / CO2_TRIPLE_PA) / (1.0 / CO2_CRITICAL_K - 1.0 / CO2_TRIPLE_K)
    return CO2_TRIPLE_PA * math.exp(slope * (1.0 / temperature_k - 1.0 / CO2_TRIPLE_K))


def co2_phase(temperature_k, pressure_pa):
    """`"liquid"`, `"supercritical"`, `"gas"` or `"solid"`: the phase carbon dioxide takes at this temperature and pressure."""
    if temperature_k >= CO2_CRITICAL_K:
        return "supercritical" if pressure_pa >= CO2_CRITICAL_PA else "gas"
    if temperature_k < CO2_TRIPLE_K:
        return "solid"
    return "liquid" if pressure_pa >= co2_saturation_pa(temperature_k) else "gas"


# --- The rule file ---------------------------------------------------------------------------

class LifeClass:
    """One class: its words, its tier, and the rules that score it."""

    __slots__ = ("id", "number", "variant", "name", "drinks", "eats", "breathes", "gives_off", "reaches", "tier", "solvent",
                 "energy", "gates", "windows", "supplies", "life", "note")

    def __init__(self, row):
        for name in self.__slots__:
            setattr(self, name, row.get(name))
        self.gates = tuple(row.get("gates") or ())
        self.windows = tuple(row.get("windows") or ())
        self.supplies = supplies_of(row)
        self.note = row.get("note") or ""

    @property
    def metabolism(self):
        """How the class makes its living (`METABOLISMS`)."""
        return self.life["metabolism"]

    @property
    def oxidant(self):
        """What the class breathes (`OXIDANTS`)."""
        return self.life["oxidant"]

    @property
    def releases_oxygen(self):
        """Whether the class splits water for its electrons and so makes free oxygen (`biogenic_oxygen`)."""
        return self.life["metabolism"] == "oxygenic"

    @property
    def light(self):
        """Whether the class eats starlight (a producer that needs a sunlit surface and can make Flora)."""
        return self.energy["source"] == "light"


def supplies_of(row):
    """The supplies of a rule row as a tuple: its `supply` is one entry, a list of them or none."""
    supply = row.get("supply")
    if not supply:
        return ()
    return tuple(supply) if isinstance(supply, list) else (supply,)


def _check_standin(owner, entry):
    if entry.get("tier") not in TIERS or entry.get("tweak") not in TWEAKS or not entry.get("basis"):
        raise LifeClassError(f"{owner}: every number needs a tier T1..T9, a tweak L0..L4 and a basis, got {entry!r}")


def _check(row):
    owner = f"class {row.get('id')}"
    if row.get("tier") not in TIERS:
        raise LifeClassError(f"{owner}: no tier")
    solvent = row.get("solvent") or {}
    if solvent.get("medium") not in SOLVENT_MEDIA or not isinstance(solvent.get("factor"), (int, float)):
        raise LifeClassError(f"{owner}: the solvent is a medium from {SOLVENT_MEDIA} and a factor")
    _check_standin(f"{owner} solvent", solvent)
    energy = row.get("energy") or {}
    if energy.get("source") not in ENERGY_SOURCES:
        raise LifeClassError(f"{owner}: the energy source is one of {ENERGY_SOURCES}")
    if energy["source"] == "light" and not (len(energy.get("band_nm") or ()) == 2 and energy["band_nm"][0] < energy["band_nm"][1]
                                           and isinstance(energy.get("efficiency"), (int, float))):
        raise LifeClassError(f"{owner}: a light class has a band and an efficiency")
    _check_standin(f"{owner} energy", energy)
    for gate in row.get("gates") or ():
        operators = [op for op in OPERATORS if op in gate]
        if "quantity" not in gate or len(operators) != 1:
            raise LifeClassError(f"{owner}: a gate is a quantity and one of {OPERATORS}, got {gate!r}")
        _check_standin(f"{owner} gate {gate['quantity']}", gate)
    for window in row.get("windows") or ():
        keys = ("quantity", "low", "high", "sigma_low", "sigma_high")
        if any(key not in window for key in keys) or window["low"] >= window["high"] or min(window["sigma_low"], window["sigma_high"]) <= 0.0:
            raise LifeClassError(f"{owner}: a window has {keys} with low under high and sigmas over zero, got {window!r}")
        _check_standin(f"{owner} window {window['quantity']}", window)
    for supply in supplies_of(row):
        if "quantity" not in supply or not 0.0 <= supply.get("placeholder", -1.0) <= 1.0 or not supply.get("ref", 0.0) > 0.0:
            raise LifeClassError(f"{owner}: a supply is a quantity, a ref over 0 and a placeholder from 0 to 1, got {supply!r}")
        if "npp_per_unit" in supply and not supply["npp_per_unit"] > 0.0:
            raise LifeClassError(f"{owner}: npp_per_unit is over 0, got {supply!r}")
        _check_standin(f"{owner} supply {supply['quantity']}", supply)
    life = row.get("life") or {}
    if life.get("metabolism") not in METABOLISMS or life.get("oxidant") not in OXIDANTS:
        raise LifeClassError(f"{owner}: life has a metabolism from {METABOLISMS} and an oxidant from {OXIDANTS}, got {life!r}")
    if (life["oxidant"] == "o2") != (life.get("oxidant_yield") is None) or not (life.get("oxidant_yield") is None or 0.0 < life["oxidant_yield"] <= 1.0):
        raise LifeClassError(f"{owner}: an oxygen breather has no oxidant_yield (the air's oxygen sets it), any other oxidant has one over 0 and at most 1")
    if not life.get("first_step_scale", 0.0) > 0.0:
        raise LifeClassError(f"{owner}: life needs a first_step_scale over 0")
    _check_standin(f"{owner} life", life)


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    classes = []
    for row in document["classes"]:
        _check(row)
        classes.append(LifeClass(row))
    _check_standin("non_water_chemistry", document.get("non_water_chemistry") or {})
    ids = [c.id for c in classes]
    if len(set(ids)) != len(ids):
        raise LifeClassError("a class id is used by two rows")
    return document, tuple(classes)


DOCUMENT, CLASSES = _load()
"""tuple: The `LifeClass` rows in file order."""

BY_ID = {c.id: c for c in CLASSES}
IDS = tuple(BY_ID)
VERSION = DOCUMENT["version"]
REFERENCE_STAR_K = DOCUMENT["reference_star_k"]
REFERENCE_BAND_NM = tuple(DOCUMENT["reference_band_nm"])


def describe(class_id):
    """The `LifeClass` row of `class_id` (`"1a"`, `"4c"`, `"6"` ...)."""
    try:
        return BY_ID[class_id]
    except KeyError:
        raise LifeClassError(f"{class_id!r} is not a life class id") from None


# --- The quantities the rules read -----------------------------------------------------------

def values_from_body(body):
    """A planet or moon (in memory) as the column-named dict `score_classes` reads: the index's inputs, the heat flows, the redox state and the
    temperature of the star (`habitability_world.values_from_body`)."""
    return habitability_world.values_from_body(body)


def _earth_band_flux():
    return band_fraction(REFERENCE_STAR_K, REFERENCE_BAND_NM)


def quantities(values):
    """The named quantities the rules read, from a body's stored `values` (a column-named mapping). A quantity the body does not
    store is absent, so a rule that reads it fails and the class names it as missing."""
    pressure = values.get("atmospheric_pressure_pa") or 0.0
    found = {
        "temperature_k": values.get("surface_temperature_k"),
        "pressure_pa": pressure,
        "hydrosphere": values.get("hydrosphere") or "dry",
        "ocean_ph": values.get("ocean_ph"),
        "water_mass_fraction": values.get("water_mass_fraction") or 0.0,
    }
    for gas in atmosphere.SPECIES:
        partial_pa = atmosphere.partial_kpa(values, gas) * 1000.0
        found[f"share_{gas}"] = min(1.0, partial_pa / pressure) if pressure > 0.0 else 0.0
    found["co2_partial_pa"] = atmosphere.partial_kpa(values, "co2") * 1000.0
    habitat = life_supply.habitat_temperature_k(values)
    if habitat is not None:
        found["habitat_temperature_k"] = habitat
    temperature = found["temperature_k"]
    if temperature:
        phase = co2_phase(temperature, pressure) if found["share_co2"] >= CO2_DENSE_SHARE else "gas"
        found["co2_phase"] = phase
    teff, flux = values.get("star_teff_k"), values.get("energy_flux_w_m2")
    if teff and flux:
        found["starlight_w_m2"] = flux / habitability_world.LIGHT_SHARE
    if teff:
        found["_teff_k"] = teff
    for name in OPTIONAL_QUANTITIES:
        if values.get(name) is not None:
            found[name] = values[name]
    return found


def holds(conditions, found):
    """`(ok, failed, missing)`: whether every gate holds for `found`, the quantities of the gates that failed, and those the body
    does not store (a gate that reads one fails)."""
    failed, missing = [], []
    for gate in conditions:
        name = gate["quantity"]
        value = found.get(name)
        if value is None:
            missing.append(name)
            failed.append(name)
            continue
        ok = ((("min" not in gate) or value >= gate["min"]) and (("lt" not in gate) or value < gate["lt"])
              and (("equals" not in gate) or value == gate["equals"]) and (("not_equals" not in gate) or value != gate["not_equals"])
              and (("in" not in gate) or value in gate["in"]) and (("not_in" not in gate) or value not in gate["not_in"]))
        if not ok:
            failed.append(name)
    return not failed, failed, list(dict.fromkeys(missing))


# --- Scoring ---------------------------------------------------------------------------------

def _window(window, value):
    """1 inside the window, a half Gaussian of the rule's sigma outside it."""
    if value is None:
        return 0.0
    if value < window["low"]:
        return math.exp(-0.5 * ((window["low"] - value) / window["sigma_low"]) ** 2)
    if value > window["high"]:
        return math.exp(-0.5 * ((value - window["high"]) / window["sigma_high"]) ** 2)
    return 1.0


def _energy_curve(flux_w_m2):
    """The index's `L_ener` for a flux, W/m2 (log scale, 0 under the minimum, 1 from the reference)."""
    if not flux_w_m2 or flux_w_m2 <= habitability.ENERGY_FLUX_MIN_W_M2:
        return 0.0
    span = math.log10(habitability.ENERGY_FLUX_REF_W_M2 / habitability.ENERGY_FLUX_MIN_W_M2)
    return min(1.0, math.log10(flux_w_m2 / habitability.ENERGY_FLUX_MIN_W_M2) / span)


def _solvent_term(row, found, world, values):
    """`(L_solv, missing quantity or None)`: the solvent's term."""
    medium = row.solvent["medium"]
    if medium in ("water", "hydrocarbon"):
        return (habitability.l_solv(world) if world.solvent == medium else 0.0), None
    if medium == "co2":
        ok = found.get("co2_phase") in ("liquid", "supercritical")
        return (row.solvent["factor"] if ok else 0.0), None
    quantity, line = ("acid_medium", 1.0) if medium == "sulfuric" else ("ammonia_water_mass_fraction", 0.1)
    stored = found.get(quantity)
    if stored is None:
        return 0.0, quantity
    return (row.solvent["factor"] if stored >= line else 0.0), None


def _chemistry_term(row, world):
    """The index's `L_chem`; a class that does not drink water has its hydrogen from the supply term and a phosphate of the index's
    non-aqueous value (`DOCUMENT["non_water_chemistry"]`)."""
    if row.solvent["medium"] == "water":
        return habitability.l_chem(world)
    floor = DOCUMENT["non_water_chemistry"]["hydrogen_inventory_floor"]
    inventories = {**world.inventories, "H": max(world.inventories.get("H", 0.0), floor)}
    return habitability.l_chem(dataclasses.replace(
        world, inventories=inventories, phosphate_umol_l=habitability_world.NON_AQUEOUS_PHOSPHATE_UMOL_L))


def band_ratio(row, found):
    """A light class's energy against Class 1a's on the Sun: its band's share of the star's photons, times its efficiency, over the Sun's
    share of 400 to 700 nm, at most 1 (a body with no star temperature is taken as Sun-like). 1 for a chemical class."""
    if row.energy["source"] != "light":
        return 1.0
    teff = found.get("_teff_k")
    if not teff:
        return 1.0
    return min(1.0, band_fraction(teff, row.energy["band_nm"]) * row.energy["efficiency"] / _earth_band_flux())


def _energy_term(row, found, world, values):
    energy = row.energy
    if energy["source"] == "chemical":
        heat = (values.get("internal_heat_flux_w_m2") or 0.0) + (values.get("tidal_flux_w_m2") or 0.0)
        return _energy_curve(habitability_world.CHEMICAL_SHARE * heat)
    sunlit = values.get("hydrosphere") not in ("ice-covered ocean",) and world.life_on_surface
    if not sunlit:
        return 0.0
    teff = found.get("_teff_k")
    return _energy_curve(world.energy_flux_w_m2 * band_ratio(row, found))


def score_class(row, found, world, values):
    """The `ClassScore` of one class for a body."""
    local = found
    if row.energy["source"] == "light" and "starlight_w_m2" in found and found.get("_teff_k"):
        earth = constants.SOLAR_LUMINOSITY / (4.0 * math.pi * constants.AU_TO_M ** 2)
        local = dict(found)
        local["band_light_vs_earth"] = (found["starlight_w_m2"] * band_fraction(found["_teff_k"], row.energy["band_nm"])
                                        / (earth * _earth_band_flux()))
    ok, failed, missing = holds(row.gates, local)
    solv, missing_solvent = _solvent_term(row, found, world, values)
    if missing_solvent:
        missing = list(dict.fromkeys(missing + [missing_solvent]))
        failed = failed + [missing_solvent]
    elif solv <= 0.0:
        failed = failed + ["solvent"]
    terms = {"solvent": solv, "chemistry": _chemistry_term(row, world), "energy": _energy_term(row, found, world, values),
             "dose": habitability.l_rad(world)}
    blend = list(terms.values())
    lowest = min(blend)
    likelihood = 0.0 if lowest <= 0.0 else (lowest ** habitability.LIEBIG_MIN_WEIGHT
                                            * math.prod(blend) ** ((1.0 - habitability.LIEBIG_MIN_WEIGHT) / len(blend)))
    if terms["energy"] <= 0.0:
        failed = failed + ["energy"]
    windows_only = 1.0
    for window in row.windows:
        factor = _window(window, found.get(window["quantity"]))
        terms[f"window:{window['quantity']}"] = factor
        windows_only *= factor
    placeholders = []
    availability = 1.0
    productivity = 1.0
    for supply in row.supplies:
        stored = values.get(supply["quantity"])
        if stored is None:
            factor = productive = supply["placeholder"]
            placeholders.append(supply["quantity"])
        else:
            factor = min(1.0, stored / supply["ref"])
            productive = min(1.0, stored * supply["npp_per_unit"]) if "npp_per_unit" in supply else factor
        terms[f"supply:{supply['quantity']}"] = factor
        availability *= factor
        productivity *= productive
    likelihood = likelihood * windows_only * availability if ok else 0.0
    npp = net_primary_production(row, found, world, values) * productivity * windows_only if likelihood > 0.0 else 0.0
    return ClassScore(row.id, likelihood, terms, list(dict.fromkeys(failed)), missing, placeholders,
                      band_ratio(row, found) if row.energy["source"] == "light" else None, npp)


def earth_light_w_m2():
    """The usable light Earth's surface gets, W/m2 (`habitability_world.LIGHT_SHARE` of the solar constant): a producer's light of 1."""
    return habitability_world.LIGHT_SHARE * constants.SOLAR_LUMINOSITY / (4.0 * math.pi * constants.AU_TO_M ** 2)


def net_primary_production(row, found, world, values):
    """
    A class's net primary production before its supplies and windows, as a share of Earth's (1 is an oxygenic biosphere under the Sun):
    a light class's is the usable light of the surface times its band's share of the star's photons over the Sun's (`band_ratio`
    without its cap at 1, itself capped at `tuning.LIFE_NPP["light_ratio_max"]`), a chemical class's is the share of Earth's primary
    production that chemosynthesis gives (`tuning.LIFE_NPP["chemosynthesis_share"]`) times its heat flow in Earth's
    (`life_supply.heat_share`). Which class is the dominant one of a world is the one with the largest (GEN.889).
    """
    if row.energy["source"] == "chemical":
        return tuning.LIFE_NPP["chemosynthesis_share"] * life_supply.heat_share(values)
    if not (values.get("hydrosphere") != "ice-covered ocean" and world.life_on_surface):
        return 0.0
    teff = found.get("_teff_k")
    ratio = 1.0 if not teff else band_fraction(teff, row.energy["band_nm"]) * row.energy["efficiency"] / _earth_band_flux()
    return world.energy_flux_w_m2 / earth_light_w_m2() * min(ratio, tuning.LIFE_NPP["light_ratio_max"])


def derived_quantities(values):
    """`{name: value}` for each `OPTIONAL_QUANTITIES` that a body's `values` give no number of, worked out from its air and bulk composition
    (the stored acid layer and ammonia percent of the surface liquid, `surface_liquid_world`, GEN.1215; else `acid_medium.level` and `strange_solvents`); a quantity with nothing to work it out from (`None`) is left out."""
    found = {
        "acid_medium": surface_liquid_world.acid_level(values),
        "ammonia_water_mass_fraction": surface_liquid_world.ammonia_water_mass_fraction(values),
        "silicon_rich_crust": strange_solvents.silicon_rich_crust(values),
    }
    return {name: value for name, value in found.items() if value is not None and values.get(name) is None}


def score_classes(body):
    """`{class id: ClassScore}` for a finished planet or moon (an object, or a column-named mapping of its values). A gas giant, or a
    body with no surface temperature, is scored nowhere: `{}`."""
    values = dict(body) if isinstance(body, Mapping) else values_from_body(body)
    if values.get("body_type") == "g" or values.get("surface_temperature_k") is None:
        return {}
    for quantity, supplied in life_supply.supplies(values).items():
        values.setdefault(quantity, supplied)
    values.update(derived_quantities(values))
    found = quantities(values)
    world = habitability_world.world_from(values)
    return {row.id: score_class(row, found, world, values) for row in CLASSES}


def best(scores):
    """The `ClassScore` with the largest likelihood among those at or above the line, ties to the lower class number; `None` if none."""
    live = [s for s in scores.values() if s.score > 0]
    return max(live, key=lambda s: (s.likelihood, -IDS.index(s.id))) if live else None


def best_npp(scores, among=None):
    """The `ClassScore`s that are supported (a score of 1 or more), largest net primary production first; those of the ids `among` only."""
    live = [s for s in scores.values() if s.score > 0 and (among is None or s.id in among)]
    return sorted(live, key=lambda s: (-s.npp, IDS.index(s.id)))


def dominant(scores, draw_stream):
    """
    The class that lives on a world (GEN.889, life-classes.md 8.2): of the supported classes, the one with the largest net primary
    production; classes within `tuning.LIFE_NPP["tie_factor"]` of it are tied, and a tie is broken by a draw weighted by net primary
    production (`draw_stream` is anything with `random()`: a `draw.Stream`, or the `draw` module). `None` when no class is supported.
    A single candidate costs no draw.
    """
    ranked = best_npp(scores)
    if not ranked or ranked[0].npp <= 0.0:
        return None
    tied = [s for s in ranked if s.npp * tuning.LIFE_NPP["tie_factor"] >= ranked[0].npp]
    if len(tied) == 1:
        return tied[0]
    pick = draw_stream.random() * sum(s.npp for s in tied)
    for candidate in tied:
        pick -= candidate.npp
        if pick < 0.0:
            return candidate
    return tied[-1]


def class_id(number, variant=None):
    """The id of a class from its stored number and variant letter (`4, "c"` gives `"4c"`; `6, None` gives `"6"`), `None` for no class."""
    if number is None:
        return None
    return f"{int(number)}{variant or ''}"


def living_class_id(number, variant, stage):
    """The id of the class that lives on a body, from its stored class number, variant letter and life stage: `None` for no class and for a class
    with no stage (GEN.1001: a world whose star is younger than the first step of its schedule suits a class but has no life yet)."""
    return None if stage is None else class_id(number, variant)
