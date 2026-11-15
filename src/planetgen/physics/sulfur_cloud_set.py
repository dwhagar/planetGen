# planetgen/physics/sulfur_cloud_set.py

"""
The sulfur-cloud chemistry a planet's air can hold, and the six gates that give it (GEN.1425, section 10 of the unseen-compounds report)
=========================================================================================================================================

The one acid cloud we have measured, Venus's, holds more than the acid: sulfur dioxide, the lower sulfur oxides, sulfur aerosol, chlorine and iron
compounds and an ultraviolet absorber nobody has identified. This module is the set (`physics/data/sulfur_cloud_set.json`, rows V1 to V10 with their
tier, tweak and evidence status) and the rule that says which of it an air gets. Nothing is stored and nothing names Venus: `chemistry(values)` reads the
values a body already has and answers from its air, its star's ultraviolet and its surface. A cold, sulfur-poor world gets none of it; a hot volcanic
carbon-dioxide world under a young active star gets V1 to V5.

The gates, in order (each later one needs the one before it):

1. **Sulfur supplied, in an oxidizing air** (V2): sulfur dioxide, an oxidiser (oxygen plus carbon dioxide) and a trace of water, the floors of `acid_medium`.
2. **An acid cloud** (V1): the air has `cloud_water_min_share` of water and, at some level of its dry adiabat, the acid vapor it holds
   (`acid_vapor_share_of_so2` of its sulfur dioxide) is above the pressure over the solution its water fixes (`sulfuric_vapor.condenses`, GEN.1424). The base is
   the warmest such level and the strength is that solution's.
3. **Reduced sulfur** (V4, V5): the star's ultraviolet-C flux at the body (`uv_c_min_w_m2`) photolyses the sulfur dioxide and carbon monoxide or carbonyl sulfide
   (`reductant_min_share`) is the reductant; V5 needs more sulfur (`sulfur_aerosol_min_so2_share`).
4. **Chlorine** (V3 hydrogen chloride; V6, V9): hydrogen chloride in the air at `chlorine_min_share` (the star's chlorine got into the air through
   `mixture_rocky`'s halogen budget, so the star is not read twice).
5. **Iron and carbon supply** (V6, V7, V8): iron from an exposed rocky surface (a terrestrial body with no global ocean over it); carbon from cosmic dust, whose influx on the body
   (`influx.py` over the belts and comets of its own system, GEN.1543: `generation/dust_supply.py` keeps it on the body) must reach `dust_min_influx_kg_m2_yr` per area (`dust_supply`; the caller of
   `chemistry` says whether it does with `carbon_dust`, default none) and V8 stays closed without it.
6. **Colour** (V10): the measured absorber band, with the acid cloud and sulfur dioxide, replaced by V8 or V9 if either is on (their published spectra
   would be used, not added).

V8 and V9 are single-source hypotheses and are off by default (`include_off_by_default`).

The functions are pure. The numbers are in `physics/data/sulfur_cloud_set.json` with their tier and tweak level.
"""

import json
import math
import os
from dataclasses import dataclass, field
from statistics import NormalDist

from planetgen.physics import acid_medium, aerosols, air_column, atmosphere, constants, elements, sulfuric_vapor

DATA_FILE = os.path.join(elements.DATA_DIRECTORY, "sulfur_cloud_set.json")

INCLUDE_HYPOTHESES = False
"""bool: Whether the chain may give a body V8 and V9, the single-source hypotheses (iron-PAH from cosmic dust, sulfur-chlorine products). Off by default (Boss's set, 11-08);
turn it on to see them. V8 also needs the body's cosmic-dust influx from its own system's belts and comets (`dust_supply`, GEN.1543)."""

GLOBAL_OCEAN = ("surface ocean", "hycean")
"""tuple: The hydrosphere values that bury a body's rock under liquid, so no iron is exposed to the air."""


class SulfurCloudSetError(ValueError):
    """A rule file that is refused."""


def _check(name, entry):
    if entry.get("tier") not in tuple(f"T{n}" for n in range(1, 10)) or entry.get("tweak") not in tuple(f"L{n}" for n in range(5)):
        raise SulfurCloudSetError(f"{name}: every entry needs a tier T1..T9 and a tweak L0..L4, got {entry!r}")


def _load():
    with open(DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    for name, entry in document["numbers"].items():
        _check(name, entry)
        if not entry.get("basis") or "value" not in entry:
            raise SulfurCloudSetError(f"{name}: a value and a basis are needed")
    seen = []
    for row in document["species"]:
        _check(row.get("id"), row)
        if row.get("default") not in ("on", "off") or row.get("source") not in ("air", "chain") or not row.get("status") or not row.get("evidence"):
            raise SulfurCloudSetError(f"{row.get('id')}: a default (on or off), a source (air or chain), a status and an evidence are needed")
        seen.append(row["id"])
    if seen != [f"V{n}" for n in range(1, 11)]:
        raise SulfurCloudSetError(f"the set is V1 to V10 in order, got {seen}")
    return document


DOCUMENT = _load()
NUMBERS = DOCUMENT["numbers"]
SPECIES = {row["id"]: row for row in DOCUMENT["species"]}
VERSION = DOCUMENT["version"]


def number(name):
    """The value of the rule file's number `name`."""
    return NUMBERS[name]["value"]


@dataclass(frozen=True)
class Cloud:
    """The acid cloud of an air: the base (the warmest level where the liquid forms), the top (the coolest the profile reaches) and the base's strength."""

    base_temperature_k: float
    base_pressure_pa: float
    top_pressure_pa: float
    weight_percent: float


@dataclass(frozen=True)
class Chemistry:
    """What `chemistry` found: the species ids that are present, the gates that opened, the acid cloud and the reason each absent species is absent."""

    species: tuple
    gates: tuple
    cloud: object = None
    absent: dict = field(default_factory=dict)

    def has(self, species_id):
        return species_id in self.species


def _shares(values):
    gas = atmosphere.gas_partials_pa(values)
    total = sum(gas.values())
    if not gas or total <= 0.0:
        return None, {}
    return gas, {formula: pressure / total for formula, pressure in gas.items()}


def acid_cloud_of_column(column, temperature_k, skin=None):
    """
    The `Cloud` of an `air_column.AirColumn` (the ground at `temperature_k`, the air turning isothermal at `skin`, else the grey factor of the ground), or `None`:
    gate 2. The water at a level is the surface mixing ratio, never more than that level's saturation (a cold trap), the acid vapor `acid_vapor_share_of_so2`
    of the sulfur dioxide at every level, and the liquid is where `sulfuric_vapor.condenses`. The base is the warmest such level.
    """
    water, acid = column.fraction("H2O"), column.fraction("SO2") * number("acid_vapor_share_of_so2")
    if water < number("cloud_water_min_share") or acid <= 0.0:
        return None
    pressure, level_k = aerosols.profile(column, temperature_k, skin_k=skin)
    window = []
    for p, t in zip(pressure, level_k):
        saturation = aerosols.saturation_pa("H2O", float(t))
        water_pa = water * float(p)
        if math.isfinite(saturation):
            water_pa = min(water_pa, saturation)
        if sulfuric_vapor.condenses(float(t), acid * float(p), water_pa):
            window.append((float(p), float(t), water_pa))
    if not window:
        return None
    base = max(window, key=lambda level: level[1])
    solution = sulfuric_vapor.composition(base[1], base[2])
    return Cloud(base[1], base[0], min(level[0] for level in window), solution.weight_percent)


def acid_cloud(values, gas=None):
    """The `Cloud` of a body's stored `values` (see `acid_cloud_of_column`; the skin is `acid_medium.skin_k`), or `None`: gate 2."""
    if gas is None:
        gas = atmosphere.gas_partials_pa(values)
    temperature = values.get("surface_temperature_k")
    gravity = (values.get("gravity_g") or 0.0) * constants.EARTH_GRAVITY
    if not gas or not temperature or gravity <= 0.0 or sum(gas.values()) <= 0.0:
        return None
    column = air_column.from_partials(air_column.partials_pa(values), gravity)
    return acid_cloud_of_column(column, temperature, acid_medium.skin_k(values, temperature))


def uv_c_flux_w_m2(values):
    """The star's ultraviolet-C flux at the body, W/m2: the energy flux it absorbs times the star's `uv_c_share`; `None` when either is unknown."""
    flux, share = values.get("energy_flux_w_m2"), values.get("uv_c_share")
    if not flux or share is None:
        return None
    return flux * share


_UNSET = object()
"""object: The default of `chemistry`'s `cloud`: the cloud is worked out from the values."""


def chemistry(values, uv_c_w_m2=None, carbon_dust=False, include_off_by_default=False, cloud=_UNSET):
    """
    The `Chemistry` of a body's stored `values` (a column-named mapping with its gases: `atmosphere.GAS_FIELD`). `uv_c_w_m2` is the star's ultraviolet-C flux at the
    body (else `uv_c_flux_w_m2(values)`), `carbon_dust` whether cosmic dust reaches it, `include_off_by_default` whether V8 and V9 may be given. `cloud` is the acid
    cloud the body's aerosols were drawn with (`None` for none), so gate 2 and the stored optical depth agree; left out, it is worked out from the values.
    """
    gas, share = _shares(values)
    absent = {}
    if gas is None:
        return Chemistry((), (), None, {species_id: "no air" for species_id in SPECIES})
    uv_c = uv_c_flux_w_m2(values) if uv_c_w_m2 is None else uv_c_w_m2
    sulfur = share.get("SO2", 0.0)
    oxidant = share.get("O2", 0.0) + share.get("CO2", 0.0)
    gate1 = sulfur >= acid_medium.number("so2_min_share") and oxidant >= acid_medium.number("oxidant_min_share") \
        and share.get("H2O", 0.0) >= acid_medium.number("water_min_share")
    if cloud is _UNSET:
        cloud = acid_cloud(values, gas) if gate1 else None
    elif not gate1:
        cloud = None
    gate2 = cloud is not None
    reductant = share.get("CO", 0.0) + share.get("OCS", 0.0)
    gate3 = gate1 and uv_c is not None and uv_c >= number("uv_c_min_w_m2") and reductant >= number("reductant_min_share")
    gate4 = share.get("HCl", 0.0) >= number("chlorine_min_share")
    iron = values.get("body_type") == "t" and values.get("hydrosphere") not in GLOBAL_OCEAN
    gate5_iron, gate5_carbon = iron, bool(carbon_dust)
    gate6 = gate2
    needs = {
        "V1": (gate2, "no acid cloud (gate 2)"),
        "V2": (gate1, "no sulfur dioxide with an oxidiser and water (gate 1)"),
        "V3": (gate4 or share.get("HF", 0.0) >= number("chlorine_min_share") or reductant >= number("reductant_min_share"), "none of hydrogen chloride, hydrogen fluoride, carbonyl sulfide or carbon monoxide above its floor in the air"),
        "V4": (gate3, "no sulfur dioxide photolysis with a reductant (gate 3)"),
        "V5": (gate3 and sulfur >= number("sulfur_aerosol_min_so2_share"), "too little sulfur or no photolysis with a reductant (gate 3)"),
        "V6": (gate1 and gate2 and gate4 and gate5_iron, "needs the acid cloud, chlorine and exposed iron (gates 2, 4, 5)"),
        "V7": (gate1 and gate2 and gate5_iron, "needs the acid cloud and exposed iron (gates 2, 5)"),
        "V8": (gate2 and gate5_iron and gate5_carbon, "needs the acid cloud, iron and a cosmic dust supply (gates 2, 5)"),
        "V9": (gate3 and gate4, "needs photolysed sulfur and chlorine (gates 3, 4)"),
        "V10": (gate6, "needs the acid cloud with sulfur dioxide (gate 6)"),
    }
    present = []
    for species_id, (open_, why) in needs.items():
        if not open_:
            absent[species_id] = why
        elif SPECIES[species_id]["default"] == "off" and not include_off_by_default:
            absent[species_id] = "off by default (single-source hypothesis)"
        else:
            present.append(species_id)
    if ("V8" in present or "V9" in present) and "V10" in present:
        present.remove("V10")
        absent["V10"] = "replaced by the published spectrum of V8 or V9, not added to it"
    gates = tuple(i + 1 for i, opened in enumerate((gate1, gate2, gate3, gate4, gate5_iron or gate5_carbon, gate6)) if opened)
    return Chemistry(tuple(present), gates, cloud, absent)


def dust_supply(influx_kg_yr, radius_m):
    """Whether the cosmic dust that reaches a body (`influx_kg_yr`, kg a year, from the belts and comets of its own system: `generation/dust_supply`) is at least `dust_min_influx_kg_m2_yr` per m2 of its surface (gate 5, V8; GEN.1543)."""
    if not influx_kg_yr or not radius_m or radius_m <= 0.0:
        return False
    return influx_kg_yr / (4.0 * math.pi * radius_m ** 2) >= number("dust_min_influx_kg_m2_yr")


COLUMN = "cloud_chemistry"
"""str: The column of `planets` and `moons` (and `rogue_planets`, which stays empty: no star, no light) that holds `encode(chemistry)` (GEN.1425)."""


def encode(found):
    """The text a body stores for a `Chemistry`: its species ids in order, comma separated (`V1,V2,V3,V10`); `None` when it has none."""
    return ",".join(found.species) if found.species else None


def species_of(stored):
    """The species ids of a stored text (`encode`), in order; an empty tuple for `None`."""
    return tuple(part for part in (stored or "").split(",") if part)


V3_GASES = (("HCl", "hydrogen chloride"), ("HF", "hydrogen fluoride"), ("OCS", "carbonyl sulfide"), ("CO", "carbon monoxide"))
"""tuple: The gases V3 names, `(formula, words)` (the V3 gate opens when any one is above its floor, so a page names only those that are, UX.309)."""


def _v3_name(air):
    """V3's name for an air (`air`: gas-table formula to partial pressure): the gases of `V3_GASES` above their floor, else the whole set's name."""
    total = sum(air.values())
    floor = number("chlorine_min_share") if total > 0.0 else None
    if floor is None:
        return None
    held = [words for formula, words in V3_GASES if air.get(formula, 0.0) / total >= (floor if formula in ("HCl", "HF") else number("reductant_min_share"))]
    if not held:
        return None
    return (held[0] if len(held) == 1 else ", ".join(held[:-1]) + " and " + held[-1]) + (" gas" if len(held) == 1 else " gases")


def sentence(stored, air=None):
    """The stored text as words for a page: the names of the species the chain gave the air, with the ones that are only candidates marked, e.g. `sulfuric acid water droplets: observed in a planet (Venus); sulfur dioxide gas: observed in a planet (Venus)`.
    Given the body's `air` (gas-table formula to Pa), V3 names only the gases that are above their floor (UX.309), not the four of the set."""
    names = []
    for species_id in species_of(stored):
        row = SPECIES.get(species_id)
        if row is not None:
            name = _v3_name(air) if species_id == "V3" and air else None
            names.append(f"{name or row['name']}: {row['status']}")
    return "; ".join(names)


def fill(body, values):
    """
    Sets `cloud_chemistry` on a rocky planet or moon from its finished air, light and ground (`values` is `habitability_world.values_from_body`, with the energy flux set): the
    species of the Venus set the chain gives it (GEN.1425). Gate 2 is the acid cloud the body's aerosols were drawn with (`sulfuric_cloud.window` at the aerosols' skin
    temperature), the ultraviolet-C flux is the star's `uv_c_share` of the absorbed flux, the cosmic-dust supply is `dust_supply` of the influx its system's belts and comets send it (`generation/dust_supply`), and V8 and V9
    (single-source hypotheses, off by default) follow `INCLUDE_HYPOTHESES`. A body with no air stores `None`.
    """
    body.cloud_chemistry = None
    from planetgen.physics import sulfuric_cloud       # imports this module
    column = air_column.column_of(body)
    temperature = getattr(body, "surface_temperature", None)
    if column is None or not temperature:
        return
    star = getattr(body, "star", None)
    share = getattr(star, "uv_c_share", None)
    flux = values.get("energy_flux_w_m2")
    uv_c = flux * share if flux and share is not None else None
    cloud = sulfuric_cloud.window(column, temperature, aerosols.skin_of(body))
    radius_km = getattr(body, "radius", None)
    dust = bool(radius_km) and dust_supply(getattr(body, "_dust_influx_kg_yr", None), radius_km * 1000.0)       # GEN.1543: set by `generation/dust_supply` once the belts and comets are in
    body.cloud_chemistry = encode(chemistry(values, uv_c_w_m2=uv_c, carbon_dust=dust, include_off_by_default=INCLUDE_HYPOTHESES, cloud=cloud))


def absorber_depth(cloud_tau_550, so2_share, wavelength_um, center_um=None, fwhm_um=None):
    """
    The absorption optical depth of V10 at `wavelength_um`: a Gaussian band (center `absorber_center_um`, full width `absorber_fwhm_um`, or the `center_um`
    and `fwhm_um` given) whose peak is `absorber_depth_per_cloud_depth` of the cloud's depth at 550 nm, scaled down in proportion to the sulfur dioxide below
    `absorber_so2_reference_share`.
    """
    center = number("absorber_center_um") if center_um is None else center_um
    fwhm = number("absorber_fwhm_um") if fwhm_um is None else fwhm_um
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    strength = number("absorber_depth_per_cloud_depth") * cloud_tau_550 * min(1.0, so2_share / number("absorber_so2_reference_share"))
    return strength * math.exp(-0.5 * ((wavelength_um - center) / sigma) ** 2)


def absorber_spread_um():
    """`(center_sigma, fwhm_sigma)`, the stated uncertainties of the band, micron, for a caller that varies the band inside them."""
    return number("absorber_center_sigma_um"), number("absorber_fwhm_sigma_um")


GOLDEN_CENTER, GOLDEN_WIDTH = 1.6180339887498949, 2.7182818284590452
"""float: The multipliers that make the band's center and width from the one stored `aerosol_draw` (fractional part of `aerosol_draw` times each), so the two
are not the same number and no column is added."""
SPREAD_LIMIT_SIGMA = 2.0
"""float: How far, in stated sigmas, the draw may move the band (the tails past two sigma are not used)."""


_SPREAD_LOW = NormalDist().cdf(-SPREAD_LIMIT_SIGMA)
"""float: The unit below which (and above one minus which) `_shift_sigmas` holds the draw."""


def _shift_sigmas(unit):
    """A normal deviate in sigmas for `unit` in [0, 1), held to `SPREAD_LIMIT_SIGMA` either way."""
    return NormalDist().inv_cdf(min(max(unit, _SPREAD_LOW), 1.0 - _SPREAD_LOW))


def absorber_band(cloud_tau_550, so2_share, wavelengths_um, draw=0.5):
    """
    `absorber_depth` at each of `wavelengths_um`, a list: the absorption optical depth V10 puts in the acid cloud on those wavelengths. The band's center and
    width move inside their stated uncertainties (`absorber_spread_um`) by the body's stored `draw` (0.5 is the published band).
    """
    center_sigma, fwhm_sigma = absorber_spread_um()
    center = number("absorber_center_um") + center_sigma * _shift_sigmas((draw * GOLDEN_CENTER) % 1.0)
    fwhm = max(0.0, number("absorber_fwhm_um") + fwhm_sigma * _shift_sigmas((draw * GOLDEN_WIDTH) % 1.0))
    return [absorber_depth(cloud_tau_550, so2_share, wavelength, center, fwhm) for wavelength in wavelengths_um]
