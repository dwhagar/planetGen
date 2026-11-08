# planetgen/physics/isotope_record.py

"""
The isotope record of a body (GEN.1377 and GEN.1420, part 1)
============================================================

All real materials exist: a planet, moon or rogue planet carries, for every element it holds, the share of each of that element's isotopes (an atom fraction, the shares of one element add up to 1).
Deuterium is a share of hydrogen, not a fictional substance. The record is a list of `ROW_FIELDS` dicts on the body (`ATTRIBUTE`), stored in the child tables `planet_isotopes`, `moon_isotopes`
and `rogue_planet_isotopes`.

**Which elements a body holds** (`elements_held`): those of the parts it is made of. The rock gives O, Mg, Si, Al, Ca and every minor oxide with a share above zero (Na, K, Ti, Cr, Mn, Ni, P), Fe
from its FeO and metal, and U and Th where the body stores a radiogenic content; the core gives Fe, Ni and each light element with a share (S, Si, O, C, H); water gives H and O; each ice (CO2, CH4, NH3, N2,
SO2, CO, HCN, the rest as organics) gives its own elements; an envelope gives H and He and, past a heavy-element share above zero, C, N and O; the air gives the elements of every gas it holds. A carbon
rock gives C and Si in place of the oxides.

**Which isotopes** of an element: every isotope of the table with a natural abundance (`elements`, NIST and IUPAC, [R2768] for the half-lives). Nothing real is dropped. The shares default to the
table's natural mix (origin `natural`) and are replaced where a physical value is known:

* **Hydrogen** (`physical`): the deuterium share is the hydrogen-atom-weighted mean of the body's reservoirs. The gas of the star's birth cloud has the primordial D/H (2.527e-5, [R5165]) reduced by the star formation that burned deuterium since: `exp(-ASTRATION * 10^[Fe/H])` (the metal-richer the star, the more of the gas has been through stars). Every ice and the water
  are enriched over the gas by `ice_factor`, a law in the formation temperature `exp(T0 (1/T - 1/T_ref))` (ion-molecule fractionation grows as the gas gets colder, [R5168: Sec. 2.1]) that is
  calibrated at two published points of the Solar System: the ocean (1.56e-4, VSMOW) at the snow-line temperature and comet 67P's water (5.3e-4, [R0088]) at 30 K. A body that formed warmer than
  the snow line took its ice from the snow line's temperature. One log-normal scatter per body.
* **Carbon** (`physical`): 13C/12C from the Galaxy's chemical evolution: 12C/13C = 6.21 D + 18.71 ([R5169: equation 5], D in kpc from the Galactic centre) read at the galactocentric distance
  whose gas metallicity is the star's (`nebula_composition.GAS_GRADIENT_DEX_PER_KPC`), and moved in time between the local interstellar value (68) and the Solar System's (89) by the star's age.
* **Nitrogen** (`physical`): 14N/15N the same way ([R5171] as quoted in [R5166: equation 2]: 21.1 D + 123.8), between the local interstellar value (about 400) and the
  protosolar one (441, solar wind, quoted in [R5166]), then raised in 15N by `ice_factor` for the solid share of the body (comets are enriched 2 to 3 times, [R5168]). An envelope's nitrogen is the gas value;
  the two are mixed by the envelope's share of the body's mass.
* **Argon** (`physical`, `radiogenic`): the primordial 36Ar/38Ar is the solar wind's 5.501 ([R5170]); 40Ar is the share of the air's argon that the body's potassium made and its interior
  let out: the 40K that has decayed (`radiogenic.formation_mix`, the body's age), times the 10.72 percent of decays that end in 40Ar, times the share of the mantle's easily degassed element
  that the interior's released heat has brought to the surface (`surface_pressure.degassed_share`), over the moles of argon in the air; it cannot be more than all but the primordial share.
* **Potassium, uranium** (`decay`): the isotopes that decay are those of the body's age, from the formation mix of `radiogenic` (today's natural mix decayed back to the Solar System's birth, a stand-in
  for every star) decayed over the star's age. 234U and the daughters (the lead, strontium and neodymium isotopes) stay at the natural mix.

Left out, and where it goes: the short-lived isotopes (Fe-60, Mn-53, I-129, Pu-244) and the cosmogenic ones (tritium, C-14, Be-10, Al-26) are the second part (GEN.1377, GEN.1421); the xenon excess
from I-129 and Pu-244 and the mass-independent oxygen anomalies (O stays at the natural mix); the isotope split of a body's reservoirs (a mantle against an ocean); belts, comets and stars.
"""

import math

from planetgen import tuning
from planetgen.physics import bulk_draw, elements, interior_heat, nebula_composition, radiogenic, surface_pressure
from planetgen.util import draw

ATTRIBUTE = "isotopes"
"""str: The body attribute that holds the record, a list of `ROW_FIELDS` dicts (the child tables `planet_isotopes`, `moon_isotopes` and `rogue_planet_isotopes`)."""

ROW_FIELDS = ("element", "isotope", "fraction", "origin")
"""tuple: The columns of an isotope row after the owner id: the element's symbol, the isotope's name (`H-2`), its atom fraction of the element, and where the value comes from."""

NATURAL, PHYSICAL, DECAY, RADIOGENIC = "natural", "physical", "decay", "radiogenic"
ORIGINS = (NATURAL, PHYSICAL, DECAY, RADIOGENIC)
"""tuple: The `origin` of a row: the table's natural mix, a physical law of the star and the formation temperature, a radioactive isotope at the body's age, a daughter the body's own decays made."""

# --- the laws (stand-ins IS1 to IS8) -------------------------------------------------------------------------------------
PRIMORDIAL_D_H = 2.527e-5
"""float: The primordial deuterium over hydrogen by number ([R5165]). [Stand-in ISO1, tier T1, tweak L0.]"""

PROTOSOLAR_D_H = 2.1e-5
"""float: The deuterium over hydrogen of the Sun's birth gas, used only to set `ASTRATION` ([R5168: Sec. 4.2], quoting the Jovian and solar-wind values). [Stand-in ISO2, tier T2, tweak L0.]"""

ASTRATION = math.log(PRIMORDIAL_D_H / PROTOSOLAR_D_H)
"""float: The fall of ln D/H from the primordial value to the Sun's birth gas, at the Sun's metallicity; the fall is taken proportional to the metals, `ASTRATION * 10^[Fe/H]` (stars burn deuterium as they make metals). [Stand-in ISO2, tier T5, tweak L1.]"""

ICE_REFERENCE_K = 182.0
"""float: The temperature of the snow line's water, the highest the enrichment of an ice is taken at (a warmer body took its ice from there). [Stand-in ISO3, tier T1, tweak L0.]"""

COLD_ICE_K = 30.0
"""float: The temperature of the cold calibration point, comet 67P's water (the CO2 and CH4 line zone, below which the forward H3+ + HD reaction wins, [R5168]). [Stand-in ISO3, tier T2, tweak L0.]"""

FLOOR_ICE_K = 20.0
"""float: The coldest temperature the law is read at: below it the enrichment would run past the calibration. [Stand-in ISO4, tier T4, tweak L3.]"""

OCEAN_D_H = 1.56e-4
"""float: The deuterium over hydrogen of ocean water (VSMOW, [R5168: Sec. 1.2]), the enrichment at `ICE_REFERENCE_K`. [Stand-in ISO3, tier T1, tweak L0.]"""

COMET_D_H = 5.3e-4
"""float: The deuterium over hydrogen of the water of a Jupiter-family comet ((5.3 +- 0.7)e-4, [R0088], quoted in [R5167]), the enrichment at `COLD_ICE_K`. [Stand-in ISO3, tier T1, tweak L0.]"""

D_H_ENRICHMENT_REFERENCE = OCEAN_D_H / PROTOSOLAR_D_H
D_H_ENRICHMENT_COLD = COMET_D_H / PROTOSOLAR_D_H
D_H_T0_K = math.log(D_H_ENRICHMENT_COLD / D_H_ENRICHMENT_REFERENCE) / (1.0 / COLD_ICE_K - 1.0 / ICE_REFERENCE_K)
"""float: The temperature scale of the deuterium enrichment of ice, K (about 44): the law is fitted to the two published points, so it is tuned to published statistics. [Stand-in ISO3, tier T5, tweak L2.]"""

D_H_SCATTER_LN = 0.3
"""float: The log-normal scatter of one body's ice D/H about the law (comets span a factor of about 2 to 3 around their mean, [R5167]). [Stand-in ISO4, tier T4, tweak L3.]"""

CARBON_GRADIENT = (6.21, 18.71)
"""tuple: 12C/13C against galactocentric distance in kpc, `(slope, intercept)` ([R5169: equation 5], fitted to CN, CO and H2CO). [Stand-in ISO5, tier T1, tweak L0.]"""

NITROGEN_GRADIENT = (21.1, 123.8)
"""tuple: 14N/15N against galactocentric distance in kpc ([R5171], equation 2 of [R5166]). [Stand-in ISO5, tier T1, tweak L0.]"""

CARBON_RATIOS = (68.0, 89.0)
"""tuple: `(local interstellar now, Solar System birth)` 12C/13C ([R5169: Abstract]). The star's age moves the ratio between them. [Stand-in ISO5, tier T1, tweak L0.]"""

NITROGEN_RATIOS = (400.0, 441.0)
"""tuple: `(local interstellar now, protosolar)` 14N/15N ([R5166: Abstract and Sec. 3.3]; the solar wind). [Stand-in ISO5, tier T1, tweak L0.]"""

SOLAR_GALACTOCENTRIC_KPC = 7.9
"""float: The Sun's distance from the Galactic centre, kpc ([R5169]), where the two gradients read the local interstellar ratios. [Stand-in ISO5, tier T1, tweak L0.]"""

GALACTOCENTRIC_RANGE_KPC = (2.0, 20.0)
"""tuple: The range the gradients are read in (the data of [R5169] reach from 2 to 12 kpc; the inner and outer Galaxy past it hold no data). [Stand-in ISO4, tier T4, tweak L3.]"""

NITROGEN_ENRICHMENT_REFERENCE = 441.0 / 274.0
NITROGEN_ENRICHMENT_COLD = 441.0 / 140.0
"""float: 15N enrichment of the solids over the protosolar gas at `ICE_REFERENCE_K` (the terrestrial 14N/15N, 274 from the table's IUPAC mix) and at `COLD_ICE_K` (about 140 in the comets, 'very constant', [R5167]). [Stand-in ISO6, tier T2, tweak L0.]"""

N_T0_K = math.log(NITROGEN_ENRICHMENT_COLD / NITROGEN_ENRICHMENT_REFERENCE) / (1.0 / COLD_ICE_K - 1.0 / ICE_REFERENCE_K)
"""float: The temperature scale of the 15N enrichment, K (about 24), fitted to the two points above. [Stand-in ISO6, tier T5, tweak L2.]"""

ARGON_36_OVER_38 = 5.501
"""float: The primordial 36Ar/38Ar, the solar wind's ([R5170: Abstract], 5.501 +- 0.005). [Stand-in ISO7, tier T1, tweak L0.]"""

ARGON_40_SHARE_MAX = 0.9996
"""float: The most of an argon inventory that 40Ar can be, leaving the primordial share a gas always holds. [Stand-in ISO8, tier T4, tweak L3.]"""

ARGON_MOLAR_KG = 0.039948
"""float: The molar mass of argon, kg/mol (`elements` for the natural mix; written as the number the column of moles needs)."""

RADIOACTIVE_PARENTS = {"K": ("K-39", "K-40", "K-41"), "U": ("U-238", "U-235")}
"""dict: The elements whose isotopes decay to the body's age, and the isotopes followed (`radiogenic.PARENTS`; 234U is left out)."""

SOLID_ICES = {"ice_co2_mass_fraction": "CO2", "ice_ch4_mass_fraction": "CH4", "ice_nh3_mass_fraction": "NH3", "ice_n2_mass_fraction": "N2", "ice_so2_mass_fraction": "SO2",
              "ice_co_mass_fraction": "CO", "ice_hcn_mass_fraction": "HCN", "ice_other_mass_fraction": "CH2O"}
"""dict: The ice columns of a body and the formula taken for each (`CH2O` for the organics and the rest, a stand-in for their H, C and O). [Stand-in ISO8, tier T4, tweak L3.]"""

ROCK_OXIDE_ELEMENTS = {"rock_na2o_fraction": "Na", "rock_k2o_fraction": "K", "rock_tio2_fraction": "Ti", "rock_cr2o3_fraction": "Cr", "rock_mno_fraction": "Mn", "rock_nio_fraction": "Ni",
                       "rock_p2o5_fraction": "P"}
"""dict: The minor-oxide column of a rock and the element it holds."""

CORE_ELEMENTS = {"core_ni_fraction": "Ni", "core_s_fraction": "S", "core_si_fraction": "Si", "core_o_fraction": "O", "core_c_fraction": "C", "core_h_fraction": "H"}
"""dict: The light-element column of a core and the element it holds."""

_EARTH_G = 9.80665
_GY_S = radiogenic.GY_S


# --- the elements a body holds ----------------------------------------------------------------------------------------
def _positive(body, name):
    return (getattr(body, name, None) or 0.0) > 0.0


def _formula_elements(formula):
    return set(elements.parse_formula(formula))


def elements_held(body):
    """The symbols of the elements `body` holds, by atomic number (see the module docstring); empty for an object with no stored composition."""
    held = set()
    rock, metal, water = (getattr(body, name, None) or 0.0 for name in ("rock_mass_fraction", "metal_mass_fraction", "water_mass_fraction"))
    envelope = getattr(body, "envelope_mass_fraction", None) or 0.0
    if rock > 0.0:
        if getattr(body, "carbon_rock", False):
            held |= {"C", "Si"}
        else:
            held |= {"O", "Mg", "Si", "Al", "Ca", "Fe"}
            held |= {symbol for name, symbol in ROCK_OXIDE_ELEMENTS.items() if _positive(body, name)}
        if radiogenic.content_of(body) is not None:
            held |= {"K", "Th", "U"}
    if metal > 0.0:
        held |= {"Fe"} | {symbol for name, symbol in CORE_ELEMENTS.items() if _positive(body, name)}
    if water > 0.0:
        held |= {"H", "O"}
    for name, formula in SOLID_ICES.items():
        if _positive(body, name):
            held |= _formula_elements(formula)
    if envelope > 0.0:
        held |= {"H", "He"}
        if _positive(body, "envelope_z_fraction"):
            held |= {"C", "N", "O"}
    for gas in (getattr(body, "gas_partials_pa", None) or {}):
        held |= _formula_elements(gas)
    return tuple(sorted(held, key=lambda symbol: elements.element(symbol).z))


# --- the star's gas ---------------------------------------------------------------------------------------------------
def _metallicity(star):
    value = getattr(star, "metallicity_fe_h", None)
    return 0.0 if value is None else float(value)


def _age_gy(star):
    value = getattr(star, "age", None)
    return tuning.ROGUE_EARTH_REFERENCE_AGE_GY if value is None else float(value)


def gas_d_h(star):
    """The deuterium over hydrogen of the gas the star formed from: the primordial value less what star formation burned, which grows with the metals (`ASTRATION`)."""
    return PRIMORDIAL_D_H * math.exp(-ASTRATION * 10.0 ** _metallicity(star))


def _galactocentric_kpc(star):
    """The distance from the Galactic centre at which the gas has the star's metallicity, kpc: the Sun's distance moved by the gas gradient, held inside `GALACTOCENTRIC_RANGE_KPC`."""
    kpc = SOLAR_GALACTOCENTRIC_KPC - _metallicity(star) / nebula_composition.GAS_GRADIENT_DEX_PER_KPC
    return min(max(kpc, GALACTOCENTRIC_RANGE_KPC[0]), GALACTOCENTRIC_RANGE_KPC[1])


def _gradient_ratio(star, gradient, anchors):
    """A galactic-evolution ratio (12C/13C or 14N/15N) for the star: the gradient at its radius over the gradient at the Sun's, times the ratio between the local interstellar value and the birth value by age."""
    slope, intercept = gradient
    local = slope * SOLAR_GALACTOCENTRIC_KPC + intercept
    here = slope * _galactocentric_kpc(star) + intercept
    now, birth = anchors
    by_age = now + (birth - now) * min(_age_gy(star) / tuning.STELLAR_SOLAR_SYSTEM_AGE_GY, 2.0)
    return by_age * here / local


def gas_c13_over_c12(star):
    """13C/12C of the star's birth gas."""
    return 1.0 / _gradient_ratio(star, CARBON_GRADIENT, CARBON_RATIOS)


def gas_n15_over_n14(star):
    """15N/14N of the star's birth gas."""
    return 1.0 / _gradient_ratio(star, NITROGEN_GRADIENT, NITROGEN_RATIOS)


def _ice_temperature_k(body):
    """The temperature an ice of the body formed at, K: the disk's at its formation distance, never above the snow line's and never below `FLOOR_ICE_K`."""
    star = getattr(body, "star", None)
    if star is None or not getattr(star, "snow_line_au", None):
        return ICE_REFERENCE_K          # a rogue planet or a body with no host line: the snow line's ice
    return min(max(bulk_draw.formation_temperature_k(body), FLOOR_ICE_K), ICE_REFERENCE_K)


def ice_factor(temperature_k, reference, t0_k):
    """How much more of the heavy isotope an ice formed at `temperature_k` holds than the gas it froze from: `reference * exp(t0 (1/T - 1/T_ref))` (the enrichment at the snow line, growing as it gets colder)."""
    return reference * math.exp(t0_k * (1.0 / temperature_k - 1.0 / ICE_REFERENCE_K))


# --- the reservoirs of a body -----------------------------------------------------------------------------------------
def _hydrogen_reservoirs(body):
    """`[(hydrogen atoms, is_ice)]` per kg of the body: the envelope's hydrogen gas, the water and each ice's hydrogen."""
    envelope = getattr(body, "envelope_mass_fraction", None) or 0.0
    helium = getattr(body, "envelope_helium_fraction", None) or 0.0
    heavy = getattr(body, "envelope_z_fraction", None) or 0.0
    reservoirs = []
    if envelope > 0.0:
        reservoirs.append((envelope * max(0.0, 1.0 - helium - heavy) / 1.00794, False))
    water = getattr(body, "water_mass_fraction", None) or 0.0
    if water > 0.0:
        reservoirs.append((water * 2.0 / 18.015, True))
    for name, formula in SOLID_ICES.items():
        share = getattr(body, name, None) or 0.0
        count = elements.parse_formula(formula).get("H", 0)
        if share > 0.0 and count:
            reservoirs.append((share * count / (elements.molar_mass(formula) * 1000.0), True))
    return reservoirs


def deuterium_fraction(body, stream):
    """The atom fraction of hydrogen that is deuterium, over the body's hydrogen reservoirs; `None` where the body holds no hydrogen reservoir."""
    reservoirs = _hydrogen_reservoirs(body)
    total = math.fsum(atoms for atoms, _ice in reservoirs)
    if total <= 0.0:
        return None
    gas = gas_d_h(getattr(body, "star", None))
    scatter = math.exp(D_H_SCATTER_LN * stream.gauss())
    ice = gas * ice_factor(_ice_temperature_k(body), D_H_ENRICHMENT_REFERENCE, D_H_T0_K) * scatter
    return math.fsum(atoms * (r / (1.0 + r)) for atoms, is_ice in reservoirs for r in ((ice if is_ice else gas),)) / total


def carbon_13_fraction(body):
    """The atom fraction of carbon that is 13C (the star's birth gas; the ices keep it)."""
    ratio = gas_c13_over_c12(getattr(body, "star", None))
    return ratio / (1.0 + ratio)


def nitrogen_15_fraction(body):
    """The atom fraction of nitrogen that is 15N: the gas value for the envelope's share of the body and the enriched value for the rest."""
    gas = gas_n15_over_n14(getattr(body, "star", None))
    solid = gas * ice_factor(_ice_temperature_k(body), NITROGEN_ENRICHMENT_REFERENCE, N_T0_K)
    weight = min(max(getattr(body, "envelope_mass_fraction", None) or 0.0, 0.0), 1.0)
    ratio = weight * gas + (1.0 - weight) * solid
    return ratio / (1.0 + ratio)


def potassium_40_argon_moles(body):
    """The moles of 40Ar the body's potassium has made and its interior has let out; 0 where the body stores no content or interior."""
    content = radiogenic.content_of(body)
    found = interior_heat.heat_of(body)
    if content is None or found is None:
        return 0.0
    mantle = interior_heat.mantle_of(body.mass, body.radius, radiogenic.rock_of(body))
    if mantle is None:
        return 0.0
    age = _age_gy(getattr(body, "star", None))
    rock = radiogenic.rock_of(body)
    weight_kg = elements.element("K").atomic_weight * 1e-3
    potassium_mol = content.ppm("K") * 1e-6 * body.mass * rock.silicate / weight_kg
    mix, _factor = radiogenic.formation_mix("K")
    atoms = {name: share / elements.isotope(name).mass_u for name, share in mix.items()}
    k40_atoms = potassium_mol * atoms["K-40"] / math.fsum(atoms.values())
    k40 = elements.isotope("K-40")
    branch = dict((mode, percent) for mode, percent, *_ in k40.modes).get("B+", 0.0) / 100.0
    decayed = -math.expm1(-k40.decay_constant * age * _GY_S)
    return branch * k40_atoms * decayed * surface_pressure.degassed_share(found.released_j, mantle.dry_kg)


def argon_40_fraction(body):
    """The atom fraction of the air's argon that is 40Ar (radiogenic), or `None` where the air holds no argon."""
    partial = (getattr(body, "gas_partials_pa", None) or {}).get("Ar")
    if not partial or partial <= 0.0:
        return None
    air_mol = partial * 4.0 * math.pi * (body.radius * 1000.0) ** 2 / (body.gravity * _EARTH_G) / ARGON_MOLAR_KG
    if air_mol <= 0.0:
        return None
    return min(potassium_40_argon_moles(body) / air_mol, ARGON_40_SHARE_MAX)


# --- the rows ---------------------------------------------------------------------------------------------------------
def _natural(symbol):
    """`{isotope: atom fraction}` of the table's natural mix of `symbol`, normalized to 1."""
    found = {iso.name: iso.abundance for iso in elements.isotopes_of(symbol) if iso.abundance}
    total = math.fsum(found.values())
    return {name: value / total for name, value in found.items()}


def _two(symbol, light, heavy, heavy_fraction):
    """The record of an element with a light and a heavy isotope that follow a physical value; every other isotope of the element keeps its natural share of what is left."""
    natural = _natural(symbol)
    rest = {name: value for name, value in natural.items() if name not in (light, heavy)}
    rest_total = math.fsum(rest.values())
    major = 1.0 - rest_total
    found = {light: major * (1.0 - heavy_fraction), heavy: major * heavy_fraction}
    found.update(rest)
    return found


def _decayed(symbol, age_gy):
    """`{isotope: atom fraction}` of `symbol`'s `RADIOACTIVE_PARENTS` at `age_gy`: the formation mix decayed."""
    mix, _factor = radiogenic.formation_mix(symbol)
    atoms = {}
    for name in RADIOACTIVE_PARENTS[symbol]:
        isotope = elements.isotope(name)
        decay = 1.0 if isotope.stable else math.exp(-isotope.decay_constant * age_gy * _GY_S)
        atoms[name] = mix[name] / isotope.mass_u * decay
    total = math.fsum(atoms.values())
    return {name: value / total for name, value in atoms.items()}


def _argon(share_40):
    primordial = 1.0 - share_40
    return {"Ar-36": primordial * ARGON_36_OVER_38 / (1.0 + ARGON_36_OVER_38), "Ar-38": primordial / (1.0 + ARGON_36_OVER_38), "Ar-40": share_40}


def stream_of(body, key):
    """The body's own stream for one draw of the record: seeded from what is final about the body, so no draw of the run's stream moves."""
    star = getattr(body, "star", None)
    return draw.Stream(f"isotopes:{key}:{getattr(body, 'is_moon', False)}:{body.mass!r}:{body.radius!r}:{body.distance!r}:{getattr(star, 'age', None)!r}")


def record_of(body):
    """The rows of `body`'s record (`ROW_FIELDS` dicts, by atomic number then mass number), or `None` for a body with no stored composition."""
    held = elements_held(body)
    if not held:
        return None
    age = _age_gy(getattr(body, "star", None))
    rows = []
    for symbol in held:
        origin, mix = NATURAL, _natural(symbol)
        if symbol == "H" and "H-2" in mix:
            deuterium = deuterium_fraction(body, stream_of(body, "D"))
            if deuterium is not None:
                origin, mix = PHYSICAL, _two("H", "H-1", "H-2", deuterium)
        elif symbol == "C":
            origin, mix = PHYSICAL, _two("C", "C-12", "C-13", carbon_13_fraction(body))
        elif symbol == "N":
            origin, mix = PHYSICAL, _two("N", "N-14", "N-15", nitrogen_15_fraction(body))
        elif symbol == "Ar":
            share = argon_40_fraction(body)
            if share is not None:
                origin, mix = RADIOGENIC, _argon(share)
        elif symbol in RADIOACTIVE_PARENTS and radiogenic.content_of(body) is not None:
            origin, mix = DECAY, _decayed(symbol, age)
        for name in sorted(mix, key=lambda isotope_name: elements.isotope(isotope_name).a):
            rows.append({"element": symbol, "isotope": name, "fraction": mix[name], "origin": origin})
    return rows


SHEET_RATIOS = (("D/H", "H-2", "H-1", "Deuterium over hydrogen"), ("13C/12C", "C-13", "C-12", "Carbon-13 over carbon-12"), ("15N/14N", "N-15", "N-14", "Nitrogen-15 over nitrogen-14"),
                ("40Ar/36Ar", "Ar-40", "Ar-36", "Argon-40 over argon-36"), ("235U/238U", "U-235", "U-238", "Uranium-235 over uranium-238"),
                ("40K/39K", "K-40", "K-39", "Potassium-40 over potassium-39"))
"""tuple: The ratios a property sheet lists, `(short name, top isotope, bottom isotope, label)`: the ones the record moves away from the natural mix."""


def sheet_rows(rows):
    """
    The lines a property sheet shows for a stored record: `[(key, label, text, ratio)]` for each of `SHEET_RATIOS` whose two isotopes the body holds, then one line counting the elements
    and isotopes of the record. `rows` are `{"isotope", "fraction"}` mappings (a stored record). Empty for a body with no record.
    """
    if not rows:
        return []
    fractions = {row["isotope"]: row["fraction"] for row in rows}
    shown = []
    for short, top, bottom, label in SHEET_RATIOS:
        if top in fractions and fractions.get(bottom):
            ratio = fractions[top] / fractions[bottom]
            shown.append((f"isotope_{short}", label, f"{short} {ratio:.4g}", ratio))
    shown.append(("isotope_count", "Isotopes recorded", f"{len(rows)} isotopes of {len({name.split('-')[0] for name in fractions})} elements", float(len(rows))))
    return shown


def generate(body):
    """Sets `ATTRIBUTE` on a planet, moon or rogue planet from its stored state (call after the air, the interior and the ground; own streams, no draw of the run's)."""
    setattr(body, ATTRIBUTE, None)
    if getattr(body, "star", None) is None and getattr(body, "rock_mass_fraction", None) is None:
        return
    setattr(body, ATTRIBUTE, record_of(body))
