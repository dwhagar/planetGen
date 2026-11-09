# planetgen/physics/isotope_record.py

"""
The isotope record of a body (GEN.1377, GEN.1420 and GEN.1421)
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

* **Hydrogen** (`physical`): the deuterium share is the hydrogen-atom-weighted mean of the body's reservoirs. The gas of the star's birth cloud has the primordial D/H (2.527e-5, [R5173]) reduced by the star formation that burned deuterium since: `exp(-ASTRATION * 10^[Fe/H])` (the metal-richer the star, the more of the gas has been through stars). Every ice and the water
  are enriched over the gas by `ice_factor`, a law in the formation temperature `exp(T0 (1/T - 1/T_ref))` (ion-molecule fractionation grows as the gas gets colder, [R5176: Sec. 2.1]) that is
  calibrated at two published points of the Solar System: the ocean (1.56e-4, VSMOW) at the snow-line temperature and comet 67P's water (5.3e-4, [R0088]) at 30 K. A body that formed warmer than
  the snow line took its ice from the snow line's temperature. One log-normal scatter per body.
* **Carbon** (`physical`): 13C/12C from the Galaxy's chemical evolution: 12C/13C = 6.21 D + 18.71 ([R5177: equation 5], D in kpc from the Galactic center) read at the galactocentric distance
  whose gas metallicity is the star's (`nebula_composition.GAS_GRADIENT_DEX_PER_KPC`), and moved in time between the local interstellar value (68) and the Solar System's (89) by the star's age.
* **Nitrogen** (`physical`): 14N/15N the same way ([R5179] as quoted in [R5174: equation 2]: 21.1 D + 123.8), between the local interstellar value (about 400) and the
  protosolar one (441, solar wind, quoted in [R5174]), then raised in 15N by `ice_factor` for the solid share of the body (comets are enriched 2 to 3 times, [R5176]). An envelope's nitrogen is the gas value;
  the two are mixed by the envelope's share of the body's mass.
* **Argon** (`physical`, `radiogenic`): the primordial 36Ar/38Ar is the solar wind's 5.501 [R5178]; 40Ar is the share of the air's argon that the body's potassium made and its interior
  let out: the 40K that has decayed (`radiogenic.formation_mix`, the body's age), times the 10.72 percent of decays that end in 40Ar, times the share of the mantle's easily degassed element
  that the interior's released heat has brought to the surface (`surface_pressure.degassed_share`), over the moles of argon in the air; it cannot be more than all but the primordial share.
* **Potassium, uranium** (`decay`): the isotopes that decay are those of the body's age, from the formation mix of `radiogenic` (today's natural mix decayed back to the Solar System's birth, a stand-in
  for every star) decayed over the star's age. 234U and the daughters (the lead, strontium and neodymium isotopes) stay at the natural mix.

**The short-lived isotopes** (`extinct`, part 2): 26Al, 60Fe, 53Mn, 129I and 244Pu are tracers on the body. Each carries the ratio of the star's birth gas (`initial_ratio`: to 27Al, 56Fe, 55Mn, 127I and, for plutonium,
238U) and `fraction`, the survivors per atom of that reference element at the star's age (nothing for a star of a hundred million years, since the longest mean life is 115 Myr). 26Al is the ratio `radiogenic` already
draws; 60Fe is a sum of decayed additions about the gamma-ray mean of the galaxy (`birth_radioisotopes`, the way 26Al is); 53Mn, 129I and 244Pu take the same sampling about a tenth of the early Solar System's ratios
([R2768: Table 2]). One stream per star and isotope, so the bodies of a star share them.

**The cosmogenic isotopes** (`cosmogenic`, part 2): tritium and radiocarbon in a body with air. The production is the Earth's [R5203, R5204] times the share of the cascade the air column catches, the targets in a gram of
the air (tritium from N and O, radiocarbon from N), and the field's shielding (a body with no dipole gets twice the Earth's, [R5204: Table 1]); the steady-state inventory is the production over the decay constant, spread
over the pool the isotope mixes into (the air, and the surface ocean for hydrogen and, at the Earth's ratio, for carbon). The Earth's own air gives 1.6 and 0.345 atoms per cm2 per second, 14C/C of 1.3e-12 and 3H/H of 1e-20.
A body with no air has none.

Left out, and where it goes: 10Be and 26Al made in the air (no source in hand for a global rate; the pool for beryllium is the crust's); cosmic-ray production in an airless surface and an ocean-covered one; the xenon excess
from I-129 and Pu-244 and the mass-independent oxygen anomalies (O stays at the natural mix); the isotope split of a body's reservoirs (a mantle against an ocean); belts, comets and stars.
"""

import math
from types import SimpleNamespace

from planetgen import tuning
from planetgen.galaxy import shape_draw
from planetgen.physics import birth_radioisotopes, bulk_draw, condensate_shares, constants, disk_condensation, elements, interior_heat, nebula_composition, radiogenic, stellar_abundances, surface_pressure
from planetgen.util import draw

ATTRIBUTE = "isotopes"
"""str: The body attribute that holds the record, a list of `ROW_FIELDS` dicts (the child tables `planet_isotopes`, `moon_isotopes` and `rogue_planet_isotopes`)."""

ROW_FIELDS = ("element", "isotope", "fraction", "initial_ratio", "origin")
"""tuple: The columns of an isotope row after the owner id: the element's symbol, the isotope's name (`H-2`), its atom fraction of the element, its atom ratio to the stable reference isotope at the body's formation (`None` where the isotope is not followed from formation), and where the value comes from."""

NATURAL, PHYSICAL, DECAY, RADIOGENIC, EXTINCT, COSMOGENIC = "natural", "physical", "decay", "radiogenic", "extinct", "cosmogenic"
ORIGINS = (NATURAL, PHYSICAL, DECAY, RADIOGENIC, EXTINCT, COSMOGENIC)
"""tuple: The `origin` of a row: the table's natural mix, a physical law of the star and the formation temperature, a radioactive isotope at the body's age, a daughter the body's own decays made, a
short-lived isotope the star's birth gas held (a tracer: its `fraction` is the survivors per atom of its reference isotope's element, `initial_ratio` the ratio at birth), an isotope the particles of space make in the air
(a tracer: its `fraction` is the steady-state atoms per atom of its element in the pool they mix into)."""

# --- the laws (stand-ins IS1 to IS8) -------------------------------------------------------------------------------------
PRIMORDIAL_D_H = 2.527e-5
"""float: The primordial deuterium over hydrogen by number [R5173]. [Stand-in ISO1, tier T1, tweak L0.]"""

PROTOSOLAR_D_H = 2.1e-5
"""float: The deuterium over hydrogen of the Sun's birth gas, used only to set `ASTRATION` ([R5176: Sec. 4.2], quoting the Jovian and solar-wind values). [Stand-in ISO2, tier T2, tweak L0.]"""

ASTRATION = math.log(PRIMORDIAL_D_H / PROTOSOLAR_D_H)
"""float: The fall of ln D/H from the primordial value to the Sun's birth gas, at the Sun's metallicity; the fall is taken proportional to the metals, `ASTRATION * 10^[Fe/H]` (stars burn deuterium as they make metals). [Stand-in ISO2, tier T5, tweak L1.]"""

ICE_REFERENCE_K = 182.0
"""float: The temperature of the snow line's water, the highest the enrichment of an ice is taken at (a warmer body took its ice from there). [Stand-in ISO3, tier T1, tweak L0.]"""

COLD_ICE_K = 30.0
"""float: The temperature of the cold calibration point, comet 67P's water (the CO2 and CH4 line zone, below which the forward H3+ + HD reaction wins, [R5176]). [Stand-in ISO3, tier T2, tweak L0.]"""

FLOOR_ICE_K = 20.0
"""float: The coldest temperature the law is read at: below it the enrichment would run past the calibration. [Stand-in ISO4, tier T4, tweak L3.]"""

OCEAN_D_H = 1.56e-4
"""float: The deuterium over hydrogen of ocean water (VSMOW, [R5176: Sec. 1.2]), the enrichment at `ICE_REFERENCE_K`. [Stand-in ISO3, tier T1, tweak L0.]"""

COMET_D_H = 5.3e-4
"""float: The deuterium over hydrogen of the water of a Jupiter-family comet ((5.3 +- 0.7)e-4, [R0088], quoted in [R5175]), the enrichment at `COLD_ICE_K`. [Stand-in ISO3, tier T1, tweak L0.]"""

D_H_ENRICHMENT_REFERENCE = OCEAN_D_H / PROTOSOLAR_D_H
D_H_ENRICHMENT_COLD = COMET_D_H / PROTOSOLAR_D_H
D_H_T0_K = math.log(D_H_ENRICHMENT_COLD / D_H_ENRICHMENT_REFERENCE) / (1.0 / COLD_ICE_K - 1.0 / ICE_REFERENCE_K)
"""float: The temperature scale of the deuterium enrichment of ice, K (about 44): the law is fitted to the two published points, so it is tuned to published statistics. [Stand-in ISO3, tier T5, tweak L2.]"""

D_H_SCATTER_LN = 0.3
"""float: The log-normal scatter of one body's ice D/H about the law (comets span a factor of about 2 to 3 around their mean, [R5175]). [Stand-in ISO4, tier T4, tweak L3.]"""

CARBON_GRADIENT = (6.21, 18.71)
"""tuple: 12C/13C against galactocentric distance in kpc, `(slope, intercept)` ([R5177: equation 5], fitted to CN, CO and H2CO). [Stand-in ISO5, tier T1, tweak L0.]"""

NITROGEN_GRADIENT = (21.1, 123.8)
"""tuple: 14N/15N against galactocentric distance in kpc ([R5179], equation 2 of [R5174]). [Stand-in ISO5, tier T1, tweak L0.]"""

CARBON_RATIOS = (68.0, 89.0)
"""tuple: `(local interstellar now, Solar System birth)` 12C/13C ([R5177: Abstract]). The star's age moves the ratio between them. [Stand-in ISO5, tier T1, tweak L0.]"""

NITROGEN_RATIOS = (400.0, 441.0)
"""tuple: `(local interstellar now, protosolar)` 14N/15N ([R5174: Abstract and Sec. 3.3]; the solar wind). [Stand-in ISO5, tier T1, tweak L0.]"""

SOLAR_GALACTOCENTRIC_KPC = 7.9
"""float: The Sun's distance from the Galactic center, kpc [R5177], where the two gradients read the local interstellar ratios. [Stand-in ISO5, tier T1, tweak L0.]"""

GALACTOCENTRIC_RANGE_KPC = (2.0, 20.0)
"""tuple: The range the gradients are read in (the data of [R5177] reach from 2 to 12 kpc; the inner and outer Galaxy past it hold no data). [Stand-in ISO4, tier T4, tweak L3.]"""

NITROGEN_ENRICHMENT_REFERENCE = 441.0 / 274.0
NITROGEN_ENRICHMENT_COLD = 441.0 / 140.0
"""float: 15N enrichment of the solids over the protosolar gas at `ICE_REFERENCE_K` (the terrestrial 14N/15N, 274 from the table's IUPAC mix) and at `COLD_ICE_K` (about 140 in the comets, 'very constant', [R5175]). [Stand-in ISO6, tier T2, tweak L0.]"""

N_T0_K = math.log(NITROGEN_ENRICHMENT_COLD / NITROGEN_ENRICHMENT_REFERENCE) / (1.0 / COLD_ICE_K - 1.0 / ICE_REFERENCE_K)
"""float: The temperature scale of the 15N enrichment, K (about 24), fitted to the two points above. [Stand-in ISO6, tier T5, tweak L2.]"""

ARGON_36_OVER_38 = 5.501
"""float: The primordial 36Ar/38Ar, the solar wind's ([R5178: Abstract], 5.501 +- 0.005). [Stand-in ISO7, tier T1, tweak L0.]"""

ARGON_40_SHARE_MAX = 0.9996
"""float: The most of an argon inventory that 40Ar can be, leaving the primordial share a gas always holds. [Stand-in ISO8, tier T4, tweak L3.]"""

ARGON_MOLAR_KG = 0.039948
"""float: The molar mass of argon, kg/mol (`elements` for the natural mix; written as the number the column of moles needs)."""

EXTINCT_ISOTOPES = {"Al-26": ("Al", "Al-27", "Al26"), "Fe-60": ("Fe", "Fe-56", "Fe60"), "Mn-53": ("Mn", "Mn-55", "Mn53"), "I-129": ("I", "I-127", "I129"), "Pu-244": ("Pu", "U-238", "Pu244")}
"""dict: The short-lived isotopes a body keeps a birth ratio of: `{isotope: (element, the stable isotope the ratio is taken to, key of birth_radioisotopes.MEAN_LIFE_YR)}`. 244Pu is taken to 238U, since plutonium has no stable isotope."""

ESS_INITIAL_RATIO = {"Mn-53": 7.0e-6, "I-129": 1.28e-4, "Pu-244": 7.0e-3}
"""dict: The ratios the early Solar System held at the formation of its first solids: 53Mn/55Mn (7 +- 1)e-6, 129I/127I (1.28 +- 0.03)e-4, 244Pu/238U (7 +- 1)e-3 ([R2768: Table 2], 2 sigma, statistical only).
[Stand-in ISO10, tier T1, tweak L0.]"""

EXTINCT_MEAN_OVER_ESS = 0.1
"""float: The mean ratio of a short-lived isotope in the gas of a star's galaxy over the early Solar System's, for 53Mn, 129I and 244Pu, which no gamma-ray line measures. The interstellar 26Al/27Al is 3 to 25 times
under the Solar System's ([R2768: Sec. 3.2]) and the birth-ratio sampling of `birth_radioisotopes` puts the Solar System's value in the upper tail of the draw for a mean a tenth of it; the same tenth is taken
for the others, whose draw then spreads by the recurrence of the events that make them. [Stand-in ISO11, tier T6, tweak L1.]"""

COMET_ICE_K = COLD_ICE_K
"""float: The temperature a comet's ice formed at, K: comets form at the same distance in units of the star's snow line for every star, so at the temperature the cold calibration point has (comet 67P's, 30 K) whatever
the star's luminosity. [Stand-in ISO17, tier T6, tweak L1.]"""

COSMOGENIC_RATE_EARTH = {"H-3": 0.345, "C-14": 1.6}
"""dict: The global mean production of each isotope in the Earth's atmosphere, atoms per cm2 per second, at a moderate solar activity (modulation potential 650 MV): tritium 0.345 ([R5204: Sec. 3, Table 1]; 0.32
in [R5204] for [R5206]), radiocarbon 1.6 ([R5203: Table 1, quoted]; 2.36e26 atoms a year at 802 MV, about 1.5, [R5203: Sec. 3.2]; 1.6 to 2, [R5205: Sec. 2], quoting [R5207]).
[Stand-in ISO12 (14C, tier T2, tweak L0) and ISO13 (3H, tier T1, tweak L0).]"""

EARTH_COLUMN_G_CM2 = 1033.0
"""float: The mass of the Earth's air over a square centimeter at sea level, g/cm2 ([R5203: Eq. 6]). [Stand-in ISO13, tier T1, tweak L0.]"""

ATTENUATION_G_CM2 = 150.0
"""float: The depth of air in which the nucleon cascade of a cosmic ray is mostly absorbed, g/cm2 (the Regener-Pfotzer maximum of tritium lies at about 40 g/cm2, [R5204: Sec. 3], and the cascade is spent
within a few hundred); a column of mass X catches the share `1 - exp(-X / 150)` of the cascade, over the Earth's. [Stand-in ISO14, tier T5, tweak L0.]"""

EARTH_TARGETS_PER_G = {"N": 3.22e22, "O": 8.67e21}
"""dict: Atoms of the two spallation targets in a gram of the Earth's air ([R5204: Eq. 4], kappa_N and kappa_O). [Stand-in ISO13, tier T1, tweak L0.]"""

EARTH_DIPOLE_A_M2 = 7.8e22
"""float: The Earth's dipole moment, A m2 ([R5203: Sec. 3.2]). [Stand-in ISO13, tier T1, tweak L0.]"""

NO_FIELD_FACTOR = 0.72 / 0.345
"""float: How many times the Earth's global production a body with no field gets: the polar caps of the Earth (60 to 90 degrees), where no cutoff stops a particle, over the global mean, for tritium at 650 MV
([R5204: Table 1], 0.72 and 0.345). [Stand-in ISO15, tier T1, tweak L0.]"""

OCEAN_CARBON_OVER_AIR = 54.0
"""float: The carbon of the air and the surface ocean together over that of the air alone, for a body with a surface ocean: the Earth's 41 000 Gt of exchangeable carbon (air, ocean, soil, biosphere) over its 760 Gt
of air, taken for every ocean world. [Stand-in ISO16, tier T6, tweak L3.]"""

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
AVOGADRO = 6.02214076e23
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
    """The distance from the Galactic center at which the gas has the star's metallicity, kpc: the Sun's distance moved by the gas gradient, held inside `GALACTOCENTRIC_RANGE_KPC`."""
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
    fixed = getattr(body, "ice_formation_k", None)
    if fixed is not None:
        return min(max(fixed, FLOOR_ICE_K), ICE_REFERENCE_K)        # a small body (GEN.1487): the temperature its ice formed at is worked out from what its row stores
    star = getattr(body, "star", None)
    distance = getattr(body, "formation_distance_au", None)
    distance = getattr(body, "distance", None) if distance is None else distance
    if star is None or not getattr(star, "snow_line_au", None) or not distance or distance <= 0.0:
        return ICE_REFERENCE_K          # a rogue planet, a body with no host line, or one with no orbit distance (a pulsar planet): the snow line's ice
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


REFERENCE_ISOTOPE = {"K": "K-39", "U": "U-238"}
"""dict: The isotope each radioactive parent's `initial_ratio` is counted against (the stable, or far longer lived, isotope of its element)."""


def _decayed(symbol, age_gy):
    """`({isotope: atom fraction at age_gy}, {isotope: atom ratio to the reference isotope at formation})` of `symbol`'s `RADIOACTIVE_PARENTS`: the formation mix, then decayed."""
    mix, _factor = radiogenic.formation_mix(symbol)
    born = {name: mix[name] / elements.isotope(name).mass_u for name in RADIOACTIVE_PARENTS[symbol]}
    atoms = {}
    for name, count in born.items():
        isotope = elements.isotope(name)
        atoms[name] = count * (1.0 if isotope.stable else math.exp(-isotope.decay_constant * age_gy * _GY_S))
    total = math.fsum(atoms.values())
    reference = born[REFERENCE_ISOTOPE[symbol]]
    return ({name: value / total for name, value in atoms.items()},
            {name: count / reference for name, count in born.items() if name != REFERENCE_ISOTOPE[symbol] and not elements.isotope(name).stable})


def _argon(share_40):
    primordial = 1.0 - share_40
    return {"Ar-36": primordial * ARGON_36_OVER_38 / (1.0 + ARGON_36_OVER_38), "Ar-38": primordial / (1.0 + ARGON_36_OVER_38), "Ar-40": share_40}


def stream_of(body, key):
    """The body's own stream for one draw of the record: seeded from what is final about the body, so no draw of the run's stream moves."""
    seed = getattr(body, "isotope_seed", None)
    if seed is not None:
        return draw.Stream(f"isotopes:{key}:{seed}")          # a small body (GEN.1487) has no mass and orbit of its own to seed from: its row names the stream
    star = getattr(body, "star", None)
    return draw.Stream(f"isotopes:{key}:{getattr(body, 'is_moon', False)}:{body.mass!r}:{body.radius!r}:{body.distance!r}:{getattr(star, 'age', None)!r}")


def fe60_mean(body):
    """The steady-state 60Fe/56Fe in the gas of the body's star's galaxy at its [Fe/H] (`birth_radioisotopes.mean_ratio`), the way `radiogenic.al26_mean` does for 26Al."""
    shape = getattr(getattr(body, "system_config", None), "GALAXY_SHAPE", None)
    shape = shape_draw.default_shape(None) if shape is None else shape
    iron = stellar_abundances.PROTO_MASS_FRACTION["Fe"] * 10.0 ** _metallicity(getattr(body, "star", None))
    return birth_radioisotopes.mean_ratio("Fe60", radiogenic._depletion_time_yr(shape), iron)


def birth_ratio(body, isotope):
    """
    The ratio of a short-lived isotope to its stable reference in the gas the body's star was born from: 26Al/27Al is the one `radiogenic` already drew (`None` for a body with none, a giant), 60Fe/56Fe a sum of decayed
    additions about the gamma-ray mean (`fe60_mean`), 53Mn, 129I and 244Pu the same about `ESS_INITIAL_RATIO` times `EXTINCT_MEAN_OVER_ESS`. One stream per star and isotope, so every body of a star shares it.
    """
    if isotope == "Al-26":
        return getattr(body, "al26_ratio", None)
    key = radiogenic._star_key(getattr(body, "star", None))
    stream = stream_of(body, isotope) if key is None else draw.Stream(f"{isotope}-birth:{key}")
    mean = fe60_mean(body) if isotope == "Fe-60" else ESS_INITIAL_RATIO[isotope] * EXTINCT_MEAN_OVER_ESS
    return birth_radioisotopes.sample_ratio(EXTINCT_ISOTOPES[isotope][2], mean, stream)


def extinct_rows(body, held, age_gy):
    """The `EXTINCT` rows of a body: each short-lived isotope its element or a rock allows, with its birth ratio and what is left of it at the star's age (nothing in a few hundred million years)."""
    rows = []
    for isotope, (element, reference, key) in EXTINCT_ISOTOPES.items():
        if (element in ("Al", "Fe") and element not in held) or "Si" not in held:
            continue
        initial = birth_ratio(body, isotope)
        if initial is None:
            continue
        left = initial * math.exp(-age_gy * 1.0e9 / birth_radioisotopes.MEAN_LIFE_YR[key])
        reference_fraction = _natural(reference.split("-")[0])[reference]
        rows.append({"element": element, "isotope": isotope, "fraction": left * reference_fraction, "initial_ratio": initial, "origin": EXTINCT})
    return rows


def _air_atoms(body):
    """`(atoms of each element in a gram of the air, the air's mean molar mass in g/mol, the air column in g/cm2)` from the gases the body holds, or `None` for a body with no air."""
    partials = getattr(body, "gas_partials_pa", None)
    pressure = getattr(body, "atmospheric_pressure", None)
    gravity = getattr(body, "gravity", None)
    if not partials or not pressure or not gravity or sum(partials.values()) <= 0.0:
        return None
    total = sum(partials.values())
    mean_molar = 0.0
    counts = {}
    for gas, partial in partials.items():
        try:
            formula, molar = elements.parse_formula(gas), elements.molar_mass(gas) * 1000.0
        except (KeyError, ValueError):
            continue
        mean_molar += partial / total * molar
        for symbol, count in formula.items():
            counts[symbol] = counts.get(symbol, 0.0) + partial / total * count
    if mean_molar <= 0.0:
        return None
    column = float(pressure) / (float(gravity) * _EARTH_G) / 10.0
    return {symbol: count / mean_molar * AVOGADRO for symbol, count in counts.items()}, mean_molar, column


def _area_cm2(body):
    return 4.0 * math.pi * (float(body.radius) * 1.0e5) ** 2


def _ocean_water_kg(body):
    """The water of the body's surface ocean, kg (its share of the surface times its depth), 0 when there is none."""
    fraction, depth = getattr(body, "ocean_fraction", None), getattr(body, "ocean_depth_km", None)
    if not fraction or not depth:
        return 0.0
    return float(fraction) * 4.0 * math.pi * (float(body.radius) * 1.0e3) ** 2 * float(depth) * 1.0e3 * 1000.0


def cosmogenic_rate(body, isotope):
    """
    The production of `isotope` (`COSMOGENIC_RATE_EARTH`) in the body's air, atoms per cm2 per second: the Earth's rate times the share of the cascade the column catches (`ATTENUATION_G_CM2`), the
    targets in a gram of the air over the Earth's (tritium from nitrogen and oxygen, radiocarbon from the nitrogen alone), and the field's shielding (`NO_FIELD_FACTOR` with no dipole, falling to 1 at the Earth's
    moment and under it above). 0 for a body with no air. The galactic flux, the star's wind and the star's own flares are the Earth's.
    """
    air = _air_atoms(body)
    if air is None:
        return 0.0
    atoms, _molar, column = air
    caught = (1.0 - math.exp(-column / ATTENUATION_G_CM2)) / (1.0 - math.exp(-EARTH_COLUMN_G_CM2 / ATTENUATION_G_CM2))
    if isotope == "H-3":
        target = (atoms.get("N", 0.0) + atoms.get("O", 0.0)) / (EARTH_TARGETS_PER_G["N"] + EARTH_TARGETS_PER_G["O"])
    else:
        target = atoms.get("N", 0.0) / EARTH_TARGETS_PER_G["N"]
    moment = float(getattr(body, "magnetic_moment_a_m2", None) or 0.0) / EARTH_DIPOLE_A_M2
    field = NO_FIELD_FACTOR / (1.0 + (NO_FIELD_FACTOR - 1.0) * moment)
    return COSMOGENIC_RATE_EARTH[isotope] * caught * target * field


def _pool_atoms(body, element):
    """The atoms of `element` in the pool the isotopes of the air mix into: the air's, and for carbon the surface ocean's (`OCEAN_CARBON_OVER_AIR`), for hydrogen the surface ocean's water."""
    atoms, _molar, column = _air_atoms(body)
    air = atoms.get(element, 0.0) * column * _area_cm2(body)                       # atoms per gram of air times grams of air
    ocean_kg = _ocean_water_kg(body)
    if element == "H":
        return air + ocean_kg * 1000.0 / 18.015 * 2.0 * AVOGADRO
    return air * OCEAN_CARBON_OVER_AIR if ocean_kg > 0.0 else air


def cosmogenic_rows(body):
    """The `COSMOGENIC` rows of a body with air: tritium and radiocarbon at their steady state, the production over the decay constant per cm2 times the body's area, over the atoms of the element in the pool (`_pool_atoms`)."""
    rows = []
    if _air_atoms(body) is None:
        return rows
    for isotope, element in (("H-3", "H"), ("C-14", "C")):
        rate = cosmogenic_rate(body, isotope)
        pool = _pool_atoms(body, element)
        if rate <= 0.0 or pool <= 0.0:
            continue
        made = rate * _area_cm2(body) / elements.isotope(isotope).decay_constant
        rows.append({"element": element, "isotope": isotope, "fraction": made / (pool + made), "initial_ratio": None, "origin": COSMOGENIC})
    return rows


SMALL_BODY_KINDS = ("belt", "belt_member", "comet", "interstellar_comet", "asteroid_field", "field_member")
"""tuple: The kinds of small body `small_body_record` takes (GEN.1487): the belt and its members, a comet, an interstellar comet, an interstellar field and its members."""

ICE_BUDGET_MARGIN_K = 0.1
"""float: How far below the water line (K) the ice budget of a small body formed at the line is read: water condenses below the line only. [Stand-in ISO19, tier T6, tweak L0.]"""

SMALL_BODY_ICE_COLUMNS = {"NH3": "ice_nh3_mass_fraction", "CO2": "ice_co2_mass_fraction", "CH4": "ice_ch4_mass_fraction", "N2": "ice_n2_mass_fraction", "CO": "ice_co_mass_fraction"}
"""dict: The ices of the disk budget (`disk_condensation.masses_percent`) a small body's record follows beside its water, and the column of a planet's ice each one is read into. The ammonium salts, H2S and the noble
gases that the budget freezes are left out of the elements a small body holds. [Stand-in ISO18, tier T6, tweak L0.]"""


def small_body_temperature_k(kind, snow_line_au=None, distance_km=None):
    """
    The temperature the ice of a small body formed at, K: a belt's (and its members') at its distance (`condensate_shares.ICE_CONDENSATION_K` water at the star's snow line, falling as the square root of the distance,
    held inside the ice law's range; the snow line when the belt stores none), a comet's at `COMET_ICE_K`, an interstellar comet's too (a planetesimal thrown out of a disk, as a comet is), an interstellar field's
    the snow line's. [Stand-in ISO20, tier T6, tweak L0.]
    """
    if kind in ("comet", "interstellar_comet"):
        return COMET_ICE_K
    if kind in ("belt", "belt_member") and snow_line_au and distance_km and distance_km > 0.0:
        temperature = condensate_shares.ICE_CONDENSATION_K["H2O"] * math.sqrt(float(snow_line_au) * constants.AU_TO_KM / float(distance_km))
        return min(max(temperature, FLOOR_ICE_K), ICE_REFERENCE_K)
    return ICE_REFERENCE_K


def galactic_gas_metallicity(galactic_radius_pc):
    """The [Fe/H] taken for a small body that has no star (an interstellar comet or field): the gas of the galaxy at the radius it is at (`nebula_composition.gas_gradient_oh`, [O/H] standing for [Fe/H]). [Stand-in ISO20, tier T6, tweak L0.]"""
    return nebula_composition.gas_gradient_oh(float(galactic_radius_pc))


def small_body_view(kind, metallicity_fe_h, age_gy, ice_mass_fraction, distance_km=None, snow_line_au=None, seed="", rock_fraction=None, metal_fraction=None, water_fraction=None, ice_fraction=None):
    """
    A stand-in `record_of` reads for a belt, a comet, an interstellar comet, an interstellar field or a member of a belt or field (GEN.1487), made from the values its row stores: no table of its own holds a record,
    since the record is a pure function of those values. The body is rock and ice: `ice_mass_fraction` of it is ice (a member gives its own `rock_fraction`, `metal_fraction` and `water_fraction`, and
    `ice_fraction` for its other ices), the ice split into water and the other ices of the budget at the temperature it formed at (`small_body_temperature_k`), and the rest rock. The draws are made on a stream
    named by `kind` and `seed` (the belt's system and distance, a comet's id), so one body always shows the same record. `age_gy` of `None` is the Earth's reference age.
    """
    temperature = small_body_temperature_k(kind, snow_line_au, distance_km)
    masses = disk_condensation.masses_percent(min(temperature, ICE_REFERENCE_K - ICE_BUDGET_MARGIN_K))    # water is frozen below its line, not on it: a body at the line itself holds the water the line leaves
    ices = {name: value for name, value in masses.items() if name in ("H2O", *SMALL_BODY_ICE_COLUMNS)}
    total = math.fsum(ices.values())
    ice = float(ice_mass_fraction or 0.0) if ice_fraction is None and water_fraction is None else float(ice_fraction or 0.0) + float(water_fraction or 0.0)
    view = SimpleNamespace(star=SimpleNamespace(metallicity_fe_h=float(metallicity_fe_h), age=None if age_gy is None else float(age_gy), snow_line_au=snow_line_au), ice_formation_k=temperature,
                           isotope_seed=f"ice:{kind}:{seed}", rock_mass_fraction=(1.0 - ice) if rock_fraction is None else float(rock_fraction), metal_mass_fraction=float(metal_fraction or 0.0),
                           water_mass_fraction=0.0, mass=0.0, radius=0.0, distance=0.0, is_moon=False)
    if total > 0.0 and ice > 0.0:
        if water_fraction is None:
            view.water_mass_fraction = ice * ices.get("H2O", 0.0) / total
            for name, column in SMALL_BODY_ICE_COLUMNS.items():
                setattr(view, column, ice * ices.get(name, 0.0) / total)
        else:
            view.water_mass_fraction = float(water_fraction)
            other = float(ice_fraction or 0.0)
            rest = total - ices.get("H2O", 0.0)
            for name, column in SMALL_BODY_ICE_COLUMNS.items():
                setattr(view, column, other * ices.get(name, 0.0) / rest if rest > 0.0 else 0.0)
    return view


def small_body_record(kind, metallicity_fe_h, age_gy, ice_mass_fraction, distance_km=None, snow_line_au=None, seed="", **fractions):
    """The rows of the isotope record of a small body (`SMALL_BODY_KINDS`; `ROW_FIELDS` dicts), worked out from the values its row stores through `record_of`, the same laws as a planet's; `None` for a body with no composition."""
    if kind not in SMALL_BODY_KINDS:
        raise ValueError(f"no small body of kind {kind!r}")
    return record_of(small_body_view(kind, metallicity_fe_h, age_gy, ice_mass_fraction, distance_km, snow_line_au, seed, **fractions))


def record_of(body):
    """The rows of `body`'s record (`ROW_FIELDS` dicts, by atomic number then mass number), or `None` for a body with no stored composition."""
    held = elements_held(body)
    if not held:
        return None
    age = _age_gy(getattr(body, "star", None))
    rows = []
    for symbol in held:
        origin, mix, initial = NATURAL, _natural(symbol), {}
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
            origin = DECAY
            mix, initial = _decayed(symbol, age)
        for name in sorted(mix, key=lambda isotope_name: elements.isotope(isotope_name).a):
            rows.append({"element": symbol, "isotope": name, "fraction": mix[name], "initial_ratio": initial.get(name), "origin": origin})
    return rows + extinct_rows(body, held, age) + cosmogenic_rows(body)


SHEET_RATIOS = (("D/H", "H-2", "H-1", "Deuterium over hydrogen"), ("13C/12C", "C-13", "C-12", "Carbon-13 over carbon-12"), ("15N/14N", "N-15", "N-14", "Nitrogen-15 over nitrogen-14"),
                ("40Ar/36Ar", "Ar-40", "Ar-36", "Argon-40 over argon-36"), ("235U/238U", "U-235", "U-238", "Uranium-235 over uranium-238"),
                ("40K/39K", "K-40", "K-39", "Potassium-40 over potassium-39"), ("3H/1H", "H-3", "H-1", "Tritium over hydrogen (of the air and ocean)"),
                ("14C/12C", "C-14", "C-12", "Carbon-14 over carbon-12 (of the air and ocean)"))
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
