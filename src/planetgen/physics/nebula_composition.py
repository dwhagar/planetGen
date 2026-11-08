# planetgen/physics/nebula_composition.py

"""
What a nebula is made of (GEN.1269, part a of GEN.1268, reports/gen1268-nebula-composition-2026-11-05.md sections 2, 3, 4 and 11): one model for every
class of `tuning.NEBULA_CLASSES`, from the stored class, the gas metallicity of the place, the density and size, and a key.

`makeup` returns the composition by mass of the gas and dust together (hydrogen, helium, the metals), how much of each metal sits in dust and
how much is left in the gas, the dust-to-gas and dust-to-metals ratios, the grain mix, the split of the hydrogen into ionized, atomic and molecular
shares, and the gas-phase abundances as 12 + log(X/H) of C, N, O, Ne, S, Ar and Fe. It stores nothing and draws only from streams of its own
keyed on the nebula (`gas_oh`), so the same nebula reads the same at any worker count. The composition item (GEN.1276) writes what it returns.

The steps (report section 2):

1. **The gas metallicity of the place** (`gas_oh`): the Milky Way gas has its own gradient (`GAS_GRADIENT_DEX_PER_KPC`, Boss 8:49 PM EST
   2026-11-04: not the stars' 0.06), zero at the galaxy's reference radius, with a scatter drawn from the key. In a satellite component the gas is
   the component's own [Fe/H] plus a young-gas offset drawn inside the widest range of the two published estimates (`SATELLITE_GAS_OFFSETS`).
2. **The total composition** (gas and dust together): the birth composition of a thin-disk star of that metallicity
   (`stellar_abundances.composition`), so a nebula and a star born in it share one chemistry.
3. **The dust** (`dust_to_gas_law`): the broken power law of [R5102] of the oxygen abundance, times the class's survival
   (section 11 of the report: the dust that lives through the ionizing radiation, the shock or the winds of the class).
4. **The grains**: the dust is shared among grain kinds (`GRAIN_MIXES`) of fixed element content (`GRAINS`); no element goes into dust past
   its cap (`MAX_IN_DUST`), so the dust of a metal-poor or a heavily depleted nebula is whatever its metals allow.
5. **The phases** (`PHASES`): the share of the hydrogen that is ionized, atomic and molecular, a table by class; for the cold neutral classes
   (F, M, N, P, Q) the molecular share of the neutral gas comes from the column density and the metallicity (`molecular_fraction`).
6. **Enrichment**: a Wolf-Rayet ring (X) mixes the wind of its host star into its shell (`WIND_SHARE_OF_RING`); a planetary nebula or a remnant
   takes the ejecta of its progenitor (`planetary_nebula_ejecta`, GEN.1270; the remnants' ejecta are GEN.1271) through the `ejecta` argument.
   Without it, such a class reads as swept-up or shed interstellar gas and says so in `open_items`.

Every number that is a choice and not a measurement is a stand-in of the registry, section 6.179 (rows NBC1 to NBC13 of
docs/design/research/gap-filling-plan.md), with its source tier and tweak level.
"""

import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.galaxy import scales
from planetgen.physics import elements as element_data
from planetgen.physics import stellar_abundances as sa
from planetgen.physics import stellar_composition as sc
from planetgen.util import draw

# --- The gas metallicity of the place -----------------------------------------------------------------------------------------------

GAS_GRADIENT_DEX_PER_KPC = 0.034
"""float: The fall of the Milky Way gas's metallicity with galactic radius, dex per kpc: the oxygen gradient of 35 H II regions from 5 to 17 kpc
[R5107] (abstract read by the research thread). Boss's answer: the gas has its own gradient, 0.03 to 0.04. [Stand-in NBC1, tier T2, tweak L0.]"""

GAS_SCATTER_DEX = 0.14
"""float: The scatter of the gas around its gradient, dex (the same 35 H II regions). [Stand-in NBC2, tier T2, tweak L0.]"""

GAS_OH_MAX_DEX = 0.7
"""float: The highest [O/H] the gas of a place can have, dex: the gradient is written against the Milky Way's, so a galaxy with a far larger reference radius (a wide disk) would put
the gas at its center several dex above solar, where the dust law and the composition stop meaning anything (the dust would hold the hydrogen). The cap is above the gradient's value at the Milky Way's center plus three times the scatter (0.28 + 0.42), a recalled bound for the inner disk, not a measurement. [Stand-in NBC16, tier T6, tweak L0.]"""

SATELLITE_GAS_OFFSETS = ((-0.9, 0.24, 0.62), (-0.4, 0.06, 0.26))
"""tuple: `([Fe/H] of a satellite, lowest, highest offset)` of the gas's [O/H] above that [Fe/H], dex, from the smaller to the larger [Fe/H]. The
anchors are the two Clouds: the LMC (stored stand-in [Fe/H] about -0.4) has gas at [O/H] -0.34 to -0.14, the SMC (about -0.9) -0.66 to -0.28,
the lower value from collisional lines and the higher from recombination lines [R5103: Tables 5 to 7 and 10]
(read by the research thread). Boss: draw inside the widest range. Between the anchors the offsets are interpolated, outside them held.
[Stand-in NBC3, tier T2, tweak L2: tuned to the two published ranges.]"""

SOLAR_OXYGEN_12LOG = 8.69
"""float: The solar 12 + log(O/H) the dust law is written against [R5102]."""


def _stream(key):
    return draw.Stream(f"nebula-gas|{key}")


def gas_gradient_oh(galactic_radius_pc):
    """The mean [O/H] of the Milky Way gas at a galactic radius, dex: zero at the reference radius, falling by `GAS_GRADIENT_DEX_PER_KPC`."""
    reference_kpc = scales.reference_radius_pc() / 1000.0
    return -GAS_GRADIENT_DEX_PER_KPC * (galactic_radius_pc / 1000.0 - reference_kpc) + 0.0


def satellite_offsets(component_fe_h):
    """`(lowest, highest)` offset of a satellite's gas [O/H] above the component's [Fe/H] (`SATELLITE_GAS_OFFSETS`)."""
    (fe_a, low_a, high_a), (fe_b, low_b, high_b) = SATELLITE_GAS_OFFSETS
    share = min(max((component_fe_h - fe_a) / (fe_b - fe_a), 0.0), 1.0)
    return low_a + share * (low_b - low_a), high_a + share * (high_b - high_a)


def gas_oh(galactic_radius_pc, key, component_fe_h=None):
    """
    [O/H] of the gas a nebula `key` formed in, dex. In the Milky Way: the gradient at `galactic_radius_pc` plus a normal draw of `GAS_SCATTER_DEX`, not above `GAS_OH_MAX_DEX`.
    In a satellite component of stored [Fe/H] `component_fe_h`: that value plus an offset drawn inside `satellite_offsets`. Both draws are always
    taken, in this order, so a nebula's number never depends on which one is used.
    """
    stream = _stream(key)
    normal, u = stream.gauss(0.0, 1.0), stream.random()
    if component_fe_h is not None:
        low, high = satellite_offsets(component_fe_h)
        return component_fe_h + low + u * (high - low)
    return min(gas_gradient_oh(galactic_radius_pc) + normal * GAS_SCATTER_DEX, GAS_OH_MAX_DEX)


# --- The total composition ------------------------------------------------------------------------------------------------------------

_WEIGHT = {symbol: element_data.element(symbol).atomic_weight for symbol in ("H", "He", *sa.METAL_SYMBOLS)}
METALS = sa.METAL_SYMBOLS
ABUNDANCE_ELEMENTS = ("C", "N", "O", "Ne", "S", "Ar", "Fe")
"""tuple: The elements the abundances (12 + log(X/H)) are given for."""


def _birth(fe_h):
    return sa.composition(fe_h, sa.abundance_offsets(fe_h, sa.BRANCH_THIN, (0.0,) * len(sa.RESIDUAL_ELEMENTS)))


def _number_ratio(elements, symbol):
    return (elements[symbol] / _WEIGHT[symbol]) / (elements["H"] / _WEIGHT["H"])


_SOLAR_O_TO_H = _number_ratio(_birth(0.0).elements, "O")


def oh_dex(elements):
    """[O/H] of a mass-fraction dict against the solar mix of `stellar_abundances`, dex."""
    return math.log10(_number_ratio(elements, "O") / _SOLAR_O_TO_H)


def fe_h_for_oh(target_oh):
    """The [Fe/H] of the thin-disk birth composition whose [O/H] is `target_oh` (the [O/Fe] of the disk is a few hundredths of a dex), found by
    bisection (the fits of [X/Fe] bend, so [O/H] rises with [Fe/H] but not in a straight line)."""
    low, high = target_oh - 1.0, target_oh + 1.0
    for _ in range(50):
        middle = 0.5 * (low + high)
        if oh_dex(_birth(middle).elements) < target_oh:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def oxygen_12log(elements):
    """12 + log(O/H) on the scale of `SOLAR_OXYGEN_12LOG`: the solar mix is 8.69."""
    return SOLAR_OXYGEN_12LOG + oh_dex(elements)


# --- The dust ----------------------------------------------------------------------------------------------------------------------------

DUST_BREAK_12LOG = 7.96
DUST_HIGH_METALLICITY = (2.21, 1.0)
DUST_LOW_METALLICITY = (0.68, 3.08)
"""The broken power law of the gas-to-dust ratio of galaxies [R5102: Table 1] (Galactic X_CO): above
`DUST_BREAK_12LOG`, log(gas/dust) = 2.21 + 1.0 (8.69 - x); at or below it, 0.68 + 3.08 (8.69 - x), with x = 12 + log(O/H). Accuracy a factor 1.6,
scatter 0.37 dex. [Stand-in NBC4, tier T2, tweak L0.]"""


def gas_to_dust_law(x_12log):
    """The gas-to-dust mass ratio of a galaxy of 12 + log(O/H) `x_12log` (the law above)."""
    intercept, slope = DUST_HIGH_METALLICITY if x_12log > DUST_BREAK_12LOG else DUST_LOW_METALLICITY
    return 10.0 ** (intercept + slope * (SOLAR_OXYGEN_12LOG - x_12log))


def dust_to_gas_law(x_12log):
    """The dust-to-gas mass ratio of the law (before the class's survival and the metal caps)."""
    return 1.0 / gas_to_dust_law(x_12log)


MAX_IN_DUST = {"Mg": 0.95, "Si": 0.95, "Fe": 0.99, "C": 0.70, "O": 0.60}
"""dict: The most of an element that can sit in dust, as a share of its total (the rest stays in the gas). Iron is 96 percent in dust in Orion
(1.4 dex below solar, [R5106: Tables 14 and 15]); the others are recalled depletion patterns ([R5105] has the coefficients, not
yet read). [Stand-in NBC5, tier T4, tweak L0.]"""

_SILICATE = {"Mg": 1.5, "Si": 1.0, "O": 3.5}   # number ratios of a pyroxene and olivine mix, without its iron (the iron is a grain of its own)


def _mass_shares(numbers):
    masses = {symbol: count * _WEIGHT[symbol] for symbol, count in numbers.items()}
    total = sum(masses.values())
    return {symbol: mass / total for symbol, mass in masses.items()}


ICE_MANTLE_PER_100_H2O = {"H2O": 100.0, "CO": 25.0, "CO2": 26.0, "CH3OH": 8.0, "NH3": 6.0, "CH4": 4.5}
"""dict: The molecules of an interstellar ice mantle, in molecules per 100 of water: the medians of the ices toward background stars behind quiescent clouds for CO (25), CO2 (26)
and CH3OH (8), and of the low-mass young stellar objects for NH3 (6) and CH4 (4.5), where the quiescent clouds give an upper limit only (< 7 and < 3) [R0461: Table 2, read from the
arXiv text]. One mantle for every class that has ice. [Stand-in NBC17, tier T2, tweak L1: the medians of two environments joined.]"""

_MOLECULES = {"H2O": {"H": 2, "O": 1}, "CO": {"C": 1, "O": 1}, "CO2": {"C": 1, "O": 2}, "CH3OH": {"C": 1, "H": 4, "O": 1}, "NH3": {"N": 1, "H": 3}, "CH4": {"C": 1, "H": 4}}
"""dict: The atoms of each ice molecule."""


def _molecule_mass(name):
    return sum(count * _WEIGHT[symbol] for symbol, count in _MOLECULES[name].items())


def _ice_mantle():
    masses = {name: count * _molecule_mass(name) for name, count in ICE_MANTLE_PER_100_H2O.items()}
    total = sum(masses.values())
    return {name: mass / total for name, mass in masses.items()}


ICE_MANTLE = _ice_mantle()
"""dict: The mass share of each molecule in the ice mantle of a grain (add to 1)."""


def _ice_elements():
    elements = {}
    for name, share in ICE_MANTLE.items():
        molecule_mass = _molecule_mass(name)
        for symbol, count in _MOLECULES[name].items():
            elements[symbol] = elements.get(symbol, 0.0) + share * count * _WEIGHT[symbol] / molecule_mass
    return elements


GRAINS = {
    "silicate": _mass_shares(_SILICATE),
    "iron": {"Fe": 1.0},
    "carbon": {"C": 1.0},
    "PAH": {"C": 0.92, "H": 0.08},
    "ice": _ice_elements(),
    "SiC": _mass_shares({"Si": 1.0, "C": 1.0}),
}
"""dict: The element mass shares of each grain kind: a magnesium silicate (Mg 1.5, Si 1, O 3.5 by number, the average of olivine and pyroxene), iron
(metal and sulfide; the iron of the silicates is here, so the dust follows the iron a metal-poor gas has), amorphous carbon and graphite together, a PAH
(92 percent carbon by mass, the rest hydrogen), an ice mantle (`ICE_MANTLE`: water with CO, CO2, methanol, ammonia and methane), and silicon carbide.
[Stand-in NBC6, tier T6, tweak L0.]"""

SIC_SI_SHARE_OF_SILICATE_SI = 0.025
"""float: The share of the silicon of the silicate grains that sits in SiC in the interstellar mixes. The absence of the 11.3 micron SiC band toward the Galactic center puts the
silicon in SiC dust at less than 5 percent of the silicon in silicates [R5144: section 5.3]; the mixes take half the limit. [Stand-in NBC18, tier T2, tweak
L1: inside a published upper limit.]"""

GRAIN_MIXES = {
    "diffuse": {"silicate": 0.42, "iron": 0.20, "carbon": 0.3065, "PAH": 0.07, "SiC": 0.0035},
    "dense": {"silicate": 0.36, "iron": 0.17, "carbon": 0.197, "PAH": 0.04, "ice": 0.23, "SiC": 0.003},
    "reflection": {"silicate": 0.39, "iron": 0.18, "carbon": 0.2368, "PAH": 0.09, "ice": 0.10, "SiC": 0.0032},
    "oxygen-rich": {"silicate": 0.50, "iron": 0.23, "carbon": 0.20, "PAH": 0.07},
    "carbon-rich": {"silicate": 0.08, "iron": 0.02, "carbon": 0.60, "PAH": 0.20, "SiC": 0.10},
    "ejecta": {"silicate": 0.45, "carbon": 0.20, "iron": 0.35},
    "iron-group": {"silicate": 0.20, "iron": 0.80},
}
"""dict: The shares of the dust mass in each grain kind, by mix (a class names one in `CLASSES`); each adds to 1. The diffuse mix is a Milky Way
one (silicates, iron, carbon with a PAH tail of a few percent, and SiC at half the interstellar limit `SIC_SI_SHARE_OF_SILICATE_SI`); the dense one gains ice mantles; the
carbon-rich one (the winds of carbon stars) has SiC at 10 percent of the dust, the share below which 37 of 44 galactic carbon stars lie [R5146: the text];
the others follow the report's table. [Stand-ins NBC6 and NBC19, tiers T6 and T2, tweaks L3 and L1: the first chosen so that a solar-metallicity diffuse nebula's dust fits inside
its metals within a fifth of the law, the second inside the published carbon-star range.]"""

assert all(abs(sum(mix.values()) - 1.0) < 1e-12 for mix in GRAIN_MIXES.values()), "a grain mix must add to 1"

OLIVINE_SHARE_OF_AMORPHOUS_SILICATE = 0.849
"""float: The mass share of amorphous olivine in the amorphous silicate (the rest is amorphous pyroxene): 84.9 percent olivine and 15.1 percent pyroxene by mass fit the 10
micron band of the diffuse interstellar medium toward the Galactic center [R5143: abstract]. A split of the mass only; the element content of the silicate kind stays that of
`GRAINS`. [Stand-in NBC20, tier T1, tweak L0.]"""

CRYSTALLINE_SHARE_OF_SILICATE = {
    "diffuse": 0.002, "dense": 0.002, "reflection": 0.002, "oxygen-rich": 0.15, "carbon-rich": 0.15, "ejecta": 0.15, "iron-group": 0.15,
}
"""dict: By mix, the mass share of the silicate that is crystalline (forsterite and enstatite). The diffuse medium holds 0.2 +- 0.2 percent [R5143: abstract]; the silicates of
stellar ejecta are 11 to 18 percent crystalline (about 10 percent around OH/IR stars, 15 to 20 around red supergiants) [R5143: section 4.1], so the mixes of stellar ejecta take 15
percent. The ejecta of a supernova and the carbon-rich winds take the same value for want of a source of their own. [Stand-in NBC21, tier T2, tweak L1.]"""

GRAPHITE_SHARE_OF_CARBON = {
    "diffuse": 0.5, "dense": 0.25, "reflection": 0.5, "oxygen-rich": 0.5, "carbon-rich": 0.1, "ejecta": 0.3, "iron-group": 0.0,
}
"""dict: By mix, the mass share of the carbon grains (not the PAHs) that is graphite; the rest is amorphous carbon. Two models of the diffuse medium disagree: the
silicate-graphite-PAH model puts all its large carbon in graphite [R5145], the astrodust model has none [R5147], so the diffuse mix takes the middle of the widest range; a
dense cloud and a carbon star make amorphous carbon (a carbon star's dust is amorphous carbon with SiC); presolar graphite is a small part of the grains of a meteorite
(1 ppm against 6 ppm of SiC) [R5144: Table 3]. [Stand-in NBC22, tier T6, tweak L3 where a mix is not the diffuse one, T3 and L1 where it is.]"""

assert set(CRYSTALLINE_SHARE_OF_SILICATE) == set(GRAPHITE_SHARE_OF_CARBON) == set(GRAIN_MIXES), "every mix needs a crystalline and a graphite share"
SURVIVES_AS = {"SiC": "silicate"}
"""dict: A grain kind with no survival share of its own is read as another's: SiC survives like the silicates. [Stand-in NBC23, tier T6, tweak L3: no source.]"""


def molecular_column_scale():
    """The mass of a hydrogen nucleus times 1.4 (helium), g: the mass per hydrogen atom of the gas."""
    return 1.4 * 1.6735575e-24


# --- The phases of the hydrogen -------------------------------------------------------------------------------------------------------

PHASES = {
    "A": (0.10, 0.90, 0.0), "B": (0.90, 0.10, 0.0),
    "C": (0.95, 0.05, 0.0), "D": (0.97, 0.03, 0.0), "E": (0.97, 0.03, 0.0), "G": (0.60, 0.20, 0.20), "X": (0.95, 0.05, 0.0),
    "F": (0.0, None, None), "M": (0.0, None, None), "N": (0.0, None, None), "P": (0.0, None, None), "Q": (0.0, None, None),
    "H": (0.80, 0.15, 0.05), "J": (0.80, 0.15, 0.05), "K": (0.80, 0.15, 0.05), "L": (0.80, 0.15, 0.05),
    "R": (1.0, 0.0, 0.0), "S": (1.0, 0.0, 0.0), "T": (1.0, 0.0, 0.0), "U": (1.0, 0.0, 0.0), "V": (0.5, 0.5, 0.0), "W": (1.0, 0.0, 0.0),
}
"""dict: `(ionized, atomic, molecular)` shares of a class's hydrogen mass; `None` for the atomic and molecular share of a cold neutral class means
the neutral gas (1 - ionized) is split by `molecular_fraction`. A and B follow the warm neutral and warm ionized medium (recalled); the H II classes are
Stromgren spheres; G has an ionized core in a dusty rim; the planetary nebulae keep a molecular envelope; hot remnants are plasma, an old
radiative remnant half neutral; a remnant of Type Ia (W) holds no hydrogen at all, which the ejecta set. [Stand-in NBC10, tier T6, tweak L3; A and
B tier T5, tweak L1.]"""

KMT_CHI_COEFFICIENT = 0.77
KMT_CHI_METALLICITY = (3.1, 0.365)
KMT_LOG_TERMS = (0.6, 0.01)
KMT_S_COEFFICIENT = 0.04
KMT_FRACTION_TERMS = (0.75, 0.25)
"""The molecular fraction of a cloud of column density Sigma (solar masses per pc2) and metallicity Z' = Z/Z0, after [R5104]: f = 1 - 0.75 s / (1 + 0.25 s) for s < 2 and 0 above, with s = ln(1 + 0.6 chi + 0.01 chi^2) / (0.04 Sigma Z') and
chi = 0.77 (1 + 3.1 Z'^0.365). The chi coefficient is checked against its eq. 7 (2.3/3 = 0.77); the f form is recalled and its half-molecular column
(34 solar masses per pc2 at solar Z') is near the 27 of the paper's summary. [Stand-in NBC11, tier T2, tweak L0.]"""

_PC_CM = 3.0856775814913673e18
_MSUN_PER_PC2_IN_G_PER_CM2 = 1.98847e33 / _PC_CM ** 2


def column_density_msun_pc2(density_cm3, radius_ly):
    """The mean surface density of a uniform sphere of hydrogen density `density_cm3` (per cm3) and radius `radius_ly`: 4/3 n R of mass 1.4 m_H per
    hydrogen atom, in solar masses per pc2. [Stand-in NBC12, tier T6, tweak L0: a uniform cloud; the stored density is the cloud's mean.]"""
    radius_cm = radius_ly * 9.4607304725808e17
    grams = (4.0 / 3.0) * density_cm3 * radius_cm * molecular_column_scale()
    return grams / _MSUN_PER_PC2_IN_G_PER_CM2


def molecular_fraction(sigma_msun_pc2, z_over_solar):
    """The molecular share of a neutral cloud's hydrogen (the law above); 0 for a column of zero or a metal-free gas."""
    if sigma_msun_pc2 <= 0.0 or z_over_solar <= 0.0:
        return 0.0
    chi = KMT_CHI_COEFFICIENT * (1.0 + KMT_CHI_METALLICITY[0] * z_over_solar ** KMT_CHI_METALLICITY[1])
    s = math.log(1.0 + KMT_LOG_TERMS[0] * chi + KMT_LOG_TERMS[1] * chi * chi) / (KMT_S_COEFFICIENT * sigma_msun_pc2 * z_over_solar)
    if s >= 2.0:
        return 0.0
    return 1.0 - KMT_FRACTION_TERMS[0] * s / (1.0 + KMT_FRACTION_TERMS[1] * s)


# --- The classes ------------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class ClassRow:
    """
    How a nebula or remnant class differs: its grain mix and the share of each grain kind that survives (1.0 keeps the full galactic dust ratio).

    Attributes:
        mix (str): A key of `GRAIN_MIXES`.
        survival (dict): Grain kind -> share of that grain's dust that survives; a kind not listed is 1.0.
        note (str): What the row stands for.
    """

    mix: str
    survival: dict
    note: str


_IONIZED_SURVIVAL = {"silicate": 0.9, "iron": 0.9, "carbon": 0.5, "PAH": 0.5}
CLASSES = {
    "A": ClassRow("diffuse", {}, "the diffuse neutral medium keeps its dust"),
    "B": ClassRow("diffuse", {k: 0.7 for k in ("silicate", "carbon", "PAH", "iron")}, "warm ionized gas loses some dust"),
    "C": ClassRow("diffuse", {**_IONIZED_SURVIVAL, "PAH": 0.2}, "a compact H II region destroys its PAHs"),
    "D": ClassRow("diffuse", _IONIZED_SURVIVAL, "iron-bearing grains survive in H II regions, the fragile carbon half does"),
    "E": ClassRow("diffuse", _IONIZED_SURVIVAL, "as D"),
    "F": ClassRow("reflection", {}, "a cold dusty cloud keeps all its dust"),
    "G": ClassRow("reflection", _IONIZED_SURVIVAL, "an ionized core in a dusty rim"),
    "X": ClassRow("diffuse", _IONIZED_SURVIVAL, "a swept-up shell with the wind of a Wolf-Rayet star mixed in"),
    "H": ClassRow("oxygen-rich", _IONIZED_SURVIVAL, "iron-rich dust survives in a planetary nebula (more than 90 percent of its iron on grains)"),
    "J": ClassRow("oxygen-rich", _IONIZED_SURVIVAL, "as H"),
    "K": ClassRow("carbon-rich", _IONIZED_SURVIVAL, "carbon-rich gas makes carbon dust and PAHs"),
    "L": ClassRow("oxygen-rich", _IONIZED_SURVIVAL, "as H"),
    "M": ClassRow("dense", {}, "grains grow in the dense gas; none is lost"),
    "N": ClassRow("dense", {}, "as M"),
    "P": ClassRow("dense", {}, "as M, with organics in the carbon"),
    "Q": ClassRow("dense", {}, "as M, with organics in the carbon"),
    "R": ClassRow("ejecta", {k: 0.2 for k in ("silicate", "carbon", "iron")}, "fresh ejecta dust, mostly destroyed by the reverse shock (the survival is a stand-in until the dust yields of the ejecta are modeled)"),
    "S": ClassRow("diffuse", {"silicate": 0.39, "iron": 0.60, "carbon": 0.5, "PAH": 0.1}, "shocked interstellar dust: 61 percent of the Si-rich and 40 percent of the Fe-rich grains destroyed"),
    "T": ClassRow("ejecta", {k: 0.2 for k in ("silicate", "carbon", "iron")}, "filament dust of the ejecta (the survival is a stand-in until the dust yields of the ejecta are modeled)"),
    "U": ClassRow("ejecta", {k: 0.2 for k in ("silicate", "carbon", "iron")}, "filament dust of the ejecta (the survival is a stand-in until the dust yields of the ejecta are modeled)"),
    "V": ClassRow("diffuse", {k: 0.5 for k in ("silicate", "carbon", "PAH", "iron")}, "a cooling shell sputters about half its grains"),
    "W": ClassRow("iron-group", {k: 0.2 for k in ("silicate", "iron")}, "iron-group dust of a Type Ia ejecta (the survival is a stand-in until the dust yields of the ejecta are modeled)"),
}
"""dict: One `ClassRow` for every class of `tuning.NEBULA_CLASSES`. The survivals come from the dust-destruction research of the report (section 11):
more than 90 percent of the iron of planetary nebulae and H II regions is on grains [R5108]; a shocked LMC remnant keeps 39 percent of the
Si-rich and 60 percent of the Fe-rich grains [R5109]; the 0.5 of the fragile carbon, the PAHs, the old remnants and the fresh-ejecta
dust are stand-ins until the destruction literature is read (the ejecta dust is GEN.1271's). [Stand-ins NBC7 (sourced), NBC8 and NBC9 (open): tiers T2, T4 and
T6, tweak L1, L0 and L3.]"""

WIND_SHARE_OF_RING = 0.3
"""float: The share of a Wolf-Rayet ring's mass that is the host star's wind (helium and nitrogen rich); the rest is swept-up interstellar gas.
[Stand-in NBC13, tier T6, tweak L3: the observed rings are mostly swept-up gas with a wind-enriched helium and nitrogen.]"""

EJECTA_CLASSES = ("H", "J", "K", "L", "R", "T", "U", "W")
"""tuple: The classes whose gas is the ejecta of a star and not interstellar: they take `ejecta`."""

assert set(CLASSES) == set(tuning.NEBULA_CLASSES) == set(PHASES), "every class of the table needs a row"


# --- The result ---------------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class NebulaMakeup:
    """
    What a nebula is made of.

    Attributes:
        nebula_class (str): The class letter.
        oh_dex (float): [O/H] of the gas it formed in, dex.
        fe_h (float): [Fe/H] of the thin-disk composition that gave it that oxygen.
        x, y, z (float): Hydrogen, helium and metal mass fractions of gas and dust together (add to 1).
        elements (dict): The mass fraction of each element, gas and dust together (add to 1).
        dust_elements (dict): The mass fraction of the total that sits in dust, by element.
        gas_elements (dict): The mass fraction of the gas phase alone, by element (add to 1).
        dust_to_gas (float): Dust mass over gas mass, after survival and the metal caps.
        dust_to_gas_law (float): The same before them (`dust_to_gas_law` of the nebula's oxygen).
        dust_to_metals (float): Dust mass over metal mass.
        grains (dict): The share of the dust mass in each grain kind (`GRAINS`).
        species (dict): The share of the dust mass in each grain species (`SPECIES`): the grain kinds with the silicate split into amorphous olivine, amorphous pyroxene and crystalline
            silicate and the carbon into amorphous carbon and graphite (add to 1 where there is dust).
        ionized, atomic, molecular (float): The shares of the hydrogen mass in each phase (add to 1).
        abundances (dict): 12 + log(X/H) of the gas phase for `ABUNDANCE_ELEMENTS`.
        open_items (tuple): What the model leaves open for this class.
    """

    nebula_class: str
    oh_dex: float
    fe_h: float
    x: float
    y: float
    z: float
    elements: dict
    dust_elements: dict
    gas_elements: dict
    dust_to_gas: float
    dust_to_gas_law: float
    dust_to_metals: float
    grains: dict
    species: dict
    ionized: float
    atomic: float
    molecular: float
    abundances: dict
    open_items: tuple

    def values(self):
        """The values the composition item stores, as a dict."""
        return {"mass_fraction_h": self.x, "mass_fraction_he": self.y, "mass_fraction_z": self.z,
                **{f"mass_fraction_{e.lower()}": self.elements[e] for e in ABUNDANCE_ELEMENTS},
                "dust_to_gas": self.dust_to_gas, "dust_to_metals": self.dust_to_metals,
                "ionized_share": self.ionized, "atomic_share": self.atomic, "molecular_share": self.molecular,
                **{f"abundance_{e.lower()}": self.abundances[e] for e in ABUNDANCE_ELEMENTS}}


def _wind_elements(subtype):
    surface = sc.wolf_rayet_surface(subtype, 0.5)
    wind = {symbol: 0.0 for symbol in ("H", "He", *METALS)}
    wind["other metals"] = 0.0
    for name, value in surface.items():
        wind["He" if name == "He and O" else name] += value
    return wind


def _mix(first, second, share):
    return {symbol: (1.0 - share) * first.get(symbol, 0.0) + share * second.get(symbol, 0.0) for symbol in set(first) | set(second)}


def _dust_of(elements, law_ratio, row):
    """`(dust by element, grain mass, requested dust mass)` for gas and dust that together hold `elements`, the full dust-to-gas `law_ratio`
    and the survival of the class: each grain kind is asked for its share of the surviving dust, then every element is held to its cap."""
    mix = GRAIN_MIXES[row.mix]
    ratio = law_ratio / (1.0 + law_ratio)    # dust mass over the total
    wanted = {kind: ratio * share * row.survival.get(kind, row.survival.get(SURVIVES_AS.get(kind), 1.0)) for kind, share in mix.items()}
    demand = {}
    for kind, mass in wanted.items():
        for symbol, part in GRAINS[kind].items():
            demand[symbol] = demand.get(symbol, 0.0) + mass * part
    scale = {}
    for symbol, limit in MAX_IN_DUST.items():
        have = elements.get(symbol, 0.0)
        if demand.get(symbol, 0.0) > limit * have:
            scale[symbol] = limit * have / demand[symbol]
    grains = {}
    for kind, mass in wanted.items():
        factor = min([scale.get(symbol, 1.0) for symbol in GRAINS[kind]] + [1.0])
        grains[kind] = mass * factor
    dust = {}
    for kind, mass in grains.items():
        for symbol, part in GRAINS[kind].items():
            dust[symbol] = dust.get(symbol, 0.0) + mass * part
    # hydrogen in PAHs and ice comes from the gas, never past what the nebula holds
    if dust.get("H", 0.0) > 0.5 * elements["H"]:
        raise ValueError("dust cannot hold half the hydrogen of a nebula")
    return dust, grains, ratio


SPECIES = ("amorphous_olivine", "amorphous_pyroxene", "crystalline_silicate", "iron", "amorphous_carbon", "graphite", "PAH", "SiC", "ice")
"""tuple: The grain species the pages and the stored columns name."""


def species_of(mix, shares):
    """The share of the dust in each of `SPECIES` for the grain kind shares `shares` of mix `mix`; a kind the dust does not hold is 0. Empty for dust-free `shares`."""
    if not shares:
        return {}
    silicate = shares.get("silicate", 0.0)
    carbon = shares.get("carbon", 0.0)
    crystalline = silicate * CRYSTALLINE_SHARE_OF_SILICATE[mix]
    graphite = carbon * GRAPHITE_SHARE_OF_CARBON[mix]
    amorphous = silicate - crystalline
    return {"amorphous_olivine": amorphous * OLIVINE_SHARE_OF_AMORPHOUS_SILICATE, "amorphous_pyroxene": amorphous * (1.0 - OLIVINE_SHARE_OF_AMORPHOUS_SILICATE),
            "crystalline_silicate": crystalline, "iron": shares.get("iron", 0.0), "amorphous_carbon": carbon - graphite, "graphite": graphite,
            "PAH": shares.get("PAH", 0.0), "SiC": shares.get("SiC", 0.0), "ice": shares.get("ice", 0.0)}


def makeup(nebula_class, gas_oh_dex, density_cm3, radius_ly, ejecta=None, wolf_rayet_subtype=None):
    """
    The `NebulaMakeup` of a nebula or remnant of class `nebula_class`, formed in gas of [O/H] `gas_oh_dex` (`gas_oh`), of mean hydrogen density `density_cm3`
    (per cm3) and radius `radius_ly`.

    Args:
        ejecta (dict, optional): Mass fractions of H, He and the metals of the gas a planetary nebula or a remnant was made of (a
            `planetary_nebula_ejecta.PlanetaryNebulaEjecta.composition.elements`); given, it replaces the interstellar composition. Without it, classes H to
            L and R to W read as interstellar gas and say so in `open_items`.
        wolf_rayet_subtype (str, optional): For class X, the subtype of the host star (`WN`, `WC`, `WO`, with or without a subclass number);
            `WN` without a number is taken as a late WN star.

    Raises:
        ValueError: An unknown class.
    """
    if nebula_class not in CLASSES:
        raise ValueError(f"unknown nebula class {nebula_class!r}; expected one of {', '.join(sorted(CLASSES))}")
    row = CLASSES[nebula_class]
    fe_h = fe_h_for_oh(gas_oh_dex)
    interstellar = dict(_birth(fe_h).elements)
    open_items = []
    if nebula_class in EJECTA_CLASSES:
        if ejecta is None:
            open_items.append("no progenitor ejecta given: read as interstellar gas (the ejecta of a planetary nebula and of a remnant are modeled separately)")
            elements = interstellar
        else:
            elements = {symbol: ejecta.get(symbol, 0.0) for symbol in interstellar}
    elif nebula_class == "X":
        wind = _wind_elements(wolf_rayet_subtype or "WN-late")
        elements = _mix(interstellar, wind, WIND_SHARE_OF_RING)
    else:
        elements = interstellar
    total = sum(elements.values())
    elements = {symbol: value / total for symbol, value in elements.items()}
    x_12log = oxygen_12log(interstellar)       # the place's oxygen, so a hydrogen-free ejecta still has a dust law to scale
    law = dust_to_gas_law(x_12log)
    dust, grains, requested = _dust_of(elements, law, row)
    dust_total = sum(dust.values())
    gas = {symbol: elements[symbol] - dust.get(symbol, 0.0) for symbol in elements}
    gas_total = sum(gas.values())
    gas_fractions = {symbol: value / gas_total for symbol, value in gas.items()}
    abundances = {symbol: (12.0 + math.log10(max(_number_ratio(gas_fractions, symbol), 1e-30)) if gas_fractions["H"] > 0 else math.nan)
                  for symbol in ABUNDANCE_ELEMENTS}
    metals = 1.0 - elements["H"] - elements["He"]
    ionized, atomic, molecular = PHASES[nebula_class]
    if atomic is None:
        z_over_solar = max(metals, 0.0) / sa.PROTO_Z
        sigma = column_density_msun_pc2(density_cm3, radius_ly)
        molecular = (1.0 - ionized) * molecular_fraction(sigma, z_over_solar)
        atomic = 1.0 - ionized - molecular
    grain_total = sum(grains.values())
    shares = {kind: mass / grain_total for kind, mass in grains.items()} if grain_total > 0 else {}
    return NebulaMakeup(
        nebula_class, gas_oh_dex, fe_h, elements["H"], elements["He"], metals, elements, dust, gas_fractions,
        dust_total / gas_total, law, dust_total / metals if metals > 0 else 0.0,
        shares, species_of(row.mix, shares), ionized, atomic, molecular, abundances, tuple(open_items))


def field_cloud_density_cm3(nebula_class, key):
    """The mean density of a field cloud (classes M and N) of cell `key`: a log-uniform draw inside the class's range, from a stream of its own, so
    every cloud differs and each reads the same every time (Boss: varied composition, keyed on its cell)."""
    low, high = tuning.NEBULA_CLASSES[nebula_class]["density_range_cm3"]
    u = draw.Stream(f"nebula-field-density|{key}").random()
    return math.exp(math.log(low) + u * (math.log(high) - math.log(low)))
