# planetgen/physics/stellar_abundances.py

"""
The stellar composition relations (GEN.439, docs/design/composition/stellar-composition.md)
=================================================================================

What a star, a brown dwarf or a remnant is made of, as functions of the values the generator already stores: [Fe/H], the
age, the stellar branch and the mass. Constants and pure functions only: no database, no random number (the draws the
relations need come in as arguments, so the derived function of GEN.440 can key them to the star), nothing copied from one
star. The Sun appears as the zero point of the reference scale and as a test target, never as a preset.

The chain:

1. **The scale.** Asplund, Amarsi and Grevesse 2021 log eps (`tuning.STELLAR_SOLAR_LOG_EPS`) moved to the proto-solar Sun give
   the mass fraction of each element at [Fe/H] = 0 (`PROTO_MASS_FRACTION`), the proto-solar X, Y, Z and the helium-to-metal
   ratio `HELIUM_PER_METAL` that carries the scale from primordial helium (Planck 2018) to the proto-solar value.
2. **The branch.** A star younger than `tuning.STELLAR_THICK_START_GY` is thin-disc; an older one is high-alpha with the
   logistic probability of [Fe/H] (`high_alpha_probability`, Bensby 2014). `BRANCH_TABLE` is that choice as a row of the
   gated possibility-table loader (GEN.472), in the order that lets one fixed uniform number decide it: it falls to
   high-alpha when the number is below the probability, so the branch flips at most once as [Fe/H] falls and an edited
   [Fe/H] never makes a star flicker.
3. **The offsets.** The [X/Fe] of C, O, Mg, Si and Ca are the branch's mean at [Fe/H] (`mean_xfe`, hinge fits of Bensby 2014
   and Nissen 2014) plus, for C, O, Mg and Si, a residual: four standard normals through the Cholesky factor of the residual
   correlation (Brewer and Fischer 2016) times the branch's scatter. S and Ar follow Ca; Ne follows O (or Mg, a switch);
   everything else is 0. Below [Fe/H] -2 a carbon-enhanced share of stars has [C/Fe] of 0.7 or more.
4. **The fractions.** `composition` turns [Fe/H] and the offsets into X, Y, Z and each element's mass fraction, and
   `mass_ratio` and `number_ratio` give Fe/Si, Mg/Si, C/O and O/Si, which the planets inherit (GEN.428).
5. **The rest.** Surface changes of the first dredge-up, the radioactive birth abundances (Eu, Th and U for GEN.459), the white
   dwarf core type, the neutron-star constants, the black hole's formulas and the brown dwarf's limits and radius.

What the downloaded data files give (GEN.548, wired in by GEN.590) is not here: the main-sequence luminosity, temperature,
radius and lifetime (MIST), the white dwarf's radius and layers (Bedard et al. 2020) and the neutron star's outer crust
(Servais and Chamel 2026) are read in `stellar_composition` from `planetgen.datasets`, which refuses to answer when the files
were not downloaded (bulk generation refuses to start in that case, `run_galaxy.require_datasets`).

`RELATIONS` lists every relation with its source and its verification state (verified, partly verified, open or stand-in),
so the open ones are found and replaced.
"""

import bisect
import functools
import json
import math
import os
from dataclasses import dataclass

from planetgen import tuning
from planetgen.datasets import baraffe_isochrones
from planetgen.generation import possibility_table as pt
from planetgen.physics import constants
from planetgen.physics import elements, unit_registry as units
from planetgen.util import draw

_kg_to_msun = units.converter("kg", "Msun")
_msun_to_kg = units.converter("Msun", "kg")

BRANCH_THIN = "thin"
BRANCH_HIGH_ALPHA = "high-alpha"
BRANCHES = (BRANCH_THIN, BRANCH_HIGH_ALPHA)
"""tuple: The two abundance sequences of the disc and halo."""

RESIDUAL_ELEMENTS = ("C", "O", "Mg", "Si", "Ca")
"""tuple: The elements that carry a correlated residual, in the order of `tuning.STELLAR_XFE_CORRELATION`."""

DATA_DIRECTORY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

VERIFIED = pt.VERIFIED
PARTLY_VERIFIED = pt.PARTLY_VERIFIED
OPEN = pt.OPEN_STATE
STAND_IN = pt.STAND_IN


class StellarAbundanceError(ValueError):
    """A request outside what a relation covers, or a data file that is incomplete."""


@dataclass(frozen=True)
class Relation:
    """
    A relation of this module with where it comes from.

    Attributes:
        name (str): What it is.
        source (str): The paper or the derivation.
        state (str): `VERIFIED`, `PARTLY_VERIFIED`, `OPEN` or `STAND_IN`.
        limit (str): For a partly verified relation the limit the check named; for an open or stand-in one what replaces it.
    """

    name: str
    source: str
    state: str
    limit: str = ""


RELATIONS = (
    Relation("reference scale: proto-solar X, Y, Z and per-element mass fractions", "Asplund, Amarsi and Grevesse 2021 Table 2 "
             "and section 5 (diffusion offsets +0.064 dex metals, +0.070 dex helium)", VERIFIED),
    Relation("primordial helium", "Planck 2018 VI eq. 72a, arXiv:1807.06209", VERIFIED),
    Relation("helium-to-metal ratio dY/dZ", "derived from the proto-solar and primordial values; Casagrande et al. 2007 observe "
             "2.1 +- 0.9", VERIFIED),
    Relation("thin-disc age gate", "Haywood et al. 2013: the thick disc formed from about 13 to about 8 Gy ago", PARTLY_VERIFIED,
             "the age scale is uncertain by 1 to 1.5 Gy; range 6.5 to 9.5 Gy"),
    Relation("high-alpha probability logistic in [Fe/H]", "Bensby, Feltzing and Oey 2014", PARTLY_VERIFIED,
             "the paper is not among the repository's papers; the fit is part A 2.6's"),
    Relation("[X/Fe] means of Mg, Si, Ca, Al by branch", "Bensby, Feltzing and Oey 2014, hinge fits (part A 2.7)",
             PARTLY_VERIFIED, "valid inside the fitted [Fe/H] ranges; the paper is not among the repository's papers"),
    Relation("[X/Fe] means of C and O by branch", "Nissen et al. 2014 NLTE (part A 2.7)", PARTLY_VERIFIED,
             "the two oxygen sets differ by up to 0.1 dex on the plateau; neon's physics is open"),
    Relation("[X/Fe] residual scatter and correlation", "Brewer and Fischer 2016 (844 dwarfs)", VERIFIED),
    Relation("carbon-enhanced metal-poor share", "Frebel and Norris 2015: 20 to 30 percent below [Fe/H] -2", VERIFIED),
    Relation("neon follows oxygen (or magnesium)", "the solids-law filing prefers oxygen; part A 3.2 step 4 prefers magnesium", OPEN,
             "neon is not measured in photospheres (set values 7.93 to 8.15); both options are kept, 8 percent apart in Z"),
    Relation("[Eu/Fe] by branch", "Battistini and Bensby 2016 (part A 2.7)", PARTLY_VERIFIED, "paper not among the papers"),
    Relation("[Th/Eu] at zero age", "Botelho et al. 2019 Table 4", PARTLY_VERIFIED,
             "Th is measured near solar [Fe/H] only (+-0.1 dex), so Th elsewhere is inferred through Eu"),
    Relation("first dredge-up", "Boothroyd and Sackmann 1999 Tables 2 and 4", PARTLY_VERIFIED,
             "1 to 2.5 solar masses, Z 0.001 to 0.02; nitrogen at low Z not extracted, held at the Z = 0.02 factor"),
    Relation("Wolf-Rayet surface fractions", "Crowther 2007 section 3.5", PARTLY_VERIFIED, "the helium and oxygen split is open"),
    Relation("main-sequence luminosity, temperature, radius and lifetime", "MIST v1.2 basic isochrones, v/vcrit 0.4 (Choi et al. "
             "2016; Dotter 2016), downloaded at setup (SC1, datasets.isochrones)", PARTLY_VERIFIED,
             "solar-scaled, 0.1 solar masses and up, 0.1 Myr to about 20 Gy; a star's present mass stands for its initial one"),
    Relation("white dwarf radius against mass and temperature", "Bedard et al. 2020 cooling sequences, downloaded at setup (SC2, "
             "datasets.white_dwarf_models)", PARTLY_VERIFIED,
             "carbon-oxygen cores, 0.20 to 1.30 solar masses; a helium-core or oxygen-neon white dwarf takes the same sequences"),
    Relation("white dwarf radius without the downloaded sequences", "Chandrasekhar equation of state integrated by "
             "scripts/build_white_dwarf_table.py; Sirius B (Joyce et al. 2018)", STAND_IN,
             "taken only by a run told to disregard missing data; no Coulomb or temperature correction: 3 to 5 percent"),
    Relation("white dwarf core type by mass", "Althaus et al. 2010 (oxygen-neon above 1.05 solar masses)", PARTLY_VERIFIED,
             "the helium-core limit of 0.45 is recalled; the carbon-oxygen split and the layer masses are open"),
    Relation("white dwarf atmosphere shares and pollution", "Xu, Rogers and Blouin 2024", VERIFIED),
    Relation("white dwarf layers", "Bedard et al. 2020 (thick: q_H 1e-4, thin: q_H 1e-10, q_He 1e-2; SC2)", PARTLY_VERIFIED,
             "the sequences assume an even C:O core, not a measured split; none for a helium-core or oxygen-neon white dwarf"),
    Relation("neutron star mass populations", "Ozel and Freire 2016; Romani 2022", VERIFIED),
    Relation("neutron star radius band", "Miller et al. 2021 (hull of the 1-sigma intervals)", VERIFIED),
    Relation("neutron star outer crust sequence", "Servais and Chamel 2026, 20 nuclear mass models, downloaded at setup (SC3, "
             "datasets.neutron_star_matter); it replaced Chamel and Haensel 2008 Table 2", PARTLY_VERIFIED,
             "one mass model is drawn per star; the unmagnetised files only"),
    Relation("neutron star outer crust without the downloaded tables", "Chamel and Haensel 2008 Table 2", STAND_IN,
             "taken only by a run told to disregard missing data; one sequence instead of 20 nuclear mass models"),
    Relation("neutron star inner crust and core", "Malik, Pais and Providencia 2024, 21 equations of state, downloaded (SC3)",
             OPEN, "a star's radius is read from the integrated mass-radius curves of the 21 tables, weighted by NICER (physics/neutron_star.py); the inner crust and core composition is not drawn from them yet"),
    Relation("black hole quantities", "Kerr-Newman formulas (Bardeen, Press and Teukolsky 1972); the no-hair theorem "
             "(Chrusciel, Costa and Heusler 2012)", VERIFIED),
    Relation("brown dwarf deuterium limit", "Spiegel, Burrows and Milsom 2011 Table 1", VERIFIED),
    Relation("brown dwarf lithium and hydrogen limits", "Burrows et al. 2001", VERIFIED, "at solar metallicity"),
    Relation("brown dwarf hydrogen limit against metallicity", "Burrows et al. 2001 (0.092 at Z = 0)", STAND_IN,
             "linear in Z between the two values; the MIST grid (SC1) starts at 0.1 solar masses and cannot give it"),
    Relation("brown dwarf radius", "Baraffe et al. 2015 (BHAC15) track at 1 Gyr, solar Z", STAND_IN,
             "one age and one metallicity; the MIST grid (SC1) starts at 0.1 solar masses, above every brown dwarf"),
)

# --- The reference scale -----------------------------------------------------------------------------------------------

def _mass_fraction_at_proto(symbol):
    """Mass fraction of `symbol` in the proto-solar mixture: A21's log eps with the metal offset, against hydrogen."""
    number_per_hydrogen = 10.0 ** (tuning.STELLAR_SOLAR_LOG_EPS[symbol] + tuning.STELLAR_PROTO_METAL_OFFSET_DEX - 12.0)
    return (tuning.STELLAR_PROTO_X * number_per_hydrogen * elements.element(symbol).atomic_weight
            / elements.element("H").atomic_weight)


METAL_SYMBOLS = ("C", "N", "O", "Ne", "Na", "Mg", "Al", "Si", "S", "Ar", "K", "Ca", "Fe", "Ni")
"""tuple: The 14 metals of the scale; together 99.7 percent of the proto-solar Z."""

PROTO_MASS_FRACTION = {symbol: _mass_fraction_at_proto(symbol) for symbol in METAL_SYMBOLS}
"""dict: Element to its proto-solar mass fraction (A21 Table 2 with the diffusion offset)."""

PROTO_X = tuning.STELLAR_PROTO_X
PROTO_Z = tuning.STELLAR_PROTO_Z
PRIMORDIAL_Y = tuning.STELLAR_PRIMORDIAL_Y
PROTO_Y = 1.0 - PROTO_X - PROTO_Z
"""float: The proto-solar helium mass fraction (0.2725)."""

HELIUM_PER_METAL = (1.0 - PRIMORDIAL_Y - PROTO_X) / PROTO_Z - 1.0
"""float: g in Y = Y_P + g Z, derived (1.674) so that the scale reproduces the proto-solar X, Y and Z exactly. Casagrande et al.
2007 observe 2.1 +- 0.9; the rounded 1.68 of the research note would reproduce X to 1e-4 only."""

UNLISTED_METAL_FRACTION = PROTO_Z - sum(PROTO_MASS_FRACTION.values())
"""float: The proto-solar Z the 14 metals leave out (4.9e-5), scaled with iron (partly verified, flagged)."""


# --- The branch --------------------------------------------------------------------------------------------------------

def high_alpha_probability(fe_h):
    """P(high-alpha) for a star old enough to be thick-disc, a logistic in [Fe/H] (Bensby, Feltzing and Oey 2014): 0.5 at
    -0.35, 0.76 at -0.5, 0.97 at -0.8 and 0.9999 at -1.6."""
    f0, width = tuning.STELLAR_HIGH_ALPHA_LOGISTIC
    return 1.0 / (1.0 + math.exp((fe_h - f0) / width))


def branch_probability(fe_h, age_gy):
    """P(high-alpha) of a star of this [Fe/H] and age: 0 below `tuning.STELLAR_THICK_START_GY` and the logistic above."""
    return 0.0 if age_gy < tuning.STELLAR_THICK_START_GY else high_alpha_probability(fe_h)


def _high_alpha_weight(subject):
    return high_alpha_probability(subject["fe_h"])


def _thin_weight(subject):
    if subject["age_gy"] < tuning.STELLAR_THICK_START_GY:
        return 1.0
    return 1.0 - high_alpha_probability(subject["fe_h"])


BRANCH_TABLE = pt.REGISTRY.register(pt.Table(
    "stellar-abundance-branch", 1,
    (pt.Row("high-alpha sequence", (pt.ANY_CLASS,), _high_alpha_weight,
            "Bensby, Feltzing and Oey 2014 (logistic in [Fe/H]); Haywood et al. 2013 (age gate)",
            state=pt.PARTLY_VERIFIED, limit="the Bensby fit is part A 2.6's; the age gate is uncertain by 1 to 1.5 Gy",
            conditions=(pt.Condition("age_gy", minimum=tuning.STELLAR_THICK_START_GY),),
            notes="alpha-rich: thick disc and halo stars"),
     pt.Row("thin sequence", (pt.ANY_CLASS,), _thin_weight,
            "Bensby, Feltzing and Oey 2014 (logistic in [Fe/H]); Haywood et al. 2013 (age gate)",
            state=pt.PARTLY_VERIFIED, limit="the Bensby fit is part A 2.6's; the age gate is uncertain by 1 to 1.5 Gy",
            notes="every star younger than the age gate, and the old ones the logistic leaves")),
    description="Which abundance sequence a star belongs to (stellar-composition.md section 3 step 2). The high-alpha row comes "
                "first so that one fixed uniform number u decides it: high-alpha when u is below that row's probability."))
"""Table: The branch as a row of the gated possibility loader. Its subject carries `fe_h` and `age_gy`."""


def branch_from_uniform(fe_h, age_gy, u):
    """The branch of a star given the fixed uniform number `u` (0 to 1) of its key: high-alpha when `u` is below
    `branch_probability`, thin otherwise. The same result as the table's probabilities with the high-alpha row first."""
    return BRANCH_HIGH_ALPHA if u < branch_probability(fe_h, age_gy) else BRANCH_THIN


@functools.lru_cache(maxsize=1 << 16)
def _branch_uniform(fe_h):
    return draw.Stream(f"branch:{fe_h!r}").random()


def branch_of(fe_h, age_gy):
    """
    The branch of a star of this birth [Fe/H] and age (GEN.1024): high-alpha when the star's own uniform number, a stream named for its
    [Fe/H] (a float no other star shares), is under `branch_probability`. The star's track (`alpha_tracks.track_fe_h`), its composition
    and its radiogenic offsets all read it, from the two numbers a stored star keeps, so they agree and a star read back from the database
    has the branch it was made with.
    """
    return branch_from_uniform(fe_h, age_gy, _branch_uniform(fe_h))


def branch_probabilities(fe_h, age_gy):
    """`{branch: probability}` from the table, for a test or a report."""
    found = BRANCH_TABLE.probabilities(pt.Subject("star", {"fe_h": fe_h, "age_gy": age_gy}))
    result = {BRANCH_THIN: 0.0, BRANCH_HIGH_ALPHA: 0.0}
    for row, probability in found:
        result[BRANCH_HIGH_ALPHA if row.name == "high-alpha sequence" else BRANCH_THIN] = probability
    return result


# --- The offsets -------------------------------------------------------------------------------------------------------

def _fit_value(fit, fe_h):
    kind = fit[0]
    if kind == "linear":
        return fit[1] + fit[2] * fe_h
    if kind == "min":
        return fit[1] + fit[2] * min(fe_h - fit[3], 0.0)
    if kind == "max":
        return fit[1] + fit[2] * max(fe_h - fit[3], 0.0)
    raise StellarAbundanceError(f"unknown fit form {kind!r}")


def mean_xfe(symbol, fe_h, branch):
    """The branch's mean [X/Fe] (dex) of C, O, Mg, Si, Ca, Al or Eu at [Fe/H] `fe_h`, with [Fe/H] held at the edge of the range
    the fit was made over (thin -0.9 to +0.5, high-alpha -1.6 to +0.1).

    Raises:
        StellarAbundanceError: The branch or the element has no fit.
    """
    try:
        fit = tuning.STELLAR_XFE_FITS[branch][symbol]
    except KeyError:
        raise StellarAbundanceError(f"no [X/Fe] fit for {symbol!r} on the {branch!r} sequence") from None
    low, high = tuning.STELLAR_FE_H_FIT_RANGE[branch]
    return _fit_value(fit, min(max(fe_h, low), high))


def _cholesky(matrix):
    size = len(matrix)
    factor = [[0.0] * size for _ in range(size)]
    for i in range(size):
        for j in range(i + 1):
            partial = sum(factor[i][k] * factor[j][k] for k in range(j))
            factor[i][j] = math.sqrt(matrix[i][i] - partial) if i == j else (matrix[i][j] - partial) / factor[j][j]
    return factor


_CHOLESKY = _cholesky(tuning.STELLAR_XFE_CORRELATION)


def residual_offsets(branch, normals):
    """The correlated residuals (dex) of C, O, Mg, Si and Ca for five independent standard normals `normals`: the normals through
    the Cholesky factor of the residual correlation, each times the branch's scatter of its element."""
    if len(normals) != len(RESIDUAL_ELEMENTS):
        raise StellarAbundanceError("the residuals take five standard normals (C, O, Mg, Si, Ca)")
    scatter = tuning.STELLAR_XFE_SIGMA[branch]
    return {symbol: sum(_CHOLESKY[i][k] * normals[k] for k in range(i + 1)) * scatter[symbol]
            for i, symbol in enumerate(RESIDUAL_ELEMENTS)}


def carbon_enhanced(fe_h, u_enhanced, u_carbon):
    """[C/Fe] of a carbon-enhanced metal-poor star, or `None` when the star is not one: below [Fe/H] -2 a star is one when
    `u_enhanced` is under 0.25 (and always below -4, which the stored [Fe/H] floor never reaches), and its [C/Fe] is 0.7 plus
    `u_carbon` times 1.3 (Frebel and Norris 2015)."""
    cemp = tuning.STELLAR_CEMP
    if fe_h >= cemp["fe_h_below"]:
        return None
    probability = 1.0 if fe_h < cemp["rich_probability_below"] else cemp["probability"]
    if u_enhanced >= probability:
        return None
    return cemp["carbon_fe_min"] + cemp["carbon_fe_span"] * u_carbon


def abundance_offsets(fe_h, branch, normals, u_enhanced=1.0, u_carbon=0.0, neon_follows=None):
    """
    The [X/Fe] (dex) of every element of the scale for one star: C, O, Mg, Si and Ca with their residuals; S and Ar
    following Ca (neither is in the Brewer and Fischer table; both are alpha elements); Ne following O (or Mg, `neon_follows`; the
    default is `tuning.STELLAR_NEON_FOLLOWS`); N, Na, Al, K, Fe and Ni at 0. `normals` are five standard normals (C, O, Mg, Si, Ca);
    `u_enhanced` and `u_carbon` the two uniform numbers of `carbon_enhanced` (the default never makes a star carbon-enhanced). A
    star with a carbon excess keeps the residuals of the other elements.
    """
    follows = tuning.STELLAR_NEON_FOLLOWS if neon_follows is None else neon_follows
    if follows not in ("O", "Mg"):
        raise StellarAbundanceError(f"neon follows O or Mg, not {follows!r}")
    residuals = residual_offsets(branch, normals)
    offsets = {symbol: 0.0 for symbol in METAL_SYMBOLS}
    for symbol in RESIDUAL_ELEMENTS:
        offsets[symbol] = mean_xfe(symbol, fe_h, branch) + residuals[symbol]
    enhanced = carbon_enhanced(fe_h, u_enhanced, u_carbon)
    if enhanced is not None:
        offsets["C"] = enhanced
    offsets["S"] = offsets["Ar"] = offsets["Ca"]
    offsets["Ne"] = offsets[follows]
    return offsets


# --- X, Y, Z and the elements ------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class StellarComposition:
    """
    A star's composition by mass.

    Attributes:
        fe_h (float): [Fe/H], dex.
        x (float): Hydrogen mass fraction.
        y (float): Helium mass fraction.
        z (float): Metal mass fraction.
        elements (dict): The mass fraction of H, He and each metal of `METAL_SYMBOLS`; the metals the scale leaves out are in
            `"other metals"`. Together they add up to 1.
        offsets (dict): The [X/Fe] used.
    """

    fe_h: float
    x: float
    y: float
    z: float
    elements: dict
    offsets: dict

    @property
    def z_over_solar(self):
        """Z against the solar scale's own Z at the same [Fe/H] (10^[Fe/H] times the proto-solar Z)."""
        return self.z / (PROTO_Z * 10.0 ** self.fe_h)


def composition(fe_h, offsets):
    """
    X, Y, Z and the element mass fractions of a star of [Fe/H] `fe_h` and [X/Fe] `offsets` (`abundance_offsets`):
    r = Z/X = (sum of the metals' proto-solar fractions times 10^(F + [X/Fe]) plus the unlisted rest times 10^F) / X0;
    X = (1 - Y_P) / (1 + (1 + g) r); Z = r X; Y = 1 - X - Z; each metal's fraction is X times its proto-solar fraction over X0
    times 10^(F + [X/Fe]).
    """
    scale = 10.0 ** fe_h
    metals = {symbol: PROTO_MASS_FRACTION[symbol] * scale * 10.0 ** offsets.get(symbol, 0.0) for symbol in METAL_SYMBOLS}
    other = UNLISTED_METAL_FRACTION * scale
    ratio = (sum(metals.values()) + other) / PROTO_X
    x = (1.0 - PRIMORDIAL_Y) / (1.0 + (1.0 + HELIUM_PER_METAL) * ratio)
    z = ratio * x
    y = 1.0 - x - z
    factor = x / PROTO_X
    fractions = {"H": x, "He": y}
    fractions.update({symbol: value * factor for symbol, value in metals.items()})
    fractions["other metals"] = other * factor
    return StellarComposition(fe_h, x, y, z, fractions, dict(offsets))


_solar_metal_fraction = composition(0.0, {}).z
"""float: Z0, the metal mass fraction of a star of [Fe/H] 0 with every [X/Fe] 0 (the proto-solar mix of `PROTO_MASS_FRACTION`)."""


def metal_ratio_to_solar(metal_mass_fraction):
    """Z / Z0: how much more (or less) metal by mass a star of metal mass fraction `metal_mass_fraction` has than the Sun-like reference."""
    return metal_mass_fraction / _solar_metal_fraction


def number_ratio(composition_or_offsets, a, b, fe_h=None):
    """The number ratio of element `a` to element `b` (`"Fe"`, `"Si"`): the A21 log eps difference plus the two [X/Fe]. Takes a
    `StellarComposition` or an offsets dict (the ratio is independent of [Fe/H], which cancels)."""
    offsets = composition_or_offsets.offsets if isinstance(composition_or_offsets, StellarComposition) \
        else composition_or_offsets
    return 10.0 ** (tuning.STELLAR_SOLAR_LOG_EPS[a] - tuning.STELLAR_SOLAR_LOG_EPS[b] + offsets.get(a, 0.0) - offsets.get(b, 0.0))


ALPHA_ELEMENTS = ("O", "Mg", "Si", "Ca")
"""tuple: The elements whose mean [X/Fe] is a star's alpha enhancement (S and Ar follow Ca, so they add nothing)."""


def alpha_enhancement(composition_or_offsets):
    """[alpha/Fe] (dex): the mean [X/Fe] of `ALPHA_ELEMENTS`, for a `StellarComposition` or an offsets dict. A thin-disc star at [Fe/H] 0 is near 0, a
    high-alpha one near +0.2 to +0.3 (the branches' fits, `tuning.STELLAR_XFE_FITS`)."""
    offsets = composition_or_offsets.offsets if isinstance(composition_or_offsets, StellarComposition) else composition_or_offsets
    return sum(offsets[symbol] for symbol in ALPHA_ELEMENTS) / len(ALPHA_ELEMENTS)


def mass_ratio(composition_or_offsets, a, b):
    """The mass ratio of element `a` to element `b`: the number ratio times their atomic weights."""
    return number_ratio(composition_or_offsets, a, b) * elements.element(a).atomic_weight / elements.element(b).atomic_weight


# --- Changes of the surface --------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class DredgeUp:
    """
    The first dredge-up at the base of the red-giant branch.

    Attributes:
        delta_y (float): The rise in the helium mass fraction (taken from hydrogen).
        carbon_factor (float): The factor on the carbon-to-metal ratio.
        nitrogen_factor (float): The factor on the nitrogen-to-metal ratio (the Z = 0.02 one at every Z).
        weight (float): Where the star's Z lies between the two tabulated metallicities, in log Z (0 at Z = 0.001, 1 at 0.02).
        mass_held (bool): The mass was outside 1 to 2.5 solar masses and the table's edge was used.
        z_held (bool): The Z was outside 0.001 to 0.02 and the table's edge was used.
    """

    delta_y: float
    carbon_factor: float
    nitrogen_factor: float
    weight: float
    mass_held: bool
    z_held: bool


def _interpolate(xs, ys, x):
    index = min(max(bisect.bisect_right(xs, x) - 1, 0), len(xs) - 2)
    span = xs[index + 1] - xs[index]
    return ys[index] + (ys[index + 1] - ys[index]) * (x - xs[index]) / span


def first_dredge_up(mass_sol, z):
    """
    The surface change of the first dredge-up for a star of `mass_sol` solar masses and metal fraction `z` (Boothroyd and
    Sackmann 1999): the rise in helium, the factors on C/Z and N/Z, linear in mass between the table's rows (1, 1.5, 2, 2.5)
    and in log Z between its two metallicities (0.001 and 0.02). Outside either range the table's edge is held (`mass_held`,
    `z_held`). Nitrogen at low Z is not tabulated, so its factor is the Z = 0.02 one at every Z (partly verified).
    """
    masses = tuning.STELLAR_DREDGE_UP_MASSES_SOL
    high, low = tuning.STELLAR_DREDGE_UP_Z_HIGH, tuning.STELLAR_DREDGE_UP_Z_LOW
    mass = min(max(mass_sol, masses[0]), masses[-1])
    # The low-Z table has no 2.5 solar mass row: use the rows it has.
    low_masses = [m for m, value in zip(masses, low["final_y"]) if value is not None]
    low_y = [value - low["initial_y"] for value in low["final_y"] if value is not None]
    low_c = [value / low["initial_c_over_z"] for value in low["final_c_over_z"] if value is not None]
    high_y = [value - high["initial_y"] for value in high["final_y"]]
    high_c = [value / high["initial_c_over_z"] for value in high["final_c_over_z"]]
    high_n = [value / high["initial_n_over_z"] for value in high["final_n_over_z"]]
    mass_low = min(max(mass_sol, low_masses[0]), low_masses[-1])
    log_low, log_high = math.log10(low["z"]), math.log10(high["z"])
    weight_raw = (math.log10(z) - log_low) / (log_high - log_low)
    weight = min(max(weight_raw, 0.0), 1.0)
    delta_y = _interpolate_pair(_interpolate(low_masses, low_y, mass_low), _interpolate(masses, high_y, mass), weight)
    carbon = _interpolate_pair(_interpolate(low_masses, low_c, mass_low), _interpolate(masses, high_c, mass), weight)
    nitrogen = _interpolate(masses, high_n, mass)
    return DredgeUp(delta_y, carbon, nitrogen, weight, mass != mass_sol, weight != weight_raw)


def _interpolate_pair(at_low, at_high, weight):
    return at_low + weight * (at_high - at_low)


def after_first_dredge_up(star, mass_sol):
    """The `StellarComposition` of a red giant of `mass_sol`: the birth composition `star` with the helium raised from the hydrogen
    and the carbon and nitrogen scaled by the dredge-up's factors. The tables of C/Z and N/Z do not conserve the carbon and
    nitrogen together (the net change is 0.2 percent of Z), so the difference is taken from the hydrogen and Z is the sum of the
    metals; every fraction still adds up to 1."""
    dredge = first_dredge_up(mass_sol, star.z)
    fractions = dict(star.elements)
    carbon = fractions["C"] * dredge.carbon_factor
    nitrogen = fractions["N"] * dredge.nitrogen_factor
    metals_change = (carbon - fractions["C"]) + (nitrogen - fractions["N"])
    fractions["C"], fractions["N"] = carbon, nitrogen
    fractions["He"] += dredge.delta_y
    fractions["H"] -= dredge.delta_y + metals_change
    return StellarComposition(star.fe_h, fractions["H"], fractions["He"], star.z + metals_change, fractions, star.offsets)


# --- Radioactive birth abundances --------------------------------------------------------------------------------------

GY_S = 3.15576e16
"""float: Seconds in a gigayear (Julian)."""


def _decay_back(isotope_name, gy):
    return math.exp(elements.isotope(isotope_name).decay_constant * gy * GY_S)


def solar_birth_th_over_eu():
    """The birth atomic Th/Eu of the solar scale: A21's log eps difference with thorium decayed back
    `tuning.STELLAR_SOLAR_SYSTEM_AGE_GY` (0.405)."""
    logs = tuning.STELLAR_SOLAR_LOG_EPS
    return 10.0 ** (logs["Th"] - logs["Eu"]) * _decay_back("Th-232", tuning.STELLAR_SOLAR_SYSTEM_AGE_GY)


def solar_birth_u_over_th():
    """The birth atomic U/Th of the solar scale (0.572): A21's meteoritic log eps with uranium-238 and uranium-235 and thorium
    decayed back (the half-lives of GEN.438, NUBASE2020). The 2003 Palme and O'Neill cross-check sits
    6 percent away (stellar-composition.md step 10)."""
    logs = tuning.STELLAR_SOLAR_LOG_EPS
    age = tuning.STELLAR_SOLAR_SYSTEM_AGE_GY
    u238, u235 = elements.isotope("U-238"), elements.isotope("U-235")
    now = 10.0 ** (logs["U"] - logs["Th"])
    birth_factor = (u238.abundance * _decay_back("U-238", age) + u235.abundance * _decay_back("U-235", age))
    return now * birth_factor / _decay_back("Th-232", age)


def birth_radiogenic(star, branch, eu_normal=0.0, th_eu_normal=0.0):
    """
    The birth abundances of Eu, Th and U of `star` (a `StellarComposition`), for GEN.459: `{"abundances": {symbol: atoms per
    hydrogen atom}, "mass_fractions": {symbol: mass fraction}, "eu_fe", "th_eu"}`. [Eu/Fe] is the branch's mean plus
    `eu_normal` times its scatter; Th is Eu times the solar birth Th/Eu times 10^[Th/Eu] with [Th/Eu] = 0.014 plus
    `th_eu_normal` times 0.045 (Botelho et al. 2019); U is Th times the solar birth U/Th. Th is measured near solar [Fe/H]
    only, so away from it Th is inferred through Eu (partly verified). Potassium is not fitted.
    """
    logs = tuning.STELLAR_SOLAR_LOG_EPS
    eu_fe = mean_xfe("Eu", star.fe_h, branch) + eu_normal * tuning.STELLAR_XFE_SIGMA[branch]["Eu"]
    th_eu = tuning.STELLAR_TH_EU["log_ratio_mean"] + th_eu_normal * tuning.STELLAR_TH_EU["log_ratio_sigma"]
    eu = 10.0 ** (logs["Eu"] + tuning.STELLAR_PROTO_METAL_OFFSET_DEX - 12.0 + star.fe_h + eu_fe)
    th = eu * solar_birth_th_over_eu() * 10.0 ** th_eu
    abundances = {"Eu": eu, "Th": th, "U": th * solar_birth_u_over_th()}
    hydrogen = elements.element("H").atomic_weight
    mass_fractions = {symbol: star.x * number * elements.element(symbol).atomic_weight / hydrogen
                      for symbol, number in abundances.items()}
    return {"abundances": abundances, "mass_fractions": mass_fractions, "eu_fe": eu_fe, "th_eu": th_eu}


# --- White dwarfs ------------------------------------------------------------------------------------------------------

@functools.cache
def _white_dwarf_table():
    path = os.path.join(DATA_DIRECTORY, "white_dwarf_mass_radius.json")
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    for field in ("sources", "retrieved", "title", "rows", "limiting_mass_solar"):
        if not data.get(field):
            raise StellarAbundanceError(f"white_dwarf_mass_radius.json states no {field}")
    rows = data["rows"]
    return tuple(row[1] for row in rows), tuple(row[2] for row in rows)


def white_dwarf_radius_solar(mass_sol):
    """The radius (solar radii) of a zero-temperature white dwarf of `mass_sol` solar masses, interpolated in the integrated table
    (Chandrasekhar equation of state, mu_e = 2): 0.0127 at 0.6, 0.0082 at 1.0, 0.0080 at Sirius B's 1.018. The stand-in a run
    that disregards missing data (GEN.599) takes in place of the downloaded sequences' radius; it has no temperature, so a real
    star is a few percent larger.

    Raises:
        StellarAbundanceError: The mass is outside the table (0.11 to 1.41 solar masses).
    """
    masses, radii = _white_dwarf_table()
    if not masses[0] <= mass_sol <= masses[-1]:
        raise StellarAbundanceError(f"a {mass_sol} solar mass white dwarf is outside the table "
                                    f"({masses[0]:.2f} to {masses[-1]:.2f})")
    return _interpolate(masses, radii, mass_sol)


def white_dwarf_core(mass_sol):
    """The core type of a white dwarf: "helium" below `tuning.STELLAR_WD_HELIUM_CORE_MAX_SOL` (recalled, open), "oxygen-neon"
    above `tuning.STELLAR_WD_ONE_CORE_MIN_SOL` (Althaus et al. 2010), "carbon-oxygen" between. The downloaded sequences are
    carbon-oxygen only (`datasets.white_dwarf_models`), so the other two cores have no layers."""
    if mass_sol < tuning.STELLAR_WD_HELIUM_CORE_MAX_SOL:
        return "helium"
    if mass_sol > tuning.STELLAR_WD_ONE_CORE_MIN_SOL:
        return "oxygen-neon"
    return "carbon-oxygen"


# --- Neutron stars -----------------------------------------------------------------------------------------------------
# The outer crust and the equations of state come from the downloaded tables (`stellar_composition`).

NEUTRON_STAR_OUTER_CRUST = (
    ("Fe-56", 8.02e6), ("Ni-62", 2.71e8), ("Ni-64", 1.33e9), ("Ni-66", 1.50e9), ("Kr-86", 3.09e9), ("Se-84", 1.06e10),
    ("Ge-82", 2.79e10), ("Zn-80", 6.07e10), ("Zn-82", 8.46e10), ("Pd-128", 9.67e10), ("Ru-126", 1.47e11), ("Mo-124", 2.11e11),
    ("Zr-122", 2.89e11), ("Sr-120", 3.97e11), ("Kr-118", 4.27e11))
"""tuple: `(nucleus, density at which it stops being the ground state, g/cm3)` from the surface inward (Chamel and Haensel 2008
Table 2, Ruster et al. data and BSk8). The stand-in a run that disregards missing data (GEN.599) takes in place of the 20
downloaded nuclear mass models; neutron drip follows at about 4e11 g/cm3."""

NEUTRON_DRIP_DENSITY_G_CM3 = 4.0e11

NUCLEAR_SATURATION_DENSITY_G_CM3 = (2.7e14, 2.8e14)


# --- Black holes -------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class BlackHole:
    """
    The whole of a black hole (no composition: the no-hair theorem).

    Attributes:
        schwarzschild_radius_km (float): 2GM/c^2.
        hawking_temperature_k (float): hbar c^3 / (8 pi G M k_B).
        horizon_radius_km (float): The Kerr outer horizon, GM/c^2 (1 + sqrt(1 - a^2)).
        isco_gm_c2 (float): The prograde innermost stable circular orbit, in GM/c^2.
        efficiency (float): The radiative efficiency 1 - sqrt(1 - 2/(3 r_isco)) of a thin disc.
    """

    schwarzschild_radius_km: float
    hawking_temperature_k: float
    horizon_radius_km: float
    isco_gm_c2: float
    efficiency: float


def black_hole(mass_sol, spin=0.0):
    """The quantities of a black hole of `mass_sol` solar masses and dimensionless `spin` (0 to 1; charge is 0)."""
    if not 0.0 <= spin < 1.0:
        raise StellarAbundanceError("the spin is from 0 up to, not including, 1")
    mass = _msun_to_kg(mass_sol)
    gravitational_radius = constants.G * mass / constants.SPEED_OF_LIGHT_M_S ** 2
    z1 = 1.0 + (1.0 - spin * spin) ** (1.0 / 3.0) * ((1.0 + spin) ** (1.0 / 3.0) + (1.0 - spin) ** (1.0 / 3.0))
    z2 = math.sqrt(3.0 * spin * spin + z1 * z1)
    isco = 3.0 + z2 - math.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))
    temperature = (constants.REDUCED_PLANCK * constants.SPEED_OF_LIGHT_M_S ** 3
                   / (8.0 * math.pi * constants.G * mass * constants.BOLTZMANN))
    return BlackHole(2.0 * gravitational_radius / 1e3, temperature,
                     gravitational_radius * (1.0 + math.sqrt(1.0 - spin * spin)) / 1e3, isco,
                     1.0 - math.sqrt(1.0 - 2.0 / (3.0 * isco)))


# --- Brown dwarfs ------------------------------------------------------------------------------------------------------

JUPITER_MASS_SOL = _kg_to_msun(1.89813e27)
"""float: One Jupiter mass in solar masses (IAU nominal values)."""

def deuterium_burning_limit_mj(z, y=None):
    """The mass (Jupiter masses) above which a brown dwarf burns deuterium: 13.55 - 0.9 log10(Z / 0.0189) - 0.2 (Y - 0.25) / 0.01
    (a linear fit of Spiegel, Burrows and Milsom 2011 Table 1, whose solar Z is 0.0189 and not the scale's). The published
    range is 11.0 to 16.3 Jupiter masses; the fit is held inside it. `y` defaults to the helium the scale gives at `z`."""
    law = tuning.STELLAR_BROWN_DWARF
    helium = PRIMORDIAL_Y + HELIUM_PER_METAL * z if y is None else y
    value = (law["deuterium_mj"] + law["deuterium_z_slope_mj_per_dex"] * math.log10(z / law["deuterium_solar_z"])
             + law["deuterium_y_slope_mj_per_0p01"] * (helium - law["deuterium_reference_y"]) / 0.01)
    low, high = law["deuterium_range_mj"]
    return min(max(value, low), high)


def lithium_burning_limit_sol():
    """The mass (solar masses) above which a brown dwarf burns lithium, 0.06 (63 Jupiter masses; Burrows et al. 2001)."""
    return tuning.STELLAR_BROWN_DWARF["lithium_sol"]


def hydrogen_burning_limit_sol(z):
    """The mass (solar masses) above which a star burns hydrogen: 0.074 at solar Z (the upper of Burrows et al. 2001's
    0.070 to 0.074) and 0.092 at Z = 0, linear in Z between and held at 0.074 above the scale's solar Z. The dependence on Z is a
    parametrisation (stand-in), not a source's law."""
    law = tuning.STELLAR_BROWN_DWARF
    solar = law["hydrogen_solar_z_sol"][1]
    share = min(max(z / PROTO_Z, 0.0), 1.0)
    return law["hydrogen_zero_z_sol"] + share * (solar - law["hydrogen_zero_z_sol"])


def hydrogen_burning_limit_for_fe_h(fe_h):
    """`hydrogen_burning_limit_sol` of a place's birth [Fe/H] `fe_h` (its Z from `composition`, every [X/Fe] 0), the one edge between a
    star and a brown dwarf that the primaries, the companions and the star's makeup all read (GEN.936)."""
    return hydrogen_burning_limit_sol(composition(fe_h, {}).z)


def brown_dwarf_radius_solar(mass_sol, age_gy=1.0):
    """The radius (solar radii) of a brown dwarf or low-mass star at an age, solar metallicity, from the Baraffe et al. 2015 table, continued
    past the age it loses the mass (`baraffe_isochrones.cooling_track`; GEN.948: it was nine points at 1 Gyr).

    Raises:
        StellarAbundanceError: The mass is outside the table (0.01 to 1.4 solar masses).
    """
    try:
        return baraffe_isochrones.cooling_track(mass_sol, age_gy)[1]
    except ValueError as error:
        raise StellarAbundanceError(str(error)) from error
