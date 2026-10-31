# planetgen/physics/rock_types.py

"""
Rock types from a bulk composition: the normative minerals of its oxides, what it melts to, and the ground type of the result (GEN.474, part 1;
docs/design/composition/surface-composition.md, "Rock types from the oxides").

`surface_composition.GROUND_TYPES` is a library of the Earth's and the Moon's crusts (Hu, Ehlmann and Seager 2012 Table 2): each ground type is a
fixed mineral list. Here the ground type is worked out from the oxides of the body's own mantle:

1. **The norm** (`cipw_norm`). The CIPW norm (Cross, Iddings, Pirsson and Washington 1902) allocates the oxides of a rock to an ideal set of minerals by mass
   balance, in the fixed order of the method: apatite and ilmenite, then orthoclase, albite and anorthite from the alumina, magnetite, then diopside and
   hypersthene from the lime, magnesia and iron; the silica left over is quartz, and a shortage converts hypersthene to olivine, albite to nepheline and
   orthoclase to leucite (and olivine to periclase). Every oxide goes to a mineral: the norm adds up to the oxides given (tested), so a composition that cannot
   be allocated is refused and not trimmed.
2. **What melts** (`partial_melt`). A mantle does not become a crust whole: a fraction F melts and the melt leaves. The easiest-melting minerals go first, and
   together (the eutectic of a rock is the mixture of its low-melting minerals: quartz and the feldspars of a granite, the feldspar and pyroxene of a basalt;
   the orthopyroxene follows and olivine is last). `MELT_STAGES` is that grouping of the norm's minerals; within a stage the minerals melt in the proportions the
   rock holds them, and the melt is the oxides of what has melted when F of the rock's mass is melt. The heat producers K, Th and U are
   incompatible: in batch melting their concentration in the melt is `1 / (D + F (1 - D))` times the source's (`batch_enrichment`, D the bulk partition
   coefficient, taken as zero).
3. **The ground type** (`ground_type_of`), from the oxides of the melt or of the rock: anorthositic when the plagioclase is 90 percent or more of the norm
   (IUGS anorthosite), ultramafic when silica is under 45 percent (ultrabasic) or magnesia over 18 percent (komatiite), then by the silica of the total alkali-silica
   diagram's bands (Le Bas et al. 1986): basic (45 to 52 percent, basaltic), intermediate (52 to 63) and acid (63 or more, felsic).

The Earth's oxide analyses are test points and never inputs: the bulk silicate Earth (McDonough and Sun 1995, `mixture_lava.OXIDES`) melted by a tenth is a
basalt; the upper continental crust (Rudnick and Gao 2014 Table 11) is felsic and the lower crust intermediate; a lunar highland anorthosite is anorthositic.

Limits. The norm has no hydrous or carbonate minerals, no sulfur and no peralkaline allocation (acmite and sodium metasilicate): a rock whose alkalis exceed its
alumina is refused. The melt is a staged allocation, not a phase-equilibrium calculation (MELTS or pMELTS would give the real liquid at a pressure; the real
first melts of a peridotite are alkaline and richer in the pyroxene's magnesium than the stage mix gives); the iron split between minerals is proportional
(one Fe/(Fe+Mg) for all the mafic minerals, so the norm's olivine is not as magnesian as a real one). The ground type has no mixed or sedimentary rocks. Wired (GEN.474): the surface draw
takes the volcanic (secondary) crust of a body from its stored mantle through `rock_oxides`, `partial_melt` and `ground_type_of`; the library of `GROUND_TYPES` stays as the
assemblages of the types and as the calibration.
"""

import math

from planetgen.physics import elements
from planetgen.util.checks import finite_domain

OXIDES = ("SiO2", "TiO2", "Al2O3", "Fe2O3", "FeO", "MgO", "CaO", "Na2O", "K2O", "P2O5")
"""tuple: The oxides the norm reads (weight percent or any one mass unit)."""

_MASS = {oxide: elements.molar_mass(oxide) * 1000.0 for oxide in OXIDES}
_EPS = 1.0e-12

# oxide moles per formula unit of each normative mineral (FM is MgO and FeO in the rock's Fe/(Fe+Mg))
MINERALS = {
    "quartz": {"SiO2": 1.0},
    "orthoclase": {"K2O": 1.0, "Al2O3": 1.0, "SiO2": 6.0},
    "albite": {"Na2O": 1.0, "Al2O3": 1.0, "SiO2": 6.0},
    "anorthite": {"CaO": 1.0, "Al2O3": 1.0, "SiO2": 2.0},
    "nepheline": {"Na2O": 1.0, "Al2O3": 1.0, "SiO2": 2.0},
    "leucite": {"K2O": 1.0, "Al2O3": 1.0, "SiO2": 4.0},
    "corundum": {"Al2O3": 1.0},
    "diopside": {"CaO": 1.0, "SiO2": 2.0, "FM": 1.0},
    "wollastonite": {"CaO": 1.0, "SiO2": 1.0},
    "hypersthene": {"SiO2": 1.0, "FM": 1.0},
    "olivine": {"SiO2": 0.5, "FM": 1.0},
    "periclase": {"FM": 1.0},
    "magnetite": {"FeO": 1.0, "Fe2O3": 1.0},
    "hematite": {"Fe2O3": 1.0},
    "ilmenite": {"FeO": 1.0, "TiO2": 1.0},
    "rutile": {"TiO2": 1.0},
    "apatite": {"CaO": 10.0 / 3.0, "P2O5": 1.0},
}
"""dict: The normative minerals as oxide moles per formula unit (the formula unit of olivine and periclase is one FM, of the feldspars one alkali or lime)."""

MELT_STAGES = (("quartz", "orthoclase", "albite", "nepheline", "leucite", "anorthite", "diopside", "wollastonite"), ("hypersthene",), ("olivine", "periclase"))
"""tuple: The minerals of the norm by the stage they melt in: the low-melting eutectic components together (in the proportions the rock holds them: the feldspars and
pyroxene of a basalt, the quartz and feldspars of a granite), then the orthopyroxene, then the refractory olivine; the oxides (magnetite, ilmenite, apatite, corundum)
stay in the residue. [Stand-in RT1, tier T5, tweak L1.]"""

PLAGIOCLASE = ("albite", "anorthite")
FELDSPATHOID_FREE = ("quartz",)

ANORTHOSITE_PLAGIOCLASE = 0.90
"""float: The share of plagioclase in the norm from which a rock is anorthositic (the IUGS modal limit of anorthosite, 90 percent). [Stand-in RT2, tier T2, tweak L0.]"""

ULTRABASIC_SILICA = 45.0
KOMATIITE_MAGNESIA = 18.0
BASIC_SILICA = 52.0
ACID_SILICA = 63.0
"""float: The silica bands of the total alkali-silica diagram (Le Bas et al. 1986: ultrabasic under 45, basic to 52, intermediate to 63, acid above) and the
18 percent magnesia of a komatiite. [Stand-in RT2, tier T2, tweak L0.]"""

INCOMPATIBLE_D = 0.0
"""float: The bulk partition coefficient of K, Th and U between a mantle's minerals and its melt, taken as zero (all are far under 0.01 in peridotite). [Stand-in RT3, tier T5, tweak L0.]"""


def _moles(oxides):
    unknown = set(oxides) - set(OXIDES)
    if unknown:
        raise ValueError(f"no norm for the oxides {sorted(unknown)}")
    moles = {oxide: oxides.get(oxide, 0.0) / _MASS[oxide] for oxide in OXIDES}
    if any(value < 0.0 for value in moles.values()) or sum(oxides.values()) <= 0.0:
        raise ValueError("the oxides must be non-negative and not all zero")
    return moles


def _units(oxides):
    """The norm as moles of each mineral's formula unit, and the Fe/(Fe+Mg) of its mafic minerals (the algorithm of the CIPW norm)."""
    m = _moles(oxides)
    units = dict.fromkeys(MINERALS, 0.0)
    cao, alumina, feo, mgo = m["CaO"], m["Al2O3"], m["FeO"], m["MgO"]
    units["apatite"] = m["P2O5"]
    cao -= 10.0 / 3.0 * m["P2O5"]
    if cao < -_EPS:
        raise ValueError("the phosphate needs more lime than the rock has")
    units["ilmenite"] = min(m["TiO2"], feo)
    units["rutile"] = m["TiO2"] - units["ilmenite"]
    feo -= units["ilmenite"]
    units["orthoclase"] = min(m["K2O"], alumina)
    alumina -= units["orthoclase"]
    units["albite"] = min(m["Na2O"], alumina)
    alumina -= units["albite"]
    if m["K2O"] - units["orthoclase"] > _EPS or m["Na2O"] - units["albite"] > _EPS:
        raise ValueError("the alkalis exceed the alumina: no peralkaline norm")
    units["anorthite"] = min(alumina, max(cao, 0.0))
    cao -= units["anorthite"]
    alumina -= units["anorthite"]
    units["corundum"] = alumina
    units["magnetite"] = min(m["Fe2O3"], feo)
    feo -= units["magnetite"]
    units["hematite"] = m["Fe2O3"] - units["magnetite"]
    fm = feo + mgo
    iron = feo / fm if fm > 0.0 else 0.0
    units["diopside"] = min(cao, fm)
    cao -= units["diopside"]
    fm -= units["diopside"]
    units["wollastonite"] = max(cao, 0.0)
    units["hypersthene"] = fm
    silica = m["SiO2"] - sum(units[name] * MINERALS[name].get("SiO2", 0.0) for name in MINERALS if name != "quartz")
    if silica >= 0.0:
        units["quartz"] = silica
        return units, iron
    deficit = -silica
    for source, target, saved in (("hypersthene", "olivine", 0.5), ("albite", "nepheline", 4.0), ("orthoclase", "leucite", 2.0), ("olivine", "periclase", 0.5)):
        converted = min(deficit / saved, units[source])
        units[source] -= converted
        units[target] += converted
        deficit -= converted * saved
        if deficit <= _EPS:
            return units, iron
    raise ValueError("the silica is too short for any norm")


def _mass(name, iron):
    """The mass of one formula unit of a normative mineral, g/mol."""
    total = 0.0
    for oxide, count in MINERALS[name].items():
        total += count * ((1.0 - iron) * _MASS["MgO"] + iron * _MASS["FeO"] if oxide == "FM" else _MASS[oxide])
    return total


@finite_domain()
def cipw_norm(oxides):
    """
    The CIPW norm of a rock: the mass of each normative mineral in the unit the oxides were given in (weight percent adds up to the oxides' sum).

    Args:
        oxides (dict): Oxide -> mass in any one unit (the keys of `OXIDES`; absent ones are zero). Iron as `FeO` and `Fe2O3` (the caller's ferric ratio).

    Raises:
        ValueError: A negative oxide, an unknown one, a peralkaline rock or one with too little silica for any norm.
    """
    units, iron = _units(oxides)
    return {name: units[name] * _mass(name, iron) for name in MINERALS if units[name] > 0.0}




@finite_domain()
def partial_melt(oxides, melt_fraction):
    """
    The oxides (weight percent, adding to 100) of the melt of a rock of which `melt_fraction` of the mass melted, and the residue's, by the stages of `MELT_STAGES`.

    A stage melts whole before the next begins; the stage that straddles the limit melts in part, each of its minerals by the same share of itself.
    """
    if not 0.0 < melt_fraction <= 1.0:
        raise ValueError("the melt fraction must be in (0, 1]")
    units, iron = _units(oxides)
    masses = {name: units[name] * _mass(name, iron) for name in MINERALS}
    total = sum(masses.values())
    want = melt_fraction * total
    taken = dict.fromkeys(MINERALS, 0.0)
    for stage in MELT_STAGES:
        stage_mass = sum(masses[name] for name in stage)
        if stage_mass <= 0.0:
            continue
        share = min(1.0, want / stage_mass)
        for name in stage:
            taken[name] = share * units[name]
        want -= share * stage_mass
        if want <= _EPS * total:
            break
    return _oxides_of(taken, iron), _oxides_of({name: units[name] - taken[name] for name in MINERALS}, iron)


def _oxides_of(units, iron):
    out = dict.fromkeys(OXIDES, 0.0)
    for name, amount in units.items():
        for oxide, count in MINERALS[name].items():
            if oxide == "FM":
                out["MgO"] += amount * count * (1.0 - iron) * _MASS["MgO"]
                out["FeO"] += amount * count * iron * _MASS["FeO"]
            else:
                out[oxide] += amount * count * _MASS[oxide]
    total = sum(out.values())
    return {oxide: 100.0 * value / total for oxide, value in out.items()} if total > 0.0 else out


@finite_domain()
def batch_enrichment(melt_fraction, partition=INCOMPATIBLE_D):
    """The concentration of a trace element in a batch melt over its concentration in the source: `1 / (D + F (1 - D))`."""
    if not 0.0 < melt_fraction <= 1.0:
        raise ValueError("the melt fraction must be in (0, 1]")
    if not 0.0 <= partition <= 1.0:
        raise ValueError("the partition coefficient must be in [0, 1]")
    return 1.0 / (partition + melt_fraction * (1.0 - partition))


@finite_domain()
def rock_oxides(feo_fraction, al2o3_fraction, mg_si_ratio, ca_over_mg):
    """
    The oxides of a rock the generator stores as a few numbers, in weight percent adding to 100: `MgO` and `SiO2` in the molar `mg_si_ratio` (Mg/Si), `CaO` at `ca_over_mg` Ca atoms
    per Mg atom (the star's), and `Al2O3` and `FeO` as the given mass fractions of the whole rock. The alkalis, titanium and phosphorus are not stored and are zero (the norm then
    has anorthite for its feldspar and no orthoclase or albite).
    """
    if not (0.0 <= feo_fraction and 0.0 <= al2o3_fraction and feo_fraction + al2o3_fraction < 1.0 and mg_si_ratio > 0.0 and ca_over_mg >= 0.0):
        raise ValueError("the oxide fractions must be non-negative and leave room for the rest of the rock, with Mg/Si positive")
    mgo, sio2, cao = _MASS["MgO"], _MASS["SiO2"] / mg_si_ratio, ca_over_mg * _MASS["CaO"]
    rock = (mgo + sio2 + cao) / (1.0 - feo_fraction - al2o3_fraction)
    parts = {"SiO2": sio2, "MgO": mgo, "CaO": cao, "Al2O3": al2o3_fraction * rock, "FeO": feo_fraction * rock}
    return {oxide: 100.0 * mass / rock for oxide, mass in parts.items()}


@finite_domain()
def ground_type_of(oxides):
    """
    The ground type the oxides make, one of `anorthositic`, `ultramafic`, `basaltic`, `intermediate` or `felsic` (the generator's ground types of the same names).
    """
    norm = cipw_norm(oxides)
    total = sum(oxides.values())
    silica, magnesia = 100.0 * oxides.get("SiO2", 0.0) / total, 100.0 * oxides.get("MgO", 0.0) / total
    if sum(norm.get(name, 0.0) for name in PLAGIOCLASE) >= ANORTHOSITE_PLAGIOCLASE * sum(norm.values()):
        return "anorthositic"
    if silica < ULTRABASIC_SILICA or magnesia > KOMATIITE_MAGNESIA:
        return "ultramafic"
    if silica < BASIC_SILICA:
        return "basaltic"
    return "intermediate" if silica < ACID_SILICA else "felsic"
