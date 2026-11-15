# planetgen/physics/rock_assemblage.py

"""
The mineral assemblage of a rock from its own oxides (GEN.474, part 2; docs/design/composition/surface-composition.md, "Rock types from the oxides")
=======================================================================================================================================================

`rock_types` gives a melt of a body's mantle its ground type; the assemblage of that ground (the minerals the surface rows and the albedo read) is still a fixed recipe per type
(`surface_composition.GROUND_TYPES`, the Earth's and the Moon's crusts of [R2024: Table 2]). This module is what replaces the recipe: the assemblage of a rock is the **CIPW norm of its
oxides** (`rock_types.cipw_norm`), each normative mineral on the row of the materials table that has its name, as a share of the norm's mass.

* Albite, anorthite, orthoclase, nepheline, leucite, quartz, corundum, diopside, wollastonite, olivine, periclase, magnetite, hematite, ilmenite, chromite, rutile and apatite
  are rows of the materials table of the same names.
* The norm's hypersthene is the orthopyroxene `(Mg, Fe)SiO3`: it is split into `enstatite` and `ferrosilite` by the iron of the divalent pool (`rock_types.divalent_pool`).
* The two metasilicates (the alkali the alumina cannot take, GEN.1362) and eskolaite have no row; their share is returned as `unmapped` and the rest is scaled to 1.

Nothing in it is a number measured on a rock: the shares follow from the oxides by mass balance. The limits are those of the norm and of the staged melt (`rock_types`): an ideal
mineral set, not the minerals a rock crystallises, and a melt that is a stage mix and not a phase-equilibrium liquid (the melt of a tenth of a mantle comes out as plagioclase and
diopside with no orthopyroxene). Mica, amphibole and every hydrous mineral are outside the norm (the biotite of the `felsic` and `intermediate` recipes is not reached). Nothing calls this
yet: the surface draw reads the recipes until the wiring step replaces them (it changes the stored ground shares, a lock item).

Use:

    shares, unmapped = rock_assemblage.assemblage(oxides)      # {materials row: share of the norm's mass}, the share the rows could not take
"""

from planetgen.physics import elements, rock_types
from planetgen.util.checks import finite_domain

DIRECT_ROWS = ("quartz", "orthoclase", "albite", "anorthite", "nepheline", "leucite", "corundum", "diopside", "wollastonite", "olivine", "periclase",
               "magnetite", "hematite", "ilmenite", "chromite", "rutile", "apatite")
"""tuple: The normative minerals whose names are rows of the materials table."""

PYROXENE_ROWS = ("enstatite", "ferrosilite")
"""tuple: The rows of the orthopyroxene of the norm (`hypersthene`), magnesian and ferrous ends."""


def magnesian_mass_share(oxides):
    """
    The share of the mass of a norm's orthopyroxene that is the magnesian end-member: the pool's magnesium as `MgSiO3` against all of the pool (iron, manganese and nickel
    as their metasilicates), by the formula mass of each. 1 for a rock with no iron in the pool, and for a rock with no pool.
    """
    pool = rock_types.divalent_pool(oxides)
    silica = elements.molar_mass("SiO2")
    masses = {oxide: share * (elements.molar_mass(oxide) + silica) for oxide, share in pool.items()}
    total = sum(masses.values())
    return masses.get("MgO", 0.0) / total if total > 0.0 else 1.0


@finite_domain()
def assemblage(oxides):
    """
    The mineral assemblage of a rock of these oxides: `({row of the materials table: share of the norm's mass}, unmapped)`, the shares adding to 1 and `unmapped` the share of the norm's
    mass that has no row (the metasilicates and eskolaite), 0 to 1. Raises `ValueError` where the norm refuses the oxides (`rock_types.cipw_norm`).
    """
    norm = rock_types.cipw_norm(oxides)
    total = sum(norm.values())
    if total <= 0.0:
        raise ValueError("a rock with no norm has no assemblage")
    shares = {}
    for name, mass in norm.items():
        if name in DIRECT_ROWS:
            shares[name] = shares.get(name, 0.0) + mass
        elif name == "hypersthene":
            magnesian = magnesian_mass_share(oxides)
            shares["enstatite"] = shares.get("enstatite", 0.0) + mass * magnesian
            shares["ferrosilite"] = shares.get("ferrosilite", 0.0) + mass * (1.0 - magnesian)
    mapped = sum(shares.values())
    if mapped <= 0.0:
        raise ValueError("no mineral of the norm has a row")
    return {name: value / mapped for name, value in shares.items() if value > 0.0}, max(1.0 - mapped / total, 0.0)
