# planetgen/physics/accessory_minerals.py

"""
The accessory minerals of the asteroid taxon mixes (GEN.1368)
=============================================================

`small_body_composition.TAXON_DATA` writes each taxon's mix as about fifteen major minerals of its meteorite analogue (olivine, pyroxene, metal, troilite, serpentine and the like).
The meteorites carry more: chromite and phosphates in the ordinary chondrites, the calcium-aluminium-rich inclusion phases (melilite, spinel, hibonite, perovskite) in the CO and CV
chondrites, the Fe-Ni sulfides and the carbonates, sulfates and phyllosilicates of the altered carbonaceous ones, the oldhamite family of the enstatite chondrites and
achondrites, nanodiamond and the presolar silicon carbide. This module gives each taxon its accessory minerals and scales them by the abundance of the element that sets them in the
belt's host star.

**The reference mixes** (`ACCESSORIES`): the mass share of each accessory mineral in the taxon's grains at the reference abundance (the Sun's, as a scale, never as content), a modal
list recalled from the meteorite groups' reviews (the research file `chondrite_groups_modal_and_minerals.csv` and `differentiated_meteorite_minerals.csv`; stand-ins ACC1 and ACC2).

**The weights from the bulk** (`scale`): a mineral that one element sets (chromite and daubreelite by chromium, ilmenite and perovskite by titanium, the phosphates by phosphorus, spinel, hibonite
and melilite by aluminium, pentlandite by nickel, alabandite by manganese, nepheline by sodium) is scaled by that element's atoms per magnesium atom in the host star over the same ratio at the
reference abundance (`tuning.STELLAR_SOLAR_LOG_EPS`, used as a scale). A belt with no host (a field, a hand-made belt) takes the reference. The minerals no tracked element sets (the
carbonates, sulfates and phyllosilicates, which the water and the heating set; the organics and nanodiamond; the oldhamite family, set by sulfur and the redox state) are not scaled (ACC3).

**The mix** (`mix`): the accessories are added to the taxon's base mix and the whole is renormalised, so the base minerals give way in proportion to what the accessories take. No draw
is made.

Limits: the lists are modal averages of a group, not a section of one body; the scaling reads one element per mineral, not a speciation of the whole bulk (a body poor in phosphorus has fewer
phosphates, but the aluminium that the spinel takes does not come out of the plagioclase); GEN.1369 reconciles the 84-name vocabulary of `tuning.ASTEROID_COMPONENTS` with these rows.
"""

import math

from planetgen import tuning

ACCESSORIES = {
    "C": {"cronstedtite": 0.100, "saponite": 0.040, "dolomite": 0.020, "magnesite": 0.010, "pentlandite": 0.012, "pyrrhotite": 0.012, "tochilinite": 0.020, "epsomite": 0.015,
          "gypsum": 0.008, "spinel": 0.003, "hibonite": 0.002, "carbon_diamond": 0.0005, "moissanite": 1.0e-5},
    "B": {"cronstedtite": 0.080, "saponite": 0.040, "dolomite": 0.025, "magnesite": 0.012, "pentlandite": 0.012, "pyrrhotite": 0.010, "tochilinite": 0.015, "epsomite": 0.020,
          "gypsum": 0.010, "carbon_diamond": 0.0005, "moissanite": 1.0e-5},
    "S": {"chromite": 0.006, "merrillite": 0.004, "chlorapatite": 0.002, "pentlandite": 0.002, "carbon_diamond": 5.0e-5, "moissanite": 1.0e-5},
    "V": {"ilmenite": 0.010, "chromite": 0.005, "merrillite": 0.003, "chlorapatite": 0.001},
    "A": {"chromite": 0.010, "merrillite": 0.004, "chlorapatite": 0.002, "pentlandite": 0.002},
    "E": {"oldhamite": 0.010, "daubreelite": 0.008, "niningerite": 0.015, "sinoite": 0.002, "carbon_graphite": 0.003},
    "M": {"chromite": 0.004, "daubreelite": 0.005, "carbon_graphite": 0.010, "cohenite": 0.010},
    "K": {"melilite": 0.040, "spinel": 0.020, "hibonite": 0.010, "perovskite": 0.004, "chromite": 0.004, "pentlandite": 0.010, "nepheline": 0.010, "carbon_diamond": 3.0e-4,
          "moissanite": 1.0e-5},
    "L": {"melilite": 0.050, "hibonite": 0.010, "perovskite": 0.006, "chromite": 0.004, "pentlandite": 0.010, "nepheline": 0.010, "carbon_diamond": 3.0e-4, "moissanite": 1.0e-5},
    "P": {"pyrrhotite": 0.020, "pentlandite": 0.005, "carbon_diamond": 0.002, "moissanite": 1.0e-5, "dolomite": 0.010},
    "D": {"pyrrhotite": 0.020, "pentlandite": 0.005, "carbon_diamond": 0.003, "moissanite": 1.0e-5, "dolomite": 0.010},
}
"""dict: Taxon letter -> {materials row: mass share of the taxon's grains at the reference abundance}. A recalled modal list, written as mass shares (stand-ins ACC1, ACC2): the research
files give ranges, the middle taken."""

CONTROL = {
    "chromite": "Cr", "daubreelite": "Cr", "ilmenite": "Ti", "perovskite": "Ti", "merrillite": "P", "chlorapatite": "P", "whitlockite": "P", "hydroxylapatite": "P",
    "spinel": "Al", "hibonite": "Al", "melilite": "Al", "pentlandite": "Ni", "alabandite": "Mn", "nepheline": "Na",
}
"""dict: Row -> the element whose abundance sets it. A row not here is not scaled."""

HOST_FIELD = {"Cr": "cr", "Ti": "ti", "P": "p", "Al": "al", "Ni": "ni", "Mn": "mn", "Na": "na"}
"""dict: Element -> the `bulk_draw.HostRatios` attribute that holds its atoms per magnesium atom."""


def reference_ratio(element):
    """The element's atoms per magnesium atom at the reference abundance (`tuning.STELLAR_SOLAR_LOG_EPS`), the scale the shares are written at."""
    logs = tuning.STELLAR_SOLAR_LOG_EPS
    return 10.0 ** (logs[element] - logs["Mg"])


def scale(element, host):
    """How much of the element's minerals the host star's bulk gives against the reference: the star's atoms per magnesium atom over the reference ratio; 1 with no host."""
    if host is None:
        return 1.0
    ratio = getattr(host, HOST_FIELD[element])
    return ratio / reference_ratio(element)


def accessories(taxon, host=None):
    """The taxon's accessory minerals as mass shares of its grains, scaled by the host star's bulk (`scale`)."""
    return {row: share * (scale(CONTROL[row], host) if row in CONTROL else 1.0) for row, share in ACCESSORIES.get(taxon, {}).items()}


def mix(taxon, base, host=None):
    """
    The taxon's mix with its accessories: `base` (row -> mass fraction) and the scaled accessories, renormalised to 1. A row in both adds.

    Args:
        taxon (str): The taxon letter.
        base (dict): The base mix, row -> mass fraction.
        host (bulk_draw.HostRatios or None): The belt's host star's ratios; None is the reference.
    """
    total = math.fsum(base.values())
    parts = {row: value / total for row, value in base.items()}
    for row, share in accessories(taxon, host).items():
        parts[row] = parts.get(row, 0.0) + share
    whole = math.fsum(parts.values())
    return {row: value / whole for row, value in parts.items()}
