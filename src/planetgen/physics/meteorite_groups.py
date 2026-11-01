# planetgen/physics/meteorite_groups.py

"""
The chondrite groups as regions of composition space, with a nearest-group label (GEN.475, part 1; docs/design/composition/small-bodies-composition.md section 15).

`small_body_composition.TAXON_DATA` takes each taxon's mix from a meteorite analogue (CI/CM, ordinary, CV and so on): the groups are inputs there. The item
(Boss's stand-in rule of 1:43 AM EDT 2026-10-17) turns them round: a body's composition comes from its own bulk, and the group is a label on the result, a region of a
space whose axes are what sets a chondrite's kind:

* the **oxidation state**, the olivine's mol percent fayalite (`fayalite_mol_percent`), from the body's FeO and MgO;
* the **metal** volume percent, from its metal fraction;
* the **refractory** inclusions' volume percent (CAIs and amoeboid olivine aggregates), from its refractory (Ca, Al) enrichment;
* the **matrix**, the volume percent of fine grains, which carries the water and organics: high in the aqueously altered CI and CM and in the ice-bearing bodies.

`CHONDRITE_GROUPS` is Table 1 of [R4658] (Meteorites and the Early Solar System II), the averages of the chondrite groups (the table is held in
`/mnt/project-files/design-drafts/gen465-second-pass/weisberg2006_table1_chondrite_groups.txt`); the CI to CK carbonaceous groups, the enstatite chondrites (EH, EL) and the
ordinary chondrites (H, L, LL) are the same table's columns. `nearest_group` returns the group with the smallest scaled distance and the distance, which the caller can compare with `NEAREST_LIMIT`:
a body far from every group (a metal-rich or ice-rich one, a differentiated one) is labelled `none`, not forced into one.

Limits. The axes are four of the table's rows; the oxygen isotopes, the chondrule sizes and the petrologic types are not in the space. The groups' matrix and CAI abundances are
volume percent of a rock, and a derived body's matrix is a stand-in for its volatile content (the wiring step maps the bulk's water and organics). A group with only a
range in the table (CO and CV olivine, <1 to 47 and <1 to 36 mol percent Fa) has no coordinate on that axis and is placed by the others. Nothing in generation
calls this yet: the wiring step labels belts and fields by it from the derived bulk and keeps the groups as calibration (`TEST.215`).
"""

import math

from planetgen.util.checks import finite_domain

# name: (fayalite mol percent or None, metal vol percent, refractory (CAI and AOA) vol percent, matrix vol percent): [R4658: Table 1]
CHONDRITE_GROUPS = {
    "CI": (None, 0.0, 0.0, 99.5),
    "CM": (None, 0.1, 5.0, 70.0),
    "CO": (None, 3.0, 13.0, 34.0),
    "CV": (None, 2.5, 10.0, 40.0),
    "CK": (31.0, 2.5, 10.0, 40.0),
    "CR": (2.0, 6.5, 0.5, 40.0),
    "CH": (2.0, 20.0, 0.1, 5.0),
    "CB": (3.0, 70.0, 0.0, 0.0),
    "H": (19.3, 8.0, 0.0, 12.5),
    "L": (25.2, 4.0, 0.0, 12.5),
    "LL": (31.3, 2.0, 0.0, 12.5),
    "EH": (0.4, 10.0, 0.0, 8.5),
    "EL": (0.4, 10.0, 0.0, 8.5),
    "R": (38.0, 0.1, 0.0, 36.0),
    "K": (2.2, 7.0, 0.0, 73.0),
}
"""dict: Group -> (fayalite mol percent, metal, refractory, matrix volume percent), the averages of the table (ranges given in it are taken at their midpoint; the carbonaceous
groups' matrix and CAI rows as printed, `<<1` as 0; CO and CV have no olivine average). [Stand-in MG1, tier T1, tweak L0.]"""

AXIS_SCALES = (10.0, 15.0, 5.0, 30.0)
"""tuple: The width of each axis (Fa, metal, refractory, matrix) that distances are scaled by: about a quarter of the span the groups cover on it. [Stand-in MG2, tier T9, tweak L4.]"""

NEAREST_LIMIT = 1.5
"""float: A body whose scaled distance from the nearest group is over this is labelled `none`. [Stand-in MG2, tier T9, tweak L4.]"""

NONE = "none"


@finite_domain()
def fayalite_mol_percent(feo_mass_fraction, mgo_mass_fraction):
    """The fayalite content of an olivine that holds the mantle's iron and magnesium: `100 Fe / (Fe + Mg)` by moles from the oxides (FeO 71.84 and MgO 40.30 g/mol)."""
    if feo_mass_fraction < 0.0 or mgo_mass_fraction < 0.0 or feo_mass_fraction + mgo_mass_fraction <= 0.0:
        raise ValueError("the oxides must be non-negative and not both zero")
    iron = feo_mass_fraction / 71.844
    magnesium = mgo_mass_fraction / 40.304
    return 100.0 * iron / (iron + magnesium)


def _distance(coordinates, group):
    """The root-mean-square scaled difference over the axes both have, and how many that is."""
    total, used = 0.0, 0
    for index, (value, centre) in enumerate(zip(coordinates, CHONDRITE_GROUPS[group])):
        if value is None or centre is None:
            continue
        total += ((value - centre) / AXIS_SCALES[index]) ** 2
        used += 1
    return (math.sqrt(total / used) if used else math.inf), used


def nearest_group(fayalite_percent, metal_vol_percent, refractory_vol_percent, matrix_vol_percent):
    """
    The chondrite group nearest a body in the space (scaled distance over the axes both have) and the distance, or `("none", d)` past `NEAREST_LIMIT`.

    Args:
        fayalite_percent (float or None): The olivine's mol percent fayalite, None when the body has no olivine to speak of.
        metal_vol_percent, refractory_vol_percent, matrix_vol_percent (float): The body's volume percent of metal, of refractory inclusions and of matrix.
    """
    coordinates = (fayalite_percent, metal_vol_percent, refractory_vol_percent, matrix_vol_percent)
    # equal distances go to the group compared on more axes (CK, with an olivine average, over CV, which has none)
    best = min(CHONDRITE_GROUPS, key=lambda name: (_distance(coordinates, name)[0], -_distance(coordinates, name)[1]))
    distance = _distance(coordinates, best)[0]
    return (best if distance <= NEAREST_LIMIT else NONE), distance
