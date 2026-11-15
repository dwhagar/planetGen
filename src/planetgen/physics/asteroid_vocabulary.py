# planetgen/physics/asteroid_vocabulary.py

"""
The 84 asteroid component names and the minerals that host them (GEN.1369)
==========================================================================

`tuning.ASTEROID_COMPONENTS` is the old flat list of 84 names a belt once drew four of (the metals, the platinum-group elements, the rare-earth and uranium minerals, the gems, the polycyclic
aromatic hydrocarbons). A belt no longer draws from it: `small_body_composition.TAXON_DATA` and `accessory_minerals.ACCESSORIES` name the minerals a body is made of. The test that every
name is a row of the materials table (GEN.427) and the reach scan of the table (TEST.685, which reads `tuning.py`) both count a name as read because the list names it, so a name no body can
carry looked read. This module settles each of the 84 one of two ways:

**Hosted** (`HOSTED`): a primitive body carries the name in a mineral its taxon mix draws, the host the ore and meteorite literature gives for the element group (the platinum-group elements
and gold in the iron-nickel metal at 0.01 to 10 ppm, the inclusion metal nuggets of the osmium-iridium alloy in the calcium-aluminium-rich inclusion phases, titanium in ilmenite and perovskite, the
polycyclic aromatic hydrocarbons in the organic matter, the major elements in the silicates, metal and sulfides). `hosted` lists the names a drawn composition carries, from the rows it holds; the
belt stores them in its makeup detail (`small_body_composition._from_taxa`). A name hosted by a row is a statement that the body can carry it, not an abundance: no mass is drawn for it.

**Vocabulary only** (`VOCABULARY_ONLY`): no taxon mix makes the mineral, because it forms in a setting a primitive body lacks (a magmatic sulfide ore, a carbonatite, a granite or pegmatite, a placer,
a supergene weathering zone, an evaporite). The reason is written beside the name; the ore processes of a differentiated body (GEN.1380 to GEN.1382) and the evaporite sequence (GEN.1366) are
where the rows are made. The host table of a primitive body is a recalled stand-in (stand-in ACV1, tier T3, tweak L0).

A key in the mixes must be a row id, not an alias (`water_ice` and `water ice` meet by alias, which is fragile): `tests/test_asteroid_vocabulary.py` checks every key the taxon mixes, the accessory lists
and the host table use.
"""

from planetgen import tuning
from planetgen.physics import materials

# Why a name is vocabulary-only: the setting that makes the mineral and the item that builds it.
MAGMATIC_SULFIDE = "forms in a magmatic sulfide ore of a differentiated body (ore processes, GEN.1380); a primitive body carries its platinum-group elements in its metal"
RARE_EARTH_ORE = "a rare-earth mineral of carbonatite, granite and pegmatite (ore processes, GEN.1382); a primitive body holds its rare earths in hibonite, perovskite, oldhamite and merrillite"
URANIUM_ORE = "a uranium or thorium ore mineral (ore processes, GEN.1381); a primitive body holds its uranium and thorium in the inclusion phases and merrillite"
SUPERGENE = "a supergene mineral of a weathered ore on a wet, oxidising body (ore processes, GEN.1381 and GEN.1382)"
IGNEOUS = "an igneous, pegmatitic or metamorphic mineral of an evolved crust (the crust and rock types of a differentiated body), not of a primitive taxon's mix"
TITANIUM_OXIDE = "a titanium oxide or silicate of igneous and placer settings (ilmenite and perovskite carry a primitive body's titanium)"
EVAPORITE = "an evaporite salt (the evaporite sequence, GEN.1366), made where a brine dries"
OXIDISED = "an oxidation product of iron on a weathered surface (a primitive body's iron oxide is the magnetite its mix draws)"

PGE_HOSTS = ("kamacite", "taenite")
"""tuple: The rows that carry the platinum-group elements in a primitive body: the iron-nickel metal, at 0.01 to 10 ppm (stand-in ACV1)."""

HOSTED = {
    # the major and minor elements: the minerals that hold each in a primitive body's mix
    "carbon_graphite": ("carbon_graphite", "carbon_diamond", "complex_organics"),
    "silicon": ("forsterite", "fayalite", "enstatite", "ferrosilite", "diopside", "albite", "anorthite"),
    "magnesium": ("forsterite", "enstatite", "spinel", "serpentine"),
    "aluminum": ("anorthite", "albite", "spinel", "hibonite"),
    "calcium": ("anorthite", "diopside", "oldhamite", "calcite"),
    "sulfur": ("troilite", "pyrrhotite", "pentlandite"),
    "phosphorus": ("schreibersite", "merrillite", "chlorapatite"),
    "iron": ("kamacite", "taenite", "troilite", "fayalite", "magnetite"),
    "nickel": ("kamacite", "taenite", "pentlandite"),
    # the platinum-group elements and gold in the metal; the inclusion nuggets in the inclusion phases
    "iridium": PGE_HOSTS,
    "palladium": PGE_HOSTS,
    "platinum": PGE_HOSTS,
    "osmium": PGE_HOSTS,
    "ruthenium": PGE_HOSTS,
    "rhodium": PGE_HOSTS,
    "gold": ("kamacite", "taenite", "troilite"),
    "osmiridium": ("hibonite", "melilite", "spinel", "perovskite"),
    "titanium": ("ilmenite", "perovskite", "chromite"),
    # the mineral groups: the members the mixes draw
    "olivine": ("forsterite", "fayalite"),
    "pyroxene": ("enstatite", "ferrosilite", "diopside"),
    "plagioclase": ("albite", "anorthite"),
    # the organic matter holds the polycyclic aromatic hydrocarbons
    "naphthalene": ("complex_organics",),
    "anthracene": ("complex_organics",),
    "phenanthrene": ("complex_organics",),
    "pyrene": ("complex_organics",),
    "coronene": ("complex_organics",),
    "fluoranthene": ("complex_organics",),
    # the names that are themselves rows of the mixes
    "kamacite": ("kamacite",), "taenite": ("taenite",), "troilite": ("troilite",), "schreibersite": ("schreibersite",), "cohenite": ("cohenite",),
    "serpentine": ("serpentine",), "magnetite": ("magnetite",), "chromite": ("chromite",), "moissanite": ("moissanite",), "ilmenite": ("ilmenite",),
    "perovskite": ("perovskite",), "spinel": ("spinel",), "calcite": ("calcite",),
}
"""dict: Row id of a vocabulary name -> the rows of a primitive body's mix that host it (`hosted`)."""

VOCABULARY_ONLY = {
    "hematite": OXIDISED,
    "sperrylite": MAGMATIC_SULFIDE, "cooperite": MAGMATIC_SULFIDE, "braggite": MAGMATIC_SULFIDE, "laurite": MAGMATIC_SULFIDE, "vysotskite": MAGMATIC_SULFIDE,
    "bastnasite": RARE_EARTH_ORE, "monazite": RARE_EARTH_ORE, "xenotime": RARE_EARTH_ORE, "cerite": RARE_EARTH_ORE, "gadolinite": RARE_EARTH_ORE,
    "samarskite": RARE_EARTH_ORE, "fergusonite": RARE_EARTH_ORE, "euxenite": RARE_EARTH_ORE,
    "uraninite": URANIUM_ORE, "thorianite": URANIUM_ORE, "coffinite": URANIUM_ORE, "brannerite": URANIUM_ORE,
    "carnotite": SUPERGENE, "autunite": SUPERGENE, "torbernite": SUPERGENE, "turquoise": SUPERGENE, "malachite": SUPERGENE,
    "rutile": TITANIUM_OXIDE, "anatase": TITANIUM_OXIDE, "brookite": TITANIUM_OXIDE, "titanite": TITANIUM_OXIDE,
    "quartz": IGNEOUS, "corundum": IGNEOUS, "beryl": IGNEOUS, "topaz": IGNEOUS, "zircon": IGNEOUS, "schorl": IGNEOUS, "orthoclase": IGNEOUS, "fluorite": IGNEOUS,
    "fluorapatite": IGNEOUS, "almandine": IGNEOUS, "pyrope": IGNEOUS, "spessartine": IGNEOUS, "grossular": IGNEOUS, "muscovite": IGNEOUS, "biotite": IGNEOUS, "chrysoberyl": IGNEOUS,
    "halite": EVAPORITE,
}
"""dict: Row id of a vocabulary name no taxon mix makes -> why (the setting that forms it and the item that builds it)."""


def canonical(name):
    """The row id a vocabulary name resolves to (`silicon dioxide` is `quartz`, `plagioclase feldspars` is `plagioclase`)."""
    return materials.material(name).id


def _check_rows():
    """Every host row and every listed name is a row of the materials table, written by its id and not an alias (a row with no entry is refused when the module loads)."""
    for name, rows in HOSTED.items():
        assert canonical(name) == name, f"{name} is written by an alias; use the row id"
        for row in rows:
            assert canonical(row) == row, f"{row} (a host of {name}) is written by an alias; use the row id"
    for name in VOCABULARY_ONLY:
        assert canonical(name) == name, f"{name} is written by an alias; use the row id"


_check_rows()


def hosted(components):
    """
    The vocabulary names a drawn composition carries, with the row that hosts each: for a composition `components` (row -> mass fraction), every hosted name whose host rows it holds, mapped to
    the largest of them. A statement that the body can carry the name, never an abundance. A name that is itself the only row of its host (`kamacite` in `kamacite`) is left out: it is in the
    composition already.
    """
    found = {}
    for name, rows in HOSTED.items():
        if rows == (name,):
            continue
        present = [row for row in rows if components.get(row, 0.0) > 0.0]
        if present:
            found[name] = max(present, key=lambda row: (components[row], row))
    return dict(sorted(found.items()))
