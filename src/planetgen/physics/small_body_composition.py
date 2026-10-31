# planetgen/physics/small_body_composition.py

"""
What belts, comets, icy bodies and rings are made of (GEN.442, docs/design/small-bodies-composition.md sections 3 to 5)
=======================================================================================================================

Every small body draws its composition first and derives its bulk density from it (decision S1): mass fractions of the rows
of the one shared materials table (GEN.427), a porosity, and from them `physics.compaction`'s mixing rule. Four kinds:

* **A belt or field of asteroids** (`draw_belt`). Its mass is shared out among the taxa of Bus and DeMeo by a Dirichlet draw
  around a kernel of the disk temperature at the belt's distance (the blackbody one for the star's luminosity) and the aluminium-26
  its star was born with (`taxon_kernel`, GEN.473): dry and igneous (the more so the more aluminium-26) where the disk is above the
  water line, hydrated carbonaceous below it, primitive and ice-rich colder still. Each taxon's mix of materials is a meteorite analogue,
  so every taxon row stays a flagged `stand-in` that GEN.475 replaces, and the Dirichlet concentration is open (Boss, 1:43 and 1:44 AM EDT on 2026-10-17). Each taxon is a mix of the shared table's
  rows (the chemistry of its meteorite analogue), a porosity from its family's range, and an albedo.
* **A comet** (`draw_comet`). A refractory-to-ice ratio (Choukroun et al. 2020), the organic share of the refractories, a
  porosity (65 to 79 percent, Choukroun), and the volatile percentages of water of the population of GEN.495
  (`comet_population.draw_volatiles`), turned into mass fractions of the ice. The mass is the density times the nucleus sphere.
* **An icy body** (`draw_icy_body`): a Kuiper belt object, a dwarf planet or a class D body. Its rock fraction is a truncated
  normal about 0.70 (Bierson and Nimmo 2019: 70 percent rock matches 14 of 17 Kuiper belt objects within 2 sigma; the width and
  limits are the repository's), or for a dwarf planet that lost its ice in a collision a second row (the weights are open); the bulk density is the
  tapered brittle law of `physics.compaction`.
* **A ring** (`draw_ring`): material of a circumplanetary disk, whose ice fraction rises with the formation distance past the
  disk's ice line (`disk_ice_fraction`, decision S8; Cilibrasi 2018 sets the bracket 0.5 to 0.85 for the moons).

What each draw returns is a `SmallBodyMakeup`: the mass fractions by component, their sums in four groups (metal, rock, water,
other ice: the shape of lane I's `BulkFractions`), the porosity, the grain and bulk density. The draws read only `util.draw` in a
fixed order, so a body is the same at any worker count; nothing here reads a Solar System body as a centre (the bodies of the
design note's section 8 are test targets).

Use:

    makeup = draw_belt(distance_au=3.1, luminosity_w=3.828e26, al26_ratio=5.2e-5)    # SmallBodyMakeup, taxon_shares in makeup.detail
    comet = draw_comet("jupiter_family")                            # refractory_ice_ratio, volatiles in comet.detail
    mass_kg = comet.mass_kg(nucleus_diameter_km=3.5)
    body = draw_icy_body("dwarf-planet")                            # then body.bulk_density_kg_m3(radius_km)
"""

import math
import statistics
from dataclasses import dataclass, field

from planetgen import tuning
from planetgen.generation import possibility_table as pt
from planetgen.physics import comet_population, compaction, condensate_shares, elements, materials, meteorite_groups, taxon_kernel
from planetgen.util import draw

# ---------------------------------------------------------------------------------------------------------------------
# The kinds the tables are gated by (never a Solar System body)

BELT = "belt"
FIELD = "asteroid-field"
COMET = "comet"
KUIPER_BELT_OBJECT = "kuiper-belt-object"
DWARF_PLANET = "dwarf-planet"
CLASS_D = "D"
ICY_KINDS = (CLASS_D, DWARF_PLANET, KUIPER_BELT_OBJECT)

DIRICHLET_CONCENTRATION = 10.0
"""float: The concentration of the Dirichlet draw of the taxon shares around the kernel. Open: no source gives the spread among
belts (the Solar System has one belt); a lower value is the wider reading (design note section 5.1)."""

# ---------------------------------------------------------------------------------------------------------------------
# The groups a component belongs to (the four-way split of lane I's BulkFractions)

METAL_ROWS = frozenset({"iron", "nickel", "kamacite", "taenite", "schreibersite", "cohenite", "iron_epsilon", "iron_liquid"})
"""frozenset: Rows that count as metal: the iron-nickel alloys, the phosphide and the carbide of a metal phase."""


def group_of(component):
    """`"metal"`, `"water"`, `"ice"` (every other ice) or `"rock"` for a row of the shared table."""
    row = materials.material(component)
    if row.id in METAL_ROWS or row.family == "metal":
        return "metal"
    if row.id in ("water_ice", "water_liquid"):
        return "water"
    if row.family in ("ice", "organic ice"):
        return "ice"
    return "rock"


@dataclass(frozen=True)
class SmallBodyMakeup:
    """
    What a small body is made of.

    Attributes:
        kind (str): `belt`, `asteroid-field`, `comet`, `ring` or an icy kind.
        components (dict): Row of the shared materials table -> mass fraction of the grains, adding up to 1.
        metal, rock, water, ice (float): The components summed by `group_of`; they add up to 1.
        porosity (float): The void share of the body's volume (a mean over the taxa for a belt).
        grain_density_kg_m3 (float): The zero-pressure mixing density of the grains.
        bulk_density_kg_m3 (float or None): `grain_density (1 - porosity)`; `None` for an icy body, whose density depends on
            its size (`SmallBodyMakeup.density_at`).
        albedo (float or None): The mean albedo where the draw sets one.
        detail (dict): What else was drawn (the taxon shares of a belt, the refractory-to-ice ratio of a comet, ...).
        densities (dict): Component -> kg/m3 for the components whose density this draw states and the table does not.
    """

    kind: str
    components: dict
    metal: float
    rock: float
    water: float
    ice: float
    porosity: float
    grain_density_kg_m3: float
    bulk_density_kg_m3: float | None
    albedo: float | None = None
    detail: dict = field(default_factory=dict)
    densities: dict = field(default_factory=dict)

    def groups(self):
        return {"metal": self.metal, "rock": self.rock, "water": self.water, "ice": self.ice}

    def to_dict(self):
        """A JSON-safe dict of everything the draw produced (`from_dict` is its inverse)."""
        return {"kind": self.kind, "components": dict(self.components), "porosity": self.porosity,
                "bulk_density_kg_m3": self.bulk_density_kg_m3, "albedo": self.albedo, "detail": _plain(self.detail),
                "densities": dict(self.densities)}

    @classmethod
    def from_dict(cls, data):
        """The makeup `to_dict` wrote: the groups and the grain density are worked out again from the components."""
        made = _makeup(data["kind"], data["components"], data["porosity"], data.get("albedo"), data.get("detail"),
                       data.get("densities"), size_dependent=data.get("bulk_density_kg_m3") is None)
        # the stored fractions and density are kept as they were drawn (a round trip changes no digit)
        return cls(made.kind, dict(data["components"]), made.metal, made.rock, made.water, made.ice, made.porosity,
                   made.grain_density_kg_m3, data.get("bulk_density_kg_m3"), made.albedo, made.detail, made.densities)

    def mass_kg(self, nucleus_diameter_km):
        """The mass of a sphere of this diameter at the bulk density (a comet nucleus)."""
        if self.bulk_density_kg_m3 is None:
            raise ValueError("this body's density depends on its size: use density_at(radius_km)")
        radius_m = nucleus_diameter_km * 500.0
        return 4.0 / 3.0 * math.pi * radius_m ** 3 * self.bulk_density_kg_m3

    def density_at(self, radius_km):
        """The bulk density of an icy body of this radius (kg/m3): the tapered brittle law on its rock fraction."""
        return compaction.compacted_density_kg_m3(radius_km, self.rock)

    def bulk_fractions(self):
        """Lane I's `BulkFractions` for this body (no envelope; the porosity set)."""
        from planetgen.physics.bulk_composition import BulkFractions
        return BulkFractions(metal=self.metal, rock=self.rock, water=self.water, ice=self.ice, envelope=0.0,
                             porosity=self.porosity)


def _plain(value):
    """`value` with every dict and tuple turned into a JSON-safe dict or list."""
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


def _makeup(kind, components, porosity, albedo=None, detail=None, densities=None, size_dependent=False):
    total = math.fsum(components.values())
    components = {name: value / total for name, value in components.items() if value > 0.0}
    sums = {"metal": 0.0, "rock": 0.0, "water": 0.0, "ice": 0.0}
    for name, value in components.items():
        sums[group_of(name)] += value
    grain = compaction.grain_density_kg_m3(components, densities)
    bulk = None if size_dependent else grain * (1.0 - porosity)
    return SmallBodyMakeup(kind, components, sums["metal"], sums["rock"], sums["water"], sums["ice"], porosity, grain, bulk,
                           albedo, detail or {}, dict(densities or {}))


def makeup_from_components(kind, components, porosity, detail=None):
    """The makeup of a body whose grains are the rows `components` of the shared materials table (mass fractions, normalised here) at the porosity `porosity`; `kind` and `detail` as stored."""
    return _makeup(kind, components, porosity, detail=detail)


# ---------------------------------------------------------------------------------------------------------------------
# The taxa of a belt or field

ZONES = (
    # name, x from, x to, percent of the zone's mass by taxon (DeMeo and Carry 2013, Table 6, read in part D 2.4)
    ("hungaria", 0.66, 0.76, {"E": 50, "C": 21, "S": 9, "A": 7, "L": 7, "K": 5}),
    ("inner belt", 0.76, 0.93, {"V": 69, "S": 21, "C": 6}),
    ("middle belt", 0.93, 1.04, {"C": 70, "B": 15, "S": 8, "P": 4}),
    ("outer belt", 1.04, 1.21, {"C": 52, "P": 15, "B": 13, "M": 10, "S": 5}),
    ("cybele", 1.21, 1.37, {"P": 84, "C": 13}),
    ("hilda", 1.37, 1.56, {"P": 71, "D": 15, "C": 14}),
    ("trojans", 1.87, 2.00, {"D": 67, "P": 26}),
)
"""tuple: The zones of the Solar System's belt in units of its snow line (2.7 AU): the only place a Solar System number sets a
boundary, and why every taxon row is a stand-in (GEN.473 replaces them)."""

TAXA = ("E", "C", "B", "S", "A", "L", "K", "V", "P", "M", "D")
"""tuple: The taxa, in the table's fixed row order."""

CARBONACEOUS = "carbonaceous"
ROCKY = "rocky"
POROSITY_RANGES = {CARBONACEOUS: (0.25, 0.55), ROCKY: (0.10, 0.35)}
"""dict: The surface porosity ranges by family (Britt and Consolmagno 2003, Carry 2012; parametrisations)."""


@dataclass(frozen=True)
class Taxon:
    """
    One taxon of a belt.

    Attributes:
        name (str): The Bus and DeMeo class letter.
        analogue (str): The meteorite group whose chemistry the mix follows (proposed links, part D 2.3).
        family (str): `carbonaceous` or `rocky`, which sets the porosity range.
        albedo (float): The mean albedo of the class (DeMeo and Carry 2013 Table 4, SDSS).
        fractions (dict): Row of the shared table -> mass fraction of the grains.
    """

    name: str
    analogue: str
    family: str
    albedo: float
    fractions: dict


TAXON_DATA = {
    # The mixes are the analogue meteorite's modal chemistry written in rows of the shared table: a stand-in at every row
    # (GEN.473), checked in the tests against the grain densities of Britt and Consolmagno 2003 (part D 2.2).
    "C": Taxon("C", "CI/CM chondrite", CARBONACEOUS, 0.083,
               {"serpentine": 0.56, "forsterite": 0.09, "fayalite": 0.07, "enstatite": 0.05, "magnetite": 0.07, "troilite": 0.06,
                "calcite": 0.03, "complex_organics": 0.06, "kamacite": 0.01}),
    "B": Taxon("B", "CI/CM chondrite", CARBONACEOUS, 0.071,
               {"serpentine": 0.50, "magnetite": 0.12, "forsterite": 0.12, "fayalite": 0.07, "enstatite": 0.06, "troilite": 0.05,
                "calcite": 0.03, "complex_organics": 0.05}),
    "S": Taxon("S", "H/L/LL ordinary chondrite", ROCKY, 0.258,
               {"forsterite": 0.28, "fayalite": 0.17, "enstatite": 0.19, "ferrosilite": 0.06, "diopside": 0.04, "albite": 0.09,
                "kamacite": 0.08, "taenite": 0.03, "troilite": 0.06}),
    "V": Taxon("V", "howardite-eucrite-diogenite", ROCKY, 0.352,
               {"enstatite": 0.20, "ferrosilite": 0.20, "diopside": 0.15, "anorthite": 0.40, "kamacite": 0.02, "troilite": 0.03}),
    "A": Taxon("A", "brachinite and pallasite mantle", ROCKY, 0.274,
               {"forsterite": 0.80, "fayalite": 0.08, "kamacite": 0.08, "troilite": 0.04}),
    "E": Taxon("E", "aubrite (enstatite chondrite)", ROCKY, 0.536,
               {"enstatite": 0.88, "troilite": 0.05, "kamacite": 0.04, "albite": 0.03}),
    "M": Taxon("M", "iron, enstatite and mesosiderite", ROCKY, 0.143,
               {"kamacite": 0.50, "taenite": 0.12, "enstatite": 0.20, "forsterite": 0.07, "troilite": 0.06, "schreibersite": 0.05}),
    "K": Taxon("K", "CO/CV chondrite", ROCKY, 0.178,
               {"forsterite": 0.30, "fayalite": 0.12, "enstatite": 0.12, "serpentine": 0.10, "anorthite": 0.08, "kamacite": 0.08,
                "magnetite": 0.08, "troilite": 0.06, "complex_organics": 0.06}),
    "L": Taxon("L", "CV chondrite with refractory inclusions", ROCKY, 0.183,
               {"forsterite": 0.30, "fayalite": 0.10, "anorthite": 0.15, "diopside": 0.10, "spinel": 0.05, "serpentine": 0.08,
                "kamacite": 0.06, "magnetite": 0.06, "troilite": 0.05, "complex_organics": 0.05}),
    "P": Taxon("P", "anhydrous silicate with ice and organics", CARBONACEOUS, 0.053,
               {"forsterite": 0.33, "enstatite": 0.31, "complex_organics": 0.10, "troilite": 0.03, "magnetite": 0.03,
                "water_ice": 0.20}),
    "D": Taxon("D", "ice-rich, organic-rich", CARBONACEOUS, 0.098,
               {"forsterite": 0.28, "enstatite": 0.27, "complex_organics": 0.15, "troilite": 0.05, "water_ice": 0.25}),
}
"""dict: The taxa by letter."""
TAXON_DATA = {name: TAXON_DATA[name] for name in TAXA}


REPLACEMENT = "the meteorite-group map (each taxon's mix of materials is still a meteorite analogue; the weights are derived, `taxon_kernel`)"
"""str: What replaces the stand-in taxon rows' mixes of materials (the meteorite-group map of GEN.475)."""


def _taxon_row(taxon):
    data = TAXON_DATA[taxon]
    return pt.Row(
        name=f"taxon-{taxon}", classes=(BELT, FIELD), weight=lambda subject, taxon=taxon: taxon_kernel.kernel(subject["temperature_k"], subject["al26_ratio"])[taxon],
        source="the belt's taxa from the disk's ice lines and the early heating by aluminium-26 (taxon_kernel); the mix of each taxon is a meteorite analogue (DeMeo and Carry 2013, arXiv:1307.2424)",
        state=pt.STAND_IN, replacement=REPLACEMENT, fractions=data.fractions, concentration=None,
        notes=f"{data.analogue}; {data.family} family; albedo {data.albedo}")


TAXON_TABLE = pt.REGISTRY.register(pt.Table(
    "small-body-taxa", 1, [_taxon_row(t) for t in TAXA], components=pt.SHARED_COMPONENTS,
    description="The taxa of a belt or field of asteroids: the weight is the share of the belt's mass at the disk temperature of its "
                "distance for a star born with a given aluminium-26 ratio, drawn as a Dirichlet around it; the fractions are the taxon's mix of the shared materials."))
"""Table: The belt and field taxa, a flagged stand-in at every row until GEN.475 (the weights are derived)."""


def draw_belt(distance_au, luminosity_w, al26_ratio, kind=BELT, concentration=DIRICHLET_CONCENTRATION):
    """
    The composition of a belt or a field of asteroids at `distance_au` from a star of luminosity `luminosity_w` that was born with the aluminium-26 ratio `al26_ratio`
    (`radiogenic.star_al26_ratio`): the disk is at its blackbody temperature there and `taxon_kernel` gives the taxa.

    The draws, in order: the Dirichlet shares of the taxa the kernel gives weight (in table order), then one porosity for each of
    those taxa. Each taxon's bulk density is its grains' mixing density and its porosity; the belt's is the mass-weighted
    harmonic mean of those, and its albedo the mass-weighted mean.

    Returns:
        SmallBodyMakeup: With `detail["taxon_shares"]` (all `TAXA`, zero where the kernel gave none), `detail["temperature_k"]`, `detail["al26_ratio"]`,
        `detail["taxon_porosity"]` and `detail["taxon_bulk_density_kg_m3"]`.
    """
    if not distance_au > 0.0 or not luminosity_w > 0.0:
        raise ValueError("the distance and the luminosity must be positive")
    if not al26_ratio >= 0.0:
        raise ValueError("the aluminium-26 ratio must be zero or more")
    temperature_k = condensate_shares.disk_temperature_k(distance_au, luminosity_w)
    subject = pt.Subject(kind, {"temperature_k": temperature_k, "al26_ratio": al26_ratio})
    probabilities = TAXON_TABLE.probabilities(subject)
    names = [row.name[len("taxon-"):] for row, _p in probabilities]
    weights = [p for _row, p in probabilities]
    return _from_taxa(kind, names, weights, concentration, {"temperature_k": temperature_k, "al26_ratio": al26_ratio})


METAL_ROWS = ("kamacite", "taenite", "schreibersite")
"""tuple: The rows of the materials table that are the iron-nickel metal of the group map's metal axis."""

REFRACTORY_ROWS = ("spinel", "anorthite")
"""tuple: The rows that stand for the refractory inclusions (calcium-aluminium-rich) of its refractory axis."""

MATRIX_ROWS = ("serpentine", "magnetite", "complex_organics", "calcite")
"""tuple: The rows of the fine-grained matrix of its matrix axis: the water-altered minerals and the organics (ice is not matrix; an icy mix is far from every group, `ICE_MAX`)."""

OLIVINE_MIN_VOLUME = 0.05
"""float: The share of the grains' volume below which a mix has no olivine to give an oxidation state. [Stand-in MG3, tier T9, tweak L1.]"""

ICE_MAX = 0.05
"""float: The ice (and water) mass share from which a mix is far from every chondrite group: the groups are rocks, and the nearest one would mislabel an icy body. [Stand-in MG3, tier T9, tweak L1.]"""


def group_axes(components):
    """
    `(fayalite mol percent or None, metal, refractory, matrix volume percent)` of a mix of grains, the axes of `meteorite_groups`: the olivine's iron against magnesium from its two
    rows, and the volume percent of the metal, refractory and matrix rows of `components` (mass fractions of rows of the materials table).
    """
    volumes = {name: share / (compaction.grain_density_kg_m3({name: 1.0})) for name, share in components.items() if share > 0.0}
    total = math.fsum(volumes.values())
    percent = lambda rows: 100.0 * math.fsum(volumes.get(name, 0.0) for name in rows) / total    # noqa: E731
    fayalite = components.get("fayalite", 0.0) / (materials.material("fayalite").molar_mass * 1000.0)
    forsterite = components.get("forsterite", 0.0) / (materials.material("forsterite").molar_mass * 1000.0)
    olivine_volume = percent(("fayalite", "forsterite"))
    iron = 100.0 * fayalite / (fayalite + forsterite) if olivine_volume >= 100.0 * OLIVINE_MIN_VOLUME else None
    return iron, percent(METAL_ROWS), percent(REFRACTORY_ROWS), percent(MATRIX_ROWS)


def chondrite_group(components):
    """The nearest chondrite group of a mix of grains (`meteorite_groups.nearest_group` of `group_axes`) and its distance: `none` for a mix far from every group or with `ICE_MAX` of ice."""
    group, distance = meteorite_groups.nearest_group(*group_axes(components))
    ice = math.fsum(share for name, share in components.items() if group_of(name) in ("water", "ice"))
    return (meteorite_groups.NONE if ice >= ICE_MAX else group), distance


def _from_taxa(kind, names, weights, concentration, detail, porosity=None):
    """
    The makeup of a mixture of the taxa `names` around the centre `weights`: the Dirichlet shares, then a porosity for each taxon
    (`porosity` states one for all, drawing none); the rest as `draw_belt` says.
    """
    drawn = draw.dirichlet([concentration * w for w in weights])
    shares = {taxon: 0.0 for taxon in TAXA}
    shares.update(zip(names, drawn))
    taxon_porosity, bulk, components = {}, {}, {}
    for taxon in names:
        data = TAXON_DATA[taxon]
        if porosity is None:
            low, high = POROSITY_RANGES[data.family]
            taxon_porosity[taxon] = draw.uniform(low, high)
        else:
            taxon_porosity[taxon] = porosity
        grain = compaction.grain_density_kg_m3(data.fractions)
        bulk[taxon] = grain * (1.0 - taxon_porosity[taxon])
        total = math.fsum(data.fractions.values())
        for name, value in data.fractions.items():
            components[name] = components.get(name, 0.0) + shares[taxon] * value / total
    mean_bulk = 1.0 / math.fsum(shares[t] / bulk[t] for t in names)
    mean_porosity = math.fsum(shares[t] * taxon_porosity[t] for t in names)
    albedo = math.fsum(shares[t] * TAXON_DATA[t].albedo for t in names)
    makeup = _makeup(kind, components, 0.0, albedo)
    group, group_distance = chondrite_group(makeup.components)
    detail = dict(detail, group=group, group_distance=group_distance)
    # the belt's own figures: the grains' mixing density of the whole belt and the mass-weighted bulk density of its taxa
    return SmallBodyMakeup(kind, makeup.components, makeup.metal, makeup.rock, makeup.water, makeup.ice, mean_porosity,
                           makeup.grain_density_kg_m3, mean_bulk, albedo,
                           dict(detail, taxon_shares=shares, taxon_porosity=taxon_porosity, taxon_bulk_density_kg_m3=bulk),
                           makeup.densities)


FIELD_FAMILIES = {
    # family: taxa it is a mixture of (the labels of the field classes of tuning.ASTEROID_FIELD_COMPOSITIONS)
    "carbonaceous": ("C", "B", "P"),
    "stony": ("S", "L", "K"),
    "metallic": ("M",),
    "icy": ("P", "D"),
    "basaltic": ("V",),
    "mixed": TAXA,
    "dust": ("C",),
    "collisional": (),
}
"""dict: The taxa of each family of a standalone field (part D 5.1): carbonaceous = C, B and P; stony = S, L and K (Q is not a taxon of
the table); metallic = M; icy = P and D; basaltic = V; mixed = all; dust = CI and CM chemistry (taxon C) at porosity 0; collisional =
the fragments of one dominant taxon, drawn uniformly. The family names are labels; the grouping is a parametrisation."""
DUST_POROSITY = 0.0

FIELD_FAMILY_TABLE = pt.REGISTRY.register(pt.Table(
    "asteroid-field-family", 1, [
        pt.Row(name=family, classes=(FIELD,), weight=tuning.ASTEROID_FIELD_COMPOSITIONS[family]["frequency"],
               source="the frequency of the field class in tuning.ASTEROID_FIELD_COMPOSITIONS (no source)", state=pt.STAND_IN,
               replacement="the composition of the system whose bodies the field was stripped from, or a derived population "
                           "(docs/design/small-bodies-composition.md section 7)",
               notes=tuning.ASTEROID_FIELD_COMPOSITIONS[family]["description"])
        for family in FIELD_FAMILIES], components=pt.SHARED_COMPONENTS,
    description="The family of a standalone asteroid field; the weights are frequencies without a source."))
"""Table: The eight families of a standalone field (GEN.477 E3), every weight a flagged stand-in."""


def field_family_text(family):
    """What a family of standalone fields is made of, in words (the class reference page)."""
    if family == "mixed":
        return "A mixture of every asteroid taxon (" + ", ".join(TAXA) + ")"
    if family == "collisional":
        return "The fragments of one parent body: one asteroid taxon, drawn when the field is made"
    names = FIELD_FAMILIES[family]
    kinds = "; ".join(f"{name}: {TAXON_DATA[name].analogue}" for name in names)
    suffix = ", with no pores" if family == "dust" else ""
    return f"Asteroid taxa {kinds}{suffix}"


def draw_field_family():
    """The family of one standalone field: a draw of `FIELD_FAMILY_TABLE`."""
    return FIELD_FAMILY_TABLE.draw(FIELD).row.name


def draw_field(family, concentration=DIRICHLET_CONCENTRATION):
    """
    The composition of a standalone field of this family: a Dirichlet mixture of the family's taxa (a collisional field first draws
    its dominant taxon uniformly; a dust field has no pores).

    Returns:
        SmallBodyMakeup: As `draw_belt` makes it, with `detail["family"]`.
    """
    if family not in FIELD_FAMILIES:
        raise ValueError(f"the family {family!r} is none of {', '.join(FIELD_FAMILIES)}")
    names = FIELD_FAMILIES[family]
    if family == "collisional":
        names = (draw.choice(TAXA),)
    weights = [1.0 / len(names)] * len(names)
    return _from_taxa(FIELD, list(names), weights, concentration, {"family": family},
                      porosity=DUST_POROSITY if family == "dust" else None)


# ---------------------------------------------------------------------------------------------------------------------
# Comets

REFRACTORY_ICE_RATIO_RANGE = (0.5, 8.0)
"""tuple: The refractory-to-ice mass ratio, log-uniform (Choukroun et al. 2020: Rosetta RSI 3 to 7, CONSERT above 3, Fulle 7.5,
COSIMA 0.2 to 3; a verified range, a parametrised shape)."""
ORGANIC_SHARE = (0.60, 0.15)
"""tuple: The organic share of the refractories, normal (mean, sd) clipped to 0 to 1: a labelled stand-in at the widest bracket
(Boss, 1:23 PM EDT on 2026-10-17: no measurement exists; COSIMA's 0.45 lies inside it)."""
COMET_POROSITY_RANGE = (0.5, 0.8)
"""tuple: A comet nucleus' porosity, uniform (Choukroun et al. 2020: 65 to 79 percent)."""
KUIPER_POROSITY_RANGE = (0.3, 0.6)
"""tuple: The porosity of a Kuiper belt object under 100 km, uniform (design note section 5.2)."""
REFRACTORY_DENSITIES = {"silicate_dust": 3500.0}
"""dict: The grain density of a comet's silicate dust, kg/m3 (the shared table holds none for the row). A labelled stand-in: Blum and
Bischoff 2022 (arXiv:2207.12731) say the density of cometary refractory material "is not known and can only be constrained within a
range", 2600 to 5000 for the minerals; 3500 is the value the design note took from the Rosetta reviews and sits in that range. The
organics take the row of the shared table (1.5 g/cm3, the middle of Blum and Bischoff's 1000 to 2000; GEN.712)."""
ICE_STAND_IN_DENSITIES = {"ch4_ice": 500.0, "nh3_ice": 850.0, "c2h6_ice": 920.0, "hcn_ice": 920.0, "formaldehyde": 920.0,
                          "so2_ice": 920.0}
"""dict: Grain densities of ices the shared table holds no number for (kg/m3): the bracket of the design note (bulk-composition.md
section 6: other ices are close to water ice; NH3 850, CH4 500), a labelled stand-in until Fray and Schmitt 2009 is read. These
ices are a few percent of a comet's mass at most, so the density moves by under 1 percent."""
COMET_DENSITIES = {**REFRACTORY_DENSITIES, **ICE_STAND_IN_DENSITIES}
VOLATILE_ROWS = {"CO": "co_ice", "CO2": "co2_ice", "CH4": "ch4_ice", "C2H6": "c2h6_ice", "NH3": "nh3_ice", "HCN": "hcn_ice",
                 "CH3OH": "ch3oh_ice", "H2CO": "formaldehyde", "H2S": "h2s_ice"}
"""dict: The species of the volatile population (GEN.495) and the row of the shared table each is."""
WATER_FORMULA = "H2O"
COMET_DYNAMICAL_CLASS = {None: "long_period"}


def draw_comet(period_class=None, interstellar=False, ice_ratios=None):
    """
    The composition of a comet nucleus.

    The draws, in order: the refractory-to-ice ratio, the organic share, the porosity, then the volatiles of
    `comet_population.draw_volatiles`. The ice is the water and the volatiles at the number ratios the population gives, turned
    into masses with the molar masses; the refractories are silicate dust (a stand-in density, `REFRACTORY_DENSITIES`) and organics (the shared table's density).

    Args:
        period_class (str or None): `jupiter_family`, `halley_type`, `long_period`; `None` (a parabolic or interstellar
            comet) takes `long_period`.
        interstellar (bool): True for an interstellar comet (the wider carbon-oxide ranges).
        ice_ratios (tuple or None): `(CO2/H2O, CO/H2O)` of the disc zone the comet formed in (`comet_ice_budget.ice_ratios`): the median
            carbon oxides come from the zone (`comet_population.draw_volatiles`); `None` keeps the population's.

    Returns:
        SmallBodyMakeup: With `detail` holding `refractory_ice_ratio`, `organic_share`, `volatiles` (percent of water by number),
        `carbon_chain_depleted` and `ice_fraction`.
    """
    low, high = REFRACTORY_ICE_RATIO_RANGE
    delta = math.exp(draw.uniform(math.log(low), math.log(high)))
    organic = min(max(draw.gauss(*ORGANIC_SHARE), 0.0), 1.0)
    porosity = draw.uniform(*COMET_POROSITY_RANGE)
    volatiles = comet_population.draw_volatiles(COMET_DYNAMICAL_CLASS.get(period_class, period_class), interstellar=interstellar,
                                              ice_ratios=ice_ratios)
    ice_fraction = 1.0 / (1.0 + delta)
    water_mass = elements.molar_mass(WATER_FORMULA)
    masses = {"water_ice": water_mass}
    for species, percent in volatiles.percent.items():
        row = materials.material(VOLATILE_ROWS[species])
        masses[row.id] = percent / 100.0 * elements.molar_mass(row.formula)
    ice_total = math.fsum(masses.values())
    components = {name: ice_fraction * mass / ice_total for name, mass in masses.items()}
    components["silicate_dust"] = (1.0 - ice_fraction) * (1.0 - organic)
    components["complex_organics"] = (1.0 - ice_fraction) * organic
    detail = {"refractory_ice_ratio": delta, "organic_share": organic, "volatiles": dict(volatiles.percent),
              "carbon_chain_depleted": volatiles.carbon_chain_depleted, "ice_fraction": ice_fraction}
    return _makeup(COMET, components, porosity, detail=detail, densities=COMET_DENSITIES)


# ---------------------------------------------------------------------------------------------------------------------
# Icy bodies and rings

ROCK_FRACTION_LAWS = {
    # row name: (centre, sd, low, high) of a truncated normal in the rock mass fraction
    "icy-intact": (0.70, 0.08, 0.45, 0.95),
    "icy-stripped": (0.88, 0.05, 0.75, 0.98),
}
"""dict: The rock fraction of an icy body by row. The centre of the first is Bierson and Nimmo's 70 percent (verified); its width
and limits are the repository's. The second is the collisional tail of Barr and Schwamb 2016 (Haumea, Eris), whose weight is open."""

ICY_TABLE = pt.REGISTRY.register(pt.Table(
    "icy-body-rock-fraction", 1, [
        pt.Row("icy-intact", ICY_KINDS, pt.OPEN, "open", state=pt.OPEN_STATE,
               notes="rock mass fraction N(0.70, 0.08) on [0.45, 0.95]: the centre is Bierson and Nimmo 2019 (Icarus 326, 10: a "
                     "constant 70 percent matches 14 of 17 Kuiper belt objects within 2 sigma), the width and limits are the "
                     "repository's; the weight against the second row is open"),
        pt.Row("icy-stripped", (DWARF_PLANET,), pt.OPEN, "open", state=pt.OPEN_STATE,
               notes="a dwarf planet whose ice a collision removed (Barr and Schwamb 2016: Haumea, Eris): rock fraction toward 0.9; "
                     "the weight is open"),
    ], description="The rock mass fraction of a Kuiper belt object, dwarf planet or class D body, gated by kind."))
"""Table: The icy-body rock fraction rows."""


def _truncated_normal(centre, sd, low, high):
    """One truncated normal draw with one random number (inverse transform)."""
    normal = statistics.NormalDist(centre, sd)
    a, b = normal.cdf(low), normal.cdf(high)
    q = min(max(a + draw.random() * (b - a), 1e-16), 1.0 - 1e-16)
    return min(max(normal.inv_cdf(q), low), high)


def draw_icy_body(kind=DWARF_PLANET):
    """
    The composition of an icy body: a row of `ICY_TABLE` for its kind, then the rock fraction of that row.

    The ice is water ice (the other volatile ices a cold body keeps are the exosphere's business, GEN.455), and the body has no
    metal (a differentiated core is the planetary bulk draw's, GEN.428). The density depends on the size: `density_at(radius_km)`.

    Returns:
        SmallBodyMakeup: With `detail["row"]`, `detail["rock_fraction"]`; `bulk_density_kg_m3` is `None`.
    """
    possibility = ICY_TABLE.draw(kind)
    rock = _truncated_normal(*ROCK_FRACTION_LAWS[possibility.name])
    components = {"forsterite": rock * 0.5, "enstatite": rock * 0.5, "water_ice": 1.0 - rock}
    low, high = KUIPER_POROSITY_RANGE
    return _makeup(kind, components, draw.uniform(low, high), detail={"row": possibility.name, "rock_fraction": rock},
                   size_dependent=True)


ICE_LINE_RADII_COEFFICIENT = 0.4
"""float: The ramp starts at this fraction of the ice line (decision S8)."""
RAMP_LENGTH = math.log(2.5)
F_INF = (0.35, 0.10, 0.1, 0.6)
"""tuple: The asymptotic ice fraction of a disk moon: normal (mean, sd) on [low, high] (decision S8)."""


def disk_ice_fraction(r_over_r_ice, f_inf):
    """The ice mass fraction of material that formed at `r / r_ice`: `f_inf clamp(ln(r / (0.4 r_ice)) / ln 2.5, 0, 1)` (decision S8)."""
    if not r_over_r_ice > 0.0:
        raise ValueError("the distance must be positive")
    return f_inf * min(max(math.log(r_over_r_ice / ICE_LINE_RADII_COEFFICIENT) / RAMP_LENGTH, 0.0), 1.0)


def draw_ring(r_over_r_ice):
    """
    The composition of ring material that formed at `r / r_ice` in its parent's circumplanetary disk: ice by the ramp of
    `disk_ice_fraction` with a drawn asymptote, the rest rock. A ring has no porosity of its own here (its particles are grains).

    Returns:
        SmallBodyMakeup: With `detail["r_over_r_ice"]` and `detail["f_inf"]`.
    """
    f_inf = _truncated_normal(*F_INF)
    ice = disk_ice_fraction(r_over_r_ice, f_inf)
    components = {"forsterite": (1.0 - ice) * 0.5, "enstatite": (1.0 - ice) * 0.5}
    if ice > 0.0:
        components["water_ice"] = ice
    return _makeup("ring", components, 0.0, detail={"r_over_r_ice": r_over_r_ice, "f_inf": f_inf})


# ---------------------------------------------------------------------------------------------------------------------
# The words a body's description uses

LEVELS = (("high", 0.20), ("moderate", 0.05), ("small", 0.01), ("trace", 0.0))
"""tuple: The four-level words and the mass fraction each starts at (above 20 percent high, 5 to 20 moderate, 1 to 5 small, below 1
trace). Stored, flagged parametrisations (GEN.477): the words come from the fractions by these thresholds and not from a draw."""
TEXT_LINES = 4
"""int: At most this many named lines (the description and the database columns still fit)."""
METAL_LABEL = "iron-nickel metal"


def display_label(component):
    """
    The name a description gives a component: the group of the shared table that holds it (olivine for forsterite and fayalite,
    pyroxene, plagioclase feldspars), `METAL_LABEL` for the iron-nickel alloys, otherwise the row's own name.
    """
    if group_of(component) == "metal":
        return METAL_LABEL
    for group in ("olivine", "pyroxene", "plagioclase"):
        row = materials.material(group)
        if materials.material(component).id in row.members:
            return row.name
    return materials.material(component).name


def level_of(fraction):
    """The four-level word of a mass fraction."""
    for word, start in LEVELS:
        if fraction >= start:
            return word
    return LEVELS[-1][0]


def describe(makeup, lines=TEXT_LINES):
    """
    The `(label, level)` pairs of a body's description, largest first, at most `lines` of them: the components summed by
    `display_label`, each with the four-level word of its mass fraction.
    """
    shares = {}
    for name, value in makeup.components.items():
        label = display_label(name)
        shares[label] = shares.get(label, 0.0) + value
    ranked = sorted(shares.items(), key=lambda item: (-item[1], item[0]))[:lines]
    return [(label, level_of(value)) for label, value in ranked]


def component_names(makeup, lines=3):
    """The names of the `lines` largest components of a comet (the shared table's own names), largest first."""
    ranked = sorted(makeup.components.items(), key=lambda item: (-item[1], item[0]))[:lines]
    return [materials.material(name).name for name, _value in ranked]
