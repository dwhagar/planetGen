# planetgen/generation/phenomena/registry.py

"""
The phenomena registry (GEN.250): one row per kind of object the generator
makes, in one place. A kind used to be spelled out in five places (the
sector roll list, the scatter list, the density table, the object ID codes
and the class pages); each of those now reads its list from here, and
`test_phenomena_registry` checks that every other place that names a kind
(the database tables, the insert functions, the class pages, the scatter
table's constraint) agrees with it. Adding a kind is one row, a builder and
a test.

This module holds data only and imports nothing from the package, so the
ID codec, the tuning tables and the generators can all import it.

Columns (`Kind`):
    kind: The kind's spelling everywhere (the density key, the scatter
        row's `kind`, the object ID kind).
    family: `star`, `compact-remnant`, `nebula`, `small-body` or
        `nucleus`.
    status: `built` (every row today; a kind that is only designed is
        added by the item that builds it).
    rate: The key of the kind's density in
        `tuning.PHENOMENON_DENSITY_PC3` (the number stays there, with its
        sources), or `None` for a kind with no per-star rate.
    gate: What decides whether a cell can hold the kind today:
        `stellar density` (the cell's star count scales the rate),
        `nebula field` (the seeded cloud field), `hot stars` (grows around
        a sector's own hot stars) or `none`.
    mode: How the kind comes about: `scatter` (a row placed after the
        plan), `fill` (made when a sector is filled), `field` (derived from
        a seed), `flag` (an attribute of an ordinary system) or
        `component` (made by its parent).
    scatter_pass: For a `scatter` kind, which pass places it:
        `band-density` (the per-band Poisson walk of
        `phenomenon_scatter`), `bright-star` (the bright-star scatter),
        `wanderer` (hypervelocity stars) or `nucleus` (the one nucleus
        row); else `None`.
    sector_roll: `True` when `generate_sector_phenomena` rolls the kind
        for a sector that was not pre-scattered (draw order is registry
        order).
    builder: A `module:attr` path to what makes the kind. Every
        `sector_roll` builder takes `(config, galactic_center_dist_ly)`.
    system_model: A key of `SYSTEM_MODELS`: the system the kind gets today.
    carries: Kinds the object holds inside it, by registry kind.
    id_code: The 6-bit type code of the object ID (`names/object_id`), or
        `None` when the kind's objects take their `phenomenon_type`'s code.
    page_class: The slug of the class page (`web/lib/classref`), or `None`.
    marker: The map marker (MAP.132 builds these), or `None`.
    phenomenon_type: The `--type` value / table type the kind is saved as
        (`tuning.PHENOMENON_TYPE_CHOICES`), or `None` for a kind that is not
        a stored phenomenon.
    evidence_note: Where the rate and the rules come from.
"""

import importlib
from dataclasses import dataclass


SYSTEM_MODELS = {
    "none": "No system: the object stands alone.",
    "star-system": "An ordinary star system: companion, planets, belts and moons by the star model.",
    "central-star": "Its own new hot white dwarf system, made with the object.",
    "rogue-moon": "One captured moon with the chance of SYSTEM_PASS_ODDS['rogue-planet'].",
    "remnant-planets": "Planets around a compact remnant at the odds of SYSTEM_PASS_ODDS, rolled by the fill's system pass.",
}
"""dict: Every system model a row may name, with what it builds."""

MODES = ("scatter", "fill", "field", "flag", "component")
"""tuple: The valid `Kind.mode` values."""

SCATTER_PASSES = ("band-density", "bright-star", "wanderer", "nucleus")
"""tuple: The valid `Kind.scatter_pass` values."""

FAMILIES = ("star", "compact-remnant", "nebula", "small-body", "nucleus")
"""tuple: The valid `Kind.family` values."""

GATES = ("stellar density", "nebula field", "hot stars", "none")
"""tuple: The valid `Kind.gate` values."""

COLUMNS = ("kind", "family", "status", "rate", "gate", "mode", "scatter_pass", "sector_roll", "builder",
           "system_model", "carries", "id_code", "page_class", "marker", "phenomenon_type", "evidence_note")
"""tuple: Every column of a row, in order; the completeness test reads it."""


@dataclass(frozen=True)
class Kind:
    """One row of the registry (see the module docstring for the columns)."""
    kind: str
    family: str
    status: str
    rate: str | None
    gate: str
    mode: str
    scatter_pass: str | None
    sector_roll: bool
    builder: str | None
    system_model: str
    carries: tuple
    id_code: int | None
    page_class: str | None
    marker: str | None
    phenomenon_type: str | None
    evidence_note: str


_P = "planetgen.generation.phenomena.builders:"
_S = "planetgen.generation.run_sector:"

KINDS = (
    Kind("star-system", "star", "built", None, "stellar density", "fill", None, False,
         "planetgen.generation.system:StarSystem", "star-system", (), None, None, None, None,
         "The ordinary system every star gets; docs/design/phenomena-catalogue-and-fill.md 4.1."),
    Kind("bright-star", "star", "built", None, "stellar density", "scatter", "bright-star", False,
         _S + "_add_preplaced_systems", "star-system", (), 9, None, None, None,
         "The bright and massive stars placed by the bright-star scatter (bright_stars.py)."),
    Kind("black-hole", "compact-remnant", "built", "black-hole", "stellar density", "scatter", "band-density",
         True, _P + "build_black_hole", "none", (), 2, "black-hole", None, "black-hole",
         "About 1e8 in the Galaxy; tuning.PHENOMENON_DENSITY_PC3."),
    Kind("neutron-star", "compact-remnant", "built", "neutron-star", "stellar density", "scatter", "band-density",
         True, _P + "build_neutron_star", "remnant-planets", (), 3, None, None, "neutron-star",
         "About 1e9 in the Galaxy; tuning.PHENOMENON_DENSITY_PC3."),
    Kind("black-hole-core", "compact-remnant", "built", None, "none", "component", None, False,
         None, "none", (), 10, None, None, None,
         "The collapsed core of a core-collapse supernova remnant (a name of its own)."),
    Kind("neutron-star-core", "compact-remnant", "built", None, "none", "component", None, False,
         None, "none", (), 11, None, None, None,
         "The collapsed core of a core-collapse supernova remnant (a name of its own)."),
    Kind("nebula", "nebula", "built", None, "hot stars", "fill", None, False,
         _P + "build_nebula", "none", (), 4, "nebula", None, "nebula",
         "Emission and reflection nebulae grow around a sector's own hot stars (add_star_hosted_nebulae)."),
    Kind("planetary-nebula", "nebula", "built", "planetary-nebula", "stellar density", "scatter", "band-density",
         True, _P + "build_planetary_nebula", "central-star", ("star-system",), None, None, None, "nebula",
         "About 20,000 in the Galaxy (Frew & Parker 2010); built around its own white dwarf."),
    Kind("molecular-cloud", "nebula", "built", "molecular-cloud", "nebula field", "field", None, True,
         _P + "build_molecular_cloud", "none", (), None, None, None, "nebula",
         "Giant molecular clouds (Kennicutt & Evans 2012), from the seeded nebula field."),
    Kind("supernova-remnant", "nebula", "built", "supernova-remnant", "stellar density", "scatter", "band-density",
         True, _P + "build_supernova_remnant", "none", ("neutron-star-core", "black-hole-core"), 5,
         "supernova-remnant", None, "supernova-remnant",
         "About 3,000 in the Galaxy (Green 2019); a core-collapse remnant may hold a compact core."),
    Kind("rogue-planet", "small-body", "built", "rogue-planet", "stellar density", "fill", None, True,
         _P + "build_rogue_planet", "rogue-moon", (), 1, "rogue-planet", None, "rogue-planet",
         "Mroz et al. 2017 limit for giants; six mass bins (tuning.ROGUE_PLANET_MASS_BINS)."),
    Kind("brown-dwarf", "small-body", "built", "brown-dwarf", "stellar density", "fill", None, True,
         _P + "build_brown_dwarf", "none", (), None, None, None, "rogue-planet",
         "0.025 to 0.035 per pc^3 (Kirkpatrick et al. 2021); saved as a rogue planet of the brown-dwarf bin."),
    Kind("comet", "small-body", "built", "comet", "stellar density", "fill", None, True,
         _P + "build_comet", "none", (), 7, "comet", None, "comet",
         "Notable interstellar comets kept as rows, a design rate."),
    Kind("asteroid-field", "small-body", "built", "asteroid-field", "none", "fill", None, True,
         _P + "build_asteroid_field", "none", (), 8, "asteroid-field", None, "asteroid-field",
         "Rate 0: free fields disperse in 1e6 to 1e7 years (Raymond et al. 2020); hand-made only."),
    Kind("quasar", "nucleus", "built", None, "none", "scatter", "nucleus", False,
         _P + "build_quasar", "none", (), 6, None, None, "quasar",
         "One per galaxy, at the nucleus, with QUASAR_ACTIVE_NUCLEUS_CHANCE."),
    Kind("runaway-star", "star", "built", "runaway-star", "stellar density", "flag", None, False,
         _S + "flag_fast_stars", "star-system", (), None, None, None, None,
         "Stars above 30 km/s, 1 to 2% (Tauris 2015); a flag on an ordinary system."),
    Kind("hypervelocity-star", "star", "built", "hypervelocity-star", "stellar density", "scatter", "wanderer",
         False, _S + "_add_hypervelocity_star", "star-system", (), None, None, None, None,
         "Stars above 500 km/s ejected by the central black hole (Brown 2015); an ordinary young system at fill."),
)
"""tuple: The registry, one `Kind` per kind, in the order the generators draw
them: the sector roll order and the `--type` order are this order."""

BY_KIND = {row.kind: row for row in KINDS}
"""dict: `kind` -> its row."""


def rows(**match):
    """The rows whose columns equal every keyword in `match`, in registry order."""
    return tuple(row for row in KINDS if all(getattr(row, column) == value for column, value in match.items()))


def kinds(**match):
    """The `kind` of each row `rows(**match)` returns."""
    return tuple(row.kind for row in rows(**match))


def phenomenon_types():
    """The stored phenomenon types (`tuning.PHENOMENON_TYPE_CHOICES`) in registry order."""
    return tuple(dict.fromkeys(row.phenomenon_type for row in KINDS if row.phenomenon_type))


def id_codes():
    """`{kind: id_code}` for each kind with an object ID code, by code."""
    return {row.kind: row.id_code for row in sorted((r for r in KINDS if r.id_code is not None),
                                                    key=lambda r: r.id_code)}


def rated_kinds():
    """The kinds with a density in `tuning.PHENOMENON_DENSITY_PC3`, in registry order."""
    return tuple(row.rate for row in KINDS if row.rate)


def resolve(path):
    """The object a `module:attr` path names (a builder), imported on demand."""
    module, _, attr = path.partition(":")
    return getattr(importlib.import_module(module), attr)
