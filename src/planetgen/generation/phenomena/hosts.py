# planetgen/generation/phenomena/hosts.py

"""
Host classes and what can orbit them (GEN.552; docs/design/architecture/unified-star-like-objects.md section 3.3)
===================================================================================================

The registry used to say what a kind holds with `system_model`, one of five fixed values, and the odds of the two things
that hang on an object beyond its stars were a dict in `tuning` (`SYSTEM_PASS_ODDS`). A kind now names a **host class**
(`registry.Kind.host`), and a host class points to a **companion table** (the other bodies at its own level: the stars of a
group) and a **satellite table** (what orbits it: planets, moons), each a gated table in the format of
`possibility_table` (GEN.472): rows with the classes they support, a weight, a source and a source state, and a zone rule
in the row's notes.

This step only moves what the code already does into the tables; no drawn value changes, so it does not bump the
generator epoch. What a table does not hold yet is said by the host class (`built_by`), and is moved in by the item
that builds it (U7 of the design note).

The host classes:

* **group-node**: the members of a multi-star group. Companion table `star-count`, the number of stars (1 to 7, the tail
  to 12 is drawn by `multistar`).
* **star-like**: a star, a white dwarf or a planetary nebula's central star. Companion table `star-count` (the same
  table: a group node is the star-like object's root). Its planets, belts and comets are made by `StarSystem`.
* **remnant**: a neutron star or a black hole. Satellite table `remnant-satellites`; the other pairs are made by
  `multistar` from `MEMBER_STATE_BY_PLACED`.
* **brown-dwarf**: holds nothing yet (GEN.331 and the disk-limited planets of U7 will add rows).
* **planet**: a planet, bound or free. Satellite table `planet-satellites` (a free planet's captured companion).
* **moon**: holds nothing until the tidal-survival gate row of GEN.579 exists (Boss's cards S5 and S6, 12:10 to
  12:13 AM EDT, 2026-10-18, as the coordinator relayed, which superseded "moons are the lowest level"). Its satellite
  table `moon-satellites` has only the `no satellites` row.
* **small-body**: a comet or an asteroid field. Its `may_host` gate is closed, so it has no satellite table (its orbiter
  table is empty by the gate).

Every class carries a `may_host` gate: whether anything may orbit it at all. The recursion of the later `populate(host)`
ends by that gate and the satellite table's rows, not by a depth of three.

A row that draws an orbiter's mass states it as a ratio to the host's mass where the source gives a ratio (a moon to its
planet, a planet to its star), so a draw cannot make an orbiter heavier than its host except through the pair builder,
where the primary-first step lives (GEN.581). No row has a mass yet; the convention is for the rows U7 adds.

`SATELLITES` and `COMPANIONS` are looked up by name at call time, so a test can swap a table (`reweighted` in `tests/test_host_tables.py`).
"""

import functools
from dataclasses import dataclass

from planetgen.generation import possibility_table as pt
from planetgen.util import draw

NO_SATELLITES = "no satellites"
"""str: The name of the row every satellite table has: the object holds nothing."""


@dataclass(frozen=True)
class HostClass:
    """
    One host class.

    Attributes:
        name (str): The class (`registry.Kind.host`).
        label (str): What it is, in words.
        companions (str or None): The name of its companion table, or `None` (it has none).
        satellites (str or None): The name of its satellite table, or `None` (it holds nothing, or the table is not
            built yet: see `built_by`).
        built_by (str): What makes the parts that no table holds yet.
        may_host (bool): The gate: whether anything may orbit an object of this class at all. A closed gate means no
            satellite table.
    """

    name: str
    label: str
    companions: str | None
    satellites: str | None
    built_by: str
    may_host: bool = True


class ClassWeight:
    """
    A row's weight that depends on the subject's class: `by_class[class]` times `scale(subject)` (1 if no scale). It
    prints its numbers, so a table's digest changes with them.
    """

    def __init__(self, by_class, scale=None):
        self.by_class = dict(by_class)
        self.scale = scale

    def __call__(self, subject):
        weight = self.by_class[subject.body_class]
        return weight * self.scale(subject) if self.scale else weight

    def __repr__(self):
        return f"ClassWeight({sorted(self.by_class.items())!r}, {getattr(self.scale, '__name__', None)})"


STAR_COUNT_PROBABILITIES = {
    "M": (0.73000, 0.23500, 0.02800, 0.00560, 0.00123, 0.00015, 0.00001),
    "K": (0.60000, 0.31000, 0.06500, 0.02002, 0.00440, 0.00053, 0.00005),
    "G": (0.54500, 0.33500, 0.08000, 0.03202, 0.00705, 0.00085, 0.00008),
    "F": (0.50000, 0.34000, 0.11000, 0.04003, 0.00881, 0.00106, 0.00011),
    "A": (0.40000, 0.40000, 0.15000, 0.03741, 0.01048, 0.00189, 0.00023),
    "B": (0.25000, 0.40000, 0.22000, 0.08962, 0.03137, 0.00784, 0.00118),
    "B-mid": (0.24000, 0.36000, 0.25143, 0.10241, 0.03585, 0.00896, 0.00135),
    "B-early": (0.16000, 0.32000, 0.32686, 0.13314, 0.04660, 0.01165, 0.00175),
    "O": (0.10000, 0.25000, 0.30000, 0.21354, 0.09610, 0.03363, 0.00673),
}
"""dict: The chance of 1 to 7 stars to a system, by the placed star's spectral class (multistar-and-compact-systems.md
3.3). The B-mid and B-early rows (GEN.449) take the single and binary fractions of `tuning.MULTIPLICITY_BY_PRIMARY_MASS_BIN`
and split the rest across 3 to 7 stars in the proportions of the B row. They are the weights of the `star-count` table."""

FIVE_OR_MORE = 5
"""int: From this many stars the weight is scaled by the subject's `prevalence`."""


def _prevalence_for_five_or_more(subject):
    return subject.get("prevalence", 1.0)


def _star_count_row(count):
    column = count - 1
    by_class = {letter: row[column] for letter, row in STAR_COUNT_PROBABILITIES.items()}
    return pt.Row(
        name=f"{count} star{'s' if count > 1 else ''}", classes=tuple(by_class),
        weight=ClassWeight(by_class, _prevalence_for_five_or_more if count >= FIVE_OR_MORE else None),
        source="docs/design/systems/multistar-and-compact-systems.md 3.3; B-mid and B-early from [R3108] "
               "(ApJS 230:15) Table 13 (the B-mid and B-early rows)",
        state=pt.PARTLY_VERIFIED,
        limit="the A to M rows are the design note's and were moved here unchanged; they were not rechecked against "
              "[R3646] or [R1175] when the tables were built",
        notes="A system of 7 stars may have more (`multistar.TAIL_RATIOS`, to 12); the five-or-more columns are scaled by "
              "the subject's `prevalence`.")


STAR_COUNT = pt.REGISTRY.register(pt.Table(
    "star-count", 1, [_star_count_row(count) for count in range(1, 8)],
    description="How many stars a system has, by the placed star's spectral class (the companion table of a group "
                "node and of a star-like object)."))

MILLISECOND_PULSAR_PLANETS = 0.007
"""float: The chance a millisecond pulsar holds planets: the 0.7% of PSR B1257+12-like systems [R4750, R4751]. The magnetars, young pulsars, intermediate-mass and supermassive holes, brown dwarfs, comets and asteroid
fields hold none; white dwarfs and stars take their planets from the star model."""

ROGUE_PLANET_MOON = 0.2
"""float: The chance a generated rogue planet kept a captured companion (a design default; GEN.466 will derive it from
ejection survival)."""


def _no_satellites_weight(table_name):
    def weight(subject):
        others = [row for row in REGISTRY_TABLES[table_name].passing(subject) if row.name != NO_SATELLITES]
        return max(0.0, 1.0 - sum(row.evaluate(subject) for row in others))
    weight.__qualname__ = f"no_satellites_weight[{table_name}]"
    return weight


REMNANT_SATELLITES = pt.REGISTRY.register(pt.Table(
    "remnant-satellites", 1, [
        pt.Row("planets", ("neutron-star:millisecond",), MILLISECOND_PULSAR_PLANETS,
               source="[R4751] (Nature 355:145); [R4750] (Acta Astronomica 62:3)", state=pt.PARTLY_VERIFIED,
               limit="the 0.7% is the share of known millisecond pulsars with planets; a recalled figure, not "
                     "rechecked against the papers when the tables were built",
               notes="Zone rule: the planet model's orbit ceiling (the object's galactic Hill sphere), "
                     "`StarSystem(compact_remnant=...)`; a hit whose model puts no planet in the system stays bare."),
        pt.Row(NO_SATELLITES, (pt.ANY_CLASS,), _no_satellites_weight("remnant-satellites"),
               source="the rest of the row above (1 minus the chance of planets)", state=pt.PARTLY_VERIFIED, limit="derived from the planets row")],
    description="Whether a compact remnant holds planets, by the object key `kind` or `kind:subtype` "
                "(`system_pass.object_key`)."))

PLANET_SATELLITES = pt.REGISTRY.register(pt.Table(
    "planet-satellites", 1, [
        pt.Row("captured moon", ("rogue-planet",), ROGUE_PLANET_MOON,
               source="design default (docs/design/systems/phenomena-catalogue-and-fill.md section 5)", state=pt.STAND_IN, replacement="a captured-moon rate derived from ejection survival (docs/design/architecture/unified-star-like-objects.md 3.5)",
               notes="Zone rule: Roche limit to 0.4895 of the Hill radius, prograde; for a free planet its history "
                     "sets the limit (derived from ejection survival). Today the moon is a flag (`has_moons`)."),
        pt.Row(NO_SATELLITES, (pt.ANY_CLASS,), _no_satellites_weight("planet-satellites"),
               source="the rest of the row above (1 minus the chance of a moon)", state=pt.STAND_IN, replacement="a captured-moon rate derived from ejection survival (docs/design/architecture/unified-star-like-objects.md 3.5)")],
    description="Whether a planet holds a moon, by the planet's object key."))

MOON_SATELLITES = pt.REGISTRY.register(pt.Table(
    "moon-satellites", 1, [
        pt.Row(NO_SATELLITES, (pt.ANY_CLASS,), 1.0,
               source="no row yet: the tidal-survival gate row of the System design (docs/design/architecture/unified-system-barycentre.md) will add the moon of a moon", state=pt.STAND_IN, replacement="the tidal-survival gate row (docs/design/architecture/unified-system-barycentre.md)",
               notes="A moon may host through the same gate as any body (Boss's cards S5 and S6); until a sourced "
                     "row passes it, no moon holds anything.")],
    description="Whether a moon holds anything: empty until the tidal-survival row is added."))

REGISTRY_TABLES = {table.name: table for table in (STAR_COUNT, REMNANT_SATELLITES, PLANET_SATELLITES, MOON_SATELLITES)}
"""dict: The host tables by name; a test swaps one with `monkeypatch.setitem` (`reweighted` in `tests/test_host_tables.py`)."""

HOSTS = (
    HostClass("group-node", "the stars of a multi-star group", "star-count", None,
              "multistar.draw_member_tree (the pair rule, periods and stability screens)"),
    HostClass("star-like", "a star, a white dwarf or a central star", "star-count", None,
              "StarSystem (planets, belts and comets by the star model; the pair rule of multistar)"),
    HostClass("remnant", "a neutron star or a black hole", None, "remnant-satellites",
              "multistar.MEMBER_STATE_BY_PLACED (its stellar companions, drawn at scatter time)"),
    HostClass("brown-dwarf", "a brown dwarf", None, None, "nothing yet: holds no companion or satellite"),
    HostClass("planet", "a planet, bound or free", None, "planet-satellites",
              "RoguePlanet and StarSystem (a bound planet's moons: generate_moons)"),
    HostClass("moon", "a moon", None, "moon-satellites",
              "nothing yet: the tidal-survival gate row adds what a moon holds"),
    HostClass("small-body", "a comet or an asteroid field", None, None, "nothing: the gate is closed",
              may_host=False),
)
"""tuple: The host classes, in the order of the design note's table."""

BY_NAME = {host.name: host for host in HOSTS}
"""dict: `name` -> its `HostClass`."""


def _table(name):
    return REGISTRY_TABLES[name]


@functools.lru_cache(maxsize=None)
def _probabilities(table, body_class, prevalence):
    return table.probabilities(pt.Subject(body_class, {"prevalence": prevalence}))


def draw_star_count(body_class, prevalence=1.0, rng=draw):
    """
    The number of stars, 1 to 7, of a system with a star of spectral class `body_class` (a key of
    `STAR_COUNT_PROBABILITIES`): one draw of the `star-count` table on `rng`, the five-or-more columns scaled by
    `prevalence`.
    """
    table = _probabilities(_table("star-count"), body_class, prevalence)
    return int(rng.choices([row for row, _p in table], weights=[p for _row, p in table])[0].name.split()[0])


def star_count_weights(body_class, prevalence=1.0):
    """The weights of 1 to 7 stars (the table's, before they are normalised), the five-or-more columns scaled by
    `prevalence`."""
    subject = pt.Subject(body_class, {"prevalence": prevalence})
    return tuple(row.evaluate(subject) for row in _table("star-count").rows)


def holds(table_name, object_key):
    """Whether the satellite table `table_name` lets an object with this key hold anything at all: some row besides
    `NO_SATELLITES` passes its gate."""
    table = _table(table_name)
    return any(row.name != NO_SATELLITES for row in table.passing(pt.Subject(object_key)))


def chance(table_name, object_key):
    """The chance that an object with this key holds something (the chance the row that is not `NO_SATELLITES` is
    drawn), 0 for a key whose table gives it nothing."""
    table = _table(table_name)
    subject = pt.Subject(object_key)
    return sum(p for row, p in table.probabilities(subject) if row.name != NO_SATELLITES)


def holds_something(table_name, object_key, rng=draw):
    """Rolls the table for an object with this key on `rng` (one draw): whether it holds something."""
    table = _table(table_name)
    rows_and_chances = table.probabilities(pt.Subject(object_key))
    chosen = rng.choices([row for row, _p in rows_and_chances], weights=[p for _row, p in rows_and_chances])[0]
    return chosen.name != NO_SATELLITES
