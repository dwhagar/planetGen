# planetgen/generation/system_pass.py

"""
The fill's system pass (GEN.251, stage C)
=========================================

A scatter row is a point, a seed and a kind: no planet or moon exists around it until its sector is filled
(docs/design/phenomena-catalogue-and-fill.md rules 1, 3 and 4). When a fill has built its scattered objects (stage B),
and the neutron stars and black holes below the scatter's mass cut (stage G), this pass gives each one that can hold
planets its roll. The stellar companions of those objects are not drawn here: they are stored with the row at scatter
time (`scatter_companions`) and built from the stored members by the fill (GEN.264).

* **What can hold what** is one dict, `tuning.SYSTEM_PASS_ODDS`, by `object_key`: the kind, and for a neutron star
  its pulsar type, for a black hole its mass class. Only a millisecond pulsar holds planets (0.7%); the other
  compact remnants and the nuclei hold none.
* **The roll is the object's own.** It is drawn from `draw.Stream(f"{seed}:system-pass")`, `seed` the scatter row's
  seed (a below-cut remnant's: its sector's below-cut seed and its place in the draw), and the system is built on
  that same stream. So the answer is the same whenever and in whatever order the sector is filled, and no other
  stage's draws move.
* **Order.** Compact remnants first, then the rest (`registry` family), each group in the order it was made.
* **A hit replaces the bare object with a system** around it, at its point (`SpaceSector.replace_phenomenon_with_system`):
  the system's galactic Hill sphere is its primary's, the sphere the object was placed with, and the planet model
  already keeps every body inside its orbit ceiling, so every system fits and none needs its outermost body dropped.
  The object itself always stays. A roll that hits but whose planet model puts no planet in the system keeps the bare
  object (counted as `:empty`).

`run` returns the systems it made; the counts per object key are in `sector.stage_counts` under the stage
`fill_stages.SYSTEM_PASS`: how many rolled (`key`), how many became systems (`key:system`) and how many were empty
(`key:empty`).
"""

from planetgen import tuning
from planetgen.generation import fill_stages
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena import registry
from planetgen.generation.system import StarSystem
from planetgen.util import draw

STREAM = "system-pass"
"""str: The suffix of an object's seed that names the stream its roll and its system are drawn from."""

_FAMILY_ORDER = {"compact-remnant": 0}


def object_key(entry):
    """
    The key of a phenomenon entry in `tuning.SYSTEM_PASS_ODDS`: its kind, a neutron star's with its pulsar type
    (`neutron-star:millisecond`), a black hole's with its mass class (`black-hole:stellar`).
    """
    kind = entry.phenomenon_type
    phenomenon = entry.phenomenon
    if kind == "neutron-star":
        return f"{kind}:{getattr(phenomenon, 'pulsar_type', None)}"
    if kind == "black-hole":
        return f"{kind}:{getattr(phenomenon, 'mass_class', None)}"
    return kind


def hit(key, seed):
    """Whether the object with `seed` and `key` rolls planets: the first draw of its own stream against the odds
    (`False` for a key without odds)."""
    odds = tuning.SYSTEM_PASS_ODDS.get(key)
    return bool(odds) and draw.Stream(f"{seed}:{STREAM}").random() < odds.get("planets", 0.0)


def _rank(entry):
    row = registry.BY_KIND.get(entry.phenomenon_type)
    return _FAMILY_ORDER.get(row.family if row is not None else None, 1)


def run(sector, args, entries, galactic_center_dist_ly=None):
    """
    Rolls the system of each phenomenon entry in `entries` that `tuning.SYSTEM_PASS_ODDS` gives planets to, and
    replaces a hit with the system.

    Args:
        sector (SpaceSector): The sector the entries are in.
        args (argparse.Namespace): Parsed arguments; only `args.markdown` is read.
        entries (iterable): `SectorPhenomenonEntry` instances with a `system_seed` (the number their roll is drawn
            from: a scatter row's seed, or a string for an object made at fill).
        galactic_center_dist_ly (float, optional): Threaded to the system, as `generate_sector` does for any.

    Returns:
        list: The `SectorSystemEntry` instances made, in the order they were rolled.
    """
    made = []
    for entry in sorted(entries, key=_rank):
        key = object_key(entry)
        odds = tuning.SYSTEM_PASS_ODDS.get(key)
        seed = getattr(entry, "system_seed", None)
        if not odds or seed is None:
            continue
        fill_stages.count(sector, fill_stages.SYSTEM_PASS, key)
        rng = draw.Stream(f"{seed}:{STREAM}")
        if rng.random() >= odds.get("planets", 0.0):
            continue
        config = SystemConfig()
        config.MARKDOWN = args.markdown
        with draw.bound(rng):
            system = StarSystem(system_config=config, galactic_center_dist_ly=galactic_center_dist_ly,
                                compact_remnant=entry.phenomenon)
        if not system.planets:
            fill_stages.count(sector, fill_stages.SYSTEM_PASS, f"{key}:empty")
            continue
        new_entry = sector.replace_phenomenon_with_system(entry, system, system_config=config)
        scatter_id = getattr(entry, "scatter_id", None)
        if scatter_id is not None:
            new_entry.scatter_id = scatter_id
        fill_stages.count(sector, fill_stages.SYSTEM_PASS, f"{key}:system")
        made.append(new_entry)
    return made
