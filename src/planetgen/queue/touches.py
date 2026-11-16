# planetgen/queue/touches.py

"""
What each command touches, and who may overlap whom (OPS.241, the independence harness)
=========================================================================================

A command that does one thing declares the tables it **reads** and the tables it **writes** (`READS`, `WRITES`, keyed by command kind, and by the name of an end command, `conclude-<kind>` and
`finish-<kind>`). `needs.py` says what a command needs made before it can run; this says what it can disturb while it runs. Two commands **conflict** when one writes a table the other reads or
writes (`conflict`); a conflict is an edge the planner must keep, and two commands that do not conflict may run together only after a test has proven the pair (`independence.json`, `proven_pair`).
A pair that is not in the file runs in plan order, so a new command is safe by default and overlap is earned one pair at a time (OPS.242).

**The declarations are checked against reality.** `tests/test_command_touches.py` runs a toy galaxy and the commands around it under the SQL tap (`sqltap.Tap`), and fails a command that touches a table
it did not declare, and a registered command that has no declaration. A command the run does not exercise is listed in the test's `UNTAPPED` with the reason, and its declaration is read from its code,
made wide where in doubt (`ANY`, every table). Where a command touches a table only some of the time (the epoch stamps), it is declared all the same: a declaration is what the command may touch.

**What is not in the tables.** The job ledger and the run record (`ENGINE_TABLES`) are written by the engine for every command and are left out of the declarations. A write that is not to a table at all
(a file, a Redis key, a cache of the process, a lock) breaks independence the same way and is declared in `OTHER`: the files and caches the tap does not see are named from the audit of the design note
(`docs/design/operations/job-graph-and-memory-broker.md`, section 6) and are not checked by the tap. Row ids taken from `id_counters` and `id_blocks` in the order the commands finish are real conflicts
and are declared as the tables they are.

**Grain.** The default grain is the table. A finer grain (a shell, a sector range) is allowed only for a pair proven at that grain (`independence.json`, its `grain`).
"""

import json
import os

ANY = "*"
"""str: A declaration that names every table: a command whose reads (or writes) are not listed table by table (a warm-up answers through the site's own pages; a reset wipes every content table)."""

ENGINE_TABLES = frozenset({"engine_runs", "engine_jobs", "engine_attempts", "engine_items", "engine_piece_results", "generation_runs", "generation_run_arguments"})
"""frozenset: The job ledger and the run record, which the engine writes for every command; they are never a conflict between two commands and are left out of `READS` and `WRITES`."""

SECTOR_READS = (
    "asteroid_belts", "asteroid_fields", "black_holes", "bright_stars", "chart_ledger", "comets", "facilities", "galaxy_block_light_parts", "galaxy_components",
    "galaxy_shape", "galaxy_shape_pins", "galaxy_shell", "id_counters", "interstellar_comets", "moon_deposit_types", "moon_gases", "moon_interior_layers",
    "moon_isotopes", "moon_life", "moon_life_class_scores", "moon_reflection_spectrum", "moons", "nebulae", "neutron_stars", "orbit_simulation_state",
    "phenomenon_scatter",
    "planet_life", "planets", "quasars", "rogue_planet_deposit_types", "rogue_planet_gases", "rogue_planet_interior_layers", "rogue_planet_isotopes",
    "rogue_planet_life_class_scores", "rogue_planet_reflection_spectrum", "rogue_planets", "scatter_members", "scatter_nodes",
    "sector_facets", "sector_name_registry", "sector_stats", "sectors", "sky_tier_draws", "species", "star_systems", "stars", "supernova_remnants", "system_configs",
    "system_interests", "system_name_registry",
)
SECTOR_WRITES = (
    "asteroid_belt_composition", "asteroid_belts", "belt_members", "body_history", "bright_stars", "chart_ledger", "comet_composition", "comets", "dark_members",
    "effect_ledger", "galaxy_block_light_parts", "galaxy_shape", "galaxy_totals", "id_blocks", "id_counters", "interstellar_comet_composition", "interstellar_comets", "irregular_moons", "irregular_sets", "moon_deposit_types",
    "moon_gases", "moon_isotopes", "moon_life", "moon_life_class_scores", "moons", "orbit_nodes", "planet_deposit_types", "planet_gases", "planet_interior_layers",
    "planet_isotopes", "planet_life", "planet_life_class_scores", "planets", "rogue_planet_deposit_types", "rogue_planet_gases", "rogue_planet_isotopes",
    "rogue_planets", "scatter_members", "scatter_nodes", "sector_facets", "sector_name_registry", "sector_stats", "sector_tree_nodes",
    "sectors", "star_systems", "stars", "system_configs", "system_interests", "system_name_registry", "system_summary",
)

SCATTER_READS = (
    "bright_stars", "chart_ledger", "galaxy_block_light_parts", "galaxy_components", "galaxy_shape", "galaxy_shape_pins", "galaxy_shell", "id_counters", "orbit_simulation_state",
    "phenomenon_scatter",
)
SCATTER_STAR_WRITES = ("bright_stars", "galaxy_block_light_parts", "id_counters", "scatter_members", "scatter_nodes")
SCATTER_PHENOMENON_WRITES = ("galaxy_block_light_parts", "id_counters", "phenomenon_scatter", "scatter_members", "scatter_nodes")
"""tuple: What a scatter stage's units read and write, from the tap on a toy galaxy dense enough that each draws rows (`tests/pair_overlap_support.py`). The scatters work on a skeleton and write
only their own scatter table, the light parts of the blocks, the members of what they drew and the ids of the rows; the rest is the stage's end commands'."""

READS = {
    "check-counts": (
        "asteroid_belts", "galaxy_totals", "moon_life", "moons", "planet_life", "planets", "sector_facets", "sector_stats", "sectors", "star_systems", "stars",
    ),
    "check-datasets": (),
    "check-ids": (
        "asteroid_belts", "asteroid_fields", "black_holes", "comets", "id_blocks", "interstellar_comets", "moons", "nebulae", "neutron_stars", "planets", "quasars",
        "rogue_planets", "sectors", "star_systems", "stars", "supernova_remnants", "system_configs",
    ),
    "check-math": (),
    "check-object-index": (
        "asteroid_belts", "asteroid_fields", "black_holes", "comets", "facilities", "interstellar_comets", "moons", "nebulae", "neutron_stars", "object_index", "planets",
        "quasars", "rogue_planets", "star_systems", "stars", "supernova_remnants",
    ),
    "check-orphans": ("*",),
    "check-revision": ("alembic_version", "schema_migrations"),
    "check-storage": ("sectors",),
    "check-systems": (
        "asteroid_belt_composition", "asteroid_belts", "belt_members", "black_holes", "body_history", "comet_composition", "comets", "dark_members", "irregular_moons",
        "irregular_sets", "moon_deposit_types", "moon_gases", "moon_interior_layers", "moon_isotopes", "moon_life", "moon_life_class_scores", "moon_reflection_spectrum", "moons",
        "neutron_stars", "orbit_nodes", "planet_deposit_types", "planet_gases", "planet_interior_layers", "planet_isotopes", "planet_life", "planet_life_class_scores",
        "planet_reflection_spectrum", "planets", "star_systems", "stars", "system_config_slots", "system_configs",
    ),
    "check-tables": (),
    "check-tables-extended": (),
    "check-values": ("moons", "planets", "sectors", "star_systems", "stars"),
    "check-version-keys": ("sectors",),
    "clean-nebula-field": ("galaxy_components", "galaxy_shape", "galaxy_shape_pins", "nebulae"),
    "clear-bright-stars": ("bright_stars", "galaxy_components", "galaxy_shape", "galaxy_shape_pins", "scatter_members", "scatter_nodes"),
    "clear-scatter": (
        "bright_stars", "galaxy_block_light_parts", "galaxy_block_lights", "galaxy_components", "galaxy_shape", "galaxy_shape_pins", "nebula_field", "phenomenon_scatter",
        "phenomenon_scatter_classes", "scatter_members", "scatter_nodes",
    ),
    "conclude-draw-nebula-field": ("galaxy_components", "galaxy_shape", "galaxy_shape_pins", "nebula_field", "nebulae"),
    "conclude-rescore-habitability": (
        "asteroid_belts", "asteroid_fields", "black_holes", "chart_ledger", "comets", "galaxy_shape", "galaxy_shape_pins", "interstellar_comets", "moons",
        "neutron_stars", "planets", "quasars", "rogue_planets", "sector_stats", "sectors", "star_systems", "stars",
    ),
    "conclude-scatter-backfill": (),
    "conclude-scatter-luminosity": (),
    "conclude-scatter-mass": (),
    "conclude-scatter-phenomena": (),
    "conclude-warm-cluster-field": ("*",),
    "draw-nebula-field": ("galaxy_components", "galaxy_shape", "galaxy_shape_pins", "galaxy_shell"),
    "drop-phenomenon-indexes": (),
    "empty-nebula-field": (),
    "empty-phenomenon-scatter": (),
    "end-map-rebuild-maintenance": (),
    "fill-map-luminosity": ("phenomenon_scatter",),
    "build-phenomenon-indexes": (
        "bright_stars", "chart_ledger", "galaxy_block_light_parts", "galaxy_components", "galaxy_shape", "galaxy_shape_pins", "galaxy_shell", "id_counters",
        "orbit_simulation_state", "phenomenon_scatter",
    ),
    "record-phenomenon-scatter": (
        "bright_stars", "chart_ledger", "galaxy_block_light_parts", "galaxy_components", "galaxy_shape", "galaxy_shape_pins", "galaxy_shell", "id_counters",
        "orbit_simulation_state", "phenomenon_scatter",
    ),
    "write-special-phenomena": (
        "bright_stars", "chart_ledger", "galaxy_block_light_parts", "galaxy_components", "galaxy_shape", "galaxy_shape_pins", "galaxy_shell", "id_counters",
        "orbit_simulation_state", "phenomenon_scatter",
    ),
    "finish-warm-galaxy-deep-nodes": ("*",),
    "finish-warm-galaxy-nodes": ("*",),
    "finish-warm-phenomena": ("*",),
    "finish-warm-sectors": ("*",),
    "finish-warm-systems": ("*",),
    "generate-neighborhood": SECTOR_READS,
    "generate-sectors": SECTOR_READS,
    "generate-start-sector": SECTOR_READS,
    "generate-system": SECTOR_READS,
    "generate-unplaced-sector": SECTOR_READS,
    "link-neighbors": (
        "asteroid_fields", "black_holes", "galaxy_block_stats", "interstellar_comets", "nearest_systems", "nebulae", "neutron_stars", "quasars", "rogue_planets",
        "sectors", "star_systems", "supernova_remnants",
    ),
    "plan-skeleton": (),
    "refill-effect-ledger": ("star_systems", "system_summary"),
    "refresh-block-lights": ("bright_stars", "galaxy_block_light_parts", "galaxy_block_lights", "phenomenon_scatter", "phenomenon_scatter_classes"),
    "refresh-run-blocks": (
        "asteroid_belts", "asteroid_fields", "black_holes", "bright_stars", "chart_ledger", "comets", "galaxy_block_light_parts", "galaxy_block_lights",
        "galaxy_block_stats", "galaxy_shape", "galaxy_shape_pins", "interstellar_comets", "moons", "nebula_block_state", "nebulae", "neutron_stars", "phenomenon_scatter",
        "planets", "quasars", "rogue_planets", "sector_stats", "sectors", "star_systems", "stars",
    ),
    "rescore-habitability": ("star_systems", "system_summary"),
    "reset-galaxy": ("*",),
    "run-population": ("moons", "planet_life", "planets", "polities", "population_state", "species", "star_systems", "system_configs"),
    "run-territories": ("polities", "sectors", "star_systems"),
    "scatter-backfill": (
        "bright_stars", "chart_ledger", "galaxy_block_light_parts", "galaxy_components", "galaxy_shape", "galaxy_shape_pins", "galaxy_shell", "id_counters",
        "orbit_simulation_state", "phenomenon_scatter", "sector_stats", "sectors", "sky_tier_draws",
    ),
    "scatter-band": SCATTER_READS,
    "scatter-luminosity": SCATTER_READS,
    "scatter-mass": SCATTER_READS,
    "scatter-phenomena": SCATTER_READS,
    "settle-sector-paths": (
        "black_holes", "interstellar_comets", "neutron_stars", "rogue_planets", "sector_paths", "sectors", "star_systems", "system_summary",
    ),
    "store-scatter-limits": (),
    "warm-cluster-field": ("*",),
    "warm-galaxy-deep-nodes": ("*",),
    "warm-galaxy-map": ("*",),
    "warm-galaxy-nodes": ("*",),
    "warm-lists": ("*",),
    "warm-object-index": ("*",),
    "warm-phenomena": ("*",),
    "warm-search": ("*",),
    "warm-sectors": ("*",),
    "warm-systems": ("*",),
}

WRITES = {
    "check-counts": (),
    "check-datasets": (),
    "check-ids": (),
    "check-math": (),
    "check-object-index": (),
    "check-orphans": (),
    "check-revision": (),
    "check-storage": (),
    "check-systems": (),
    "check-tables": (),
    "check-tables-extended": (),
    "check-values": (),
    "check-version-keys": (),
    "clean-nebula-field": (),
    "clear-bright-stars": (
        "bright_stars", "galaxy_block_light_parts", "galaxy_shape", "id_counters", "scatter_epochs", "scatter_members", "scatter_nodes", "sector_stats", "sky_tier_draws",
    ),
    "clear-scatter": (
        "bright_stars", "galaxy_block_light_parts", "galaxy_block_lights", "galaxy_shape", "id_counters", "nebula_field", "phenomenon_scatter",
        "phenomenon_scatter_classes", "scatter_epochs", "scatter_members", "scatter_nodes", "sector_stats", "sky_tier_draws",
    ),
    "conclude-draw-nebula-field": (
        "generator_epoch", "nebula_block_classes", "nebula_block_state", "nebula_block_stats", "nebula_field", "nebulae_through_id", "scatter_epochs",
    ),
    "conclude-rescore-habitability": ("galaxy_block_stats", "region_stats", "sector_stats"),
    "conclude-scatter-backfill": ("galaxy_shape", "generator_epoch", "scatter_epochs"),
    "conclude-scatter-luminosity": ("galaxy_shape", "generator_epoch", "scatter_epochs"),
    "conclude-scatter-mass": ("galaxy_shape", "generator_epoch", "scatter_epochs"),
    "conclude-scatter-phenomena": ("phenomenon_scatter_classes",),
    "conclude-warm-cluster-field": (),
    "draw-nebula-field": ("nebula_field",),
    "drop-phenomenon-indexes": ("phenomenon_scatter",),
    "empty-nebula-field": ("nebula_field",),
    "empty-phenomenon-scatter": ("phenomenon_scatter", "phenomenon_scatter_classes"),
    "end-map-rebuild-maintenance": (),
    "fill-map-luminosity": ("phenomenon_scatter",),
    "build-phenomenon-indexes": (
        "galaxy_block_light_parts", "galaxy_shape", "generator_epoch", "id_counters", "phenomenon_scatter", "phenomenon_scatter_classes", "scatter_epochs",
    ),
    "record-phenomenon-scatter": (
        "galaxy_block_light_parts", "galaxy_shape", "generator_epoch", "id_counters", "phenomenon_scatter", "phenomenon_scatter_classes", "scatter_epochs",
    ),
    "write-special-phenomena": (
        "galaxy_block_light_parts", "galaxy_shape", "generator_epoch", "id_counters", "phenomenon_scatter", "phenomenon_scatter_classes", "scatter_epochs",
    ),
    "finish-warm-galaxy-deep-nodes": (),
    "finish-warm-galaxy-nodes": (),
    "finish-warm-phenomena": (),
    "finish-warm-sectors": (),
    "finish-warm-systems": (),
    "generate-neighborhood": SECTOR_WRITES,
    "generate-sectors": SECTOR_WRITES,
    "generate-start-sector": SECTOR_WRITES,
    "generate-system": SECTOR_WRITES,
    "generate-unplaced-sector": SECTOR_WRITES,
    "link-neighbors": ("nearest_systems", "rogue_planets"),
    "plan-skeleton": ("galaxy_components", "galaxy_shape", "galaxy_shape_pins", "galaxy_shell"),
    "refill-effect-ledger": (),
    "refresh-block-lights": ("galaxy_block_lights",),
    "refresh-run-blocks": ("galaxy_block_light_parts", "galaxy_block_lights", "galaxy_block_stats", "region_stats", "sector_stats"),
    "rescore-habitability": (),
    "reset-galaxy": ("*",),
    "run-population": ("population_state", "species", "system_interests"),
    "run-territories": ("system_owners",),
    "scatter-backfill": (
        "bright_stars", "galaxy_block_light_parts", "id_counters", "scatter_members", "scatter_nodes", "sector_stats",
    ),
    "scatter-band": SCATTER_STAR_WRITES,
    "scatter-luminosity": SCATTER_STAR_WRITES,
    "scatter-mass": SCATTER_STAR_WRITES,
    "scatter-phenomena": SCATTER_PHENOMENON_WRITES,
    "settle-sector-paths": ("sector_path_knots", "sector_paths"),
    "store-scatter-limits": ("galaxy_shape",),
    "warm-cluster-field": (),
    "warm-galaxy-deep-nodes": (),
    "warm-galaxy-map": (),
    "warm-galaxy-nodes": (),
    "warm-lists": (),
    "warm-object-index": (),
    "warm-phenomena": (),
    "warm-search": (),
    "warm-sectors": (),
    "warm-systems": (),
}


OTHER = {
    "scatter-backfill": ("file:cluster-field",),
    "generate-start-sector": ("file:cluster-field",),
    "generate-neighborhood": ("file:cluster-field",),
    "generate-sectors": ("file:cluster-field",),
    "warm-cluster-field": ("file:cluster-field",),
    "conclude-warm-cluster-field": ("file:cluster-field",),
}
"""dict: The writes and reads that are not to a table, by command kind: a file, a Redis key or a cache that two commands could share (design note section 6). Named, not tapped."""

INDEPENDENCE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "independence.json")
"""str: The pairs of commands a test has run overlapped, in both orders and at 1 and 15 workers, with the same content hashes as the plan-order run (OPS.242 adds a pair with its proof)."""

PAIR_KEYS = ("a", "b", "grain", "shared", "proof")


def declared(kind):
    """Whether `kind` has a reads row and a writes row."""
    return kind in READS and kind in WRITES


def reads_of(kind):
    """The tables `kind` reads (`ANY` for all); `()` for a command with no declaration."""
    return READS.get(kind, ())


def writes_of(kind):
    """The tables `kind` writes (`ANY` for all)."""
    return WRITES.get(kind, ())


def _hits(first, second):
    return ANY in first or ANY in second or bool(set(first) & set(second))


def conflicting_names(kind_a, kind_b):
    """
    What two commands share: the tables one writes that the other reads or writes, and the files and caches both name in `OTHER`, sorted. `ANY` for a command that touches every table (the
    answer is then `["*"]`); `["?"]` for a command with no declaration (nothing is known of it). `ENGINE_TABLES` do not count.
    """
    if not (declared(kind_a) and declared(kind_b)):
        return ["?"]
    reads_a, writes_a, reads_b, writes_b = reads_of(kind_a), writes_of(kind_a), reads_of(kind_b), writes_of(kind_b)
    if (ANY in writes_a and (writes_b or reads_b)) or (ANY in writes_b and (writes_a or reads_a)) or (ANY in reads_a and writes_b) or (ANY in reads_b and writes_a):
        return [ANY]
    found = (set(writes_a) & (set(reads_b) | set(writes_b))) | (set(writes_b) & set(reads_a)) | (set(OTHER.get(kind_a, ())) & set(OTHER.get(kind_b, ())))
    return sorted(found)


def conflict(kind_a, kind_b):
    """Whether two commands conflict (`conflicting_names` names something): one writes what the other reads or writes, or both touch the same file or cache named in `OTHER`."""
    return bool(conflicting_names(kind_a, kind_b))


def group_of(kind):
    """A stage's command and the end commands it has (`conclude-<kind>`, `finish-<kind>`), which run when the stage's units are done and are part of what a pair overlaps."""
    from planetgen.queue import commands, planrows
    return [kind] + [name for name in (commands.end_name(kind, role) for role in (planrows.CONCLUDE, planrows.FINISH)) if name in READS and name != kind]


def shared_between(kind_a, kind_b):
    """What the stages `kind_a` and `kind_b` share, their end commands counted (`conflicting_names` over `group_of` of each), sorted; empty for two that share nothing."""
    found = set()
    for first in group_of(kind_a):
        for second in group_of(kind_b):
            found |= set(conflicting_names(first, second))
    return sorted(found)


def load_pairs(path=None):
    """
    The pairs in `independence.json` as `{frozenset({a, b}): entry}`.

    Raises:
        ValueError: For a file that is not the format below, a pair that names itself or an unregistered kind, or a pair listed twice.
    """
    from planetgen.queue import commands
    with open(path or INDEPENDENCE_FILE, encoding="utf-8") as handle:
        data = json.load(handle)
    if data.get("format") != 1 or not isinstance(data.get("pairs"), list):
        raise ValueError("independence.json: the format is {\"format\": 1, \"pairs\": [...]}")
    found = {}
    for entry in data["pairs"]:
        missing = [key for key in PAIR_KEYS if key not in entry]
        if missing:
            raise ValueError(f"independence.json: a pair lacks {', '.join(missing)}: {entry}")
        pair = frozenset((entry["a"], entry["b"]))
        if len(pair) != 2:
            raise ValueError(f"independence.json: {entry['a']!r} cannot be independent of itself")
        for kind in pair:
            if kind not in commands.REGISTRY and kind not in {commands.end_name(k, role) for role, k in commands.END_COMMANDS}:
                raise ValueError(f"independence.json: {kind!r} is not a registered command")
        if pair in found:
            raise ValueError(f"independence.json: the pair {sorted(pair)} is listed twice")
        found[pair] = entry
    return found


def proven_pair(kind_a, kind_b, pairs=None):
    """Whether `kind_a` and `kind_b` are listed in `independence.json`: the only pairs that may overlap."""
    return frozenset((kind_a, kind_b)) in (load_pairs() if pairs is None else pairs)


def may_overlap(kind_a, kind_b, pairs=None):
    """Whether two commands may run together: they are a proven pair (`proven_pair`). A conflict does not forbid a pair proven at a finer grain, but a pair that is not listed never overlaps."""
    return proven_pair(kind_a, kind_b, pairs)


def must_wait_for(kinds, index, pairs=None):
    """
    The positions before `index` in a plan (`kinds`, the command kind of each entry in plan order) that entry `index` must wait for: every earlier entry that it may not overlap with (`may_overlap`),
    so a pair that is not listed in `independence.json` keeps its plan order.
    """
    pairs = load_pairs() if pairs is None else pairs
    return [place for place in range(index) if not may_overlap(kinds[place], kinds[index], pairs)]


def plan_edges(kinds, pairs=None):
    """
    The edges of a plan: for each entry, the earlier entries it waits for, without the ones it already waits for through another (the transitive reduction of `must_wait_for`). With no pair listed in
    `independence.json` each entry waits for the one before it only, which is the chain the planner wrote before the harness; a proven pair adds no edge between its two commands and the edges the others need.

    Returns:
        list: One list of positions for each entry of `kinds`.
    """
    pairs = load_pairs() if pairs is None else pairs
    reach = []
    edges = []
    for index in range(len(kinds)):
        waits = must_wait_for(kinds, index, pairs)
        covered = set()
        for place in waits:
            covered |= reach[place]
        edges.append([place for place in waits if place not in covered])
        reach.append(set(waits) | covered)
    return edges
