# planetgen/queue/commands.py

"""
Commands, and the step list the wrapper builds from a request (PERF.118, docs/design/operations/redis-and-the-work-queue.md
5.12 and docs/design/architecture/command-inventory.md).

**A command is a kind in `REGISTRY`.** It does one thing and nothing else: it does not know what ran before it, what
runs after it or where it sits in the order. A command with a `unit` is queued as a stage of units; one with an
`expand` is cut into those units when its stage starts, from its `params` (plain JSON) and the database alone. A
command with neither is *inline*: it runs in the process that started the run. PERF.106 to PERF.108 and PERF.121 to
PERF.128 turned the stages into units one at a time, and none is inline now; the inline path stays for a command
that has no units to give. Every command the wrapper can name is registered here, so a list can be printed and
checked without running anything.

**A step is a command and its parameters.** `build_steps` is the wrapper's pure function: the arguments of
`planetgen galaxy` or `planetgen plan` in, the ordered `Step`s that will run out. A stage the request does not
select is not in the list (it is never "skipped"), and the same arguments give the same list every time. The order
lives only in that list, which is data; no command reads it.

**A unit is a chunk of a command.** `run_alone` runs one unit from its JSON arguments with nothing else around it
(`planetgen command run`), for a test or a rerun.
"""

import importlib
from dataclasses import dataclass, field


@dataclass(frozen=True)
class CommandSpec:
    """
    A registered command. `unit` and `expand` are `(module, function)` pairs (what Redis holds can pick a command
    but never names a function): `unit(args, report)` does one unit and returns a JSON-able result; `expand(params)`
    returns the `Expansion` of the command's stage; `finish(params, plan, results, seconds)`, when given, is called
    once by the unit that completes the stage (PERF.106) with the step's parameters, the stage's plan, every unit's
    result in unit order and the seconds the stage took, records what the command's speed statistics need and returns
    the metrics the stage's row is stored with. `warm(unit_args)`, when given, is called by the worker before it forks
    a horse for a unit of this command (`units.warm_unit`, PERF.146): it works out in the worker what every unit
    would otherwise work out for itself. A command with no `unit` is inline.

    The work a stage does that is not a unit (PERF.172) is a hook of its command, run by whichever process starts or
    ends the stage, and a failure in one fails the run: `prepare(params)` is called once when the stage opens, before
    it is cut (empties the rows the stage is about to rewrite, say), and `conclude(params, plan, results, seconds)`
    once when its last unit is done, before the next stage is cut (stores what the stage drew as a whole); it returns a
    JSON-able dict the process following the run reads (`wrapper.run`'s `on_stage_done`). The unit that completes
    the stage calls `conclude`; a cut that gives no unit at all is concluded by the process that cut it
    (`wrapper._conclude_empty_stage`).

    The unit-length ceiling (PERF.155, `queue.unit_ceiling_seconds`): `whole` is set for a command whose stage is one
    unit by nature, which the ceiling does not hold; any other queued command is cut into units and declares in
    `unit_seconds` how long its longest unit runs at the default galaxy size, which a test holds to the ceiling. A
    unit that runs longer than the ceiling anyway is written to the run's log.
    """

    kind: str
    summary: str
    unit: tuple = None
    expand: tuple = None
    finish: tuple = None
    warm: tuple = None
    prepare: tuple = None
    conclude: tuple = None
    whole: bool = False
    unit_seconds: float = None

    @property
    def inline(self):
        return self.unit is None


@dataclass
class Expansion:
    """
    What a command's `expand` gives: the `units` as `(args, cost)` pairs, in unit order (JSON-able arguments and the
    predicted cost the stage's bar and ETA use), and `plan`, a JSON-able summary of the stage that any process
    following the run reads the units' results by.
    """

    units: list = field(default_factory=list)
    plan: dict = field(default_factory=dict)


@dataclass(frozen=True)
class Step:
    """
    One step of a run: stage `key` (the name the bars, the stats and the Generate page know it by), the command
    `kind` that does it, its `label` and its `params` (a flat JSON dict).
    """

    key: str
    kind: str
    label: str
    params: dict = field(default_factory=dict)

    def as_dict(self):
        return {"key": self.key, "kind": self.kind, "label": self.label, "params": dict(self.params)}

    @classmethod
    def from_dict(cls, data):
        return cls(data["key"], data["kind"], data["label"], dict(data.get("params") or {}))


def _commands(*specs):
    return {spec.kind: spec for spec in specs}


REGISTRY = _commands(
    CommandSpec("scatter-mass", "Scatter the massive stars: one unit per run of shells",
                unit=("planetgen.generation.run_plan", "scatter_shells_unit"),
                expand=("planetgen.generation.run_plan", "expand_scatter_mass"),
                finish=("planetgen.generation.run_plan", "finish_scatter"),
                prepare=("planetgen.generation.run_plan", "prepare_scatter_pass"),
                conclude=("planetgen.generation.run_plan", "conclude_scatter_pass"),
                warm=("planetgen.generation.run_plan", "scatter_shells_warm"),
                unit_seconds=10),
    CommandSpec("plan-skeleton", "Plan the galaxy's density skeleton: one unit",
                unit=("planetgen.generation.run_plan", "plan_skeleton_unit"),
                expand=("planetgen.generation.run_plan", "expand_plan_skeleton"), whole=True),
    CommandSpec("scatter-phenomena", "Scatter the phenomena: one unit per run of shells",
                unit=("planetgen.generation.run_plan", "scatter_phenomena_unit"),
                expand=("planetgen.generation.run_plan", "expand_scatter_phenomena"),
                finish=("planetgen.generation.run_plan", "finish_scatter_phenomena"),
                prepare=("planetgen.generation.run_plan", "prepare_scatter_phenomena"),
                conclude=("planetgen.generation.run_plan", "conclude_scatter_phenomena"),
                warm=("planetgen.queue.units", "warm_shape"),
                unit_seconds=6),
    CommandSpec("finish-phenomena", "End the phenomena scatter: the special phenomena, the indexes and the record, one unit",
                unit=("planetgen.generation.run_plan", "finish_phenomena_unit"),
                expand=("planetgen.generation.run_plan", "expand_finish_phenomena"), whole=True),
    CommandSpec("scatter-luminosity", "Scatter the bright stars: one unit per run of shells",
                unit=("planetgen.generation.run_plan", "scatter_shells_unit"),
                expand=("planetgen.generation.run_plan", "expand_scatter_luminosity"),
                finish=("planetgen.generation.run_plan", "finish_scatter"),
                prepare=("planetgen.generation.run_plan", "prepare_scatter_pass"),
                conclude=("planetgen.generation.run_plan", "conclude_scatter_pass"),
                warm=("planetgen.generation.run_plan", "scatter_shells_warm"),
                unit_seconds=5),
    CommandSpec("scatter-backfill", "Scatter the massive stars from the neighborhood: one unit per run of sectors",
                unit=("planetgen.generation.run_galaxy", "backfill_bright_stars_unit"),
                expand=("planetgen.generation.run_galaxy", "expand_backfill_bright_stars"),
                warm=("planetgen.queue.units", "warm_shape"),
                unit_seconds=12),
    CommandSpec("scatter-band", "Add the dimmer band of bright stars: one unit per run of shells",
                unit=("planetgen.generation.run_plan", "scatter_shells_unit"),
                expand=("planetgen.generation.run_plan", "expand_scatter_band"),
                finish=("planetgen.generation.run_plan", "finish_scatter"),
                warm=("planetgen.generation.run_plan", "scatter_shells_warm"),
                unit_seconds=30),
    CommandSpec("refresh-block-lights", "Make the Galaxy Map's block light totals again: one unit",
                unit=("planetgen.generation.run_plan", "refresh_block_lights_unit"),
                expand=("planetgen.generation.run_plan", "expand_refresh_block_lights"), whole=True),
    CommandSpec("refresh-run-blocks", "Bring the Galaxy Map's block totals up to what a galaxy run stored: one unit",
                unit=("planetgen.generation.run_galaxy", "refresh_run_blocks_unit"),
                expand=("planetgen.generation.run_galaxy", "expand_refresh_run_blocks"), whole=True),
    CommandSpec("fill-map-luminosity", "Fill the map luminosity: one unit per run of ids",
                unit=("planetgen.generation.run_plan", "fill_map_luminosity_unit"),
                expand=("planetgen.generation.run_plan", "expand_fill_map_luminosity"),
                unit_seconds=2),
    CommandSpec("rescore-habitability", "Re-score the habitability: one unit per run of ids",
                unit=("planetgen.generation.run_plan", "rescore_habitability_unit"),
                expand=("planetgen.generation.run_plan", "expand_rescore_habitability"),
                conclude=("planetgen.generation.run_plan", "conclude_rescore_habitability"),
                unit_seconds=25),
    CommandSpec("clean-nebula-field", "Clean up the nebula field: one unit",
                unit=("planetgen.generation.run_plan", "clean_nebula_field_unit"),
                expand=("planetgen.generation.run_plan", "expand_clean_nebula_field"), whole=True),
    CommandSpec("clear-scatter", "Clear one scatter's stored rows: one unit",
                unit=("planetgen.generation.run_plan", "clear_scatter_unit"),
                expand=("planetgen.generation.run_plan", "expand_clear_scatter"), whole=True),
    CommandSpec("draw-nebula-field", "Draw the nebula field: one unit per slab of cells",
                unit=("planetgen.generation.run_plan", "draw_nebula_field_unit"),
                expand=("planetgen.generation.run_plan", "expand_draw_nebula_field"),
                prepare=("planetgen.generation.run_plan", "prepare_draw_nebula_field"),
                conclude=("planetgen.generation.run_plan", "conclude_draw_nebula_field"),
                warm=("planetgen.queue.units", "warm_shape"),
                unit_seconds=6),
    CommandSpec("generate-start-sector", "Generate the starting sector: one unit per run of sectors",
                unit=("planetgen.generation.run_galaxy", "generate_sectors_unit"),
                expand=("planetgen.generation.run_galaxy", "expand_generate_sectors"),
                warm=("planetgen.queue.units", "warm_shape"),
                unit_seconds=20),
    CommandSpec("generate-neighborhood", "Generate the neighborhood: one unit per run of sectors",
                unit=("planetgen.generation.run_galaxy", "generate_sectors_unit"),
                expand=("planetgen.generation.run_galaxy", "expand_generate_sectors"),
                warm=("planetgen.queue.units", "warm_shape"),
                unit_seconds=20),
    CommandSpec("generate-sectors", "Generate the sectors: one unit per run of sectors",
                unit=("planetgen.generation.run_galaxy", "generate_sectors_unit"),
                expand=("planetgen.generation.run_galaxy", "expand_generate_sectors"),
                warm=("planetgen.queue.units", "warm_shape"),
                unit_seconds=20),
    CommandSpec("link-neighbours", "Link the neighbours: one unit per run of sectors",
                unit=("planetgen.generation.run_galaxy", "link_neighbours_unit"),
                expand=("planetgen.generation.run_galaxy", "expand_link_neighbours"),
                unit_seconds=90),
    CommandSpec("generate-unplaced-sector", "Generate the unplaced sectors of a `sector` run: one unit per run of sectors",
                unit=("planetgen.generation.run_sector", "unplaced_sectors_unit"),
                expand=("planetgen.generation.run_sector", "expand_unplaced_sectors"),
                unit_seconds=20),
    CommandSpec("generate-system", "Generate one system (`planetgen system`): one unit",
                unit=("planetgen.generation.run_system", "system_unit"),
                expand=("planetgen.generation.run_system", "expand_system"), whole=True),
    CommandSpec("run-population", "Run the population pass: species, civilizations and polities, one unit",
                unit=("planetgen.generation.run_population", "population_unit"),
                expand=("planetgen.generation.run_population", "expand_population"), whole=True),
    CommandSpec("run-territories", "Recompute the territories: which polity owns which system, one unit",
                unit=("planetgen.generation.run_population", "territories_unit"),
                expand=("planetgen.generation.run_population", "expand_territories"), whole=True),
    CommandSpec("reset-galaxy", "Reset the galaxy: wipe every content table of the database, one unit",
                unit=("planetgen.cli.reset", "reset_galaxy_unit"),
                expand=("planetgen.cli.reset", "expand_reset_galaxy"), whole=True),
    CommandSpec("check-math", "Check the math against known answers, one unit",
                unit=("planetgen.physics.mathcheck", "check_math_unit"),
                expand=("planetgen.physics.mathcheck", "expand_check_math"), whole=True),
    CommandSpec("check-datasets", "Check where each dataset stands, fetching nothing: the data check, one unit",
                unit=("planetgen.cli.datasets", "check_datasets_unit"),
                expand=("planetgen.cli.datasets", "expand_check_datasets"), whole=True),
    CommandSpec("settle-sector-paths", "Save the sector paths: one unit per run of sectors",
                unit=("planetgen.generation.run_galaxy", "save_sector_paths_unit"),
                expand=("planetgen.generation.run_galaxy", "expand_save_sector_paths"),
                unit_seconds=20),
    CommandSpec("warm-cluster-field", "Make the galaxy's cluster field: the mean system mass of each age in small parts, "
                "then the field as a file",
                unit=("planetgen.galaxy.field_job", "field_unit"), expand=("planetgen.galaxy.field_job", "expand_field"),
                conclude=("planetgen.galaxy.field_job", "conclude_field"), unit_seconds=5),
    CommandSpec("warm-object-index", "Build or refresh the object index: one row per object of every kind the site lists",
                unit=("planetgen.web.warmup", "warm_object_index_unit"),
                expand=("planetgen.web.warmup", "expand_warm_object_index"), whole=True),
    CommandSpec("warm-galaxy-map", "Warm the Galaxy Map's opening view, its shape, sectors and phenomena",
                unit=("planetgen.web.warmup", "warm_galaxy_map_unit"),
                expand=("planetgen.web.warmup", "expand_warm_galaxy_map"), whole=True),
    CommandSpec("warm-galaxy-nodes", "Warm the Galaxy Map's sector-block nodes: every populated block down to level 243",
                unit=("planetgen.web.warmup", "warm_galaxy_nodes_unit"),
                expand=("planetgen.web.warmup", "expand_warm_galaxy_nodes"),
                finish=("planetgen.web.warmup", "finish_warm_galaxy_nodes"),
                unit_seconds=64),
    CommandSpec("warm-galaxy-deep-nodes", "Warm the Galaxy Map's deep nodes: the zoom run in to the centre",
                unit=("planetgen.web.warmup", "warm_galaxy_deep_nodes_unit"),
                expand=("planetgen.web.warmup", "expand_warm_galaxy_deep_nodes"),
                finish=("planetgen.web.warmup", "finish_warm_galaxy_deep_nodes"),
                unit_seconds=64),
    CommandSpec("warm-lists", "Warm the first page of every list, and its counts",
                unit=("planetgen.web.warmup", "warm_lists_unit"), expand=("planetgen.web.warmup", "expand_warm_lists"),
                whole=True),
    CommandSpec("warm-sectors", "Warm every sector page: one unit per slice of the listing",
                unit=("planetgen.web.warmup", "warm_sectors_unit"),
                expand=("planetgen.web.warmup", "expand_warm_sectors"),
                finish=("planetgen.web.warmup", "finish_warm_sectors"),
                unit_seconds=75),
    CommandSpec("warm-systems", "Warm every system page: one unit per slice of the listing",
                unit=("planetgen.web.warmup", "warm_systems_unit"),
                expand=("planetgen.web.warmup", "expand_warm_systems"),
                finish=("planetgen.web.warmup", "finish_warm_systems"),
                unit_seconds=75),
    CommandSpec("warm-phenomena", "Warm every phenomenon page: one unit per slice of the listing",
                unit=("planetgen.web.warmup", "warm_phenomena_unit"),
                expand=("planetgen.web.warmup", "expand_warm_phenomena"),
                finish=("planetgen.web.warmup", "finish_warm_phenomena"),
                unit_seconds=75),
    CommandSpec("warm-search", "Warm the search page",
                unit=("planetgen.web.warmup", "warm_search_unit"), expand=("planetgen.web.warmup", "expand_warm_search"),
                whole=True),
    CommandSpec("check-revision", "Check the schema revision records and that the tables match the models",
                unit=("planetgen.db.check", "check_revision_unit"), expand=("planetgen.db.check", "expand_check"),
                whole=True),
    CommandSpec("check-tables", "Check every table (CHECK TABLE ... QUICK)",
                unit=("planetgen.db.check", "check_tables_unit"), expand=("planetgen.db.check", "expand_check"),
                whole=True),
    CommandSpec("check-tables-extended", "Check every table thoroughly (CHECK TABLE ... EXTENDED, slow)",
                unit=("planetgen.db.check", "check_tables_extended_unit"),
                expand=("planetgen.db.check", "expand_check"), whole=True),
    CommandSpec("check-orphans", "Check that every row has its parent",
                unit=("planetgen.db.check", "check_orphans_unit"), expand=("planetgen.db.check", "expand_check"),
                whole=True),
    CommandSpec("check-ids", "Check that the id blocks are ahead of every id",
                unit=("planetgen.db.check", "check_ids_unit"), expand=("planetgen.db.check", "expand_check"),
                whole=True),
    CommandSpec("check-values", "Check that columns hold values a generator could write",
                unit=("planetgen.db.check", "check_values_unit"), expand=("planetgen.db.check", "expand_check"),
                whole=True),
    CommandSpec("check-counts", "Check that each sector's stats count its systems",
                unit=("planetgen.db.check", "check_counts_unit"), expand=("planetgen.db.check", "expand_check"),
                whole=True),
    CommandSpec("check-version-keys", "Check that every sector has a well formed version key",
                unit=("planetgen.db.check", "check_version_keys_unit"),
                expand=("planetgen.db.check", "expand_check"), whole=True),
    CommandSpec("check-storage", "Check that every sealed sector still hashes to its storage checksum",
                unit=("planetgen.db.check", "check_storage_unit"), expand=("planetgen.db.check", "expand_check"),
                whole=True),
    CommandSpec("check-systems", "Check every star system against the validation rules (slow)",
                unit=("planetgen.db.check", "check_systems_unit"),
                expand=("planetgen.db.check", "expand_check_systems"), whole=True),
    CommandSpec("check-object-index", "Check that the object index holds every object as the tables do",
                unit=("planetgen.db.check", "check_object_index_unit"), expand=("planetgen.db.check", "expand_check"),
                whole=True),
    CommandSpec("end-map-rebuild-maintenance", "Open the site again after the sector-grid switch closed it",
                unit=("planetgen.web.warmup", "end_map_rebuild_unit"),
                expand=("planetgen.web.warmup", "expand_end_map_rebuild"), whole=True),
)
"""dict: Every command the wrapper can name, by kind."""

KIND_OF_STAGE = {
    "mass": "scatter-mass", "skeleton": "plan-skeleton", "phenomena": "scatter-phenomena",
    "luminosity": "scatter-luminosity", "backfill": "scatter-backfill", "band": "scatter-band",
    "map-luminosity": "fill-map-luminosity", "habitability": "rescore-habitability",
    "nebula-field": "clean-nebula-field", "nebula-draw": "draw-nebula-field", "start": "generate-start-sector",
    "neighborhood": "generate-neighborhood", "clear-mass": "clear-scatter", "clear-luminosity": "clear-scatter",
    "clear-phenomena": "clear-scatter", "clear-nebula-field": "clear-scatter",
    "sectors": "generate-sectors", "link": "link-neighbours", "population": "run-population",
    "territories": "run-territories", "settle": "settle-sector-paths", "phenomena-index": "finish-phenomena",
    "block-lights": "refresh-block-lights", "run-blocks": "refresh-run-blocks",
    "unplaced": "generate-unplaced-sector", "system": "generate-system",
}
"""dict: The command that does each stage of `generation/stages.py`."""


KIND_OF_WARM_FAMILY = {family: f"warm-{family}" for family in
                       ("cluster-field", "object-index", "galaxy-map", "lists", "galaxy-nodes", "galaxy-deep-nodes", "sectors", "systems",
                        "phenomena", "search")}
"""dict: The command that warms each cache family (`web/warmup.FAMILIES`), the steps of `planetgen warm`."""

KIND_OF_CHECK = {name: f"check-{name}" for name in (
    "revision", "tables", "tables-extended", "orphans", "ids", "values", "counts", "version-keys", "storage", "systems",
    "object-index")}
"""dict: The command that does each database check (`db/check.CHECKS`), the steps of `planetgen check-db`, in order."""

DEFAULT_CHECKS = ("revision", "tables", "orphans", "ids", "values", "counts", "version-keys", "storage")
"""tuple: The checks `planetgen check-db` runs when none is named: every check but the two slow ones."""

DEEP_CHECKS = ("revision", "tables-extended", "orphans", "ids", "values", "counts", "version-keys", "storage", "systems")
"""tuple: The list of the slow checks: the default one with the thorough table check and the star systems."""

END_MAP_REBUILD = "end-map-rebuild-maintenance"
"""str: The command that opens the site after a full warm-up; the last step of a list that warms every family."""


def spec(kind):
    """The registered command `kind`.

    Raises:
        KeyError: No such command.
    """
    return REGISTRY[kind]


def _function(pair):
    module, name = pair
    return getattr(importlib.import_module(module), name)


def unit_function(kind):
    """The function that runs one unit of `kind`.

    Raises:
        KeyError: No such command.
        ValueError: The command is inline, so it has no units.
    """
    command = REGISTRY[kind]
    if command.unit is None:
        raise ValueError(f"command {kind!r} is inline: it has no units yet")
    return _function(command.unit)


def warm_function(kind):
    """The function a worker calls with a unit's arguments before it forks the horse that runs it, or `None`
    (`CommandSpec.warm`).

    Raises:
        KeyError: No such command.
    """
    command = REGISTRY[kind]
    return None if command.warm is None else _function(command.warm)


def finish_function(kind):
    """The function the unit that completes a stage of `kind` calls with its results, or `None` for a command that
    has none (`CommandSpec.finish`).

    Raises:
        KeyError: No such command.
    """
    command = REGISTRY[kind]
    return None if command.finish is None else _function(command.finish)


def prepare_function(kind):
    """The function called once when a stage of `kind` opens, before it is cut, or `None` (`CommandSpec.prepare`).

    Raises:
        KeyError: No such command.
    """
    command = REGISTRY[kind]
    return None if command.prepare is None else _function(command.prepare)


def conclude_function(kind):
    """The function called once when the last unit of a stage of `kind` is done, or `None` (`CommandSpec.conclude`).

    Raises:
        KeyError: No such command.
    """
    command = REGISTRY[kind]
    return None if command.conclude is None else _function(command.conclude)


def expand_function(kind):
    """The function that cuts a stage of `kind` into units (`Expansion`).

    Raises:
        KeyError: No such command.
        ValueError: The command has no expansion (it is inline, or its stage is one unit given directly).
    """
    command = REGISTRY[kind]
    if command.expand is None:
        raise ValueError(f"command {kind!r} cannot be cut into units")
    return _function(command.expand)


STAGED_COMMANDS = ("galaxy", "sectors", "plan", "scatter", "clear-scatter", "fill-map-luminosity", "rescore-habitability",
                   "clean-nebula-field", "link", "backfill", "settle", "population", "territories")
"""tuple: The commands whose stages (`generation/stages.py`) the wrapper turns into steps."""


def build_steps(command, args):
    """
    The steps a `planetgen galaxy`, `plan`, `warm`, `check-db`, `sector` or `system` run with the parsed `args` will take, in order:
    exactly the commands the request selects (the stages of `generation.stages`, the cache families of
    `web/warmup`, or the checks of `db/check`), each with the settings that shape it. A pure function of the
    arguments.

    Raises:
        ValueError: For any other command, a cache family or a check that does not exist, or a setting for a
            command that is not in the list.
    """
    if command == "warm":
        return _warm_steps(args)
    if command == "check-db":
        return _check_steps(args)
    if command in ONE_COMMAND_RUNS:
        return [_one_command_step(command, args)]
    from planetgen.generation import stages
    if command not in STAGED_COMMANDS:
        raise ValueError(f"the wrapper builds steps for {', '.join(STAGED_COMMANDS + ONE_COMMAND_RUNS)}, warm and check-db, not {command!r}")
    return [Step(stage.key, KIND_OF_STAGE[stage.key], stage.label, stages.stage_settings(stage.key, args))
            for stage in stages.for_command(command, args) if stage.skip is None]


ONE_COMMAND_RUNS = ("sector", "system")
"""tuple: The commands that are a run of one command (PERF.119): `planetgen sector` and `planetgen system`."""


def _one_command_step(command, args):
    """
    The one step of a `planetgen sector` or `planetgen system` run: the generation arguments as plain JSON
    (`run_galaxy._sector_args_json`), the database (`args.database`, `wrapper.database_args`) and the seed the run's
    draws come from (`args.run_seed`, drawn by the starter, which keeps this function pure).
    """
    from planetgen.generation import run_galaxy
    params = {"database": args.database, "args": run_galaxy._sector_args_json(args), "run_seed": args.run_seed}
    if command == "sector":
        count = getattr(args, "num_sectors", None) or 1
        params["count"] = count
        return Step("unplaced", "generate-unplaced-sector", f"Sectors ({count} unplaced)", params)
    return Step("system", "generate-system", "Generate a system", params)


def _check_steps(args):
    """
    The steps of `planetgen check-db`: one command for each check named (`args.checks`, `DEFAULT_CHECKS` when none),
    in the order of `KIND_OF_CHECK`, each with the database (`args.database`, `wrapper.database_args`) and only the
    settings it takes: the scope (`args.shell`, `args.sector`) for the value, count, version-key and systems checks,
    and `wide_binaries` for the systems check.

    Raises:
        ValueError: A check that does not exist, or `--wide-binaries` with no systems check in the list.
    """
    asked = list(getattr(args, "checks", None) or DEFAULT_CHECKS)
    unknown = [name for name in asked if name not in KIND_OF_CHECK]
    if unknown:
        raise ValueError(f"no such check {unknown[0]!r}; one of {', '.join(KIND_OF_CHECK)}")
    scope = {"shells": list(args.shell or ()), "addresses": [list(address) for address in args.sector or ()]}
    steps = []
    for name in KIND_OF_CHECK:
        if name not in asked:
            continue
        params = {"database": args.database}
        if name in ("values", "counts", "version-keys", "storage", "systems"):
            params.update(scope)
        if name == "systems":
            params["wide_binaries"] = bool(args.wide_binaries)
        steps.append(Step(name, KIND_OF_CHECK[name], f"Check {name}", params))
    if getattr(args, "wide_binaries", False) and "systems" not in asked:
        raise ValueError("--wide-binaries is for the systems check")
    return steps


def _warm_steps(args):
    """
    The steps of `planetgen warm`: one command for each cache family asked for (`args.families`, every family when
    none), in the order of `web/warmup.FAMILIES`, each with the site's database (`args.database`) and, for the
    pages of objects, the most to warm (`args.cap`). A list of every family ends with the command that opens the site
    after the sector-grid switch: it is the list, not any family, that knows the maps are then rebuilt.
    """
    from planetgen.web import warmup
    asked = list(getattr(args, "families", None) or ())
    families = warmup.expand_phases(asked) if asked else list(warmup.FAMILIES)
    settings = {"database": args.database, "cap": getattr(args, "cap", None)}
    if getattr(args, "index_mode", None):
        settings["index_mode"] = args.index_mode
    steps = [Step(family, KIND_OF_WARM_FAMILY[family], f"Warm {family}: {warmup.FAMILIES[family]}", dict(settings))
             for family in families]
    if len(families) == len(warmup.FAMILIES):
        steps.append(Step("reopen", END_MAP_REBUILD, "Open the site after the map rebuild", {}))
    return steps


def run_alone(kind, args, report=None):
    """
    Runs one unit of `kind` from its JSON `args` and nothing else, as a worker would: a test or a rerun of one chunk
    of work (`planetgen command run`).

    Args:
        report (callable, optional): Called with `(done, estimate)` as the unit works.

    Returns:
        The unit's JSON-able result.
    """
    return unit_function(kind)(args, report or (lambda done, estimate: None))
