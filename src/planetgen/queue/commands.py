# planetgen/queue/commands.py

"""
Commands, and the step list the wrapper builds from a request (PERF.118, docs/design/redis-and-the-work-queue.md
5.12 and docs/design/command-inventory.md).

**A command is a kind in `REGISTRY`.** It does one thing and nothing else: it does not know what ran before it, what
runs after it or where it sits in the order. A command with a `unit` is queued as a stage of units; one with an
`expand` is cut into those units when its stage starts, from its `params` (plain JSON) and the database alone. A
command with neither is *inline*: it still runs in the process that started the run (the stages PERF.106 to PERF.108
and the commands of PERF.121 to PERF.128 turn into units, one at a time). Every command the wrapper can name is
registered here, so a list can be printed and checked without running anything.

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
    returns the `Expansion` of the command's stage. A command with no `unit` is inline.
    """

    kind: str
    summary: str
    unit: tuple = None
    expand: tuple = None

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
                expand=("planetgen.generation.run_plan", "expand_scatter_mass")),
    CommandSpec("plan-skeleton", "Plan the galaxy's density skeleton"),
    CommandSpec("scatter-phenomena", "Scatter the phenomena"),
    CommandSpec("scatter-luminosity", "Scatter the bright stars"),
    CommandSpec("scatter-backfill", "Scatter the massive stars from the neighborhood"),
    CommandSpec("scatter-band", "Add the dimmer band of bright stars"),
    CommandSpec("fill-map-luminosity", "Fill the map luminosity"),
    CommandSpec("rescore-habitability", "Re-score the habitability"),
    CommandSpec("clean-nebula-field", "Clean up the nebula field"),
    CommandSpec("generate-start-sector", "Generate the starting sector"),
    CommandSpec("generate-neighborhood", "Generate the neighborhood"),
    CommandSpec("generate-sectors", "Generate the sectors"),
    CommandSpec("link-neighbours", "Link the neighbours"),
    CommandSpec("run-population", "Run the population pass"),
    CommandSpec("settle-sector-paths", "Save the sector paths"),
    CommandSpec("warm-galaxy-map", "Warm the Galaxy Map's opening view, its shape, sectors and phenomena",
                unit=("planetgen.web.warmup", "warm_galaxy_map_unit"),
                expand=("planetgen.web.warmup", "expand_warm_galaxy_map")),
    CommandSpec("warm-lists", "Warm the first page of every list, and its counts",
                unit=("planetgen.web.warmup", "warm_lists_unit"), expand=("planetgen.web.warmup", "expand_warm_lists")),
    CommandSpec("warm-sectors", "Warm every sector page",
                unit=("planetgen.web.warmup", "warm_sectors_unit"),
                expand=("planetgen.web.warmup", "expand_warm_sectors")),
    CommandSpec("warm-systems", "Warm every system page",
                unit=("planetgen.web.warmup", "warm_systems_unit"),
                expand=("planetgen.web.warmup", "expand_warm_systems")),
    CommandSpec("warm-phenomena", "Warm every phenomenon page",
                unit=("planetgen.web.warmup", "warm_phenomena_unit"),
                expand=("planetgen.web.warmup", "expand_warm_phenomena")),
    CommandSpec("warm-search", "Warm the search page",
                unit=("planetgen.web.warmup", "warm_search_unit"), expand=("planetgen.web.warmup", "expand_warm_search")),
    CommandSpec("end-map-rebuild-maintenance", "Open the site again after the sector-grid switch closed it",
                unit=("planetgen.web.warmup", "end_map_rebuild_unit"),
                expand=("planetgen.web.warmup", "expand_end_map_rebuild")),
)
"""dict: Every command the wrapper can name, by kind."""

KIND_OF_STAGE = {
    "mass": "scatter-mass", "skeleton": "plan-skeleton", "phenomena": "scatter-phenomena",
    "luminosity": "scatter-luminosity", "backfill": "scatter-backfill", "band": "scatter-band",
    "map-luminosity": "fill-map-luminosity", "habitability": "rescore-habitability",
    "nebula-field": "clean-nebula-field", "start": "generate-start-sector", "neighborhood": "generate-neighborhood",
    "sectors": "generate-sectors", "link": "link-neighbours", "population": "run-population",
    "settle": "settle-sector-paths",
}
"""dict: The command that does each stage of `generation/stages.py`."""


KIND_OF_WARM_FAMILY = {family: f"warm-{family}" for family in
                       ("galaxy-map", "lists", "sectors", "systems", "phenomena", "search")}
"""dict: The command that warms each cache family (`web/warmup.FAMILIES`), the steps of `planetgen warm`."""

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


def build_steps(command, args):
    """
    The steps a `planetgen galaxy`, `planetgen plan` or `planetgen warm` run with the parsed `args` will take, in
    order: exactly the commands the request selects (the stages of `generation.stages`, or the cache families of
    `web/warmup`), each with the settings that shape it. A pure function of the arguments.

    Raises:
        ValueError: For any other command, or a cache family that does not exist.
    """
    if command == "warm":
        return _warm_steps(args)
    from planetgen.generation import stages
    if command not in ("galaxy", "plan"):
        raise ValueError(f"the wrapper builds steps for galaxy, plan and warm, not {command!r}")
    return [Step(stage.key, KIND_OF_STAGE[stage.key], stage.label, stages.stage_settings(stage.key, args))
            for stage in stages.for_command(command, args) if stage.skip is None]


def _warm_steps(args):
    """
    The steps of `planetgen warm`: one command for each cache family asked for (`args.families`, every family when
    none), in the order of `web/warmup.FAMILIES`, each with the site's database (`args.database`) and, for the
    pages of objects, the most to warm (`args.cap`). A list of every family ends with the command that opens the site
    after the sector-grid switch: it is the list, not any family, that knows the maps are then rebuilt.
    """
    from planetgen.web import warmup
    asked = list(getattr(args, "families", None) or ())
    unknown = [family for family in asked if family not in warmup.FAMILIES]
    if unknown:
        raise ValueError(f"no such cache family {unknown[0]!r}; one of {', '.join(warmup.FAMILIES)}")
    families = [family for family in warmup.FAMILIES if not asked or family in asked]
    settings = {"database": args.database, "cap": getattr(args, "cap", None)}
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
