# planetgen/generation/job.py

"""
One run for a Generate job (PERF.203, docs/design/operations/one-run-per-job.md): `planetgen run-job SPEC`.

A job spec is a JSON list of `{"command", "argv"}`, the commands the Generate page runs one after another today and the options
each takes (what follows the command's name on its command line). `run-job` puts the followers of every one of them in a single
`followed.run_followed`, so the job is one run in the queue with one stage list and one estimate, and what a command checks before
it starts is checked when its own stage opens, because an earlier stage of the same run may make the thing true (a plan, a
clear, the sectors a link works on).

The followers of a command are the ones the verb runs: the verbs of `generation/run_plan.py` and `run_galaxy.py` build their
runs from the same functions, so the single verb and the job cannot differ.
"""

import json

from planetgen.generation import followed, run_common, run_galaxy, run_plan, stages
from planetgen.queue import units, wrapper
from planetgen.util import log

JOB_COMMANDS = ("plan", "sectors", "link", "backfill", "settle", "scatter", "clear-scatter")
"""tuple: The commands a job can hold (the ones the Generate page's jobs run), in the words of `planetgen <command>`."""

SCATTER_KINDS = ("mass", "luminosity", "phenomena")
"""tuple: The kinds of `scatter` a job can hold; the band and the nebula field do work outside a stage."""

SHARED = ("mysql_host", "mysql_port", "mysql_user", "mysql_password", "mysql_database", "workers", "debug", "debug_areas",
          "record_file", "record_json", "quiet", "disregard_missing_data")
"""tuple: The options of `run-job` itself that every command of the job takes from it."""


class Part:
    """
    One command of a job, ready to run.

    Attributes:
        command (str): The command's name.
        args (argparse.Namespace): Its parsed and validated arguments.
        stages (list[stages.Stage]): Its stages, in order (`stages.for_command`).
        followers (list[followed.Follower]): One for each of its stages.
        after (callable or None): Says what the command made, once the run has ended.
    """

    def __init__(self, command, args, stages_, followers, after=None):
        self.command, self.args, self.stages, self.followers, self.after = command, args, stages_, followers, after


def parse_spec(text):
    """
    The entries of a job spec: `text` is JSON, or `@FILE` to read it from a file.

    Returns:
        list[tuple]: `(command, argv)` for each entry.

    Raises:
        ValueError: Not a JSON list of `{"command", "argv"}`, or a command a job cannot hold.
    """
    if text.startswith("@"):
        with open(text[1:], encoding="utf-8") as handle:
            text = handle.read()
    try:
        spec = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"the job spec is not JSON: {exc}") from exc
    if not isinstance(spec, list) or not spec:
        raise ValueError("the job spec is a list of {command, argv}, with at least one entry")
    entries = []
    for number, entry in enumerate(spec, 1):
        if not (isinstance(entry, dict) and isinstance(entry.get("command"), str) and isinstance(entry.get("argv", []), list)
                and all(isinstance(word, str) for word in entry.get("argv", []))):
            raise ValueError(f"job entry {number} is not {{\"command\": text, \"argv\": [text, ...]}}")
        if entry["command"] not in JOB_COMMANDS:
            raise ValueError(f"job entry {number}: a job holds {', '.join(JOB_COMMANDS)}, not {entry['command']!r}")
        entries.append((entry["command"], entry.get("argv", [])))
    return entries


def _unsupported(command, args):
    """Why a job cannot hold this command with these arguments, or `None`."""
    if command == "scatter" and stages.scatter_kind(args) not in SCATTER_KINDS:
        return f"a job holds scatter --{', --'.join(SCATTER_KINDS)}; the band and the nebula field are not stages of a run"
    if command == "settle" and getattr(args, "sector_id", None):
        return "a job settles the sectors it saved; --sector-id names others"
    if command == "sectors" and getattr(args, "region", False):
        return "--region moves the run's center by reading the galaxy, which a job has not made yet"
    return None


def parse_entries(text, shared, command_parsers, validate):
    """
    The commands of a job spec as `(command, args)`: each entry's options parsed and validated by its own command's parser (what
    `planetgen <command> ...` does), the options of `run-job` itself (`SHARED`) taken over.

    Args:
        text (str): The spec (`parse_spec`).
        shared (argparse.Namespace): The arguments of `run-job`.
        command_parsers (dict): Each command's parser (`cli.generate.build_parser`).
        validate (callable): `validate(args, parser)` for an entry's parsed arguments (`cli.generate.validate_command`).

    Returns:
        list[tuple]: `(command, args)` for each entry.

    Raises:
        ValueError: A bad spec, or an entry that is not valid (the parser has said why).
    """
    entries = []
    for number, (command, argv) in enumerate(parse_spec(text), 1):
        parser = command_parsers[command]
        try:
            args = parser.parse_args(argv)
            args.command = command
            validate(args, parser)
        except SystemExit as exc:
            raise ValueError(f"job entry {number} ({command} {' '.join(argv)}) is not valid") from exc
        reason = _unsupported(command, args)
        if reason:
            raise ValueError(f"job entry {number} ({command}): {reason}")
        for name in SHARED:
            if hasattr(shared, name):
                setattr(args, name, getattr(shared, name))
        entries.append((command, args))
    return entries


def _plan(args, _index, _since):
    follower = run_plan.SkeletonFollower(args)
    return Part("plan", args, stages.plan_stages(args), [follower],
                lambda: run_plan.log_skeleton_summary(args, follower.summary))


def _sectors(args, _index, _since):
    args.link_later = True              # the run links its sectors after them, as `planetgen galaxy` does
    followers = run_galaxy._deferred_sectors(args, run_common._edge_pc())
    followers[0].extra["check_plan"] = True         # the plan may be made by an earlier stage of the job: found when this one opens
    return Part("sectors", args, stages.sectors_stages(args), followers)


def _after_sectors(command):
    """The part of `link`, `backfill` or `settle`: the sectors saved since the job's sectors stage opened, else since `--since` (or
    the newest `sectors` run's start) when the job holds no sectors stage of its own."""
    follower_of = {"link": run_galaxy.LinkFollower, "backfill": run_galaxy.BackfillFollower, "settle": run_galaxy.SettleFollower}[command]

    def build(args, _index, since):
        if since is None:
            run_galaxy._resolve_since(args)
            since = args.since
        elif getattr(args, "since", None) is not None:
            since = args.since
        return Part(command, args, stages.for_command(command, args), [follower_of(args, since)])

    return build


def _scatter(args, index, _since):
    kind = stages.scatter_kind(args)
    tag = f"#{index}"
    if kind == "phenomena":
        planned = run_plan.planned_phenomena(args, kind, tag)
    else:
        planned = run_plan.planned_bright_stars(args, (kind,), kind, tag)
    return Part("scatter", args, stages.scatter_stages(args), planned.followers)


def _clear_scatter(args, _index, _since):
    return Part("clear-scatter", args, stages.clear_stages(args), [run_plan.ClearScatterFollower(args)])


BUILDERS = {"plan": _plan, "sectors": _sectors, "link": _after_sectors("link"), "backfill": _after_sectors("backfill"),
            "settle": _after_sectors("settle"), "scatter": _scatter, "clear-scatter": _clear_scatter}
"""dict: Each command's `build(args, index, since)`, which gives its `Part`. `since` is the deferred clock of the job's first sectors
stage, `None` when the job has none."""


def build_parts(entries, run_args):
    """
    The `Part` of each entry (`entries`: `(command, args)`), built before anything runs: no builder reads the database for what an
    earlier stage of the job makes (the settings a stage works from are deferred parameters, resolved when it opens).

    Raises:
        SystemExit: A `link`, `backfill` or `settle` with no sectors stage before it and no run to take its start from.
    """
    since, parts = None, []
    for index, (command, args) in enumerate(entries):
        args.generation_run_id = getattr(run_args, "generation_run_id", None)
        part = BUILDERS[command](args, index, since)
        if command == "sectors" and since is None:
            since = wrapper.opened_by(part.followers[0].key)
        parts.append(part)
    return parts


def stage_list(parts):
    """The stages of the job, in order: each part's own, the stages `followed.dedupe` drops left out."""
    kept = {id(f) for f in followed.dedupe([f for part in parts for f in part.followers])}
    listed = []
    for part in parts:
        for stage, follower in zip(part.stages, part.followers):
            if id(follower) in kept:
                listed.append(stage)
    return listed


def stage_estimates(conn, entries):
    """
    The expected seconds of each stage of the job `entries` (`parse_entries`) make, in the order of `stage_list`: each command's own
    (`stages.stage_estimates`, at the settings that command was given), the stage every scatter ends with only at its last place.
    """
    rows = []
    for command, args in entries:
        rows += [(stage.key, seconds)
                 for stage, seconds in zip(stages.for_command(command, args), stages.stage_estimates(conn, command, args))]
    last = max((place for place, (key, _seconds) in enumerate(rows) if key == stages.BLOCK_LIGHTS_STAGE.key), default=None)
    return [seconds for place, (key, seconds) in enumerate(rows) if key != stages.BLOCK_LIGHTS_STAGE.key or place == last]


def run_job(args):
    """
    `planetgen run-job SPEC`: every command of the spec as the stages of one run (`followed.run_followed`), numbered as one
    (`stages.begin`), each stage's time stored with the settings its own command was given. A command's refusal (a stage that opens
    and finds what it needs is missing) ends the run with the words the verb says.

    Args:
        args (argparse.Namespace): Validated arguments (`command == "run-job"`), with `job` the entries `parse_entries` made.
    """
    parts = build_parts(args.job, args)
    for part in parts:
        assert [f.key for f in part.followers] == [stage.key for stage in part.stages], f"the followers of {part.command} are its stages"
    followers = followed.dedupe([f for part in parts for f in part.followers])
    listed = stage_list(parts)
    key_args = {follower.key: part.args for part in parts for follower in part.followers}
    scatter = next((part.args for part in parts if part.command == "scatter"), args)
    key_args[stages.JOB_KEY] = scatter
    stages.begin(listed, args, "run-job", key_args)
    settings = {key: stages.stage_settings(key, key_args[key]) for key in key_args if key != stages.JOB_KEY}
    try:
        followed.run_followed(args, followers, "Generate job", settings=settings)
    except (units.RunStopped, RuntimeError) as exc:
        refusal = _refusal(exc)
        if refusal is None:
            raise
        log.error(refusal)
        raise SystemExit(1) from exc
    stages.finish()
    for part in parts:
        if part.after is not None:
            part.after()


REFUSALS = ("RuntimeError: ", "ValueError: ")
"""tuple: What a stage's error starts with when it refused to run (a worker raised it in words)."""


def _refusal(exc):
    """The words a stage refused with (`REFUSALS`), or `None` for any other failure."""
    text = str(exc)
    if not isinstance(exc, units.RunStopped):       # a stage that opened in this process raised it as it was
        return text
    for marker in REFUSALS:
        if marker in text:
            return text.split(marker, 1)[1]
    return None
