# planetgen/generation/stages.py

"""
The stages of a generation command, numbered and announced as the run reaches
them (UX.89).

A command that does several things in turn (`galaxy`) declares its
stages up front with `begin`; each is entered with `enter(key)` as the run
reaches it. A stage the run will not do (no shells to scatter, no sector made to link) is still
numbered and announced, as "skipped" with its reason, so the count of stages is
always the whole list. Each announcement is a log line ("Stage 4 of 9: ...")
and goes into the progress file (`progress_file.set_stage`), where the Generate
page reads which stage a step is in.

`for_command` gives the same list from a command line's arguments without
running anything, so the Generate page can number every stage of a job before
it starts.
"""

import json
import os
import time

from planetgen.galaxy import dials
from planetgen.queue import progress_file
from planetgen.util import log, logrecord

SECTOR_MODES = ("block", "span", "course", "around", "shell", "center_sector")

STAGE_OFFSET_ENV = "PLANETGEN_STAGE_OFFSET"
"""str: Set by a web job's runner to the number of stages the job's earlier steps hold, so this run's stages
are numbered across the whole job ("Stage 4 of 12") instead of from 1."""

STAGE_TOTAL_ENV = "PLANETGEN_STAGE_TOTAL"
"""str: Set with `STAGE_OFFSET_ENV` to the number of stages in the whole job."""

JOB_KEY = "job"
"""str: The `stage_key` the whole-job row is stored under (`finish`)."""

SCATTER_KEY = "scatter"
"""str: The `stage_key` of the row that adds up every scatter stage of a run (`finish`)."""

SCATTER_KEYS = ("mass", "luminosity", "phenomena", "backfill")
"""tuple[str]: The stage keys that are scatters: the massive stars, the bright stars, the phenomena and the
backfill, each also stored on its own."""

_state = {"stages": [], "next": 0, "args": None, "command": None, "open": None,
          "job_shells": 0, "job_seconds": 0.0, "overall": None, "ran": 0, "scatters": [], "key_args": {}}
"""The run's stage list, the next one to enter, its arguments (for the settings each stage is recorded
with), its command, the stage open now (`_open`), the shells and seconds of the stages finished so far (the
whole-job stat), the command-line overall bar's clock and estimate (`overall`), how many stages have
finished running (`ran`, skipped ones not counted) and the arguments of the stages that have their own
(`key_args`, PERF.203: a run of several commands)."""


class Stage:
    """One stage: `key`, `label` and `skip` (the reason it will not run, or `None`)."""

    def __init__(self, key, label, skip=None):
        self.key, self.label, self.skip = key, label, skip

    def as_dict(self):
        return {"label": self.label, "skipped": self.skip}


def _random_start(args):
    return not any(getattr(args, name, None) is not None and getattr(args, name) is not False for name in SECTOR_MODES)


def sectors_stages(args):
    """The stages of one `planetgen sectors` run, in order: the sectors, and nothing else (PERF.123)."""
    if _random_start(args):
        return [Stage("start", "Generate the starting sector"), Stage("neighborhood", "Generate the neighborhood")]
    if getattr(args, "wedge", None) is not None and getattr(args, "radius_pc", None) is not None and not any(
            getattr(args, name, None) is not None for name in ("span", "course", "block", "around")):
        return [Stage("sectors", "Generate the sector"), Stage("neighborhood", "Generate the neighborhood")]
    return [Stage("sectors", "Generate the sectors")]


PHENOMENA_STAGES = (Stage("phenomena", "Scatter the phenomena"), Stage("phenomena-index", "Index the phenomena"),
                    Stage("nebula-draw", "Draw the nebula field"))
"""tuple[Stage]: The stages of the phenomena scatter (PERF.172): the shells, the special phenomena and the indexes, and the
nebula field."""

BLOCK_LIGHTS_STAGE = Stage("block-lights", "Make the block light totals")
"""Stage: The stage that makes the Galaxy Map's block light totals again, the last of a scatter (PERF.172)."""


def galaxy_stages(args):
    """The stages of one `planetgen galaxy` run, in order: all of them are stages of the one run (PERF.172)."""
    return sectors_stages(args) + [
        Stage("link", "Link the neighbours"),
        Stage("mass", "Scatter the massive stars"),
        Stage("luminosity", "Scatter the bright stars"),
        *PHENOMENA_STAGES,
        BLOCK_LIGHTS_STAGE,
        Stage("backfill", "Scatter the massive stars from the neighborhood"),
        Stage("population", "Run the population pass"),
        Stage("territories", "Recompute the territories"),
        Stage("settle", "Save the sector paths"),
        Stage("run-blocks", "Update the map's block totals"),
    ]


def plan_stages(args):
    """The stages of one `planetgen plan` run: the skeleton, and nothing else (PERF.121)."""
    return [Stage("skeleton", "Plan the galaxy")]


SCATTER_KINDS = ("mass", "luminosity", "phenomena", "nebula-field", "band")
"""tuple[str]: The kinds of `planetgen scatter`, one per run (PERF.121): the options `--mass`, `--luminosity`,
`--phenomena`, `--nebula-field` and `--band`."""

SCATTER_STAGE = {"mass": Stage("mass", "Scatter the massive stars"),
                 "luminosity": Stage("luminosity", "Scatter the bright stars"),
                 "phenomena": Stage("phenomena", "Scatter the phenomena"),
                 "nebula-field": Stage("nebula-draw", "Draw the nebula field"),
                 "band": Stage("band", "Add the dimmer band of bright stars")}
"""dict: The stage of each kind of scatter."""


def scatter_kind(args):
    """The one kind a `planetgen scatter` run draws (`args.mass`, `args.luminosity`, ...), `None` when none is set."""
    chosen = [kind for kind in SCATTER_KINDS if getattr(args, kind.replace("-", "_"), False)]
    return chosen[0] if len(chosen) == 1 else None


def scatter_stages(args):
    """The stages of one `planetgen scatter` run: its kind's stage (PERF.121) and what follows it (PERF.172): the
    phenomena scatter is its shells, its special phenomena and indexes and the nebula field, and every scatter that
    changes the stars or phenomena ends with the block light totals made again."""
    kind = scatter_kind(args)
    if kind is None:
        return []
    if kind == "phenomena":
        return [*PHENOMENA_STAGES, BLOCK_LIGHTS_STAGE]
    if kind in ("mass", "luminosity"):
        return [SCATTER_STAGE[kind], BLOCK_LIGHTS_STAGE]
    return [SCATTER_STAGE[kind]]


CLEAR_KINDS = ("mass", "luminosity", "phenomena", "nebula-field")
"""tuple[str]: The kinds of `planetgen clear-scatter` (PERF.121, Boss's F7): the scatters that can be cleared, one per
run. The band is part of the bright stars and goes with them."""

CLEAR_STAGE = {"mass": Stage("clear-mass", "Clear the massive stars"),
               "luminosity": Stage("clear-luminosity", "Clear the bright stars"),
               "phenomena": Stage("clear-phenomena", "Clear the phenomena"),
               "nebula-field": Stage("clear-nebula-field", "Clear the nebula field")}
"""dict: The stage of each kind of clear."""


CLEAR_KEYS = tuple(stage.key for stage in CLEAR_STAGE.values())
"""tuple[str]: The keys of the clear stages."""


def clear_kind(args):
    """The one kind a `planetgen clear-scatter` run clears, `None` when none is set."""
    chosen = [kind for kind in CLEAR_KINDS if getattr(args, kind.replace("-", "_"), False)]
    return chosen[0] if len(chosen) == 1 else None


def clear_stages(args):
    """The stages of one `planetgen clear-scatter` run: its kind's stage alone."""
    kind = clear_kind(args)
    return [] if kind is None else [CLEAR_STAGE[kind]]


SINGLE_STAGE_COMMANDS = {"fill-map-luminosity": Stage("map-luminosity", "Fill the map luminosity"),
                         "rescore-habitability": Stage("habitability", "Re-score the habitability"),
                         "fill-nebula-blocks": Stage("nebula-blocks", "Make the nebula block totals"),
                         "fill-block-totals": Stage("block-totals", "Make the block totals"),
                         "clean-nebula-field": Stage("nebula-field", "Clean up the nebula field"),
                         "link": Stage("link", "Link the neighbours"),
                         "backfill": Stage("backfill", "Scatter the massive stars from the neighborhood"),
                         "settle": Stage("settle", "Save the sector paths"),
                         "population": Stage("population", "Run the population pass"),
                         "territories": Stage("territories", "Recompute the territories")}
"""dict: The commands that do one stage and nothing else (PERF.122, PERF.123), by command name."""


def for_command(command, args):
    """The stages of `command` (`galaxy`, `plan`, `scatter`, `clear-scatter`, or a one-stage command of
    `SINGLE_STAGE_COMMANDS`) for
    parsed `args`, `[]` for any other."""
    if command == "galaxy":
        return galaxy_stages(args)
    if command == "sectors":
        return sectors_stages(args)
    if command == "plan":
        return plan_stages(args)
    if command == "scatter":
        return scatter_stages(args)
    if command == "clear-scatter":
        return clear_stages(args)
    if command == "fill-map-luminosity":
        return [SINGLE_STAGE_COMMANDS[command], BLOCK_LIGHTS_STAGE]
    if command in SINGLE_STAGE_COMMANDS:
        return [SINGLE_STAGE_COMMANDS[command]]
    return []


def stage_settings(key, args):
    """
    The arguments that shape stage `key` of a run with `args` (PERF.56), a
    flat dict of numbers, strings and booleans: what a later estimate matches
    runs on. Defaults are resolved (`--phenomenon-min-mass` unset is the
    default mass limit), so runs at the default and at the same explicit value
    match.
    """
    from planetgen import tuning
    workers = getattr(args, "workers", None) or 1
    mass = getattr(args, "phenomenon_min_mass", None)
    mass = tuning.PHENOMENON_MIN_MASS_SOLAR if mass is None else float(mass)
    floor = getattr(args, "bright_star_min_luminosity", None)
    floor = tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL if floor is None else float(floor)
    settings = {"workers": int(workers)}
    if key in ("start", "neighborhood", "sectors"):
        for name in ("radius_pc", "neighborhoods", "max_shell", "min_start_density"):
            value = getattr(args, name, None)
            if value is not None:
                settings[name] = value
    elif key == JOB_KEY:
        settings["mass_limit_sol"], settings["luminosity_floor_sol"] = mass, floor
    elif key == "mass":
        settings["mass_limit_sol"] = mass
    elif key == "luminosity":
        settings["mass_limit_sol"], settings["luminosity_floor_sol"] = mass, floor
    elif key == "phenomena":
        compact = getattr(args, "compact_min_mass", None)
        settings["compact_limit_sol"] = mass if compact in (None, tuning.COMPACT_MIN_MASS_STAR) else float(compact)
    elif key in CLEAR_KEYS:
        settings = {"part": key[len("clear-"):]}
    elif key == "skeleton":
        for name in ("edge", "seed"):
            value = getattr(args, name, None)
            if value is not None:
                settings[name] = dials.lever(name).format(value) if name == "edge" else str(value)
    return settings


def note(**metrics):
    """Adds what the stage open now did (shells visited, shells that placed something, objects ...)."""
    if _state["open"] is not None:
        _state["open"]["metrics"].update(metrics)


def recording_enabled():
    """Whether stage rows and speed stats are stored at all (`PLANETGEN_GENERATION_STATS`, off in the test suite)."""
    from planetgen.generation import stats as generation_stats
    return os.environ.get(generation_stats.STATS_ENV_VAR, "1").strip().lower() not in ("0", "off", "no", "false")


def _recording_on():
    return _state["args"] is not None and recording_enabled()


def _write(record):
    """Stores one stage's row in the control database; a failure only warns in the debug log."""
    if not _recording_on():
        return
    from planetgen.db import store
    write_row(store.mysql_config_from_args(_state["args"]), _state["command"] or "", record)


def write_row(config, command, record):
    """Stores one stage's row for the database `config` names (`record`: `key`, `n`, `total`, `label`, `skipped`,
    `settings`, `metrics` and `seconds`) in the control database; a failure only warns in the debug log."""
    try:
        from planetgen.db import store
        from planetgen.galaxy.version_key import version_key
        conn = store.get_control_connection(store.control_mysql_config(config))
        try:
            conn.execute(
                "INSERT INTO generation_stage_runs (database_name, command, stage_key, stage_n, stage_total, label,"
                " skipped, skip_reason, started_at, finished_at, seconds, workers, settings, metrics, version_key)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, DATE_SUB(NOW(6), INTERVAL ? MICROSECOND), NOW(6), ?, ?, ?, ?, ?)",
                (config.database, command, record["key"], record["n"], record["total"],
                 record["label"][:80], 1 if record["skipped"] else 0,
                 (record["skipped"] or None) and record["skipped"][:200],
                 int(record["seconds"] * 1e6), record["seconds"], record["settings"].get("workers", 1),
                 json.dumps(record["settings"], sort_keys=True), json.dumps(record["metrics"], sort_keys=True),
                 version_key()))
            conn.commit()
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- statistics are a nicety, never a reason to fail a run
        log.debug(f"Stage statistics: not recorded ({exc}).")


def _with_rate(record):
    """Adds the stage's seconds per shell to its metrics when it visited shells; returns the shells."""
    shells = int(record["metrics"].get("shells") or 0)
    if shells > 0 and record["seconds"] > 0:
        record["metrics"]["seconds_per_shell"] = record["seconds"] / shells
    return shells


def open_key():
    """The key of the stage open now, or `None` when no stage is open (an estimate, or a caller outside a run)."""
    record = _state["open"]
    return None if record is None else record["key"]


def hand_over_open_stage():
    """
    Hands the stage open now to the queue (PERF.106): its units run elsewhere, so the unit that completes the stage
    writes its row (`write_queued_row`, from what the run record holds) and this process does not. This process still
    counts the stage in the whole job's totals.

    Returns:
        dict or None: The row's `command`, `n`, `total`, `label` and `settings` for the run record to carry; `None`
            when no row is stored (no stage open, or statistics off).
    """
    record = _state["open"]
    if record is None or not _recording_on():
        return None
    record["queued"] = True
    return {"command": _state["command"] or "", "n": record["n"], "total": record["total"], "label": record["label"],
            "settings": dict(record["settings"])}


def row_for(key):
    """
    The row of stage `key` of the run begun now, before the stage is entered (PERF.172): a run that queues every
    stage up front hands each its row at the start, for the unit that completes it to write (`write_queued_row`).
    `None` when statistics are off or the run's list has no such stage.
    """
    keys = [stage.key for stage in _state["stages"]]
    if not _recording_on() or key not in keys:
        return None
    index = keys.index(key)
    return {"command": _state["command"] or "", "n": index + 1, "total": len(keys),
            "label": _state["stages"][index].label, "settings": stage_settings(key, args_of(key))}


def write_skipped_row(config, key, row, reason):
    """The row of a stage that had nothing to do (`row` as `hand_over_open_stage` gives it, `reason` the words for it):
    stored as skipped, so no estimate counts it as a run of the stage."""
    write_row(config, row["command"], {"key": key, "n": row["n"], "total": row["total"], "label": row["label"],
                                       "skipped": reason, "settings": row["settings"], "metrics": {}, "seconds": 0.0})


def write_queued_row(config, key, row, seconds, metrics):
    """The row of a stage that ran on the queue (`hand_over_open_stage`'s `row`), written by the unit that completed
    it: `seconds` from the first unit's start to the last one's end and the `metrics` its command found."""
    record = {"key": key, "n": row["n"], "total": row["total"], "label": row["label"], "skipped": None,
              "settings": row["settings"], "metrics": dict(metrics), "seconds": seconds}
    _with_rate(record)
    write_row(config, row["command"], record)


def _say_end(record, seconds, reason=None):
    """
    The stage's finish line (`seconds` after its start) or, with a `reason`, its failure line, on the command line and in a record file,
    as the start line is (UX.122). A step of one stage under a web job's banner says neither: the job runner writes both for it.
    """
    if _under_banner():
        return
    offset, job_total = _numbering()
    if reason is None:
        mark = logrecord.stage_record(logrecord.STAGE_FINISH, "stages", offset + record["n"], job_total, seconds=seconds)
        log.normal(mark.message, extra={"fields": mark.fields})
    else:
        mark = logrecord.stage_record(logrecord.STAGE_FAILED, "stages", offset + record["n"], job_total, seconds=seconds,
                                      reason=reason)
        log.error(mark.message, extra={"fields": mark.fields})


def _close_open():
    """Records the stage open now, with the time it took, and says it finished."""
    record = _state["open"]
    _state["open"] = None
    if record is not None:
        _state["ran"] += 1
        record["seconds"] = time.perf_counter() - record["t0"]
        _say_end(record, record["seconds"])
        _state["job_seconds"] += record["seconds"]
        _state["job_shells"] += _with_rate(record)
        if record["key"] in SCATTER_KEYS:
            _state["scatters"].append(record)
        if not record.get("queued"):
            _write(record)


def _write_job():
    """
    Stores the totals of the run: the whole process under `JOB_KEY` (its seconds, the stages it ran, and the
    shells per second of the stages that visited shells, PERF.55) and, when scatter stages ran, one row under
    `SCATTER_KEY` that adds them all up (operations, shells, objects, seconds, seconds per shell and per
    operation), so a later estimate has the scatter kinds' common figure when one kind has no record of its
    own. Both carry the settings the run had; nothing when no stage ran.
    """
    shells, seconds = _state["job_shells"], _state["job_seconds"]
    if _state["ran"] < 1 or not _state["stages"]:
        return
    total = len(_state["stages"])
    args = _state["args"]
    settings = stage_settings(JOB_KEY, args_of(JOB_KEY)) if args is not None else {}
    metrics = {"stages_run": _state["ran"], "shells": shells}
    if shells > 0 and seconds > 0:
        metrics["shells_per_second"] = shells / seconds
    _write({"key": JOB_KEY, "n": total, "total": total, "label": "Whole job", "skipped": None,
            "settings": settings, "metrics": metrics, "seconds": seconds})
    scatters = _state["scatters"]
    if scatters:
        scatter_seconds = sum(record["seconds"] for record in scatters)
        scatter_shells = sum(int(record["metrics"].get("shells") or 0) for record in scatters)
        scatter_metrics = {"operations": len(scatters), "shells": scatter_shells,
                           "objects": sum(int(record["metrics"].get("objects") or 0) for record in scatters),
                           "seconds_per_operation": scatter_seconds / len(scatters)}
        if scatter_shells > 0 and scatter_seconds > 0:
            scatter_metrics["seconds_per_shell"] = scatter_seconds / scatter_shells
        _write({"key": SCATTER_KEY, "n": total, "total": total, "label": "All scatters", "skipped": None,
                "settings": settings, "metrics": scatter_metrics, "seconds": scatter_seconds})


# ---------------------------------------------------------------------
# The command line's overall bar (PERF.55)
# ---------------------------------------------------------------------

OVERALL_MARGIN = 1.05
"""float: A run that outlasts its estimate keeps its bar at 1/this (95%) until it ends, so the bar never
claims to be done early."""


def _start_overall(stages, args, command):
    """The overall bar's clock and estimate for a run that does more than one stage, `None` for any other or
    for a run the web interface started (its job page draws the overall bar). The estimate is the stored time
    of each stage that will run (`estimate_seconds`), `None` when one has no record."""
    running = [stage for stage in stages if not stage.skip]
    if len(running) < 2 or os.environ.get(progress_file.ENV_VAR):
        return None
    estimate = None
    if args is not None and _recording_on():
        try:
            from planetgen.db import store
            config = store.mysql_config_from_args(args)
            conn = store.get_control_connection(store.control_mysql_config(config))
            try:
                estimate = estimate_seconds(conn, command or "", args)
            finally:
                conn.close()
        except Exception as exc:  # noqa: BLE001 -- an estimate is a nicety
            log.debug(f"Overall bar: no estimate ({exc}).")
    return {"t0": time.perf_counter(), "estimate": estimate}


def overall_progress(now=None):
    """
    `(description, completed, total)` of the command line's overall bar, in seconds: the time since the run's
    stages began against the stored time of all of them (`total` is `None` without one, and grows to keep the
    bar under full when the run outlasts its estimate), or `None` when there is no overall bar.
    """
    overall = _state["overall"]
    if overall is None:
        return None
    elapsed = (time.perf_counter() if now is None else now) - overall["t0"]
    description = f"Whole job (stage {min(max(_state['next'], 1), len(_state['stages']))} of {len(_state['stages'])})"
    estimate = overall["estimate"]
    if estimate is None:
        estimate = _average_estimate(elapsed)
    if estimate is None:
        return description, elapsed, None
    return description, elapsed, max(estimate, elapsed * OVERALL_MARGIN)


def _average_estimate(elapsed):
    """
    The whole run's seconds when no stored time covers every stage: the stages finished so far took
    `job_seconds` between them, so each stage still to run is taken to need their average, and the stage running
    now at least that (or what it has taken, if longer). So the bar's time left is longer than the running
    stage's own except in the last stage. `None` before any stage has finished.
    """
    ran = _state["ran"]
    if ran < 1:
        return None
    to_run = sum(1 for stage in _state["stages"] if not stage.skip)
    average = _state["job_seconds"] / ran
    current = max(elapsed - _state["job_seconds"], 0.0)
    return _state["job_seconds"] + max(average, current) + average * max(to_run - ran - 1, 0)


def attach_overall(progress):
    """
    Adds the overall bar to a started generation display (`run_common._ReportingProgress.start`) and keeps it
    moving once a second until the display stops (`progress.overall_stop`); does nothing without an overall bar
    or on a display that draws nothing. Returns the task id or `None`.
    """
    view = overall_progress()
    if view is None or getattr(progress, "disable", False):
        return None
    import threading
    description, completed, total = view
    task = progress.add_task(description, total=total, completed=completed, percent=True, prior=1.0)
    stop = threading.Event()

    def tick():
        while not stop.wait(1.0):
            now = overall_progress()
            if now is not None:
                progress.update(task, description=now[0], completed=now[1], total=now[2])

    threading.Thread(target=tick, name="overall-bar", daemon=True).start()
    progress.overall_stop = stop
    return task


def _numbering():
    """`(offset, total)` the log numbers this run's stages by: a web job's runner sets both, so a step's stages
    carry on from its earlier steps' and the total is the whole job's; a run on its own counts from 1 to its own
    last stage."""
    try:
        offset = int(os.environ.get(STAGE_OFFSET_ENV, ""))
        total = int(os.environ.get(STAGE_TOTAL_ENV, ""))
    except ValueError:
        return 0, len(_state["stages"])
    return (offset, total) if offset >= 0 and total >= offset + len(_state["stages"]) else (0, len(_state["stages"]))


def begin(stages, args=None, command=None, key_args=None):
    """Starts a run's stage list. With `args`, each stage's time is stored with the settings it ran with (PERF.56).
    `key_args` gives the arguments of the stages that have their own, by key (a run of several commands, PERF.203: the
    settings of a stage are the ones its own command was given). It logs nothing: each stage is logged once, when it is
    entered (`_announce`)."""
    _state.update(stages=list(stages), next=0, args=args, command=command, open=None, job_shells=0,
                  job_seconds=0.0, overall=_start_overall(stages, args, command), ran=0, scatters=[],
                  key_args=dict(key_args or {}))


def args_of(key):
    """The arguments the settings of stage `key` come from: its own (`begin`'s `key_args`), else the run's."""
    return _state["key_args"].get(key, _state["args"])


def _under_banner():
    """Whether the job runner's banner ("=== Stage 4 of 13: label ===") already names the one
    stage this run has: a step of a web job holds the stages of its command, and the banner says which; a step of
    one stage is that stage."""
    return len(_state["stages"]) == 1 and os.environ.get(STAGE_OFFSET_ENV, "") != ""


def _announce(index, skipped=None, queued=False):
    _close_open()
    stage = _state["stages"][index]
    total = len(_state["stages"])
    record = {"key": stage.key, "n": index + 1, "total": total, "label": stage.label, "skipped": skipped,
              "settings": stage_settings(stage.key, args_of(stage.key)) if _state["args"] is not None else {},
              "metrics": {}, "t0": time.perf_counter()}
    if queued:
        record["queued"] = True
    if skipped:
        record["seconds"] = 0.0
        _write(record)
    else:
        _state["open"] = record
    offset, job_total = _numbering()
    if skipped:
        mark = logrecord.stage_record(logrecord.STAGE_SKIPPED, "stages", offset + index + 1, job_total, stage.label,
                                      reason=skipped)
    elif not _under_banner():
        mark = logrecord.stage_record(logrecord.STAGE_START, "stages", offset + index + 1, job_total, stage.label)
    else:
        mark = None
    if mark is not None:        # the console draws it as the one divider, the record file holds it with its fields
        log.normal(mark.message, extra={"fields": mark.fields})
    progress_file.set_stage(index + 1, total, stage.label, skipped)


def enter(key, queued=False):
    """
    Enters the stage `key`, first announcing as skipped every stage before it
    that was not entered. A stage already entered, or a key the run's list
    does not have (another caller of the same function), does nothing. `queued` is for a stage whose row the
    queue stores (`row_for`): the stage still counts in the run's totals, but this process writes no row for it.
    """
    keys = [stage.key for stage in _state["stages"]]
    if key not in keys or keys.index(key) < _state["next"]:
        return
    index = keys.index(key)
    for earlier in range(_state["next"], index):
        _announce(earlier, _state["stages"][earlier].skip or "not needed")
    _announce(index, queued=queued)
    _state["next"] = index + 1


def skip(key, reason):
    """Enters the stage `key` as skipped for `reason` found while running (nothing to do)."""
    keys = [stage.key for stage in _state["stages"]]
    if key not in keys or keys.index(key) < _state["next"]:
        return
    index = keys.index(key)
    for earlier in range(_state["next"], index):
        _announce(earlier, _state["stages"][earlier].skip or "not needed")
    _announce(index, reason)
    _state["next"] = index + 1


def finish(reason=None):
    """Announces every stage the run never entered as skipped (`reason`, else its own), and ends the list."""
    for index in range(_state["next"], len(_state["stages"])):
        _announce(index, _state["stages"][index].skip or reason or "not needed")
    _close_open()
    _write_job()
    _state.update(stages=[], next=0, args=None, command=None, open=None, job_shells=0, job_seconds=0.0,
                  overall=None, ran=0, scatters=[], key_args={})


def fail(reason):
    """
    Ends a run that stopped on an error or a cancel: the stage open now says it failed after the time it took, with `reason` (the same
    record as a finished stage's, `logrecord.STAGE_FAILED`), no time is stored for it, and the list is dropped. Does nothing when no run
    is open.
    """
    record = _state["open"]
    if record is not None:
        _say_end(record, time.perf_counter() - record["t0"], reason=str(reason) or "failed")
    _state.update(stages=[], next=0, args=None, command=None, open=None, job_shells=0, job_seconds=0.0,
                  overall=None, ran=0, scatters=[], key_args={})


def stage_estimates(conn, command, args):
    """
    The expected seconds of each stage of a `galaxy` or `plan` run with `args`, in the stages' order: the stored
    time of the stage run with the settings it runs with (`stage_settings`, `generation.stats.stage_seconds`,
    PERF.56), `0.0` for a stage that will be skipped, `None` for one with no stored time. `[]` for any other
    command.
    """
    from planetgen.generation import stats

    return [0.0 if stage.skip else stats.stage_seconds(conn, stage.key, stage_settings(stage.key, args))
            for stage in for_command(command, args)]


def estimate_seconds(conn, command, args):
    """
    How long a `galaxy` or `plan` run with `args` is expected to take: the sum of `stage_estimates`, `None`
    unless every stage that will run has a stored time, so the overall bar says "about" only when it knows.
    """
    found = stage_estimates(conn, command, args)
    if not found or any(seconds is None for seconds in found):
        return None
    return sum(found)
