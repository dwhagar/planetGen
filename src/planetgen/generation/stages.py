# planetgen/generation/stages.py

"""
The stages of a generation command, numbered and announced as the run reaches
them (UX.89).

A command that does several things in turn (`galaxy`, `plan`) declares its
stages up front with `begin`; each is entered with `enter(key)` as the run
reaches it. A stage the run will not do (the scatter without `--then-scatter`,
the backfill with `--backfill-from none`, no layers to scatter) is still
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

from planetgen.queue import progress_file
from planetgen.util import log

SECTOR_MODES = ("block", "span", "column", "shell", "slot", "ring", "center_sector")

JOB_KEY = "job"
"""str: The `stage_key` the whole-job row is stored under (`finish`)."""

_state = {"stages": [], "next": 0, "args": None, "command": None, "open": None,
          "job_layers": 0, "job_seconds": 0.0, "overall": None}
"""The run's stage list, the next one to enter, its arguments (for the settings each stage is recorded
with), its command, the stage open now (`_open`), the layers and seconds of the stages finished so far (the
whole-job stat) and the command-line overall bar's clock and estimate (`overall`)."""


class Stage:
    """One stage: `key`, `label` and `skip` (the reason it will not run, or `None`)."""

    def __init__(self, key, label, skip=None):
        self.key, self.label, self.skip = key, label, skip

    def as_dict(self):
        return {"label": self.label, "skipped": self.skip}


def _random_start(args):
    return not any(getattr(args, name, None) not in (None, False) for name in SECTOR_MODES)


def galaxy_stages(args):
    """The stages of one `planetgen galaxy` run, in order."""
    then_scatter = getattr(args, "then_scatter", False)
    no_scatter = None if then_scatter else "the star and phenomena scatter was not asked for"
    if _random_start(args):
        sectors = [Stage("start", "Generate the starting sector"), Stage("neighborhood", "Generate the neighborhood")]
    else:
        sectors = [Stage("sectors", "Generate the sectors")]
    return sectors + [
        Stage("link", "Link the neighbours"),
        Stage("mass", "Scatter the massive stars", no_scatter),
        Stage("luminosity", "Scatter the bright stars", no_scatter),
        Stage("phenomena", "Scatter the phenomena", no_scatter),
        Stage("backfill", "Scatter the massive stars from the neighborhood",
              "backfill was turned off (--backfill-from none)"
              if (getattr(args, "backfill_from", "edge") or "edge") == "none" else None),
        Stage("population", "Run the population pass",
              None if getattr(args, "population", False) else "the population pass was not asked for"),
        Stage("settle", "Save the sector paths",
              "turned off with --no-settle" if getattr(args, "no_settle", False) else None),
    ]


def plan_stages(args):
    """The stages of one `planetgen plan` run, in order."""
    if getattr(args, "bright_stars_down_to", None) is not None:
        return [Stage("band", "Add the dimmer band of bright stars")]
    redo = getattr(args, "redo_scatters", None)
    if redo:
        not_chosen = "not chosen to be redone"
        new_limit = getattr(args, "phenomenon_min_mass", None) is not None
        luminosity = None if "luminosity" in redo else (
            None if "mass" in redo and new_limit else not_chosen + (
                " (it is redone if the mass limit changes)" if "mass" in redo else ""))
        return [Stage("phenomena", "Scatter the phenomena", None if "phenomena" in redo else not_chosen),
                Stage("mass", "Scatter the massive stars", None if "mass" in redo else not_chosen),
                Stage("luminosity", "Scatter the bright stars", luminosity)]
    if getattr(args, "phenomena_only", False):
        return [Stage("phenomena", "Scatter the phenomena")]
    if getattr(args, "bright_stars_only", False):
        return [Stage("mass", "Scatter the massive stars"), Stage("luminosity", "Scatter the bright stars")]
    stages = [Stage("skeleton", "Plan the galaxy")]
    if not getattr(args, "no_bright_stars", False):
        stages += [Stage("phenomena", "Scatter the phenomena"), Stage("mass", "Scatter the massive stars"),
                   Stage("luminosity", "Scatter the bright stars")]
    return stages


def for_command(command, args):
    """The stages of `command` (`"galaxy"` or `"plan"`) for parsed `args`, `[]` for any other."""
    if command == "galaxy":
        return galaxy_stages(args)
    if command == "plan":
        return plan_stages(args)
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
        for name in ("radius_pc", "neighborhoods", "max_ring", "min_start_density"):
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
    elif key == "backfill":
        settings["backfill_from"] = getattr(args, "backfill_from", "edge") or "edge"
    elif key == "skeleton":
        for name in ("max_ring", "seed"):
            value = getattr(args, name, None)
            if value is not None:
                settings[name] = str(value) if name == "seed" else value
    return settings


def note(**metrics):
    """Adds what the stage open now did (layers visited, layers that placed something, objects ...)."""
    if _state["open"] is not None:
        _state["open"]["metrics"].update(metrics)


def _recording_on():
    from planetgen.generation import stats as generation_stats
    return (_state["args"] is not None and
            os.environ.get(generation_stats.STATS_ENV_VAR, "1").strip().lower() not in ("0", "off", "no", "false"))


def _write(record):
    """Stores one stage's row in the control database; a failure only warns in the debug log."""
    if not _recording_on():
        return
    try:
        from planetgen.db import store
        from planetgen.galaxy.version_key import version_key
        config = store.mysql_config_from_args(_state["args"])
        conn = store.get_control_connection(store.control_mysql_config(config))
        try:
            conn.execute(
                "INSERT INTO generation_stage_runs (database_name, command, stage_key, stage_n, stage_total, label,"
                " skipped, skip_reason, started_at, finished_at, seconds, workers, settings, metrics, version_key)"
                " VALUES (?, ?, ?, ?, ?, ?, ?, ?, DATE_SUB(NOW(6), INTERVAL ? MICROSECOND), NOW(6), ?, ?, ?, ?, ?)",
                (config.database, _state["command"] or "", record["key"], record["n"], record["total"],
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


def _close_open():
    """Records the stage open now, with the time it took."""
    record = _state["open"]
    _state["open"] = None
    if record is not None:
        record["seconds"] = time.perf_counter() - record["t0"]
        _state["job_seconds"] += record["seconds"]
        _state["job_layers"] += int(record["metrics"].get("layers") or 0)
        _write(record)


def _write_job():
    """
    Stores the whole job's layers per second (PERF.55): every layer the scatter stages visited, the ones that
    placed nothing too, over the seconds of all the stages the run did. One row under `JOB_KEY`, with the
    settings the run had, for `stats.job_layers_per_second`. Nothing when no stage visited a layer.
    """
    layers, seconds = _state["job_layers"], _state["job_seconds"]
    if layers <= 0 or seconds <= 0 or not _state["stages"]:
        return
    total = len(_state["stages"])
    settings = stage_settings(JOB_KEY, _state["args"]) if _state["args"] is not None else {}
    _write({"key": JOB_KEY, "n": total, "total": total, "label": "Whole job", "skipped": None, "settings": settings,
            "metrics": {"layers": layers, "layers_per_second": layers / seconds}, "seconds": seconds})


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
    if overall["estimate"] is None:
        return description, elapsed, None
    total = max(overall["estimate"], elapsed * OVERALL_MARGIN)
    return description, elapsed, total


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


def begin(stages, args=None, command=None):
    """Starts a run's stage list and logs it (skipped ones with their reasons). With `args`, each stage's time
    is stored with the settings it ran with (PERF.56)."""
    _state.update(stages=list(stages), next=0, args=args, command=command, open=None, job_layers=0,
                  job_seconds=0.0, overall=_start_overall(stages, args, command))
    if not stages:
        return
    parts = [f"{n}. {s.label}" + (f" (skipped: {s.skip})" if s.skip else "") for n, s in enumerate(stages, start=1)]
    log.normal(f"{len(stages)} stage{'s' if len(stages) != 1 else ''}: " + "; ".join(parts) + ".")


def _announce(index, skipped=None):
    _close_open()
    stage = _state["stages"][index]
    total = len(_state["stages"])
    record = {"key": stage.key, "n": index + 1, "total": total, "label": stage.label, "skipped": skipped,
              "settings": stage_settings(stage.key, _state["args"]) if _state["args"] is not None else {},
              "metrics": {}, "t0": time.perf_counter()}
    if skipped:
        record["seconds"] = 0.0
        _write(record)
    else:
        _state["open"] = record
    if skipped:
        log.normal(f"Stage {index + 1} of {total}: {stage.label} -- skipped: {skipped}.")
    else:
        log.normal(f"Stage {index + 1} of {total}: {stage.label}.")
    progress_file.set_stage(index + 1, total, stage.label, skipped)


def enter(key):
    """
    Enters the stage `key`, first announcing as skipped every stage before it
    that was not entered. A stage already entered, or a key the run's list
    does not have (another caller of the same function), does nothing.
    """
    keys = [stage.key for stage in _state["stages"]]
    if key not in keys or keys.index(key) < _state["next"]:
        return
    index = keys.index(key)
    for earlier in range(_state["next"], index):
        _announce(earlier, _state["stages"][earlier].skip or "not needed")
    _announce(index)
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
    _state.update(stages=[], next=0, args=None, command=None, open=None, job_layers=0, job_seconds=0.0,
                  overall=None)


def estimate_seconds(conn, command, args):
    """
    How long a `galaxy` or `plan` run with `args` is expected to take: the
    stored time of each of its stages that will run, looked up by the settings
    that stage runs with (`stage_settings`, `generation.stats.stage_seconds`,
    PERF.56), skipped stages counting nothing. `None` unless every stage that
    will run has a stored time, so the overall bar says "at least" rather than
    quoting a guess.
    """
    from planetgen.generation import stats

    total = 0.0
    found = for_command(command, args)
    if not found:
        return None
    for stage in found:
        if stage.skip:
            continue
        seconds = stats.stage_seconds(conn, stage.key, stage_settings(stage.key, args))
        if seconds is None:
            return None
        total += seconds
    return total
