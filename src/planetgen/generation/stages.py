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

_state = {"stages": [], "next": 0, "args": None, "command": None, "open": None}
"""The run's stage list, the next one to enter, its arguments (for the settings each stage is recorded
with), its command, and the stage open now (`_open`)."""


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
    elif key == "mass":
        settings["mass_limit_sol"] = mass
    elif key == "luminosity":
        settings["mass_limit_sol"], settings["luminosity_floor_sol"] = mass, floor
    elif key == "phenomena":
        settings["mass_limit_sol"] = mass
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
        _write(record)


def begin(stages, args=None, command=None):
    """Starts a run's stage list and logs it (skipped ones with their reasons). With `args`, each stage's time
    is stored with the settings it ran with (PERF.56)."""
    _state.update(stages=list(stages), next=0, args=args, command=command, open=None)
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
    _state.update(stages=[], next=0, args=None, command=None, open=None)


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
