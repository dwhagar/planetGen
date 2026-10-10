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

from planetgen.queue import progress_file
from planetgen.util import log

SECTOR_MODES = ("block", "span", "column", "shell", "slot", "ring", "center_sector")

_state = {"stages": [], "next": 0}


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


def begin(stages):
    """Starts a run's stage list and logs it (skipped ones with their reasons)."""
    _state["stages"], _state["next"] = list(stages), 0
    if not stages:
        return
    parts = [f"{n}. {s.label}" + (f" (skipped: {s.skip})" if s.skip else "") for n, s in enumerate(stages, start=1)]
    log.normal(f"{len(stages)} stage{'s' if len(stages) != 1 else ''}: " + "; ".join(parts) + ".")


def _announce(index, skipped=None):
    stage = _state["stages"][index]
    total = len(_state["stages"])
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
    _state["stages"], _state["next"] = [], 0
