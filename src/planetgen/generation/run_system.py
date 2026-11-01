# planetgen/generation/run_system.py

"""
Run System
==========

The `system` command: one star system, built from `SystemConfig` (the
tri-state `+name`/`-name` options and `--system-file`). The sector
command builds each of its systems through `build_system_config`.
"""

import json
import sys

from planetgen.db import store
from planetgen import tuning as program_constants
from planetgen.util import log
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen import describe
from planetgen.db import render
from planetgen.db.render import render_star_system
from planetgen.generation import followed, run_common
from planetgen.queue import commands, redisqueue, runrecord, units, wrapper
from planetgen.util import draw


DOCUMENT_FORMAT = describe.DOCUMENT_FORMAT

OUTPUT_FORMATS = (*render.FORMATS, DOCUMENT_FORMAT)
"""tuple: The formats `--output` writes: the text formats of `db.render`, then the Document itself."""

# Tri-state options settable with `+name` (force True) / `-name` (force False).
# Each maps directly onto a `SystemConfig` attribute of the same name.
TRISTATE_OPTIONS = [
    ("habitable_world", "HABITABLE_WORLD", "a habitable world"),
    ("asteroid_belt", "ASTEROID_BELT", "an asteroid belt"),
    ("comets", "COMETS", "a star-bound comet"),
    ("large_star", "LARGE_STAR", "a large, massive star"),
    ("moons", "MOONS", "moons on every planet"),
    ("max_planets", "MAX_PLANETS", "the maximum number of orbital objects"),
    ("intelligent_life", "INTELLIGENT_LIFE", "intelligent life on a planet"),
    ("binary_system", "BINARY_SYSTEM", "a binary star system"),
    ("wide_binary", "WIDE_BINARY", "an S-type (wide) binary instead of a P-type (close) one"),
    ("planets", "PLANETS", "at least one planet or asteroid belt"),
]


ODDS_OPTIONS = ("habitable_world", "intelligent_life")
"""tuple: The options that raise the odds of their outcome and force nothing (GEN.1130): the system may come out without it and says so."""


def load_system_file(path):
    """
    Loads a system generation specification from a JSON file.

    The JSON file may contain any of the following keys, each corresponding
    to a `SystemConfig` attribute of the same name (see that class's
    docstrings for details): `star_type`, `name`, `age`, `num_orbits`,
    `slots`, `habitable_world`, `asteroid_belt`, `comets`, `large_star`, `moons`,
    `max_planets`, `intelligent_life`, `binary_system`, `wide_binary`,
    `stars`, `max_stars`, `planets`, `flavor_chance_system`, `flavor_chance_planet`,
    `max_planet_flavor`, `output`.

    `slots`, if present, is a list whose entries are either `null` or an
    object with a required `"type"` ("planet" or "asteroid_belt") and, for
    planets, optional `"planet_class"` (e.g. "M") and `"moons"` (an exact
    moon count) keys, e.g.:

        {
          "star_type": "G2V",
          "num_orbits": 5,
          "habitable_world": true,
          "slots": [
            {"type": "planet", "planet_class": "M", "moons": 1},
            {"type": "asteroid_belt"},
            null,
            {"type": "planet", "planet_class": "J", "moons": 4},
            null
          ]
        }

    Args:
        path (str): Path to the JSON system specification file.

    Returns:
        dict: The parsed JSON content.
    """
    import json
    with open(path, 'r') as f:
        return json.load(f)


def apply_system_file(system_config, data):
    """
    Applies a system specification (as loaded by `load_system_file`) onto a
    `SystemConfig` instance.

    Args:
        system_config (SystemConfig): The config object to update in place.
        data (dict): The parsed JSON system specification.
    """
    simple_keys = [
        "star_type", "name", "age", "num_orbits", "slots",
        "habitable_world", "asteroid_belt", "comets", "large_star", "moons",
        "max_planets", "intelligent_life", "binary_system", "wide_binary",
        "stars", "max_stars", "planets",
    ]
    key_to_attr = {key: key.upper() for key in simple_keys}

    for key, attr in key_to_attr.items():
        if key in data:
            setattr(system_config, attr, data[key])

    if "flavor_chance_system" in data and data["flavor_chance_system"] is not None:
        program_constants.FLAVOR_CHANCE_SYSTEM = data["flavor_chance_system"]

    if "flavor_chance_planet" in data and data["flavor_chance_planet"] is not None:
        program_constants.FLAVOR_CHANCE_PLANET = data["flavor_chance_planet"]

    if data.get("max_planet_flavor"):
        program_constants.MAX_FLAVOR_TOTAL = 99
        program_constants.FLAVOR_CHANCE_PLANET = 1


def build_system_config(args):
    """
    Builds a fully-configured `SystemConfig` from parsed command-line
    arguments: a `--system-file` JSON spec (if given) sets the baseline,
    then any explicitly-given tri-state or value option on the command
    line overrides it, then the cross-option normalization rules
    (intelligent life implies a habitable world; a habitable world plus
    an asteroid belt implies a large star) are applied last.

    Shared by the `system` and `sector` subcommands (`iter_sector_configs`
    calls this once per system in the sector), so the two can never
    silently drift apart on what a given set of options actually means.

    Args:
        args (argparse.Namespace): Parsed arguments (from the `system`
            subcommand, or a namespace shaped the same way, as
            `iter_sector_configs` builds for each system in a sector).

    Returns:
        SystemConfig: The resolved configuration.

    Raises:
        SystemExit: If a habitable world and an asteroid belt are both
                   forced while a large star is explicitly forbidden (see
                   TRISTATE_OPTIONS's `large_star`) -- that combination has
                   no room to be satisfied.
    """
    system_config = SystemConfig()

    if args.system_file:
        file_data = load_system_file(args.system_file)
        apply_system_file(system_config, file_data)

    # Command-line tri-state options override anything set by --system-file.
    for name, attr, _description in TRISTATE_OPTIONS:
        value = getattr(args, name)
        if value is not None:
            setattr(system_config, attr, value)

    # Command-line value options override anything set by --system-file.
    if args.star_type is not None:
        system_config.STAR_TYPE = args.star_type
    if args.num_orbits is not None:
        system_config.NUM_ORBITS = args.num_orbits
    if args.name is not None:
        system_config.NAME = args.name
    if args.age is not None:
        system_config.AGE = args.age
    # `--stars` is for a single system only (sector and galaxy parsers refuse it), `--max-stars` is shared (GEN.267).
    if getattr(args, "stars", None) is not None:
        system_config.STARS = args.stars
    if getattr(args, "max_stars", None) is not None:
        system_config.MAX_STARS = args.max_stars

    # A sector or galaxy run's --prevalence (GEN.52), the last one given
    # for a feature winning.
    for feature, percent in getattr(args, "prevalence", None) or ():
        system_config.PREVALENCE[feature] = percent

    if args.flavor_chance_system is not None:
        program_constants.FLAVOR_CHANCE_SYSTEM = args.flavor_chance_system

    if args.flavor_chance_planet is not None:
        program_constants.FLAVOR_CHANCE_PLANET = args.flavor_chance_planet

    if args.max_planet_flavor:
        program_constants.MAX_FLAVOR_TOTAL = 99
        program_constants.FLAVOR_CHANCE_PLANET = 1

    # Intelligent life raises the odds of a habitable world as well; -intelligent_life leaves them alone.
    if system_config.INTELLIGENT_LIFE is not None:
        system_config.HABITABLE_WORLD = True

    if system_config.HABITABLE_WORLD is True and system_config.ASTEROID_BELT is True:
        if system_config.LARGE_STAR is False:
            log.error("Error: forcing both a habitable world and an asteroid belt requires a large "
                      "star; -large_star cannot be combined with +habitable_world and +asteroid_belt.")
            raise SystemExit(1)
        system_config.LARGE_STAR = True

    return system_config


def max_stars_problem(max_stars):
    """
    Why `max_stars` (`--max-stars`, the system file's `max_stars`) is not allowed, or `None`: it is a whole number from 1
    to `tuning.MAX_STARS_CEILING`. `None` (left out) is fine.
    """
    if max_stars is None:
        return None
    if isinstance(max_stars, bool) or not isinstance(max_stars, int) or not 1 <= max_stars <= program_constants.MAX_STARS_CEILING:
        return f"max_stars {max_stars!r} must be a whole number from 1 to {program_constants.MAX_STARS_CEILING}."
    return None


def stars_problem(stars, max_stars, binary_system, wide_binary):
    """
    Why a forced star count (`--stars N`, GEN.267) cannot be met, or `None`: `stars` is not a whole number of at least 1,
    is more than the cap in force (`max_stars`, else `tuning.MAX_STARS_DEFAULT`) or than any system may hold
    (`tuning.MAX_STARS_CEILING`), or contradicts the pair options. One star is a single star, so `+binary_system` and
    `+wide_binary` (a companion of the wide kind) have none to apply to; two or more stars are not a single star, so
    `-binary_system` contradicts them. `+binary_system` with three or more stars stands: the system is a binary only in
    having a companion, and the member-tree host (GEN.630) decides what the rest of the tree is. `-wide_binary` with one
    star stands: no companion is of the wide kind.

    Args:
        stars (int or None): The forced count, `None` when left to chance.
        max_stars (int or None): The cap, `None` for the default.
        binary_system, wide_binary (bool or None): The tri-state options.

    Returns:
        str or None: The message.
    """
    if stars is None:
        return None
    if isinstance(stars, bool) or not isinstance(stars, int) or stars < 1:
        return f"stars {stars!r} must be a whole number of at least 1."
    cap = program_constants.MAX_STARS_DEFAULT if max_stars is None else max_stars
    if stars > program_constants.MAX_STARS_CEILING:
        return f"stars {stars} is more than any system may hold ({program_constants.MAX_STARS_CEILING})."
    if stars > cap:
        return f"stars {stars} is more than the cap of {cap} stars; raise the cap with --max-stars (at most {program_constants.MAX_STARS_CEILING})."
    if stars == 1 and binary_system is True:
        return "stars 1 is a single star, which +binary_system contradicts."
    if stars == 1 and wide_binary is True:
        return "stars 1 is a single star, which has no companion for +wide_binary."
    if stars > 1 and binary_system is False:
        return f"stars {stars} is more than one star, which -binary_system contradicts."
    return None


def single_system_conflict(system_config):
    """
    The contradiction in a single system's resolved options, or `None`:
    `-planets` (no planets or belts at all) with anything that needs a
    planet or a belt (GEN.50), and a forced star count the cap or the pair
    options rule out (`stars_problem`). `validate_system_args` already
    rejects these on the command line; this also catches them from a
    `--system-file`.

    Args:
        system_config (SystemConfig): The resolved configuration.

    Returns:
        str or None: The error message.
    """
    problem = max_stars_problem(system_config.MAX_STARS) or stars_problem(
        system_config.STARS, system_config.MAX_STARS, system_config.BINARY_SYSTEM, system_config.WIDE_BINARY)
    if problem:
        return problem
    if system_config.PLANETS is not False:
        return None
    forced = [flag for attr, flag in (("MOONS", "+moons"), ("MAX_PLANETS", "+max_planets"),
                                      ("HABITABLE_WORLD", "+habitable_world"),
                                      ("ASTEROID_BELT", "+asteroid_belt"))
              if getattr(system_config, attr) is True]
    if system_config.NUM_ORBITS:
        forced.append("num_orbits")
    if system_config.SLOTS:
        forced.append("slots")
    if not forced:
        return None
    return f"-planets cannot be combined with {', '.join(forced)}."


class SystemFollower(followed.Follower):
    """
    The one system of a `system` run as its one stage (PERF.119, PERF.202): the `generate-system` command's one unit
    generates it, and this side keeps what the unit says (`result`) for the caller to log or write.
    """

    key = "system"
    label = "Generate a system"

    def __init__(self, args):
        super().__init__(args)
        self.result = None

    def params(self):
        return commands.build_steps("system", self.args)[0].params

    def on_result(self, unit, result):
        self.result = result


def run_system(args):
    """
    Generates one star system and saves it to the database, or, with
    `--output`, writes its page to a file (or stdout) and touches no
    database at all -- the admin site's one-off system page runs it that
    way (`src/planetgen/web/system_page.py`).

    The run is one command through the wrapper (PERF.119, `generate-system`, one unit): a worker generates the system
    and this process says what it made. With `--output` nothing is saved, so the run is a maintenance one and does not
    take the generation slot.

    A forced option the system can't meet is never saved silently
    (GEN.49): contradictory options are refused up front, and a system
    still missing a forced body after `SINGLE_SYSTEM_GENERATION_ATTEMPTS`
    whole systems ends the run with an error, saving nothing.

    Args:
        args (argparse.Namespace): Validated arguments (`command ==
            "system"`).

    Raises:
        SystemExit: On contradictory options, or when no attempt met every
            forced option.
    """
    system_config = build_system_config(args)
    conflict = single_system_conflict(system_config)
    if conflict:
        log.error(f"Error: {conflict}")
        raise SystemExit(1)

    mysql_config = store.mysql_config_from_args(args)
    args.database = wrapper.database_args(mysql_config)
    args.run_seed = draw.getrandbits(128)
    follower = SystemFollower(args)
    try:
        followed.run_followed(args, [follower], follower.label, workers=1, galaxy=False,
                              kind=runrecord.MAINTENANCE if args.output else runrecord.GENERATION)
    except (redisqueue.Unavailable, runrecord.Busy, units.RunStopped) as exc:
        log.error(f"error: the system could not be generated: {exc}")
        raise SystemExit(1)
    result = follower.result
    if result["warning"]:
        run_common._refuse_or_warn(args, result["warning"])
    if args.output:
        text = result["text"]
        if args.output == "-":
            sys.stdout.write(text if text.endswith("\n") else text + "\n")
        else:
            with open(args.output, "w", encoding="utf-8") as f:
                f.write(text)
            log.normal(f"Wrote system '{result['name']}' to {args.output} (not saved to the database).")
        return
    run_common.RUN_COUNTS["systems"] += 1
    log.normal(f"Saved system '{result['name']}' to the database (star_system_id={result['star_system_id']}, "
               f"{mysql_config.database}@{mysql_config.host}:{mysql_config.port}).")


def expand_system(params):
    """The cut of a `system` run (PERF.119; `generate-system`): one unit of cost 1."""
    return commands.Expansion([({"database": params["database"], "args": params["args"],
                                 "run_seed": params["run_seed"]}, 1.0)], {"systems": 1})


def system_unit(unit_args, report):
    """
    The unit of a `system` run: generates the system on a stream seeded from the run's (so it is the same at any
    worker count), and renders it (`--output`) or saves it.

    Returns:
        dict: `warning` (the message of a system that kept no forced body, which the starter refuses or logs by
            `--strict`, else `None`), `name`, and either `text` (the page, with `--output`) or `star_system_id`.
    """
    from planetgen.generation import run_galaxy
    args = run_galaxy._sector_args_of(unit_args)
    with draw.bound(unit_args["run_seed"]):
        system_config = build_system_config(args)
        attempts = program_constants.SINGLE_SYSTEM_GENERATION_ATTEMPTS
        warning = None
        for attempt in range(1, attempts + 1):
            system = StarSystem(system_config=system_config)
            if not system.unmet_requirements:
                break
        else:
            star = f"--star-type {system_config.STAR_TYPE}" if system_config.STAR_TYPE else "the drawn star"
            # GEN.81: the last try is kept, without the body it had no room for.
            warning = (f"No system with {' and '.join(system.unmet_requirements)} came out of {attempts} tries for "
                       f"{star}; the last one is kept without it. Hot, short-lived and very large stars often have "
                       f"no room for one.")
        report(1, 1)
        if warning and run_common._strict(args):
            return {"warning": warning, "name": system.name, "text": None, "star_system_id": None}
        if args.output:
            text = (json.dumps(describe.describe(system)) if args.format == DOCUMENT_FORMAT
                    else render_star_system(system, args.format))
            return {"warning": warning, "name": system.name, "text": text, "star_system_id": None}
        star_system_id = store.save_system(system, system_config, config=store.mysql_config_from_args(args))
        return {"warning": warning, "name": system.name, "text": None, "star_system_id": star_system_id}
