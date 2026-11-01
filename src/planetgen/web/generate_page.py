# planetgen/web/generate_page.py

"""
`/admin/generate`: galaxy generation from the browser, so an admin never
needs a terminal on the server. Four actions, each a background job
(`web/jobs.py`) running the same command-line tools an admin would type:

- New galaxy: wipe the database (`python3 -m planetgen.cli.reset --yes`), build the
  density skeleton (`planetgen plan`), then generate
  a first neighborhood around a random start and only then scatter the
  bright stars, leaving those sectors out (`planetgen sectors`, `link`, the scatters, `backfill` and `settle` as steps of
  one job, GEN.30, PERF.123).
- Plan: rebuild the density skeleton, then scatter the bright stars.
  The scatter is its own step so the job shows its progress bar and the
  count it placed; a checkbox on either form skips it.
- Rebuild bright stars: the scatter alone, on the stored plan (filled
  sectors are always left out).
- Add a dimmer step: keep the bright stars already placed and add only
  those from a lower level up to the current one
  (`planetgen scatter --band --down-to N`).
- Generate sectors: `planetgen sectors` in any of its modes (random
  start, a whole shell, around a sector, one address, ...). "Around a
  sector" names its center as a filled sector (searched by name or
  picked from `/admin/generate/sectors`'s list), a sector address, or a
  galaxy-frame position in pc. The bright-star backfill runs once the
  sectors are done, out from the run's edge (GEN.30, GEN.98).
- Reset: wipe the database only.

New galaxy and Reset delete every generated sector and system, so both
ask for the database name to be typed back, the same confirmation
`planetgen.cli.reset` asks for in a terminal.

The page shows the running job (step, progress bar, elapsed time and the
output in a terminal), kept live by `static/generatejobs.js` over a
Server-Sent Events stream (`/admin/generate/jobs/<id>/stream`, ADM.22;
`/admin/generate/status` is the polling fallback), plus the last few
jobs, each with its page at `/admin/generate/jobs/<id>` and its full
output to download at `/admin/generate/jobs/<id>/log`.

The Generate page is a Start screen with a card per task; each task has
its own page and a Review step before anything runs (`generate_flow.py`).

Admins only: a visitor who isn't logged in is sent to the login page,
and any POST or status request without an admin session gets a 403.
Every form carries `csrf_field()` (checked app-wide by `csrf.protect`).
"""

import contextlib
import contextvars
import json
import os
import threading
import time
from typing import Annotated
from urllib.parse import urlsplit

from flask import (Response, abort, current_app, jsonify, make_response, redirect, request, send_file,
                   stream_with_context, url_for)

from pydantic import Field, TypeAdapter, ValidationError

from planetgen.admin import cachebuild
from planetgen.queue import api_jobs, progress_rate, redisqueue, roughtime, runline, runrecord, runview
from planetgen.queue.commands import DEEP_CHECKS, STAGED_COMMANDS
from planetgen.web import limit_sliders
from planetgen.web.lib import apiclient, estimatecache
from planetgen.web.lib.fmt import utc_time_html
from planetgen.web.lib.pagination import PAGE_SIZE, clamp_page, page_offset, parse_page
from planetgen.admin import activity_log, updateguard
from planetgen.generation import finish_run, luminosity_floor, plan_outline, prevalence, stages, stats
from planetgen import tuning
from planetgen.db import dataversion, store
from planetgen.util import generator_epoch, log, stagemark
from planetgen.galaxy import bring_to_epoch, dials, drill, mass_chain, objectref
from planetgen.galaxy.spherical import sector_address_at
from planetgen.galaxy.span import Span, SpanError, parse_range
from planetgen.physics.units import ly_to_pc, pc_to_ly
from planetgen.util.format import counted, format_number
from planetgen.generation.limits import (
    MAX_GENERATE_LIMIT, MAX_GENERATE_NEIGHBORHOODS, MAX_GENERATE_RADIUS_PC, MAX_GENERATE_SHELL,
)

from . import bp, jobs
from .admin_pages import came_from, safe_next
from .helpers import crumb, current_admin, db_name, page_url, pager, render_page, trusted_html

# ---------------------------------------------------------------------
# Form fields
# ---------------------------------------------------------------------

PLAN_FIELDS = dials.LEVERS
"""tuple: The `plan` options the page offers: the rows of the lever table (`planetgen.galaxy.dials`), which also builds
`planetgen plan`'s flags, so the page and the command line cannot differ. Every row is Auto (blank) until given."""

MASS_LIMIT_LABEL = "Mass limit (solar masses)"
COMPACT_LIMIT_LABEL = "Neutron star and black hole limit (solar masses)"
COMPACT_STORAGE_WARNINGS = ((1.0, "about 1.2 billion rows (161 GB)"), (2.0, "about 390 million rows (54 GB)"),
                            (4.0, "about 220 million rows (30 GB)"), (6.0, "about 200 million rows (28 GB)"))
"""tuple: What each compact-object limit costs in the database (GEN.195, Boss 2026-10-10), beside the control."""
"""str: The mass-limit slider's label (GEN.183)."""


def _backfill_text():
    """The backfill steps in words ("down to 1 solar mass in the step of
    sectors around the generated ones, 2 in the next, ... and 8 in the
    fourth")."""
    masses = tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL
    ordinals = ("first", "second", "third", "fourth", "fifth", "sixth")
    parts = [f"{mass:g} in the {ordinals[index]}" for index, mass in enumerate(masses[1:], start=1)]
    first = f"down to {masses[0]:g} solar mass{'' if masses[0] == 1 else 'es'} in the first step out from the generated sectors"
    return ", ".join([first, *parts[:-1]]) + (f" and {parts[-1]} step" if parts else "")


BACKFILL_TEXT = _backfill_text()
"""str: How far down the bright-star backfill around the generated sectors
goes, step by step (GEN.187, `tuning.BRIGHT_STAR_BACKFILL_RING_MASSES_SOL`)."""

PLACE_TABS = (
    ("random", "A random populated spot", "Picks a spot at random where stars are dense enough and generates around it."),
    ("sector", "A named sector", "Find a filled sector by name, pick it from the list, or type its ID."),
    ("address", "A sector address", "Shell, band and wedge, as the Galaxy Map shows them. "
                                    "The sector there is generated first if it isn't yet."),
    ("position", "A position (pc)", "Galaxy-frame x, y and z in parsecs (the center is 0, 0, 0; "
                                    "z is height above the plane). The sector holding that point is the center."),
)
"""tuple: The tabs of the Add sectors place picker (`place`), `(value, label, hint)`: the random start and the three
ways to name the center of "Around a sector" (`center_by`)."""

GRID_CHOICES = (
    ("shell", "A whole shell", "Every not-yet-generated sector in one shell, optionally only one band of it."),
    ("span", "A range", "A range of shells, a range of bands, or one arc of wedges. Leave a box blank for all of it."),
)
"""tuple: The two ways to pick sectors by their grid address (`mode` with `what=grid`), `(value, label, hint)`."""


def where(form):
    """
    The `(mode, center_by)` the Add sectors form means. The page posts `what` (`place` or `grid`) and, for a place,
    the picked `place` tab: `random` is the random start, any other tab is "Around a sector" naming its center that way.
    A grid choice keeps `mode` (`shell` or `span`). A post without `what` (the Galaxy Map's dialog, the sector page's
    button) names `mode` and `center_by` itself.
    """
    place = form.get("place")
    if place and form.get("what") != "grid":
        return ("random", None) if place == "random" else ("center", place)
    return form.get("mode") or "random", form.get("center_by") or "sector"


FINDER_PAGE_SIZE = 50
"""int: Sectors per page in the "Around a sector" list (every web list
pages 50 rows)."""

CONFIRM_ACTIONS = frozenset({"new_galaxy", "new_plan", "reset"})

BAND_LABEL = "Add a dimmer step of bright stars"

MASS_SCATTER_LABEL = "Scatter the massive stars"

SCATTER_LABEL = "Scatter the bright stars"

PHENOMENA_LABEL = "Scatter the phenomena"

BRIGHT_THRESHOLD_LABEL = "Bright stars from (solar luminosities)"
"""str: The galaxy-wide scatter threshold field (GEN.30), on New galaxy,
Plan and Rebuild the bright stars."""

"""str: The step that pre-places every bright star
(`tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL` and up) on the plan."""

DISREGARD_FLAG = "--disregard-missing-data"
DEBUG_FLAG = "--debug"
"""str: Added to every step of a job whose form has the Debug box ticked: each command logs every decision to stdout,
which is the job log. Debug mode in `config.json` does the same for every job."""
"""str: The option a job's generating commands and its data check get when the form's "Disregard missing data files" box is ticked (GEN.599)."""

DATA_CHECK_LABEL = "Check the data"
"""str: Every generating job's first step (GEN.590): `python -m planetgen.cli.datasets
status`, which exits 1 while a downloaded data file is missing. A failure stops the job
there, before the reset or anything else is written."""


ESTIMATE_TIMEOUT_S = 120
"""int: How long the page waits for `planetgen sectors --estimate-only`."""

ESTIMATE_PREFIX = "ESTIMATE "
"""str: `generate.ESTIMATE_PREFIX`: the line the estimate is read from."""

ESTIMATE_CONFIRM_FIELD = "estimate_ok"
"""str: The form field a confirmed estimate's "Generate" button adds."""

GENERATE_ANYWAY_FIELD = "generate_anyway"
"""str: The form field the "Generate anyway" button adds to a run the
database disk has no room for (ADM.33). The job then starts like any
confirmed one (`planetgen galaxy` warns and goes ahead without `--strict`),
and the override goes in the activity log."""


class FormError(ValueError):
    """A form value the page rejects; the message is shown to the admin."""


def run_estimate(argv, env):
    """
    PERF.3: the size and time a Generate sectors job would take, from
    `planetgen sectors ... --estimate-only` (run here, while the admin
    waits, before anything is written).

    Args:
        argv (list): The job step's command line.
        env (dict): `jobs.mysql_env`.

    Returns:
        dict: `stats.Estimate.as_dict()` plus `what`.

    Raises:
        FormError: The estimate couldn't be made; the message says why
            (`planetgen`'s own last lines, e.g. a shell too large to
            generate without a limit).
    """
    done = api_jobs.command_and_wait(argv + ["--estimate-only"], {**os.environ, **env}, jobs.REPO_DIR,
                                     ESTIMATE_TIMEOUT_S)
    if done["timed_out"]:
        raise FormError("Working out the size and time took too long. Nothing was generated.")
    if done["error"]:
        raise FormError(f"The size and time couldn't be worked out: {done['error']}")
    for line in reversed(done["stdout"].splitlines()):
        if line.startswith(ESTIMATE_PREFIX):
            try:
                return json.loads(line[len(ESTIMATE_PREFIX):])
            except ValueError:
                break
    tail = [line.strip() for line in (done["stderr"] + "\n" + done["stdout"]).splitlines() if line.strip()]
    message = " ".join(tail[-3:]) if tail else f"the generator exited with {done['returncode']}"
    raise FormError(f"Nothing was generated: {message}")


def _estimate_versions(database):
    """
    What the cached estimate (`estimatecache`, PERF.109) was made under: the generator epoch, the galaxy's
    data version (every data change bumps it) and the speed statistics' fingerprint.

    Raises:
        Exception: Redis or the control database can't be asked (the estimate then runs uncached).
    """
    incarnation, found = dataversion.versions(database, [dataversion.COUNT_FIELD, dataversion.TAG_EPOCH])
    conn = store.get_control_connection(current_app.config["MYSQL_CONFIG"])
    try:
        speeds = stats.fingerprint(conn)
    finally:
        conn.close()
    return {"generator_epoch": generator_epoch.generator_epoch(), "data": [incarnation, found], "speeds": speeds}


_problems = contextvars.ContextVar("generate_form_problems", default=None)
"""ContextVar: The list `collecting_problems` is filling, or `None`."""


@contextlib.contextmanager
def collecting_problems():
    """
    Collect every bad field of a form instead of stopping at the first:
    `_number` records its message and returns `None`, and leaving the
    block raises one `FormError` listing them all.
    """
    found = []
    token = _problems.set(found)
    try:
        yield found
    finally:
        _problems.reset(token)
    if found:
        raise FormError(" ".join(found))


def _problem(message):
    found = _problems.get()
    if found is None:
        raise FormError(message)
    if message not in found:
        found.append(message)


def _number(form, name, label, kind, required=False, minimum=None, maximum=None, exclusive_max=False):
    raw = (form.get(name) or "").strip()
    if not raw:
        if required:
            _problem(f"{label} is required.")
        return None
    try:
        number = kind(raw)
    except ValueError:
        _problem(f"{label} must be {'a whole number' if kind is int else 'a number'}.")
        return None
    if number != number or number in (float("inf"), float("-inf")):
        _problem(f"{label} must be a finite number.")
        return None
    bounds = {}
    if minimum is not None:
        bounds["ge"] = minimum
    if maximum is not None:
        bounds["lt" if exclusive_max else "le"] = maximum
    try:
        return TypeAdapter(Annotated[kind, Field(**bounds)]).validate_python(number)
    except ValidationError:
        if minimum is not None and number < minimum:
            _problem(f"{label} must be at least {minimum:g}.")
        else:
            _problem(f"{label} must be {'less than' if exclusive_max else 'at most'} {maximum:g}.")
        return None


def compact_limit_argv(form):
    """`--compact-min-mass` from the form's neutron star and black hole limit
    (GEN.195): "star" when "Use the star mass setting" is ticked, else one of
    `tuning.COMPACT_MIN_MASS_PRESETS`; nothing when the form has neither.

    Raises:
        FormError: A value that isn't one of the presets.
    """
    if form.get("compact_use_star"):
        return ["--compact-min-mass", tuning.COMPACT_MIN_MASS_STAR]
    value = _number(form, "compact_min_mass", COMPACT_LIMIT_LABEL, float)
    if value is None:
        return []
    if not any(abs(value - preset) < 1e-9 for preset in tuning.COMPACT_MIN_MASS_PRESETS):
        _problem(f"{COMPACT_LIMIT_LABEL} must be one of "
                 + ", ".join(f"{preset:g}" for preset in tuning.COMPACT_MIN_MASS_PRESETS) + ".")
        return []
    return ["--compact-min-mass", f"{value:g}"]


LIMIT_OPTIONS = (
    (limit_sliders.CORE_FLOOR, "--core-min-luminosity"),
    (limit_sliders.CORE_MASS, "--core-mass-limit"),
    (limit_sliders.HALO_FLOOR, "--halo-min-luminosity"),
    (limit_sliders.CELL_CEILING, "--cell-star-ceiling"),
    (limit_sliders.FILL_CAP, "--fill-luminosity-cap"),
)
"""tuple: The limits a plan stores with the galaxy (UX.116): each slider and the `planetgen plan` option it sets."""


def limits_argv(form):
    """
    The plan options for the core, halo and dense-cell limits (UX.116): a number inside its slider's range, blank
    leaving the default. The core's two follow the galaxy-wide limits while "Follow the galaxy-wide limits" is
    ticked, and then give no option.

    Raises:
        FormError: A value that isn't a number inside the slider's range.
    """
    argv = []
    for slider, flag in LIMIT_OPTIONS:
        if slider.name.startswith("core_") and form.get("core_follow"):
            continue
        value = _number(form, slider.name, slider.label, float, minimum=slider.low, maximum=slider.high)
        if value is not None:
            argv += [flag, limit_sliders.plain(value)]
    return argv


def plan_argv(form):
    """`planetgen plan` arguments from the Plan fields (blank means
    the default) and the limits stored with the galaxy (`limits_argv`)."""
    return limits_argv(form) + lever_argv(form)


def lever_argv(form):
    """The lever table's flags from the Plan fields (blank means Auto), each checked against its range."""
    argv = []
    for row in PLAN_FIELDS:
        value = _lever_text(form, row)
        if value is not None:
            argv += [row.flag, value]
    return argv


def _lever_text(form, row):
    """The text of a lever's flag value from the form (blank is Auto, so `None`), checked against the row's range;
    a length is typed with its unit (`2600pc`)."""
    raw = (form.get(row.name) or "").strip()
    if not raw:
        return None
    try:
        value = row.convert(raw)
    except ValueError as exc:
        _problem(f"{row.label}: {exc}.")
        return None
    message = row.problem(value)
    if message:
        _problem(f"{row.label}: {message}.")
        return None
    return raw if row.kind == "length" else str(value)


def random_start_argv(form):
    """`planetgen sectors` arguments for a random start."""
    argv = []
    radius = _number(form, "radius_pc", "Radius (pc)", float, minimum=0.1, maximum=MAX_GENERATE_RADIUS_PC)
    if radius is not None:
        argv += ["--radius-pc", str(radius)]
    max_shell = _number(form, "max_shell", "Highest shell", int, minimum=0, maximum=MAX_GENERATE_SHELL)
    if max_shell is not None:
        argv += ["--max-shell", str(max_shell)]
    density = _number(form, "min_start_density", "Minimum start density", float, minimum=0.001)
    if density is not None:
        argv += ["--min-start-density", str(density)]
    if form.get("avoid_filled_space"):
        argv += ["--avoid-filled-space"]
    count = _number(form, "neighborhoods", "Neighborhoods", int, minimum=1, maximum=MAX_GENERATE_NEIGHBORHOODS)
    if count is not None and count > 1:
        argv += ["--neighborhoods", str(count)]
        gamma = _number(form, "neighborhood_gamma", "Density bias", float, minimum=0, maximum=10)
        if gamma:
            argv += ["--neighborhood-gamma", f"{gamma:g}"]
    return argv


MIN_NEIGHBORHOOD_RADIUS_LY = 13
"""int: The smallest neighborhood radius offered, in light-years (about
one sector edge)."""

MAX_GENERATE_RADIUS_LY = int(pc_to_ly(MAX_GENERATE_RADIUS_PC))
"""int: `MAX_GENERATE_RADIUS_PC` in whole light-years (about 652)."""


def _neighborhood_radius_pc(form):
    """
    The "One address" neighborhood radius in parsecs: `address_radius_ly`
    (light-years, 13 to about 652, as the page and the Galaxy Map's
    dialog send it), converted with `ly_to_pc` and rounded to 0.1 pc, or
    the older `address_radius_pc`. `None` for just the one sector.
    """
    radius_ly = _number(form, "address_radius_ly", "Neighborhood radius (ly)", float,
                        minimum=MIN_NEIGHBORHOOD_RADIUS_LY, maximum=MAX_GENERATE_RADIUS_LY)
    if radius_ly is not None:
        return min(round(ly_to_pc(radius_ly), 1), MAX_GENERATE_RADIUS_PC)
    return _number(form, "address_radius_pc", "Radius (pc)", float, minimum=0.1, maximum=MAX_GENERATE_RADIUS_PC)


def center_argv(form, edge_pc=None):
    """
    "Around a sector": `--center-sector` for a filled sector, or the
    One address arguments with a radius (`--shell --band --wedge
    --radius-pc`, which generates that sector first when it's empty) for
    an address or a galaxy-frame position, the position turned into the
    address of the cell holding it with the plan's sector edge
    (`edge_pc`, the standard 4 pc when there is no plan). With "This
    region" ticked (`center_region`, GEN.826) the radius is dropped and
    `--region` snaps the run to the anchor of the region holding the place.

    Returns:
        tuple: `(argv, description)`.
    """
    by = where(form)[1]
    if by == "sector":
        sector_id = _number(form, "center_sector", "Sector ID", int, required=True, minimum=1)
        description = f"around sector {sector_id}"
        argv = ["--center-sector", str(sector_id)]
    elif by == "address":
        shell = _number(form, "center_shell", "Shell", int, required=True, minimum=0, maximum=MAX_GENERATE_SHELL)
        band = _number(form, "center_band", "Band", int) or 0
        wedge = _number(form, "center_wedge", "Wedge", int, required=True, minimum=0)
        description = f"around shell {shell} band {band} wedge {wedge}"
        argv = ["--shell", str(shell), "--band", str(band), "--wedge", str(wedge)]
    elif by == "position":
        point = (_number(form, "center_x_pc", "x (pc)", float, required=True),
                 _number(form, "center_y_pc", "y (pc)", float, required=True),
                 _number(form, "center_z_pc", "z (pc)", float) or 0.0)
        if None in point:
            return [], ""
        edge = float(edge_pc or tuning.DEFAULT_SECTOR_EDGE_PC)
        if sum(v * v for v in point) ** 0.5 > (MAX_GENERATE_SHELL + 0.5) * edge:
            raise FormError("That position is outside the galaxy grid.")
        shell, band, wedge = sector_address_at(point, edge)
        if shell > MAX_GENERATE_SHELL:
            raise FormError("That position is outside the galaxy grid.")
        description = (f"around {', '.join(f'{v:g}' for v in point)} pc "
                       f"(shell {shell} band {band} wedge {wedge})")
        argv = ["--shell", str(shell), "--band", str(band), "--wedge", str(wedge)]
    else:
        raise FormError("Choose how to name the center sector.")
    if form.get("center_region"):
        # GEN.826: the whole region the place is in, as the fill of its anchor's standard sphere; the radius is not used
        return argv + ["--region"], description.replace("around ", "in the region of ", 1)
    radius = _number(form, "center_radius_pc", "Radius (pc)", float, required=True, minimum=0.1,
                     maximum=MAX_GENERATE_RADIUS_PC)
    return argv + ["--radius-pc", str(radius)], description


def galaxy_argv(form, edge_pc=None):
    """
    `planetgen sectors` arguments for the chosen mode. `edge_pc` is the
    plan's sector edge, which only "Around a sector" by position needs.

    Returns:
        tuple: `(argv, description)`.
    """
    mode = where(form)[0]
    if mode == "random":
        return random_start_argv(form), "around a random start"
    if mode == "center":
        return center_argv(form, edge_pc)
    if mode == "address":
        shell = _number(form, "address_shell", "Shell", int, required=True, minimum=0, maximum=MAX_GENERATE_SHELL)
        band = _number(form, "address_band", "Band", int) or 0
        wedge = _number(form, "address_wedge", "Wedge", int, required=True, minimum=0)
        argv = ["--shell", str(shell), "--band", str(band), "--wedge", str(wedge)]
        description = f"at shell {shell} band {band} wedge {wedge}"
        radius = _neighborhood_radius_pc(form)
        if radius is not None:
            argv += ["--radius-pc", str(radius)]
            description = f"around shell {shell} band {band} wedge {wedge}"
        return argv, description
    if mode == "span":
        return span_argv(form)
    if mode == "course":
        return course_argv(form)
    if mode == "shell":
        shell = _number(form, "shell_index", "Shell", int, required=True, minimum=0, maximum=MAX_GENERATE_SHELL)
        argv = ["--shell", str(shell)]
        description = f"in shell {shell}"
        band = _number(form, "shell_band", "Band", int)
        if band is not None:
            argv += ["--band", str(band)]
            description += f" band {band}"
        limit = _number(form, "shell_limit", "Limit", int, minimum=1, maximum=MAX_GENERATE_LIMIT)
        if limit is not None:
            argv += ["--limit", str(limit)]
        elif form.get("whole_shell"):
            argv.append("--yes")
        return argv, description
    if mode == "block":
        return block_argv(form)
    raise FormError("Choose what to generate.")


def block_argv(form):
    """
    `planetgen sectors --block` arguments for a block of the Galaxy Map's ladder (UX.127, GEN.310): `block_key`
    (`m.shell.band.wedge` of a stage level, as the map's stage links write it), then `--limit` (`block_limit`) or, for a
    whole-block run the admin has confirmed (`whole_block`), `--yes`.

    Returns:
        tuple: `(argv, description)`.

    Raises:
        FormError: A key that names no block of a stage level, or a limit out of range.
    """
    text = (form.get("block_key") or "").strip()
    if not text:
        raise FormError("Block is required.")
    try:
        key = drill.format_drill_key(drill.parse_stage_key(text))
    except ValueError as exc:
        raise FormError(f"Block: {exc}") from None
    argv = ["--block", key]
    limit = _number(form, "block_limit", "Limit", int, minimum=1, maximum=MAX_GENERATE_LIMIT)
    if limit is not None:
        argv += ["--limit", str(limit)]
    elif form.get("whole_block"):
        argv.append("--yes")
    return argv, f"in block {key}"


def course_argv(form):
    """
    `planetgen sectors` arguments for charting a course (NAV.48, from the NAV chart page): the two ends
    (object references, as the NAV page writes them), `--course-border` when ticked, and `--yes` once the
    admin has confirmed a chart past `tuning.NAV_CHART_CONFIRM_SECTORS` sectors.

    Returns:
        tuple: `(argv, description)`.

    Raises:
        FormError: An end that is not an object reference, or is a sector.
    """
    ends = []
    for field, label in (("course_from", "Start"), ("course_to", "Destination")):
        text = (form.get(field) or "").strip()
        try:
            kind, printed = objectref.parse_public(text)
        except ValueError:
            kind = printed = None
        if kind is None or kind == "sector":
            raise FormError(f"{label} must be the object ID of a system, a body in one, or a phenomenon.")
        ends.append(objectref.format_public(kind, printed))
    argv = ["--course", *ends]
    description = f"the uncharted sectors on the course {ends[0]} to {ends[1]}"
    if form.get("course_border"):
        argv.append("--course-border")
        description += " and a border of one cell"
    if form.get("course_confirm"):
        argv.append("--yes")
    return argv, description


def span_argv(form):
    """
    `planetgen sectors` arguments for the span mode (ADM.29): the Shells,
    Bands and Wedges ranges, each blank for all, and the limit.

    Returns:
        tuple: `(argv, description)`.

    Raises:
        FormError: A range that can't be read, nothing asked for, or a
            wedge arc without exactly one shell and band.
    """
    argv, ranges = [], {}
    for key, label in (("shells", "Shells"), ("bands", "Bands"), ("wedges", "Wedges")):
        text = (form.get(f"span_{key}") or "").strip()
        if not text:
            continue
        try:
            ranges[key] = parse_range(text, label)
            argv.append(f"--{key}={text}")
        except SpanError as exc:
            raise FormError(str(exc)) from None
    if not argv:
        raise FormError("Give at least one range: shells, bands or wedges.")
    try:
        span = Span(**ranges)
    except SpanError as exc:
        raise FormError(str(exc)) from None
    limit = _number(form, "span_limit", "Limit", int, minimum=1, maximum=MAX_GENERATE_LIMIT)
    if limit is not None:
        argv += ["--limit", str(limit)]
    elif form.get("whole_span"):
        argv.append("--yes")
    return argv, f"in {span.describe()}"


def luminosity_floor_argv(form):
    """
    The form's galaxy-wide threshold as `--bright-star-min-luminosity` (GEN.30; blank means
    `tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL`, and gives no option).

    Raises:
        FormError: A threshold that isn't a number between the lowest and
            highest luminosity floor (GEN.184).
    """
    threshold = _number(form, "bright_min_luminosity", BRIGHT_THRESHOLD_LABEL, float,
                        minimum=tuning.BRIGHT_STAR_FLOOR_MIN_SOL, maximum=tuning.BRIGHT_STAR_FLOOR_MAX_SOL)
    return [] if threshold is None else ["--bright-star-min-luminosity", f"{threshold:.0f}"]


def mass_value_argv(form):
    """The form's mass limit as `--phenomenon-min-mass` alone (GEN.183: one of `tuning.PHENOMENON_MIN_MASS_PRESETS`;
    blank leaves the galaxy's stored limit and gives no option).

    Raises:
        FormError: A value that isn't one of the presets.
    """
    value = _number(form, "phenomenon_min_mass", MASS_LIMIT_LABEL, float)
    if value is None:
        return []
    if not any(abs(value - preset) < 1e-9 for preset in tuning.PHENOMENON_MIN_MASS_PRESETS):
        _problem(f"{MASS_LIMIT_LABEL} must be one of "
                 + ", ".join(f"{preset:g}" for preset in tuning.PHENOMENON_MIN_MASS_PRESETS) + ".")
        return []
    return ["--phenomenon-min-mass", f"{value:g}"]


SCATTER_STEP_LABELS = {"mass": MASS_SCATTER_LABEL, "luminosity": SCATTER_LABEL, "phenomena": PHENOMENA_LABEL}
"""dict: The label of each scatter step (`scatter_steps`)."""


def scatter_steps(generate, form, kinds=("mass", "luminosity", "phenomena")):
    """
    One job step per scatter in `kinds`, in the order the scatters have to run (PERF.121: `planetgen scatter` draws one
    kind per run): the massive stars at the form's mass limit and threshold, the bright stars at its threshold (the
    mass limit stored by the step before), and the phenomena at its neutron star and black hole limit, which follows
    a mass limit the form gave when the form names none of its own.

    Raises:
        FormError: A mass limit, a threshold or a neutron star limit that is not allowed.
    """
    mass = mass_value_argv(form) if "mass" in kinds else []
    floor = luminosity_floor_argv(form) if {"mass", "luminosity"} & set(kinds) else []
    compact = []
    if "phenomena" in kinds:
        compact = compact_limit_argv(form) or (["--compact-min-mass", tuning.COMPACT_MIN_MASS_STAR] if mass else [])
    argv = {"mass": ["scatter", "--mass", *mass, *floor], "luminosity": ["scatter", "--luminosity", *floor],
            "phenomena": ["scatter", "--phenomena", *compact]}
    return [{"label": SCATTER_STEP_LABELS[kind], "argv": generate + argv[kind]}
            for kind in ("mass", "luminosity", "phenomena") if kind in kinds]


PREVALENCE_LABELS = {
    "habitable_world": "Habitable worlds (% of systems)",
    "asteroid_belt": "Asteroid belts (% of systems)",
    "comets": "Comets (% of systems)",
    "large_star": "Large stars (% of systems)",
    "moons": "Moons (% of planets with any)",
    "max_planets": "Most planets a star can hold (% of stars)",
    "intelligent_life": "Intelligent life (% of worlds that reach a technological age)",
    "binary_system": "Binary systems (% of systems)",
    "wide_binary": "Wide pairs (% of binaries)",
    "planets": "Any planets (% of systems)",
    "multi_star": "Five or more stars (% of systems)",
}
"""dict: The Generate page's label for each prevalence feature (ADM.16),
naming what its share is of."""


def _usual_percent(feature):
    """`feature`'s usual share as the page shows it: a percentage, to
    three significant figures."""
    return f"{prevalence.USUAL_SHARES[feature] * 100:.3g}"


STAR_MIX_FEATURES = ("binary_system", "wide_binary")
"""tuple[str]: The features the star mix replaces on the page (ADM.45)."""

PREVALENCE_FIELDS = tuple((f"prevalence_{feature}", PREVALENCE_LABELS[feature], _usual_percent(feature))
                          for feature in prevalence.FEATURES if feature not in STAR_MIX_FEATURES)
"""tuple: `(field name, label, usual share in percent)` for each
prevalence field, in `prevalence.FEATURES` order. Each field starts at
its feature's usual share (Boss, 2026-10-08: "the page needs to have
meaningful information at all times")."""


STAR_MIX_LABELS = {
    "single": "Single stars (% of systems)",
    "close_pair": "Close binaries (% of systems)",
    "wide_pair": "Wide pairs (% of systems)",
}
"""dict: The Generate page's label for each star-mix share."""

STAR_MIX_FIELDS = tuple((f"star_mix_{kind}", STAR_MIX_LABELS[kind], f"{prevalence.usual_star_mix()[kind] * 100:.3g}")
                        for kind in prevalence.STAR_MIX)
"""tuple: `(field name, label, usual share in percent)` for the three
shares of systems by stars (ADM.45), which must total 100%."""


def star_mix_argv(form):
    """
    `--star-mix SINGLE CLOSE WIDE` when the star mix differs from the usual
    one (ADM.45). The three shares must total 100%; the error says which to
    raise or lower and by how much.

    Raises:
        FormError: A share that isn't a number from 0 to 100, or a total
            that isn't 100%.
    """
    shares = []
    for name, label, usual in STAR_MIX_FIELDS:
        share = _number(form, name, label, float, minimum=0, maximum=100)
        shares.append(float(usual) if share is None else share)
    if shares == [float(usual) for _, _, usual in STAR_MIX_FIELDS]:
        return []
    try:
        prevalence.star_mix_prevalences(*shares)
    except ValueError as exc:
        raise FormError(str(exc)[0].upper() + str(exc)[1:] + ".")
    return ["--star-mix", *(f"{share:.6g}" for share in shares)]


def prevalence_argv(form):
    """
    `--prevalence FEATURE=PERCENT` for each prevalence field changed from
    its usual share (ADM.16, GEN.52): the field is the share wanted, and
    the run gets the percentage that moves the usual share there
    (`prevalence.percent_for_share`). Blank or the usual share leaves the
    feature alone.

    Raises:
        FormError: A share that isn't a number from 0 to 100.
    """
    argv = []
    for name, label, usual in PREVALENCE_FIELDS:
        feature = name[len("prevalence_"):]
        share = _number(form, name, label, float, minimum=0, maximum=100)
        if share is None or share == float(usual):
            continue
        percent = prevalence.percent_for_share(feature, share / 100.0)
        argv += ["--prevalence", f"{feature}={percent:.6g}"]
    return argv + star_mix_argv(form) + max_stars_argv(form)


MAX_STARS_LABEL = "Most stars in a system"
"""str: The Generate page's label for the star cap (`--max-stars`, GEN.267)."""


def max_stars_argv(form):
    """
    `--max-stars N` when the Most stars field differs from the default cap (`tuning.MAX_STARS_DEFAULT`); blank or the
    default leaves it alone.

    Raises:
        FormError: A value that isn't a whole number from 1 to `tuning.MAX_STARS_CEILING`.
    """
    count = _number(form, "max_stars", MAX_STARS_LABEL, int, minimum=1, maximum=tuning.MAX_STARS_CEILING)
    if count is None or count == tuning.MAX_STARS_DEFAULT:
        return []
    return ["--max-stars", str(count)]


DIRECTIVE_FIELDS = (
    ("directive_systems", "Systems", "systems"),
    ("directive_habitable", "Habitable worlds", "habitable"),
    ("directive_type_g", "G-type stars", "type:G"),
    ("directive_type_k", "K-type stars", "type:K"),
    ("directive_type_m", "M-type stars", "type:M"),
)


def directive_argv(form):
    """
    `--directive KEY>=N` for each minimum the Override section asks for
    (GEN.96); blank or 0 asks for nothing.

    Raises:
        FormError: A count that isn't a whole number from 0 to 1,000.
    """
    argv = []
    for name, label, key in DIRECTIVE_FIELDS:
        count = _number(form, name, label, int, minimum=0, maximum=1000)
        if count:
            argv += ["--directive", f"{key}>={count}"]
    return argv


REDO_CHOICES = (("mass", "redo_mass"), ("luminosity", "redo_luminosity"), ("phenomena", "redo_phenomena"))
"""tuple: The scatter (a kind of `planetgen scatter`) and the form's checkbox for it (GEN.196)."""


CLEAR_STEP_LABELS = {"luminosity": "Clear the bright stars", "mass": "Clear the massive stars",
                     "phenomena": "Clear the phenomena"}
"""dict: The label of each clear step (`redo_scatters_steps`), in the order the clears have to run: the bright stars
before the massive stars they are split from."""


def redo_scatters_steps(generate, form):
    """
    The steps that redo the scatters whose boxes are ticked, each at its setting (GEN.196, Boss's F7): `planetgen
    clear-scatter` for each ticked scatter, then `planetgen scatter` for each (`scatter_steps`), at the stellar mass
    limit for the massive stars, the luminosity floor for the bright stars, the neutron star and black hole limit for
    the phenomena. A scatter left unticked is left as it was, except that the bright stars are split from the massive
    stars at the mass limit, so redoing the massive stars redoes the bright stars with them.

    Raises:
        FormError: Nothing ticked, or a setting that isn't one of its presets.
    """
    ticked = [word for word, field in REDO_CHOICES if form.get(field)]
    if not ticked:
        _problem("Tick at least one scatter to redo.")
        return []
    if "mass" in ticked and "luminosity" not in ticked:
        ticked.append("luminosity")
    unwanted = ({"phenomenon_min_mass"} if "mass" not in ticked else set()) | (
        {"bright_min_luminosity"} if "luminosity" not in ticked else set())
    clears = [{"label": CLEAR_STEP_LABELS[kind], "argv": generate + ["clear-scatter", f"--{kind}"]}
              for kind in CLEAR_STEP_LABELS if kind in ticked]
    return clears + scatter_steps(generate, {key: value for key, value in form.items() if key not in unwanted}, ticked)


def plan_steps(generate, form):
    """
    The plan step and, unless the form's "Skip all scatters" box is ticked, the scatters as the steps after
    it (`scatter_steps`: the massive stars, the bright stars, the phenomena), at the form's limits.

    Returns:
        list[dict]: Job steps.
    """
    steps = [{"label": "Plan the galaxy", "argv": generate + ["plan"] + plan_argv(form)}]
    if not form.get("skip_scatters"):
        steps += scatter_steps(generate, form)
    return steps


def after_sectors_steps(generate, scatters=()):
    """
    The steps after a `planetgen sectors` step (PERF.123): `link`, the `scatters` (`scatter_steps`), the bright-star
    `backfill` and `settle`, each for the sectors saved since the start of the newest `sectors` run (their default).
    """
    return [{"label": "Link the neighbours", "argv": generate + ["link"]}, *scatters,
            {"label": "Scatter the massive stars from the neighborhood", "argv": generate + ["backfill"]},
            {"label": "Save the sector paths", "argv": generate + ["settle"]}]


def build_job(action, form, database, edge_pc=None):
    """
    The job a form asks for (`edge_pc`: see `galaxy_argv`). Every job
    that generates starts with the data check (`DATA_CHECK_LABEL`); a plain reset and a database check don't.

    Returns:
        tuple: `(kind, title, steps)` for `jobs.start_job`.

    Raises:
        FormError: A bad value, or a missing/wrong typed confirmation.
    """
    if action in CONFIRM_ACTIONS and (form.get("confirm") or "").strip() != database:
        raise FormError(f"Type the database name ({database}) to confirm. Nothing was changed.")
    if action == "bring_to_epoch":
        _require_bringable(database)
    with collecting_problems():
        kind, title, steps = _build_job_steps(action, form, edge_pc, database)
    if kind not in ("reset", "check_db", "finish_run"):
        python = jobs.python_executable()
        disregard = bool(form.get("disregard_missing_data"))
        if disregard:       # GEN.599: every generating command of the job is told, the data check only reports
            generating = [python, *jobs.GENERATE_COMMAND]
            for step in steps:
                if step["argv"][:len(generating)] == generating:
                    step["argv"] = step["argv"] + [DISREGARD_FLAG]
        steps = [{"label": DATA_CHECK_LABEL,
                  "argv": [python, *jobs.DATASETS_COMMAND, "status", *([DISREGARD_FLAG] if disregard else [])]}, *steps]
    if form.get("debug"):
        for step in steps:
            step["argv"] = step["argv"] + [DEBUG_FLAG]
    for step in steps:
        step["stages"] = step_stages(step)
    return kind, title, steps


def step_stages(step):
    """
    The stages one job step goes through (UX.89), `{"label", "skipped"}` each:
    the stages its `sectors`, `galaxy` or `plan` command declares (`generation.stages`); any other step is one stage
    of its own label.
    """
    argv = step["argv"]
    start = len(jobs.GENERATE_COMMAND) + 1
    command = argv[start] if len(argv) > start and argv[:start] == [argv[0], *jobs.GENERATE_COMMAND] else None
    if command in STAGED_COMMANDS:
        try:
            from planetgen.cli import generate as generate_cli
            _parser, parsers = generate_cli.build_parser()
            args = parsers[command].parse_args(argv[start + 1:])
            found = stages.for_command(command, args)
        except (SystemExit, Exception):  # noqa: BLE001 -- the step is then one stage
            found = []
        if found:
            return [{"label": stage.label, "skipped": stage.skip} for stage in found]
    return [{"label": step["label"], "skipped": None}]


def stage_view(job):
    """
    What the page shows of a job's stages (UX.89): `stage_text` ("Stage 7 of
    12: Scatter the phenomena", with the reason when it is being skipped) and
    `stage_list`, one `{n, label, state, reason}` per stage with `state` one of
    done, running, skipped or waiting.
    """
    stages_ = job.get("stages") or []
    step = int(job.get("step") or 0)
    reported = (job.get("progress") or {}).get("stage") or {}
    current = _current_stage_number(job, stages_) if 0 < step and any(s["step"] == step for s in stages_) else None
    listing = []
    for stage in stages_:
        if job.get("finished") and current is None:
            state = "skipped" if stage["skipped"] else "waiting"
        elif current is not None and stage["n"] < current:
            state = "skipped" if stage["skipped"] else "done"
        elif current is not None and stage["n"] == current:
            state = "skipped" if stage["skipped"] or reported.get("skipped") else "running"
        else:
            state = "skipped" if stage["skipped"] else "waiting"
        reason = stage["skipped"] or (reported.get("skipped") if stage["n"] == current else None)
        listing.append({"n": stage["n"], "label": stage["label"], "state": state, "reason": reason})
    text = ""
    if current is not None and not job.get("finished"):
        entry = listing[current - 1]
        text = f"{stagemark.marker(current, len(stages_))}: {entry['label']}"
        if entry["state"] == "skipped":
            text += f" (skipped: {entry['reason']})"
    return {"stage_text": text, "stage_list": listing}


def _build_job_steps(action, form, edge_pc=None, database=None):
    """`build_job`'s `(kind, title, steps)` before the data check step."""
    python = jobs.python_executable()
    generate = [python, *jobs.GENERATE_COMMAND]
    reset_step = {"label": "Reset the galaxy", "argv": [python, *jobs.RESET_COMMAND, "--yes"]}
    if action == "new_galaxy":
        # GEN.30: the scatters run after the sectors, so they leave out every sector just filled.
        plan = {"label": "Plan the galaxy", "argv": generate + ["plan"] + plan_argv(form)}
        sectors = {"label": "Generate sectors around a random start",
                   "argv": generate + ["sectors"] + random_start_argv(form) + prevalence_argv(form)}
        # The plan step doesn't store the mass limit, so the scatters take it from here.
        scatters = [] if form.get("skip_scatters") else scatter_steps(generate, form)
        return "new_galaxy", "New galaxy", [reset_step, plan, sectors, *after_sectors_steps(generate, scatters)]
    if action == "new_plan":
        # UX.147: the plan alone, with no sectors and no scatter, to add sectors by hand on the Galaxy Map.
        plan = {"label": "Plan the galaxy", "argv": generate + ["plan"] + plan_argv(form)}
        return "new_plan", "Plan a new galaxy", [reset_step, plan]
    if action == "plan":
        return "plan", "Plan the galaxy", plan_steps(generate, form)
    if action == "redo_scatters":
        return "redo_scatters", "Redo scatters", redo_scatters_steps(generate, form)
    if action == "bright_band":
        down_to = _number(form, "down_to", "Go down to (solar luminosities)", float, required=True, minimum=1.0)
        if down_to is None:
            return "bright_band", "", []
        argv = generate + ["scatter", "--band", "--down-to", f"{down_to:g}"]
        return "bright_band", f"Bright stars down to {down_to:g} L\u2609", [{"label": BAND_LABEL, "argv": argv}]
    if action == "galaxy":
        argv, description = galaxy_argv(form, edge_pc)
        label = f"Generate sectors {description}"
        sectors = {"label": label, "argv": generate + ["sectors"] + argv + prevalence_argv(form) + directive_argv(form)}
        return "galaxy", label, [sectors, *after_sectors_steps(generate)]
    if action == "reset":
        return "reset", "Reset the galaxy", [reset_step]
    if action == "bring_to_epoch":
        label = "Bring the galaxy to the running generator epoch"
        return "bring_to_epoch", label, [{"label": label, "argv": generate + ["bring-to-epoch"]}]
    if action == "finish_run":
        # PERF.182: the command carries on the run the page just read. A job is one run (PERF.203), so what never started
        # is stages of that run, and the run's own record carries them on.
        found = finish_status(database)
        if found["error"]:
            raise FormError(found["error"] + " Nothing was changed.")
        label = f"Finish the last run: {found['title']}"
        finish = {"label": label, "argv": generate + ["finish-run", "--run", found["run"]]}
        return "finish_run", label, [finish]
    if action == "check_db":
        deep = bool(form.get("deep"))
        argv = generate + ["check-db"] + (list(DEEP_CHECKS) + ["--yes"] if deep else [])
        title = "Check the database (deep)" if deep else "Check the database"
        return "check_db", title, [{"label": title, "argv": argv}]
    raise FormError("Unknown action.")


# ---------------------------------------------------------------------
# Views
# ---------------------------------------------------------------------

@bp.after_request
def _allow_inline_styles_for_the_terminal(response):
    """The two pages with the job terminal get `JOB_LOG_CONTENT_SECURITY_POLICY`
    (see there); the page's other responses (redirects, errors, JSON, downloads) are untouched."""
    if (request.endpoint in ("web.generate", "web.generate_job") and response.status_code == 200
            and response.mimetype == "text/html"):
        from planetgen.web import JOB_LOG_CONTENT_SECURITY_POLICY
        response.headers["Content-Security-Policy"] = JOB_LOG_CONTENT_SECURITY_POLICY
    return response


def _no_store(response):
    response = make_response(response)
    response.headers["Cache-Control"] = "no-store"
    return response


def _admin_or_redirect():
    """`(admin, None)`, or `(None, redirect)` to log in (or to change the
    seeded default credentials first)."""
    admin = current_admin()
    here = request.full_path.rstrip("?")
    if admin is None:
        return None, _no_store(redirect(page_url("login", next=here), code=302))
    if admin.get("must_change_credentials"):
        return None, _no_store(redirect(page_url("account", next=here), code=302))
    return admin, None


def _admin_or_403():
    admin = current_admin()
    if admin is None or admin.get("must_change_credentials"):
        activity_log.event("AUTHZ", "admin.required", user=admin.get("username") if admin else None,
                          path=request.path)
        abort(403, description="Only a logged-in admin can do that.")
    return admin


def _galaxy_summary(database):
    """Whether the galaxy is planned and how many sectors it holds.
    Fails open: the page still works when the database doesn't answer."""
    summary = {"shape": None, "pins_text": None, "sectors": None, "bright": None, "shells": None,
               "version_warning": None, "error": None}
    try:
        summary["shape"] = apiclient.get_galaxy_shape(database)
        if summary["shape"]:
            summary["pins_text"] = dials.describe_pins(summary["shape"].get("pins"))
        summary["bright"] = apiclient.get_bright_star_status(database)
        summary["shells"] = apiclient.get_shell_specs(database)
        summary["version_warning"] = apiclient.get_version_warning(database)
        summary["sectors"] = apiclient.get_sectors(database, limit=1, offset=0)["total"]
    except (apiclient.ApiError, apiclient.NotFoundError, cachebuild.CacheRebuilding) as exc:
        # CacheRebuilding: a run or rebuild is on and the shape's cache is not ready (admin/cachebuild.py); the console
        # must stay usable then, since it is where the run is watched.
        summary["error"] = str(exc)
    return summary


def _edge_pc(action, form, database):
    """The plan's sector edge in parsecs, asked for only when "Around a
    sector" was given a position (the one form that needs it); `None`
    otherwise, or when the database doesn't answer (the standard edge is
    then used)."""
    if action != "galaxy" or where(form) != ("center", "position"):
        return None
    try:
        shape = apiclient.get_galaxy_shape(database)
    except (apiclient.ApiError, apiclient.NotFoundError, cachebuild.CacheRebuilding) as exc:
        log.error(f"Generate page: no galaxy shape for a center position: {exc}")
        return None
    return (shape or {}).get("edge_pc")


MAKES_SECTORS_KINDS = ("galaxy", "new_galaxy")
"""tuple: Job kinds whose finished page offers "Show on Galaxy Map" (ADM.31)."""


def _job_view(job):
    """A `jobs.get_job` dict plus what the page and the script show."""
    if job is None:
        return None
    view = dict(job)
    view["url"] = url_for("web.generate_job", job_id=job["id"])
    progress = job.get("progress") or {}
    view["progress_completed"] = progress.get("completed")
    view["progress_total"] = progress.get("total")
    view["progress_description"] = progress.get("description")
    # Plain UTC text for the JSON status; the page shows `created_html`,
    # which static/localtime.js turns into the viewer's own zone.
    view["created_text"] = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime(job["created_at"])) \
        if job.get("created_at") else ""
    view["created_html"] = trusted_html(utc_time_html(job.get("created_at")))
    view["elapsed_text"] = format_elapsed(job.get("elapsed_s"))
    view["rerun_url"] = generate_flow.rerun_url(job)
    # ADM.31: a finished run that makes sectors offers them on the Galaxy Map.
    view["made_url"] = None
    if job.get("finished") and job.get("kind") in MAKES_SECTORS_KINDS and job.get("started_at"):
        view["made_url"] = url_for("web.galaxy", made=f"{int(job['started_at'])},{int(job.get('finished_at') or time.time()) + 1}")
    # The run's own decaying-average ETA (PERF.7), counted down from when
    # it was written, so the page and the terminal agree.
    eta_s = progress.get("eta_s")
    remaining = None
    if eta_s is not None and not job.get("finished"):
        updated_at = progress.get("updated_at") or time.time()
        remaining = max(0.0, float(eta_s) - max(0.0, time.time() - float(updated_at)))
    view["remaining_text"] = format_elapsed(remaining)
    view["remaining_label"] = remaining_label(
        remaining, running=not job.get("finished") and bool(progress.get("description")))
    view["progress_text"] = _progress_text(progress)
    view.update(overall_view(job, remaining))
    view.update(stage_view(job))
    # PERF.4: the second bar (the bright-star shells being drawn while
    # shells are slow), or blank.
    detail = progress.get("detail") if not job.get("finished") else None
    view["progress_detail_text"] = ""
    if detail:
        text = _progress_text(detail)
        if detail.get("eta_s") is not None:
            text += f", about {format_elapsed(detail['eta_s'])} left"
        view["progress_detail_text"] = text
    return view


def job_runs(job_id):
    """
    The runs the job started (PERF.99), read from their records in Redis, in the order they started: each as the
    Queue page shows a run (`queue_page.run_view`), with its page and its log. `[]` for a job that started none, or
    when no Redis answers.
    """
    from .queue_page import run_view

    try:
        connection = redisqueue.connect()
        now = time.time()
        found = [runview.read(connection, run_id, now) for run_id in runrecord.job_runs(connection, job_id)]
    except Exception as exc:  # noqa: BLE001 -- the job page is the job's output; the runs are an extra
        log.debug(f"Generate page: the runs of job {job_id} were not read ({exc}).")
        return []
    return [run_view(run) for run in found if run is not None]


def _back_to_page():
    """The admin page the strip's button was pressed on (the referrer, when it is this site), else Generate."""
    parts = urlsplit(request.referrer or "")
    if parts.netloc == request.host and parts.path.startswith("/admin") and not parts.path.startswith("//"):
        return parts.path + (f"?{parts.query}" if parts.query else "")
    return url_for("web.generate")


def strip_view():
    """
    What the admin shell's job strip shows (ADM.119, ADM.123, PERF.154), or `None` when nothing runs or waits:
    `{"job", "line"}`. `job` is the running job's `id`, `title`, `status`, `url`, `text` (the stage, else the
    progress), `elapsed_text` and `value` (0 to 1, or `None` when it does not know how far it is); `line` is the waiting
    list in the order it arrived (`line_entries`). `static/adminjobstrip.js` keeps it current from
    `/admin/generate/status`.

    A waiting job whose predecessor has ended is started here, so it never waits for the runner's last act.
    """
    try:
        root = jobs.jobs_dir()
        jobs.promote_queued(root)
        job = jobs.active_job(root)
        line = line_entries(root)
    except OSError:
        return None
    if job is None and not line:
        return None
    view = _job_view(job) if job else None
    value = None
    if view:
        value = view.get("overall_value") if view.get("overall_shown") else None
        if value is None and view.get("progress_total"):
            value = min(1.0, (view.get("progress_completed") or 0) / view["progress_total"])
    return {
        "job": {
            "id": view["id"], "title": view.get("title") or "Job", "status": view.get("status") or "",
            "url": view["url"], "text": view.get("stage_text") or view.get("progress_text") or "",
            "elapsed_text": view.get("elapsed_text") or "", "value": value,
        } if view else None,
        "line": line,
    }


def line_entries(root=None):
    """
    The waiting list as the strip and the Generate page show it, in the order it arrived: `{"ticket", "type" (`job` or
    `run`), "title", "who", "held", "position" (in its own order, 1 next), "note" (the words), "startable" (Start now
    is offered), "discardable"}`. A job waits for the running job, a run for the machine.
    """
    found = jobs.line_view(root)
    return [{**entry, "note": _line_note(entry), "startable": entry["type"] == "job" and entry["position"] == 1 and bool(entry["held"]),
             "discardable": entry["type"] == "job"} for entry in found]


def _line_note(entry):
    if entry["held"]:
        return entry["held"]
    if entry["type"] == "run":
        return f"Waiting for the machine, place {entry['position']}."
    if entry["position"] > 1:
        return f"Starts after the {entry['position'] - 1} job{'s' if entry['position'] > 2 else ''} ahead of it."
    return "Starts when the running job ends."


def overall_view(job, step_remaining, now=None):
    """
    The one bar across all of a job's stages (PERF.55), for a job with more
    than one: `overall_shown`, `overall_value` (0 to 1, `None` while the
    time left is unknown), `overall_text` and `overall_remaining_label`.

    The time left is the running stage's own (`step_remaining`, the bar's
    ETA counted down) or else what earlier runs of that stage took less its
    time so far, plus an estimate for every stage still to run: what earlier
    runs of it took (`stage_estimates`), else the average of the stages this
    job has finished (PERF.66), else, in the job's first stage, what that
    stage is expected to take in all. So the bar's time left is longer than
    the running stage's own until the last stage. The bar is the job's
    elapsed time over elapsed plus left, so it moves forward as each stage
    reports and is corrected when one runs long.
    """
    units = job.get("stages") or [{"n": i, "step": i, "skipped": None} for i, _label in enumerate(job.get("steps") or [], 1)]
    view = {"overall_shown": len(units) > 1, "overall_value": None, "overall_text": "",
            "overall_remaining_label": ""}
    if len(units) < 2:
        return view
    if job.get("finished"):
        done = job.get("status") == "succeeded"
        view["overall_value"] = 1.0 if done else None
        view["overall_text"] = f"All {len(units)} stages" if done else ""
        return view
    now = time.time() if now is None else now
    view["overall_text"] = "Whole job"
    step = int(job.get("step") or 0)
    if step == 0:
        return view
    current = _current_stage_number(job, units)
    estimates = _stage_estimate_list(job, units)
    stage_started = ((job.get("progress") or {}).get("stage") or {}).get("started_at") or job.get("step_started_at")
    in_stage = max(0.0, now - float(stage_started)) if stage_started else 0.0
    left = step_remaining
    if left is None and estimates[current - 1] is not None and stage_started:
        left = max(0.0, estimates[current - 1] - in_stage)
    rough = False
    if left is None and stage_started:       # PERF.184: a server with no timings guesses from its hardware
        guess = roughtime.stage_seconds()
        left, rough = max(guess - in_stage, roughtime.OVERRUN_KEEP * guess), True
    if left is None:
        view["overall_remaining_label"] = "estimating the time left"
        return view
    elapsed = job.get("elapsed_s") or 0.0
    finished = [u for u in units if u["n"] < current and not u["skipped"]]
    average = (max(elapsed - in_stage, 0.0) / len(finished)) if finished else in_stage + left
    remaining = left
    for unit in units:
        if unit["n"] > current and not unit["skipped"]:
            known = estimates[unit["n"] - 1]
            remaining += known if known is not None else average
    view["overall_value"] = elapsed / (elapsed + remaining) if elapsed + remaining > 0 else None
    view["overall_remaining_label"] = remaining_label(remaining, rough=rough)
    return view


def _current_stage_number(job, units):
    """The number (1-based, across the whole job) of the stage running now: the running step's stage the run
    last reported, else the step's first."""
    step = int(job.get("step") or 0)
    inside = [unit for unit in units if unit["step"] == step]
    if not inside:
        return max(1, min(step, len(units)))
    reported = int((((job.get("progress") or {}).get("stage")) or {}).get("index") or 1)
    return inside[min(max(reported, 1), len(inside)) - 1]["n"]


def _stage_estimate_list(job, units):
    """The expected seconds of each stage, one entry per stage (`None` where nothing is known): the runner's
    per-stage list, else each step's own estimate for a step of a single stage."""
    stored = list(job.get("stage_estimates") or [])
    if len(stored) == len(units):
        return stored
    per_step = list(job.get("step_estimates") or [])
    sizes = {}
    for unit in units:
        sizes[unit["step"]] = sizes.get(unit["step"], 0) + 1
    return [per_step[unit["step"] - 1] if sizes[unit["step"]] == 1 and unit["step"] <= len(per_step) else None
            for unit in units]


def remaining_label(remaining, running=True, rough=False):
    """What the page says about the time left (PERF.33): "about 3 m 05 s to 4 m 00 s left" as a range around the
    estimate, "estimating the time left" while a running job has none yet, never dashes; `""` once it finished. A
    `rough` figure (an estimate from the hardware, PERF.184) has a wide range and says it is not timed yet (PERF.199)."""
    if remaining is None:
        return "estimating the time left" if running else ""
    if rough:
        low, high = roughtime.rough_range(remaining)
        return f"roughly {format_elapsed(low)} to {format_elapsed(high)} left (not timed yet)"
    low, high = progress_rate.eta_range(remaining)
    low_text, high_text = format_elapsed(low), format_elapsed(high)
    if low_text == high_text:
        return f"about {high_text} left"
    return f"about {low_text} to {high_text} left"


def _progress_text(progress):
    """`"Sectors: 12 of 40"`, or `"Bright stars (3 of 81 shells): 12%"` for
    a bar written as a share (PERF.9), or `""` before the first report."""
    description = progress.get("description")
    if not description:
        return ""
    completed, total = progress.get("completed") or 0, progress.get("total")
    if progress.get("percent") and total:
        return f"{description}: {min(100 * float(completed) / float(total), 100.0):.0f}%"
    text = f"{description}: {format_number(completed)}"
    if total:
        text += f" of {format_number(total)}"
    unit = progress.get("unit")
    if unit and not description.lower().endswith(counted(unit, 2)):
        text += f" {counted(unit, total or completed)}"
    return text


def format_elapsed(seconds):
    """`"1 h 02 m"`, `"3 m 05 s"`, `"12 s"`, or `""` for `None`."""
    if seconds is None:
        return ""
    seconds = int(seconds)
    hours, rest = divmod(seconds, 3600)
    minutes, secs = divmod(rest, 60)
    if hours:
        return f"{hours} h {minutes:02d} m"
    if minutes:
        return f"{minutes} m {secs:02d} s"
    return f"{secs} s"


def page_context(admin, form=None, error=None, status=200):
    """
    What every Generate screen's template is drawn from (`web/generate_flow.py`): the galaxy summary, the running job
    and the last one, the recent jobs, and the form constants (field lists, limits, presets). `form` is the draft being
    shown again.
    """
    database = db_name()
    try:
        root = jobs.jobs_dir()
        active = jobs.active_job(root)
        recent_total = jobs.count_jobs(root)
        recent_page = clamp_page(parse_page(request.args.get("jobs_page")), recent_total)
        recent = jobs.list_jobs(limit=PAGE_SIZE, root=root, offset=page_offset(recent_page))
        jobs_error = None
    except OSError as exc:
        active, recent, jobs_error = None, [], str(exc)
        recent_total, recent_page = 0, 1
        log.error(f"Generate page: {exc}")
    return dict(
        admin=admin,
        database=database,
        summary=_galaxy_summary(database),
        active=_job_view(active),
        active_log=jobs.log_tail(active["id"], max_bytes=8 * 1024) if active else None,
        recent=[_job_view(job) for job in recent],
        recent_pager=pager("jobs_page", recent_page, recent_total, anchor="recent-jobs", label="Recent jobs pages"),
        jobs_error=jobs_error,
        plan_fields=PLAN_FIELDS,
        halo_min_stellar_mass=mass_chain.HALO_MIN_STELLAR_MASS_SOLAR,
        full_structure_min_mass=dials.FULL_STRUCTURE_MIN_STELLAR_MASS_SOLAR,
        prevalence_fields=PREVALENCE_FIELDS,
        max_stars_label=MAX_STARS_LABEL,
        max_stars_default=tuning.MAX_STARS_DEFAULT,
        max_stars_ceiling=tuning.MAX_STARS_CEILING,
        star_mix_fields=STAR_MIX_FIELDS,
        directive_fields=DIRECTIVE_FIELDS,
        sliders=limit_sliders.SLIDERS,
        core_factor=tuning.CORE_TIER_FACTOR,
        core_factor_percent=round(tuning.CORE_TIER_FACTOR * 100),
        core_floor_default=round(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL * tuning.CORE_TIER_FACTOR, 6),
        core_mass_default=round(tuning.PHENOMENON_MIN_MASS_SOLAR * tuning.CORE_TIER_FACTOR, 6),
        seconds_per_system=tuning.CELL_SECONDS_PER_SYSTEM,
        bytes_per_system=tuning.CELL_BYTES_PER_SYSTEM,
        bright_min_luminosity=tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL,
        bright_floor_presets=luminosity_floor.PRESETS,
        bright_threshold_label=BRIGHT_THRESHOLD_LABEL,
        mass_limit_label=MASS_LIMIT_LABEL,
        mass_limit_presets=tuning.PHENOMENON_MIN_MASS_PRESETS,
        mass_limit_default=tuning.PHENOMENON_MIN_MASS_SOLAR,
        compact_label=COMPACT_LIMIT_LABEL,
        compact_presets=tuning.COMPACT_MIN_MASS_PRESETS,
        compact_warnings=COMPACT_STORAGE_WARNINGS,
        backfill_text=BACKFILL_TEXT,
        place_tabs=PLACE_TABS,
        grid_choices=GRID_CHOICES,
        max_radius_pc=MAX_GENERATE_RADIUS_PC,
        min_radius_ly=MIN_NEIGHBORHOOD_RADIUS_LY,
        max_radius_ly=MAX_GENERATE_RADIUS_LY,
        max_shell=MAX_GENERATE_SHELL,
        max_limit=MAX_GENERATE_LIMIT,
        max_neighborhoods=MAX_GENERATE_NEIGHBORHOODS,
        advanced_open=bool(form) and any(name.startswith(ADVANCED_PREFIXES) for name in form),
        error=error,
        form=form or {},
        estimate_confirm_field=ESTIMATE_CONFIRM_FIELD,
        generate_anyway_field=GENERATE_ANYWAY_FIELD,
        status=status,
        last_job=_job_view(recent[0]) if recent and not active else None,
        line=_line(),
        line_full=_line_full(),
    )


ADVANCED_PREFIXES = ("prevalence_", "star_mix_", "directive_", "max_stars")
"""tuple: Form fields that live in the Advanced block; a draft shown again that holds one opens it."""


def _line():
    try:
        return line_entries()
    except OSError:
        return []


def _line_full():
    """Whether the waiting list holds `queue.max_waiting` entries, so a job could not join it."""
    try:
        return len(jobs.line_view()) >= runline.max_waiting()
    except OSError:
        return False


@bp.route("/admin/generate", methods=["GET", "POST"])
def generate():
    """The Generate page (GET) and its forms (POST, answered with a 303
    back to the page, or the page again with an error when nothing
    started)."""
    if request.method in ("GET", "HEAD"):
        admin, response = _admin_or_redirect()
        return response or generate_flow.start_page(admin)

    admin = _admin_or_403()
    action = request.form.get("action") or ""
    if action == "cancel":
        job_id = request.form.get("job") or ""
        try:
            if jobs.cancel_job(job_id, who=admin.get("username")):
                activity_log.event("GEN", "job.cancel", user=admin.get("username"), job=job_id)
        except OSError as exc:
            log.error(f"Could not cancel job {job_id}: {exc}")
        # Cancel on a job's own page returns to the page the admin came from (`next`, set by `generate_job`).
        target = safe_next(request.form.get("next"), default=url_for("web.generate", _anchor="current-job"))
        return _no_store(redirect(target, code=303))

    if action in ("discard_queued", "start_queued"):
        ticket = request.form.get("ticket") or None
        if action == "discard_queued":
            if jobs.discard_queued(ticket=ticket):
                activity_log.event("GEN", "job.unqueue", user=admin.get("username"))
        else:
            started = jobs.promote_queued(force=True)
            if started:
                activity_log.event("GEN", "job.start", user=admin.get("username"), job=started, queued=True)
        return _no_store(redirect(_back_to_page(), code=303))

    database = db_name()
    wants_json = _wants_json()

    def refuse(message, status):
        """The refusal as JSON for a script, else the task's page with the draft and the message."""
        if wants_json:
            return _no_store(make_response(jsonify({"error": message}), status))
        return generate_flow.task_page_for_action(admin, action, form=request.form, error=message, status=status)

    try:
        kind, title, steps = build_job(action, request.form, database, edge_pc=_edge_pc(action, request.form, database))
    except FormError as exc:
        return refuse(str(exc), 400)
    env = jobs.mysql_env(current_app.config["MYSQL_CONFIG"], database)
    if estimate_needed(kind, request.form) and not request.form.get(ESTIMATE_CONFIRM_FIELD):
        # PERF.3: show the size and time first; the admin confirms (or the disk refuses it) before the job starts.
        try:
            estimate = estimate_for(kind, steps, database, env)
        except FormError as exc:
            return refuse(str(exc), 400)
        if wants_json:
            return _no_store(make_response(jsonify({
                "error": estimate.get("refusal") or f"Confirm first: {estimate.get('summary') or estimate.get('what')}",
                "estimate": estimate, "confirm_field": ESTIMATE_CONFIRM_FIELD,
                "generate_anyway_field": GENERATE_ANYWAY_FIELD,
            }), 409))
        return generate_flow.review_page(admin, action, kind, title, steps, request.form, estimate=_estimate_view(estimate))
    try:
        submitted = jobs.submit(kind, title, steps, env=env, admin=admin.get("username"), database=database,
                                draft=generate_flow.draft_of(action, request.form))
    except runline.LineFull as exc:
        return refuse(str(exc), 409)
    except updateguard.UpdateRunning as exc:
        return refuse(str(exc), 503)
    except OSError as exc:
        log.exception(f"Could not start a {kind} job: {exc}")
        return refuse(f"The job could not be started: {exc}", 500)
    job_id = submitted.job_id
    if job_id is None:         # another job is running: this one waits in the list and starts by itself (PERF.154)
        activity_log.event("GEN", "job.queue", user=admin.get("username"), kind=kind, db=database, title=title)
        if wants_json:
            return _no_store(make_response(jsonify({
                "queued": True, "ticket": submitted.ticket, "position": submitted.position, "waiting_for": submitted.running,
                "status_url": url_for("web.generate_status"),
            }), 202))
        return _no_store(redirect(url_for("web.generate"), code=303))
    activity_log.event("GEN", "job.start", user=admin.get("username"), job=job_id, kind=kind, db=database, title=title)
    if kind == "galaxy" and request.form.get(GENERATE_ANYWAY_FIELD):
        activity_log.event("GEN", "job.generate_anyway", user=admin.get("username"), job=job_id, db=database,
                           title=title)
    if wants_json:
        return _no_store(make_response(jsonify({
            "job": job_id,
            "url": url_for("web.generate_job", job_id=job_id),
            "status_url": url_for("web.generate_status", job=job_id),
        }), 202))
    return _no_store(redirect(url_for("web.generate", _anchor="current-job"), code=303))


def estimate_needed(kind, form):
    """Whether a job asks for its size and time first (a galaxy run, or a deep database check)."""
    return kind in ("galaxy", "check_db") and (kind == "galaxy" or bool(form.get("deep")))


def estimate_for(kind, steps, database, env):
    """The size and time of the job's `planetgen sectors` step, or its last step when it has none (`run_estimate`,
    cached by `estimatecache`).

    Raises:
        FormError: The estimate refused the arguments.
    """
    first = 1 + len(jobs.GENERATE_COMMAND)
    sectors = [step for step in steps if step["argv"][first:first + 1] == ["sectors"]]
    argv = (sectors[0] if sectors else steps[-1])["argv"]
    return estimatecache.estimate(database, kind, argv, lambda: _estimate_versions(database),
                                  lambda: run_estimate(argv, env))


def _estimate_view(estimate):
    """`run_estimate`'s dict plus the text the page shows."""
    view = dict(estimate)
    view["time_text"] = stats.time_text(estimate.get("seconds"), estimate.get("measured", True))
    if estimate.get("deep_check"):      # DB.21: a time and a count of systems, no size or disk
        return view
    view["size_text"] = stats.format_bytes(estimate.get("bytes"))
    disk = estimate.get("disk")
    view["disk_text"] = (f"{stats.format_bytes(disk['free_bytes'])} free of "
                         f"{stats.format_bytes(disk['total_bytes'])} on {disk.get('where') or disk['path']}") if disk else None
    note = estimate.get("disk_note")
    view["disk_note_text"] = f"Not measured: {note}. Nothing is refused for space." if note else None
    return view


def _wants_json():
    """True when the form was sent by script (the Galaxy Map's Generate
    buttons) asking for JSON: `Accept: application/json`. Such a caller
    gets `{"job", "url", "status_url"}` (202) or `{"error"}` instead of a
    redirect or the page."""
    best = request.accept_mimetypes.best_match(["application/json", "text/html"])
    return best == "application/json" and request.accept_mimetypes[best] > request.accept_mimetypes["text/html"]


@bp.route("/admin/generate/status")
def generate_status():
    """
    JSON for `static/generatejobs.js`: the running job (or the one named
    by `?job=`, to see it finish) with the tail of its output.
    """
    _admin_or_403()
    root = jobs.jobs_dir()
    job_id = request.args.get("job")
    job = jobs.get_job(job_id, root) if job_id else jobs.active_job(root)
    body = {"job": _job_view(job), "strip": strip_view(), "log": jobs.log_tail(job["id"], max_bytes=8 * 1024, root=root) if job else None}
    return _no_store(jsonify(body))


generate_status.json_only = True  # not a page: tests/test_web_a11y.py skips it


def epoch_status(database):
    """
    What `planetgen bring-to-epoch` would do to the stored galaxy (GEN.764, step 4): `bring_to_epoch.status`'s dict, or
    `{"error": ...}` when the database does not answer. It changes nothing.
    """
    config = current_app.config["MYSQL_CONFIG"]
    galaxy = store.MySQLConfig(host=config.host, port=config.port, user=config.user, password=config.password,
                               database=database)
    try:
        return bring_to_epoch.status(galaxy, current_app.config.get("CONTROL_MYSQL_CONFIG"))
    except Exception as exc:  # noqa: BLE001 -- the page then says why instead of failing
        log.error(f"Generate page: the epoch of the galaxy could not be read: {exc}")
        return {"error": f"The galaxy could not be read: {exc}"}


def finish_status(database):
    """
    What "Finish the last run" would do (PERF.182), read from the run record: `{"error": None, "run", "title",
    "stopped_in", "carried", "again"}` (the id and title of the run that stopped, the stage it stopped in, the stages it
    finished and keeps, and the stages that run now: those the run never started among them), or `{"error": ...}` saying
    why there is nothing to carry on. It changes nothing but the record of a
    run whose process is gone, which it marks stopped. The build is checked by the command, not here.
    """
    try:
        seed = _stored_seed(database)
        found = finish_run.find(redisqueue.connect(), None if seed is None else seed.hex(), None, database)
    except finish_run.CannotFinish as exc:
        return {"error": str(exc)}
    except Exception as exc:  # noqa: BLE001 -- the page then says why instead of failing
        log.error(f"Generate page: the last run could not be read: {exc}")
        return {"error": f"The last run could not be read: {exc}"}
    run, labels = found.run, [step.label for step in found.steps]
    return {"error": None, "run": run["id"], "title": run["title"], "stopped_in": labels[found.stopped],
            "carried": labels[:found.first], "again": labels[found.first:]}


def _require_bringable(database):
    """Raises `FormError` when `planetgen bring-to-epoch` has nothing it could do now: the galaxy could not be read, a
    stage is too old to bring up in part (the message says what to do instead), or it is at the running epoch."""
    found = epoch_status(database)
    if found["error"]:
        raise FormError(found["error"])
    if found["refuse"]:
        raise FormError(" ".join(found["refuse"]) + " Nothing was changed.")
    if found["nothing"]:
        raise FormError(found["lines"][0])


def _stored_seed(database):
    """The seed of the planned galaxy (bytes), or `None` before a plan or when the database doesn't answer."""
    config = current_app.config["MYSQL_CONFIG"]
    try:
        conn = store.get_connection(store.MySQLConfig(host=config.host, port=config.port, user=config.user,
                                                      password=config.password, database=database),
                                    ensure_schema=False)
        try:
            return store.get_galaxy_seed(conn)
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- the line then draws the Auto dials from the sample seed
        log.error(f"Generate page: no galaxy seed for the shape report: {exc}")
        return None


SHAPE_ANSWERS_KEPT = 32
"""int: How many answers of the shape report the page keeps."""

_shape_answers = {}
_shape_lock = threading.Lock()


def _shape_preview(argv, seed):
    """`plan_outline.preview(argv, galaxy_seed=seed)`, kept (TEST.219): it takes 3 to 10 s warm and a minute or more cold on a
    busy machine, the page asks again at every change of a dial (400 ms after the last), and a page that is closed leaves
    its asks running. An ask for dials already worked out is answered at once, one still being worked out is waited for
    rather than worked out a second time, and a refusal (`ValueError`) is raised afresh, not kept."""
    key = (tuple(argv), seed)
    with _shape_lock:
        found = _shape_answers.get(key)
        if found is None:
            found = plan_outline.preview(argv, galaxy_seed=seed)
            while len(_shape_answers) >= SHAPE_ANSWERS_KEPT:
                del _shape_answers[next(iter(_shape_answers))]
            _shape_answers[key] = found
        return found


@bp.route("/admin/generate/shape-report")
def generate_shape_report():
    """
    JSON for `static/generatelevers.js` (GEN.617): the line under the dials. `?mass=...&flatness=...` are the Plan
    fields; the answer is what `planetgen plan` would print for them (`generation.plan_outline`): `line`, the
    `warnings`, and `refusal` when the outline would be empty (the plan stops with that reason). The Auto dials are
    drawn from the planned galaxy's seed, else from a sample seed (`sample`).
    """
    _admin_or_403()
    try:
        with collecting_problems():
            argv = lever_argv(request.args)
    except FormError as exc:
        return _no_store(make_response(jsonify({"error": str(exc)}), 400))
    seed = _stored_seed(db_name())
    try:
        made, refusal = _shape_preview(argv, seed)
    except ValueError as exc:
        return _no_store(make_response(jsonify({"error": str(exc)}), 400))
    return _no_store(jsonify({"line": made.line(), "warnings": made.warnings, "refusal": refusal,
                              "sample": seed is None}))


generate_shape_report.json_only = True  # not a page: tests/test_web_a11y.py skips it


STREAM_SECONDS = 40
"""int: How long one `.../stream` response stays open. Under Apache a
request that outlives `WSGIDaemonProcess request-timeout` (60 s in
`examples/apache/planetgen.conf.example`) restarts the whole daemon, so
the stream ends well before and the browser's `EventSource` reconnects,
resuming from the `Last-Event-ID` (the byte offset reached)."""

STREAM_POLL_SECONDS = 0.5
STREAM_KEEPALIVE_SECONDS = 15
STREAM_REPLAY_BYTES = 256 * 1024
"""int: A first connection (no offset) replays at most this much of the log's end."""


def _sse(event, data, event_id=None):
    lines = [f"event: {event}"]
    if event_id is not None:
        lines.append(f"id: {event_id}")
    lines.append(f"data: {data}")
    return "\n".join(lines) + "\n\n"


def _stream_start_offset(job_id, root):
    """Where a stream starts reading: the `Last-Event-ID` header (a reconnect) or
    `?offset=`, else the last STREAM_REPLAY_BYTES, from the start of a line."""
    raw = request.headers.get("Last-Event-ID") or request.args.get("offset")
    try:
        return max(0, int(raw))
    except (TypeError, ValueError):
        pass
    size = jobs.log_size(job_id, root) or 0
    if size <= STREAM_REPLAY_BYTES:
        return 0
    start = size - STREAM_REPLAY_BYTES
    chunk = jobs.read_log(job_id, start, 4096, root)
    if chunk and "\n" in chunk[0]:
        return start + len(chunk[0].split("\n", 1)[0].encode("utf-8")) + 1
    return start


@bp.route("/admin/generate/jobs/<job_id>/stream")
def generate_job_stream(job_id):
    """
    Server-Sent Events for one job (`static/generatejobs.js`): `log` events
    carrying new output (`{"text": ...}`, with the byte offset reached as
    the event id, so a reconnect resumes), `state` events carrying the job
    as `_job_view` shows it whenever it changes, and `done` once the job
    has finished and its output is all sent. The response ends after
    STREAM_SECONDS; the browser reconnects by itself unless it saw `done`.
    """
    _admin_or_403()
    root = jobs.jobs_dir()
    if jobs.get_job(job_id, root) is None:
        abort(404, description="No such job.")
    offset = _stream_start_offset(job_id, root)

    @stream_with_context
    def events():
        nonlocal offset
        yield "retry: 2000\n\n"
        deadline = time.monotonic() + STREAM_SECONDS
        last_state = None
        last_sent = time.monotonic()
        while True:
            job = jobs.get_job(job_id, root)
            if job is None:
                return
            # Finished is read before the log, so output written just before it is still sent.
            finished = bool(job.get("finished"))
            sent = False
            while True:
                chunk = jobs.read_log(job_id, offset, root=root)
                if not chunk or not chunk[0]:
                    break
                text, offset = chunk
                yield _sse("log", json.dumps({"text": text}), offset)
                sent = True
            state = current_app.json.dumps(_job_view(job))
            if state != last_state:
                last_state = state
                yield _sse("state", state, offset)  # with an id, so a reconnect never replays the tail
                sent = True
            if finished:
                yield _sse("done", "{}")
                return
            now = time.monotonic()
            if sent:
                last_sent = now
            elif now - last_sent >= STREAM_KEEPALIVE_SECONDS:
                yield ": keep-alive\n\n"
                last_sent = now
            if now >= deadline:
                return
            time.sleep(STREAM_POLL_SECONDS)

    response = Response(events(), mimetype="text/event-stream")
    response.headers["Cache-Control"] = "no-store"
    response.headers["X-Accel-Buffering"] = "no"  # nginx: pass events through
    return response


generate_job_stream.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/admin/generate/jobs/<job_id>/log")
def generate_job_log(job_id):
    """A job's full output as a plain-text download; `?form=records` gives the same job as records instead (UX.122): one
    log line per line with its time in UTC, severity and component (`jobs.RECORDS_NAME`)."""
    _admin_or_403()
    root = jobs.jobs_dir()
    if jobs.get_job(job_id, root) is None:
        abort(404, description="No such job.")
    records = request.args.get("form") == "records"
    path = os.path.join(jobs._job_dir(root, job_id), jobs.RECORDS_NAME if records else "output.log")
    if not os.path.isfile(path):
        abort(404, description="This job has no output yet.")
    name = f"job-{job_id}.records.log" if records else f"job-{job_id}.log"
    return _no_store(send_file(path, mimetype="text/plain", as_attachment=True, download_name=name))


generate_job_log.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/admin/generate/sectors")
def generate_sectors():
    """
    JSON for the named-sector tab's list of filled sectors
    (`static/generatefinder.js`): `?page=N`, 50 a page, nearest the
    galactic core first (`apiclient.get_sectors`). Sectors without an
    address (placed before the cylindrical grid) are listed too: the
    center only needs their ID.
    """
    _admin_or_403()
    try:
        page = max(1, int(request.args.get("page") or 1))
    except ValueError:
        page = 1
    try:
        found = apiclient.get_sectors(db_name(), limit=FINDER_PAGE_SIZE, offset=(page - 1) * FINDER_PAGE_SIZE)
    except (apiclient.ApiError, apiclient.NotFoundError) as exc:
        log.error(f"Generate page: sector list failed: {exc}")
        return _no_store(make_response(jsonify({"error": "The sector list could not be loaded."}), 502))
    items = [{
        "id": row["id"], "name": row["name"], "systems": row.get("system_count"),
        "shell": row.get("shell_index"), "band": row.get("band_index"), "wedge": row.get("wedge_index"),
    } for row in found["items"]]
    return _no_store(jsonify({"items": items, "total": found["total"], "page": page,
                              "pages": max(1, -(-found["total"] // FINDER_PAGE_SIZE))}))


generate_sectors.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/admin/generate/jobs/<job_id>")
def generate_job(job_id):
    """One job and the last 256 KB of its output."""
    admin, response = _admin_or_redirect()
    if response:
        return response
    root = jobs.jobs_dir()
    job = jobs.get_job(job_id, root)
    if job is None:
        abort(404, description="No such job.")
    return _no_store(render_page(
        "generate_job.html",
        title=job["title"] or "Job",
        section="admin_generate",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Generate", "generate"), crumb(job["title"] or "Job")],
        admin=admin,
        job=_job_view(job),
        back_url=safe_next(request.args.get("next"), default="") or came_from(request.path),
        runs=job_runs(job_id),
        output=jobs.log_tail(job_id, max_bytes=256 * 1024, root=root),
    ))


from . import generate_flow  # noqa: E402,F401 -- the Generate screens (they use this module's logic)
