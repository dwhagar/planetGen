# html/web/generate_page.py

"""
`/admin/generate`: galaxy generation from the browser, so an admin never
needs a terminal on the server. Four actions, each a background job
(`web/jobs.py`) running the same command-line tools an admin would type:

- New galaxy: wipe the database (`src/resetDb.py --yes`), build the
  density skeleton (`generate.py plan --no-bright-stars`), scatter the
  bright stars (`generate.py plan --bright-stars-only`), then generate a
  first neighborhood around a random start (`generate.py galaxy`).
- Plan: rebuild the density skeleton, then scatter the bright stars.
  The scatter is its own step so the job shows its progress bar and the
  count it placed; a checkbox on either form skips it.
- Rebuild bright stars: the scatter alone, on the stored plan
  (`--force` leaves already filled sectors out instead of refusing).
- Add a dimmer layer: keep the bright stars already placed and add only
  those from a lower level up to the current one
  (`generate.py plan --bright-stars-down-to N`).
- Generate sectors: `generate.py galaxy` in any of its modes (random
  start, a whole ring, around a sector, one address, ...). "Around a
  sector" names its center as a filled sector (searched by name or
  picked from `/admin/generate/sectors`'s list), a sector address, or a
  galaxy-frame position in pc.
- Reset: wipe the database only.

New galaxy and Reset delete every generated sector and system, so both
ask for the database name to be typed back, the same confirmation
`resetDb.py` asks for in a terminal.

The page shows the running job (step, progress bar, elapsed time and the
tail of its output), refreshed every few seconds by
`static/generatejobs.js` from `/admin/generate/status`, plus the last few
jobs, each with its full output at `/admin/generate/jobs/<id>`.

Every section folds (ADM.4; `static/generatefolds.js` remembers which
are open per browser).

Admins only: a visitor who isn't logged in is sent to the login page,
and any POST or status request without an admin session gets a 403.
Every form carries `csrf_field()` (checked app-wide by `csrf.protect`).
"""

import json
import os
import subprocess
import time

from flask import abort, current_app, jsonify, make_response, redirect, request, url_for

import apiclient
from fmt import utc_time_html
from stellarObjects import activitylog, generationStats, log, program_constants
from stellarObjects.galaxyDrill import format_drill_key, parse_drill_key
from stellarObjects.galaxyGeometry import sector_address_at
from stellarObjects.utils import format_number, ly_to_pc, pc_to_ly
from stellarObjects.generationLimits import (
    MAX_GENERATE_LIMIT, MAX_GENERATE_RADIUS_PC, MAX_GENERATE_RING,
)

from . import bp, jobs
from .helpers import crumb, current_admin, db_name, page_url, render_page, trusted_html

# ---------------------------------------------------------------------
# Form fields
# ---------------------------------------------------------------------

PLAN_FIELDS = (
    # (form name, generate.py flag, label, type, default, minimum, exclusive maximum)
    ("disk_scale_length_pc", "--disk-scale-length-pc", "Disk scale length (pc)", float, 2800.0, 1.0, None),
    ("disk_scale_height_pc", "--disk-scale-height-pc", "Disk scale height (pc)", float, 350.0, 1.0, None),
    ("bulge_scale_radius_pc", "--bulge-scale-radius-pc", "Bulge scale radius (pc)", float, 200.0, 1.0, None),
    ("bulge_amplitude", "--bulge-amplitude", "Bulge amplitude", float, 1.0, 0.0, None),
    ("arm_count", "--arm-count", "Spiral arms", int, 2, 0, None),
    ("pitch_angle_deg", "--pitch-angle-deg", "Arm pitch angle (degrees)", float, 15.0, 1.0, 90.0),
    ("arm_amplitude", "--arm-amplitude", "Arm contrast (0 to 1)", float, 0.4, 0.0, 1.0),
)
"""tuple: The `plan` options the page offers, with `generate.py`'s own
defaults (a test checks they still match `generate.py plan`'s parser)."""

GALAXY_MODES = (
    ("random", "Around a random start",
     "Picks a random populated spot and generates the sectors within a radius of it "
     "(12 pc, about 39 ly, when the radius is left blank; the 100 ly around every generated "
     "sector gets its bright stars either way)."),
    ("ring", "A whole ring",
     "Every not-yet-generated sector in one ring at one height layer (0 is the galactic plane), "
     "or only the first few with a limit."),
    ("center", "Around a sector",
     "Every not-yet-generated sector within a radius of a center: an existing filled sector "
     "(found by name or picked from the list), a sector address, or a position in the galaxy."),
    ("slot", "One address",
     "Exactly one sector, by ring, layer and slot (the address the Galaxy Map shows). "
     "With a radius, its neighborhood too."),
    ("column", "A column",
     "Every sector at one ring and slot, through every layer the galaxy reaches there."),
    ("block", "A Galaxy Map block",
     "Every sector the galaxy allows inside one drill-down block, by its key (size.ring.wedge.slab, "
     "as in the Galaxy Map's links), or only one of its layers. A size-3 block holds about 9 sectors "
     "a layer; bigger ones need a limit or the confirmation box past 2,000."),
    ("shell", "A shell (not recommended)",
     "Every sector of one ring through every layer: a whole cylinder, usually thousands of sectors, "
     "so it needs a limit or the confirmation box."),
)

CENTER_CHOICES = (
    ("sector", "A filled sector", "Find it by name, pick it from the list, or type its ID."),
    ("address", "A sector address", "Ring, layer and slot, as the Galaxy Map shows them. "
                                    "The sector there is generated first if it isn't yet."),
    ("position", "A position (pc)", "Galaxy-frame x, y and z in parsecs (the center is 0, 0, 0; "
                                    "z is height above the plane). The sector holding that point is the center."),
)
"""tuple: How "Around a sector" names its center (`center_by`):
`(value, label, hint)`."""

SECTION_FOR_ACTION = {
    "new_galaxy": "new-galaxy", "galaxy": "generate-sectors", "plan": "plan",
    "bright_stars": "bright-stars", "bright_band": "bright-band", "reset": "reset",
}
"""dict: The page section each form lives in, so a form shown again
with an error is open whatever the browser remembered."""

FINDER_PAGE_SIZE = 50
"""int: Sectors per page in the "Around a sector" list (every web list
pages 50 rows)."""

CONFIRM_ACTIONS = frozenset({"new_galaxy", "reset"})

BAND_LABEL = "Add a dimmer layer of bright stars"

SCATTER_LABEL = "Scatter the bright stars"
"""str: The step that pre-places every bright star
(`program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL` and up) on the plan."""


ESTIMATE_TIMEOUT_S = 120
"""int: How long the page waits for `generate.py galaxy --estimate-only`."""

ESTIMATE_PREFIX = "ESTIMATE "
"""str: `generate.ESTIMATE_PREFIX`: the line the estimate is read from."""

ESTIMATE_CONFIRM_FIELD = "estimate_ok"
"""str: The form field a confirmed estimate's "Generate" button adds."""


class FormError(ValueError):
    """A form value the page rejects; the message is shown to the admin."""


def run_estimate(argv, env):
    """
    PERF.3: the size and time a Generate sectors job would take, from
    `generate.py galaxy ... --estimate-only` (run here, while the admin
    waits, before anything is written).

    Args:
        argv (list): The job step's command line.
        env (dict): `jobs.mysql_env`.

    Returns:
        dict: `generationStats.Estimate.as_dict()` plus `what`.

    Raises:
        FormError: The estimate couldn't be made; the message says why
            (`generate.py`'s own last lines, e.g. a ring too large to
            generate without a limit).
    """
    try:
        done = subprocess.run(argv + ["--estimate-only"], env={**os.environ, **env}, capture_output=True,
                              text=True, timeout=ESTIMATE_TIMEOUT_S, check=False)
    except subprocess.TimeoutExpired:
        raise FormError("Working out the size and time took too long. Nothing was generated.") from None
    except OSError as exc:
        raise FormError(f"The size and time couldn't be worked out: {exc}") from None
    for line in reversed(done.stdout.splitlines()):
        if line.startswith(ESTIMATE_PREFIX):
            try:
                return json.loads(line[len(ESTIMATE_PREFIX):])
            except ValueError:
                break
    tail = [line.strip() for line in (done.stderr + "\n" + done.stdout).splitlines() if line.strip()]
    message = " ".join(tail[-3:]) if tail else f"generate.py exited with {done.returncode}"
    raise FormError(f"Nothing was generated: {message}")


def _number(form, name, label, kind, required=False, minimum=None, maximum=None, exclusive_max=False):
    raw = (form.get(name) or "").strip()
    if not raw:
        if required:
            raise FormError(f"{label} is required.")
        return None
    try:
        value = kind(raw)
    except ValueError:
        raise FormError(f"{label} must be {'a whole number' if kind is int else 'a number'}.") from None
    if value != value or value in (float("inf"), float("-inf")):
        raise FormError(f"{label} must be a finite number.")
    if minimum is not None and value < minimum:
        raise FormError(f"{label} must be at least {minimum:g}.")
    if maximum is not None and (value >= maximum if exclusive_max else value > maximum):
        raise FormError(f"{label} must be {'less than' if exclusive_max else 'at most'} {maximum:g}.")
    return value


def plan_argv(form):
    """`generate.py plan` arguments from the Plan fields (blank means
    the default)."""
    argv = []
    for name, flag, label, kind, _default, minimum, maximum in PLAN_FIELDS:
        value = _number(form, name, label, kind, minimum=minimum, maximum=maximum, exclusive_max=True)
        if value is not None:
            argv += [flag, str(value)]
    return argv


def random_start_argv(form):
    """`generate.py galaxy` arguments for a random start."""
    argv = []
    radius = _number(form, "radius_pc", "Radius (pc)", float, minimum=0.1, maximum=MAX_GENERATE_RADIUS_PC)
    if radius is not None:
        argv += ["--radius-pc", str(radius)]
    max_ring = _number(form, "max_ring", "Highest ring", int, minimum=0, maximum=MAX_GENERATE_RING)
    if max_ring is not None:
        argv += ["--max-ring", str(max_ring)]
    density = _number(form, "min_start_density", "Minimum start density", float, minimum=0.001)
    if density is not None:
        argv += ["--min-start-density", str(density)]
    return argv


MIN_NEIGHBORHOOD_RADIUS_LY = 13
"""int: The smallest neighborhood radius offered, in light-years (about
one sector edge)."""

MAX_GENERATE_RADIUS_LY = int(pc_to_ly(MAX_GENERATE_RADIUS_PC))
"""int: `MAX_GENERATE_RADIUS_PC` in whole light-years (about 652)."""


def _neighborhood_radius_pc(form):
    """
    The "One address" neighborhood radius in parsecs: `slot_radius_ly`
    (light-years, 13 to about 652, as the page and the Galaxy Map's
    dialog send it), converted with `ly_to_pc` and rounded to 0.1 pc, or
    the older `slot_radius_pc`. `None` for just the one sector.
    """
    radius_ly = _number(form, "slot_radius_ly", "Neighborhood radius (ly)", float,
                        minimum=MIN_NEIGHBORHOOD_RADIUS_LY, maximum=MAX_GENERATE_RADIUS_LY)
    if radius_ly is not None:
        return min(round(ly_to_pc(radius_ly), 1), MAX_GENERATE_RADIUS_PC)
    return _number(form, "slot_radius_pc", "Radius (pc)", float, minimum=0.1, maximum=MAX_GENERATE_RADIUS_PC)


def center_argv(form, edge_pc=None):
    """
    "Around a sector": `--center-sector` for a filled sector, or the
    One address arguments with a radius (`--ring --layer --slot
    --radius-pc`, which generates that sector first when it's empty) for
    an address or a galaxy-frame position, the position turned into the
    address of the cell holding it with the plan's sector edge
    (`edge_pc`, the standard 4 pc when there is no plan).

    Returns:
        tuple: `(argv, description)`.
    """
    by = form.get("center_by") or "sector"
    if by == "sector":
        sector_id = _number(form, "center_sector", "Sector ID", int, required=True, minimum=1)
        description = f"around sector {sector_id}"
        argv = ["--center-sector", str(sector_id)]
    elif by == "address":
        ring = _number(form, "center_ring", "Ring", int, required=True, minimum=0, maximum=MAX_GENERATE_RING)
        layer = _number(form, "center_layer", "Layer", int) or 0
        slot = _number(form, "center_slot", "Slot", int, required=True, minimum=0)
        description = f"around ring {ring} layer {layer} slot {slot}"
        argv = ["--ring", str(ring), "--layer", str(layer), "--slot", str(slot)]
    elif by == "position":
        point = (_number(form, "center_x_pc", "x (pc)", float, required=True),
                 _number(form, "center_y_pc", "y (pc)", float, required=True),
                 _number(form, "center_z_pc", "z (pc)", float) or 0.0)
        edge = float(edge_pc or program_constants.DEFAULT_SECTOR_EDGE_PC)
        if max(abs(v) for v in point) > (MAX_GENERATE_RING + 1) * edge:
            raise FormError("That position is outside the galaxy grid.")
        ring, layer, slot = sector_address_at(point, edge)
        if ring > MAX_GENERATE_RING:
            raise FormError("That position is outside the galaxy grid.")
        description = (f"around {', '.join(f'{v:g}' for v in point)} pc "
                       f"(ring {ring} layer {layer} slot {slot})")
        argv = ["--ring", str(ring), "--layer", str(layer), "--slot", str(slot)]
    else:
        raise FormError("Choose how to name the center sector.")
    radius = _number(form, "center_radius_pc", "Radius (pc)", float, required=True, minimum=0.1,
                     maximum=MAX_GENERATE_RADIUS_PC)
    return argv + ["--radius-pc", str(radius)], description


def galaxy_argv(form, edge_pc=None):
    """
    `generate.py galaxy` arguments for the chosen mode. `edge_pc` is the
    plan's sector edge, which only "Around a sector" by position needs.

    Returns:
        tuple: `(argv, description)`.
    """
    mode = form.get("mode") or "random"
    if mode == "random":
        return random_start_argv(form), "around a random start"
    if mode == "ring":
        ring = _number(form, "ring", "Ring", int, required=True, minimum=0, maximum=MAX_GENERATE_RING)
        layer = _number(form, "layer", "Layer", int)
        argv = ["--ring", str(ring), "--layer", str(layer or 0)]
        limit = _number(form, "limit", "Limit", int, minimum=1, maximum=MAX_GENERATE_LIMIT)
        if limit is not None:
            argv += ["--limit", str(limit)]
        elif form.get("whole_ring"):
            argv.append("--yes")
        return argv, f"in ring {ring} layer {layer or 0}"
    if mode == "center":
        return center_argv(form, edge_pc)
    if mode == "slot":
        ring = _number(form, "slot_ring", "Ring", int, required=True, minimum=0, maximum=MAX_GENERATE_RING)
        layer = _number(form, "slot_layer", "Layer", int) or 0
        slot = _number(form, "slot", "Slot", int, required=True, minimum=0)
        argv = ["--ring", str(ring), "--layer", str(layer), "--slot", str(slot)]
        description = f"at ring {ring} layer {layer} slot {slot}"
        radius = _neighborhood_radius_pc(form)
        if radius is not None:
            argv += ["--radius-pc", str(radius)]
            description = f"around ring {ring} layer {layer} slot {slot}"
        return argv, description
    if mode == "column":
        ring = _number(form, "column_ring", "Ring", int, required=True, minimum=0, maximum=MAX_GENERATE_RING)
        slot = _number(form, "column_slot", "Slot", int, required=True, minimum=0)
        return (["--ring", str(ring), "--slot", str(slot), "--column"],
                f"in the column at ring {ring} slot {slot}")
    if mode == "block":
        key = (form.get("block") or "").strip()
        if not key:
            raise FormError("Block is required.")
        try:
            block = parse_drill_key(key)
        except ValueError:
            raise FormError("Block must be a Galaxy Map block key, size.ring.wedge.slab (e.g. 3.40.7.0).") from None
        if block.m == 1:
            raise FormError("That is a single sector; use One address instead.")
        argv = ["--block", format_drill_key(block)]
        description = f"in block {format_drill_key(block)}"
        layer = _number(form, "block_layer", "Layer", int)
        if layer is not None:
            argv += ["--block-layer", str(layer)]
            description += f" layer {layer}"
        limit = _number(form, "block_limit", "Limit", int, minimum=1, maximum=MAX_GENERATE_LIMIT)
        if limit is not None:
            argv += ["--limit", str(limit)]
        elif form.get("whole_block"):
            argv.append("--yes")
        return argv, description
    if mode == "shell":
        ring = _number(form, "shell_ring", "Ring", int, required=True, minimum=0, maximum=MAX_GENERATE_RING)
        argv = ["--ring", str(ring), "--shell"]
        limit = _number(form, "shell_limit", "Limit", int, minimum=1, maximum=MAX_GENERATE_LIMIT)
        if limit is not None:
            argv += ["--limit", str(limit)]
        elif form.get("whole_shell"):
            argv.append("--yes")
        return argv, f"in the shell at ring {ring}"
    raise FormError("Choose what to generate.")


def plan_steps(generate, form):
    """
    The plan step and, unless the form's "skip the bright-star scatter"
    box is ticked, the scatter as a second step (`generate.py plan` would
    otherwise run both in one command, with no separate label).

    Returns:
        list[dict]: Job steps.
    """
    steps = [{"label": "Plan the galaxy", "argv": generate + ["plan"] + plan_argv(form) + ["--no-bright-stars"]}]
    if not form.get("skip_bright_stars"):
        steps.append({"label": SCATTER_LABEL, "argv": generate + ["plan", "--bright-stars-only"]})
    return steps


def build_job(action, form, database, edge_pc=None):
    """
    The job a form asks for (`edge_pc`: see `galaxy_argv`).

    Returns:
        tuple: `(kind, title, steps)` for `jobs.start_job`.

    Raises:
        FormError: A bad value, or a missing/wrong typed confirmation.
    """
    if action in CONFIRM_ACTIONS and (form.get("confirm") or "").strip() != database:
        raise FormError(f"Type the database name ({database}) to confirm. Nothing was changed.")
    python = jobs.python_executable()
    generate = [python, jobs.GENERATE_SCRIPT]
    reset_step = {"label": "Reset the galaxy", "argv": [python, jobs.RESET_SCRIPT, "--yes"]}
    if action == "new_galaxy":
        plan = plan_steps(generate, form)
        start = random_start_argv(form)
        return "new_galaxy", "New galaxy", [
            reset_step,
            *plan,
            {"label": "Generate sectors around a random start", "argv": generate + ["galaxy"] + start},
        ]
    if action == "plan":
        return "plan", "Plan the galaxy", plan_steps(generate, form)
    if action == "bright_stars":
        argv = generate + ["plan", "--bright-stars-only"]
        if form.get("bright_force"):
            argv.append("--force")
        return "bright_stars", "Rebuild the bright stars", [{"label": SCATTER_LABEL, "argv": argv}]
    if action == "bright_band":
        down_to = _number(form, "down_to", "Go down to (solar luminosities)", float, required=True, minimum=1.0)
        argv = generate + ["plan", "--bright-stars-down-to", f"{down_to:g}"]
        return "bright_band", f"Bright stars down to {down_to:g} L\u2609", [{"label": BAND_LABEL, "argv": argv}]
    if action == "galaxy":
        argv, description = galaxy_argv(form, edge_pc)
        label = f"Generate sectors {description}"
        return "galaxy", label, [{"label": label, "argv": generate + ["galaxy"] + argv}]
    if action == "reset":
        return "reset", "Reset the galaxy", [reset_step]
    raise FormError("Unknown action.")


# ---------------------------------------------------------------------
# Views
# ---------------------------------------------------------------------

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
        activitylog.event("AUTHZ", "admin.required", user=admin.get("username") if admin else None,
                          path=request.path)
        abort(403, description="Only a logged-in admin can do that.")
    return admin


def _galaxy_summary(database):
    """Whether the galaxy is planned and how many sectors it holds.
    Fails open: the page still works when the database doesn't answer."""
    summary = {"shape": None, "sectors": None, "bright": None, "error": None}
    try:
        summary["shape"] = apiclient.get_galaxy_shape(database)
        summary["bright"] = apiclient.get_bright_star_status(database)
        summary["sectors"] = apiclient.get_sectors(database, limit=1, offset=0)["total"]
    except (apiclient.ApiError, apiclient.NotFoundError) as exc:
        summary["error"] = str(exc)
    return summary


def _edge_pc(action, form, database):
    """The plan's sector edge in parsecs, asked for only when "Around a
    sector" was given a position (the one form that needs it); `None`
    otherwise, or when the database doesn't answer (the standard edge is
    then used)."""
    if action != "galaxy" or form.get("mode") != "center" or form.get("center_by") != "position":
        return None
    try:
        shape = apiclient.get_galaxy_shape(database)
    except (apiclient.ApiError, apiclient.NotFoundError) as exc:
        log.error(f"Generate page: no galaxy shape for a center position: {exc}")
        return None
    return (shape or {}).get("edge_pc")


def _kept_section(form):
    """
    The section of a form shown again (with its error or its estimate),
    which the server renders open and `static/generatefolds.js` leaves
    open whatever the browser remembered; `None` on a plain visit.
    """
    return SECTION_FOR_ACTION.get(form.get("action")) if form else None


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
    # The run's own decaying-average ETA (PERF.7), counted down from when
    # it was written, so the page and the terminal agree.
    eta_s = progress.get("eta_s")
    remaining = None
    if eta_s is not None and not job.get("finished"):
        updated_at = progress.get("updated_at") or time.time()
        remaining = max(0.0, float(eta_s) - max(0.0, time.time() - float(updated_at)))
    view["remaining_text"] = format_elapsed(remaining)
    view["progress_text"] = _progress_text(progress)
    # PERF.4: the second bar (the bright-star layers being drawn while
    # layers are slow), or blank.
    detail = progress.get("detail") if not job.get("finished") else None
    view["progress_detail_text"] = ""
    if detail:
        text = _progress_text(detail)
        if detail.get("eta_s") is not None:
            text += f", about {format_elapsed(detail['eta_s'])} left"
        view["progress_detail_text"] = text
    return view


def _progress_text(progress):
    """`"Sectors: 12 of 40"`, or `"Bright stars (3 of 81 layers): 12%"` for
    a bar written as a share (PERF.9), or `""` before the first report."""
    description = progress.get("description")
    if not description:
        return ""
    completed, total = progress.get("completed") or 0, progress.get("total")
    if progress.get("percent") and total:
        return f"{description}: {100 * float(completed) / float(total):.0f}%"
    text = f"{description}: {format_number(completed)}"
    if total:
        text += f" of {format_number(total)}"
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


def _page(admin, error=None, status=200, form=None, estimate=None, estimate_title=None, estimate_fields=()):
    database = db_name()
    try:
        root = jobs.jobs_dir()
        active = jobs.active_job(root)
        recent = jobs.list_jobs(limit=10, root=root)
        jobs_error = None
    except OSError as exc:
        active, recent, jobs_error = None, [], str(exc)
        log.error(f"Generate page: {exc}")
    return _no_store(render_page(
        "generate.html",
        title="Generate",
        section="admin_generate",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Generate")],
        description="Generate, plan or reset this site's galaxy.",
        admin=admin,
        database=database,
        summary=_galaxy_summary(database),
        active=_job_view(active),
        active_log=jobs.log_tail(active["id"], max_bytes=8 * 1024) if active else None,
        recent=[_job_view(job) for job in recent],
        jobs_error=jobs_error,
        plan_fields=PLAN_FIELDS,
        bright_min_luminosity=program_constants.BRIGHT_STAR_MIN_LUMINOSITY_SOL,
        galaxy_modes=GALAXY_MODES,
        max_radius_pc=MAX_GENERATE_RADIUS_PC,
        min_radius_ly=MIN_NEIGHBORHOOD_RADIUS_LY,
        max_radius_ly=MAX_GENERATE_RADIUS_LY,
        max_ring=MAX_GENERATE_RING,
        max_limit=MAX_GENERATE_LIMIT,
        center_choices=CENTER_CHOICES,
        kept_section=_kept_section(form),
        error=error,
        form=form or {},
        estimate=estimate,
        estimate_title=estimate_title,
        estimate_fields=estimate_fields,
        estimate_confirm_field=ESTIMATE_CONFIRM_FIELD,
        status=status,
    ))


@bp.route("/admin/generate", methods=["GET", "POST"])
def generate():
    """The Generate page (GET) and its forms (POST, answered with a 303
    back to the page, or the page again with an error when nothing
    started)."""
    if request.method == "GET":
        admin, response = _admin_or_redirect()
        return response or _page(admin)

    admin = _admin_or_403()
    action = request.form.get("action") or ""
    if action == "cancel":
        job_id = request.form.get("job") or ""
        try:
            if jobs.cancel_job(job_id):
                activitylog.event("GEN", "job.cancel", user=admin.get("username"), job=job_id)
        except OSError as exc:
            log.error(f"Could not cancel job {job_id}: {exc}")
        return _no_store(redirect(url_for("web.generate", _anchor="current-job"), code=303))

    database = db_name()
    wants_json = _wants_json()
    try:
        kind, title, steps = build_job(action, request.form, database, edge_pc=_edge_pc(action, request.form, database))
    except FormError as exc:
        if wants_json:
            return _no_store(make_response(jsonify({"error": str(exc)}), 400))
        return _page(admin, error=str(exc), status=400, form=request.form)
    env = jobs.mysql_env(current_app.config["MYSQL_CONFIG"], database)
    if kind == "galaxy" and not request.form.get(ESTIMATE_CONFIRM_FIELD):
        # PERF.3: show the size and time first; the admin confirms (or the
        # disk refuses it) before the job starts.
        try:
            estimate = run_estimate(steps[0]["argv"], env)
        except FormError as exc:
            if wants_json:
                return _no_store(make_response(jsonify({"error": str(exc)}), 400))
            return _page(admin, error=str(exc), status=400, form=request.form)
        if wants_json:
            return _no_store(make_response(jsonify({
                "error": estimate["refusal"] or f"Confirm first: {estimate['summary']}",
                "estimate": estimate, "confirm_field": ESTIMATE_CONFIRM_FIELD,
            }), 409))
        return _page(admin, form=request.form, estimate=_estimate_view(estimate), estimate_title=title,
                     estimate_fields=[(name, value) for name, value in request.form.items(multi=True)
                                      if name not in ("csrf_token", ESTIMATE_CONFIRM_FIELD)])
    try:
        job_id = jobs.start_job(kind, title, steps, env=env, admin=admin.get("username"), database=database)
        activitylog.event("GEN", "job.start", user=admin.get("username"), job=job_id, kind=kind, db=database,
                          title=title)
    except jobs.JobBusy as exc:
        running = exc.job["title"] if exc.job else "Another job"
        message = f"{running} is still running. Wait for it to finish, or cancel it first."
        if wants_json:
            return _no_store(make_response(jsonify({"error": message}), 409))
        return _page(admin, error=message, status=409, form=request.form)
    except OSError as exc:
        log.exception(f"Could not start a {kind} job: {exc}")
        if wants_json:
            return _no_store(make_response(jsonify({"error": f"The job could not be started: {exc}"}), 500))
        return _page(admin, error=f"The job could not be started: {exc}", status=500, form=request.form)
    if wants_json:
        return _no_store(make_response(jsonify({
            "job": job_id,
            "url": url_for("web.generate_job", job_id=job_id),
            "status_url": url_for("web.generate_status", job=job_id),
        }), 202))
    return _no_store(redirect(url_for("web.generate", _anchor="current-job"), code=303))


def _estimate_view(estimate):
    """`run_estimate`'s dict plus the text the page shows."""
    view = dict(estimate)
    view["size_text"] = generationStats.format_bytes(estimate.get("bytes"))
    view["time_text"] = generationStats.format_duration(estimate.get("seconds"))
    disk = estimate.get("disk")
    view["disk_text"] = (f"{generationStats.format_bytes(disk['free_bytes'])} free of "
                         f"{generationStats.format_bytes(disk['total_bytes'])}") if disk else None
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
    body = {"job": _job_view(job), "log": jobs.log_tail(job["id"], max_bytes=8 * 1024, root=root) if job else None}
    return _no_store(jsonify(body))


generate_status.json_only = True  # not a page: tests/test_web_a11y.py skips it


@bp.route("/admin/generate/sectors")
def generate_sectors():
    """
    JSON for "Around a sector"'s list of filled sectors
    (`static/generatefolds.js`): `?page=N`, 50 a page, nearest the
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
        "ring": row.get("ring_index"), "layer": row.get("layer_index"), "slot": row.get("ring_slot_index"),
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
        output=jobs.log_tail(job_id, max_bytes=256 * 1024, root=root),
    ))
