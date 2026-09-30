# html/web/generate_page.py

"""
`/admin/generate`: galaxy generation from the browser, so an admin never
needs a terminal on the server. Four actions, each a background job
(`web/jobs.py`) running the same command-line tools an admin would type:

- New galaxy: wipe the database (`src/resetDb.py --yes`), build the
  density skeleton (`generate.py plan`), then generate a first
  neighborhood around a random start (`generate.py galaxy`).
- Plan: rebuild the density skeleton only.
- Generate sectors: `generate.py galaxy` in any of its four modes
  (random start, a whole ring, around a sector, one address).
- Reset: wipe the database only.

New galaxy and Reset delete every generated sector and system, so both
ask for the database name to be typed back, the same confirmation
`resetDb.py` asks for in a terminal.

The page shows the running job (step, progress bar, elapsed time and the
tail of its output), refreshed every few seconds by
`static/generatejobs.js` from `/admin/generate/status`, plus the last few
jobs, each with its full output at `/admin/generate/jobs/<id>`.

Admins only: a visitor who isn't logged in is sent to the login page,
and any POST or status request without an admin session gets a 403.
Every form carries `csrf_field()` (checked app-wide by `csrf.protect`).
"""

import time

from flask import abort, current_app, jsonify, make_response, redirect, request, url_for

import apiclient
from stellarObjects import log

from . import bp, jobs
from .helpers import crumb, current_admin, db_name, page_url, render_page

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
    ("edge_ly", "--edge-ly", "Sector edge (ly)", float, 11.5, 1.0, None),
)
"""tuple: The `plan` options the page offers, with `generate.py`'s own
defaults (a test checks they still match `generate.py plan`'s parser)."""

GALAXY_MODES = (
    ("random", "Around a random start",
     "Picks a random populated spot and generates the sectors within a radius of it "
     "(100 ly when the radius is left blank)."),
    ("ring", "A whole ring",
     "Every not-yet-generated sector in one ring at one height layer (0 is the galactic plane), "
     "or only the first few with a limit."),
    ("center", "Around a sector",
     "Every not-yet-generated sector within a radius of an existing sector."),
    ("slot", "One address",
     "Exactly one sector, by ring, layer and slot (the address the Galaxy Map shows)."),
)

CONFIRM_ACTIONS = frozenset({"new_galaxy", "reset"})


class FormError(ValueError):
    """A form value the page rejects; the message is shown to the admin."""


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
    radius = _number(form, "radius_pc", "Radius (pc)", float, minimum=0.1)
    if radius is not None:
        argv += ["--radius-pc", str(radius)]
    max_ring = _number(form, "max_ring", "Highest ring", int, minimum=0)
    if max_ring is not None:
        argv += ["--max-ring", str(max_ring)]
    density = _number(form, "min_start_density", "Minimum start density", float, minimum=0.001)
    if density is not None:
        argv += ["--min-start-density", str(density)]
    return argv


def galaxy_argv(form):
    """
    `generate.py galaxy` arguments for the chosen mode.

    Returns:
        tuple: `(argv, description)`.
    """
    mode = form.get("mode") or "random"
    if mode == "random":
        return random_start_argv(form), "around a random start"
    if mode == "ring":
        ring = _number(form, "ring", "Ring", int, required=True, minimum=0)
        layer = _number(form, "layer", "Layer", int)
        argv = ["--ring", str(ring), "--layer", str(layer or 0)]
        limit = _number(form, "limit", "Limit", int, minimum=1)
        if limit is not None:
            argv += ["--limit", str(limit)]
        elif form.get("whole_ring"):
            argv.append("--yes")
        return argv, f"in ring {ring} layer {layer or 0}"
    if mode == "center":
        sector_id = _number(form, "center_sector", "Sector ID", int, required=True, minimum=1)
        radius = _number(form, "center_radius_pc", "Radius (pc)", float, required=True, minimum=0.1)
        return ["--center-sector", str(sector_id), "--radius-pc", str(radius)], f"around sector {sector_id}"
    if mode == "slot":
        ring = _number(form, "slot_ring", "Ring", int, required=True, minimum=0)
        layer = _number(form, "slot_layer", "Layer", int) or 0
        slot = _number(form, "slot", "Slot", int, required=True, minimum=0)
        return (["--ring", str(ring), "--layer", str(layer), "--slot", str(slot)],
                f"at ring {ring} layer {layer} slot {slot}")
    raise FormError("Choose what to generate.")


def build_job(action, form, database):
    """
    The job a form asks for.

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
        plan = plan_argv(form)
        start = random_start_argv(form)
        return "new_galaxy", "New galaxy", [
            reset_step,
            {"label": "Plan the galaxy", "argv": generate + ["plan"] + plan},
            {"label": "Generate sectors around a random start", "argv": generate + ["galaxy"] + start},
        ]
    if action == "plan":
        return "plan", "Plan the galaxy", [
            {"label": "Plan the galaxy", "argv": generate + ["plan"] + plan_argv(form)},
        ]
    if action == "galaxy":
        argv, description = galaxy_argv(form)
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
        abort(403, description="Only a logged-in admin can do that.")
    return admin


def _galaxy_summary(database):
    """Whether the galaxy is planned and how many sectors it holds.
    Fails open: the page still works when the database doesn't answer."""
    summary = {"shape": None, "sectors": None, "error": None}
    try:
        summary["shape"] = apiclient.get_galaxy_shape(database)
        summary["sectors"] = apiclient.get_sectors(database, limit=1, offset=0)["total"]
    except (apiclient.ApiError, apiclient.NotFoundError) as exc:
        summary["error"] = str(exc)
    return summary


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
    view["created_text"] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(job["created_at"])) \
        if job.get("created_at") else ""
    view["elapsed_text"] = format_elapsed(job.get("elapsed_s"))
    return view


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


def _page(admin, error=None, status=200, form=None):
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
        galaxy_modes=GALAXY_MODES,
        error=error,
        form=form or {},
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
            jobs.cancel_job(job_id)
        except OSError as exc:
            log.error(f"Could not cancel job {job_id}: {exc}")
        return _no_store(redirect(url_for("web.generate", _anchor="current-job"), code=303))

    database = db_name()
    try:
        kind, title, steps = build_job(action, request.form, database)
    except FormError as exc:
        return _page(admin, error=str(exc), status=400, form=request.form)
    env = jobs.mysql_env(current_app.config["MYSQL_CONFIG"], database)
    try:
        jobs.start_job(kind, title, steps, env=env, admin=admin.get("username"), database=database)
    except jobs.JobBusy as exc:
        running = exc.job["title"] if exc.job else "Another job"
        return _page(admin, error=f"{running} is still running. Wait for it to finish, or cancel it first.",
                     status=409, form=request.form)
    except OSError as exc:
        log.exception(f"Could not start a {kind} job: {exc}")
        return _page(admin, error=f"The job could not be started: {exc}", status=500, form=request.form)
    return _no_store(redirect(url_for("web.generate", _anchor="current-job"), code=303))


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
