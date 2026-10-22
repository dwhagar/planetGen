"""
The Generate screens (ADM.108): the Start screen with one card per task, a page for each task, the Review step every
task ends in, and the History tab. `web/generate_page.py` keeps the logic they share (the argv builders, the job
views, the estimate and the `POST /admin/generate` that starts a job); this module is how an admin gets there.

The flow, for every task:

1. `GET /admin/generate` lists the tasks; `GET /admin/generate/task/<key>` is a task's form.
2. The form posts to its own page. A valid draft answers with the Review step (what will run, the command lines, the
   size and time of a galaxy run), a bad one with the form again and the message.
3. The Review step's Run posts the same fields to `POST /admin/generate`, which starts the job (or queues it behind
   the running one, `jobs.queue_next`). Back posts them to the task's page, which fills the form again. The draft
   travels in hidden fields, so the page works without script.
"""

import shlex
from dataclasses import dataclass
from urllib.parse import urlencode

from flask import abort, current_app, jsonify, request

from . import adminregistry, bp, generate_commands, jobs
from .generate_page import (
    ESTIMATE_CONFIRM_FIELD, GENERATE_ANYWAY_FIELD, FormError, _admin_or_403, _admin_or_redirect, _edge_pc,
    _estimate_view, _no_store, build_job, estimate_for, estimate_needed, page_context,
)
from .helpers import crumb, db_name, page_url, render_page


@dataclass(frozen=True)
class Task:
    """One card of the Start screen. `actions` are the `POST /admin/generate` actions its page offers (one: the page
    posts it itself; several: the page has a choice)."""

    key: str
    label: str
    summary: str
    danger: str
    actions: tuple
    endpoint: str = "generate_task"

    @property
    def params(self):
        return {"key": self.key} if self.endpoint == "generate_task" else {}


TASKS = (
    Task("add-sectors", "Add sectors to the galaxy",
         "Generate more sectors, around a place or by shells, bands and wedges.", "safe", ("galaxy",)),
    Task("new-galaxy", "Start a new galaxy",
         "Delete everything, plan a new galaxy and generate a first neighborhood.", "destroys", ("new_galaxy",)),
    Task("scatters", "Change what is scattered across the galaxy",
         "Plan again, redo some scatters, or add a dimmer step of bright stars. Sectors already generated are left alone.",
         "rewrites", ("plan", "redo_scatters", "bright_band")),
    Task("one-off-system", "Make one system (not saved)",
         "Generate a single star system with any command-line option and get it as Markdown or wikitext.", "safe", (),
         endpoint="generate_system"),
    Task("check", "Check the database", "Look for damage. It changes nothing.", "safe", ("check_db",)),
    Task("wipe", "Wipe the galaxy",
         "Delete every sector, system, phenomenon and the plan. Admin logins and API keys are kept.", "destroys",
         ("reset",)),
)
"""tuple: The tasks, in card order (the wipe card is drawn last, apart from the rest)."""

DESTRUCTIVE_ACTIONS = frozenset({"new_galaxy", "reset"})

RUN_LABELS = {
    "new_galaxy": "Delete everything and start a new galaxy",
    "reset": "Delete everything",
    "check_db": "Check",
}

DRAFT_SKIP = frozenset({"csrf_token", ESTIMATE_CONFIRM_FIELD, GENERATE_ANYWAY_FIELD, "back", "review"})
"""frozenset: Form fields that are not part of the draft carried between the task page and its Review step."""


NEVER_KEPT = frozenset({"confirm"})
"""frozenset: Fields a job does not keep for Run again: the typed database name that confirms a deletion."""


def draft_of(action, form):
    """The fields of a started job's form that Run again fills a task's page with (`[name, value]` pairs), or `[]`
    for one that cannot be run again (it destroys data and asks for a typed confirmation each time)."""
    if action in DESTRUCTIVE_ACTIONS:
        return []
    return [[name, value] for name, value in form.items(multi=True) if name not in DRAFT_SKIP | NEVER_KEPT]


def rerun_url(job):
    """Where "Run again with the same settings" goes for a `jobs.get_job` dict: the task's Review step filled with the
    job's draft (nothing starts until Run is pressed there), or `None` for a job with no kept draft."""
    draft = job.get("draft") or []
    action = next((value for name, value in draft if name == "action"), None)
    task = task_for_action(action) if action else None
    if task is None or action in DESTRUCTIVE_ACTIONS:
        return None
    return page_url(task.endpoint, **task.params) + "?" + urlencode([*draft, ("review", "1")])


def task_for_key(key):
    for task in TASKS:
        if task.key == key and task.actions:
            return task
    return None


def task_for_action(action):
    for task in TASKS:
        if action in task.actions:
            return task
    return None


def task_entries():
    """The tasks as registry entries (`adminregistry.tasks`), so search and the palette open them."""
    return [
        adminregistry.Entry(
            key=f"task:{task.key}", kind="task", label=task.label, summary=task.summary, destination="generate",
            endpoint=task.endpoint, params=task.params, keywords=("generate", "run", "make"),
            danger={"safe": "safe", "rewrites": "rewrites", "destroys": "destroys"}[task.danger],
        )
        for task in TASKS
    ]


def render_screen(template, admin, title, crumbs, form=None, error=None, status=200, tab="start", description=None,
                  **context):
    """`template` drawn with the shared Generate context; `crumbs` follow Admin and Generate (none: Generate is the
    current page)."""
    context = {**page_context(admin, form=form, error=error, status=status), **context}
    trail = [crumb("Admin", "admin"), crumb("Generate", "generate") if crumbs else crumb("Generate"), *crumbs]
    return _no_store(render_page(
        template, title=title, section="admin_generate", tab=tab, breadcrumbs=trail,
        description=description or "Generate, plan or reset this site's galaxy.", tasks=TASKS, **context))


def start_page(admin, error=None, status=200):
    """The Start screen: a card per task."""
    return render_screen("generate.html", admin, "Generate", [], error=error, status=status)


def task_page(admin, task, form=None, error=None, status=200):
    """A task's form, filled from `form` when it is shown again."""
    return render_screen("generate_task.html", admin, task.label, [crumb(task.label)], form=form, error=error, status=status,
                   task=task)


def task_page_for_action(admin, action, form=None, error=None, status=200):
    """The page of the task that owns `action`, else the Start screen (an unknown action)."""
    task = task_for_action(action)
    if task is None:
        return start_page(admin, error=error, status=status)
    return task_page(admin, task, form=form, error=error, status=status)


def _command(argv):
    """A step's command as an admin would type it: the interpreter's own path shortened."""
    words = [shlex.quote(str(word)) for word in argv]
    if words and words[0].rsplit("/", 1)[-1].startswith("python"):
        words[0] = "python3"
    return " ".join(words)


def review_page(admin, action, kind, title, steps, form, estimate=None):
    """The Review step for a built job: its steps with their command lines, the estimate if there is one, and Run."""
    task = task_for_action(action)
    if task is None:
        abort(400)
    review = {
        "title": title,
        "steps": [{"label": step["label"], "command": _command(step["argv"])} for step in steps],
        "estimate": estimate,
        "destructive": action in DESTRUCTIVE_ACTIONS,
        "run_label": RUN_LABELS.get(action) or (
            "Generate anyway" if estimate and estimate.get("refused") else "Run"),
        "fields": [(name, value) for name, value in form.items(multi=True) if name not in DRAFT_SKIP],
    }
    return render_screen("generate_review.html", admin, f"Review: {title}", [crumb(task.label, "generate_task", key=task.key),
                                                                         crumb("Review")],
                   form=form, task=task, review=review)


@bp.route("/admin/generate/task/<key>", methods=["GET", "POST"])
def generate_task(key):
    """
    A task's form (GET), and its Review step (POST; with `back` it is the form again, filled). A GET with query
    fields fills the form with them, and with `review=1` shows the Review step for them, which is where a typed command
    lands (`generate_commands`); nothing runs until Run is pressed on that step.
    """
    task = task_for_key(key)
    if task is None:
        abort(404)
    if request.method in ("GET", "HEAD"):
        admin, response = _admin_or_redirect()
        if response:
            return response
        form = request.args
        if not form:
            return task_page(admin, task)
        if form.get("review") and not form.get("back"):
            return _review_or_form(admin, task, form)
        return task_page(admin, task, form=form)
    admin = _admin_or_403()
    form = request.form
    if form.get("back"):
        return task_page(admin, task, form=form)
    return _review_or_form(admin, task, form)


def _review_or_form(admin, task, form):
    """The Review step for the draft `form`, or the task's form again with what is wrong."""
    action = form.get("action") or (task.actions[0] if len(task.actions) == 1 else "")
    if action not in task.actions:
        return task_page(admin, task, form=form, error="Choose what to do first.", status=400)
    database = db_name()
    try:
        kind, title, steps = build_job(action, form, database, edge_pc=_edge_pc(action, form, database))
        estimate = None
        if estimate_needed(kind, form):
            env = jobs.mysql_env(current_app.config["MYSQL_CONFIG"], database)
            estimate = _estimate_view(estimate_for(kind, steps, database, env))
    except FormError as exc:
        return task_page(admin, task, form=form, error=str(exc), status=400)
    return review_page(admin, action, kind, title, steps, form, estimate=estimate)


@bp.route("/admin/generate/command")
def generate_command():
    """
    What a typed command line means (`generate_commands.parse`), for the palette: JSON with `command` (its label, summary
    and the address it opens) or `error`, or `command: null` for text that is not a command.
    """
    admin, response = _admin_or_redirect()
    if response:
        return jsonify({"error": "sign in"}), 403
    try:
        found = command_link(request.args.get("q", ""))
    except generate_commands.CommandError as exc:
        return _no_store(jsonify({"command": None, "error": str(exc)}))
    return _no_store(jsonify({"command": found}))


generate_command.json_only = True  # not a page: tests/test_web_a11y.py skips it


def command_link(text):
    """`{"label", "summary", "url"}` for a typed line that is a command, `None` for one that is not.

    Raises:
        generate_commands.CommandError: A line that starts like a command but cannot be read.
    """
    command = generate_commands.parse(text)
    if command is None:
        return None
    task = next(t for t in TASKS if t.key == command.task)
    query = {**command.fields, **({"review": "1"} if command.review else {})}
    url = page_url(task.endpoint, **task.params)
    return {"label": f"{'Review' if command.review else 'Open'}: {task.label}", "summary": command.summary,
            "url": url + ("?" + urlencode(query) if query else "")}


@bp.route("/admin/generate/history")
def generate_history():
    """Every recent job, a page at a time."""
    admin, response = _admin_or_redirect()
    return response or render_screen("generate_history.html", admin, "Generate history", [crumb("History")], tab="history")
