"""
The admin registry (ADM.118): one list of everything an admin can open or do, which the shell's rail, the Ctrl+K
palette (ADM.120) and the plain `/admin/search` page all read.

Three kinds of entry:

* `page`: an admin page (`endpoint` is its route name).
* `task`: a job the Generate page starts; opening it lands on that task's card with its Review step.
* `setting`: one option of the settings model (`util/settings_view.describe`), opening the Settings page at that
  option. These are built from the model, never listed here, so a new option is findable the day it exists.

Each entry says how risky it is (`danger`: `"safe"` changes nothing, `"rewrites"` changes data or settings,
`"destroys"` removes data) and which rail destination holds it. A page or task is declared once, here; the test
`tests/test_web_adminregistry.py` fails when an `/admin` page route is in neither `PAGES` nor `EXEMPT`, so a new
page cannot be left out of search.
"""

from dataclasses import dataclass, field
from urllib.parse import urlencode

DESTINATIONS = (
    ("overview", "Overview", "admin"),
    ("generate", "Generate", "generate"),
    ("settings", "Settings", "admin_settings"),
    ("jobs", "Jobs", "admin_queue"),
    ("data", "Data", "admin_stats"),
    ("account", "Account", "account"),
)
"""tuple: The rail, in order: `(key, label, endpoint)`. The rail link of a destination is its endpoint's page."""

SECTION_DESTINATION = {
    "admin": "overview",
    "admin_generate": "generate",
    "admin_settings": "settings",
    "admin_queue": "jobs",
    "admin_stats": "data",
    "admin_facility_types": "data",
    "account": "account",
}
"""dict: A page's `render_page(section=...)` -> the rail destination it lights. Only these sections wear the shell."""

DANGER = ("safe", "rewrites", "destroys")


@dataclass(frozen=True)
class Entry:
    """One thing an admin can open or do. `endpoint` is a route name (no `web.` prefix); `params` are its URL
    parameters and `query` its query string, `anchor` the element to land on."""

    key: str
    kind: str
    label: str
    summary: str
    destination: str
    endpoint: str
    keywords: tuple = ()
    params: dict = field(default_factory=dict)
    query: dict = field(default_factory=dict)
    anchor: str = ""
    danger: str = "safe"

    def text(self):
        """Lower-case text the search matches against."""
        return " ".join((self.label, self.summary, " ".join(self.keywords), self.kind)).lower()


def _page(key, label, summary, destination, endpoint, *keywords, danger="safe", anchor="", query=None):
    return Entry(key, "page", label, summary, destination, endpoint, tuple(keywords), anchor=anchor,
                 query=query or {}, danger=danger)


PAGES = (
    _page("overview", "Overview", "The admin home: site status, maintenance mode, naming key and links.", "overview",
          "admin", "home", "maintenance", "naming key"),
    _page("generate", "Generate", "Make or change the galaxy: sectors, scatters and redos.", "generate",
          "generate", "make", "create", "galaxy", "sectors", "run", danger="rewrites"),
    _page("generate-system", "Generate a system", "Make one system from the settings you choose.", "generate",
          "generate_system", "system", "star", "planet", danger="safe"),
    _page("generate-history", "Generate history", "Every recent job, with its status and time.", "generate",
          "generate_history", "jobs", "recent", "log"),
    _page("settings", "Settings", "Every site setting you can change, in plain groups.", "settings",
          "admin_settings", "options", "config", "configuration"),
    _page("settings-backups", "Settings backups", "Earlier copies of the settings file, with restore.", "settings",
          "admin_backups", "restore", "backup", "undo", danger="rewrites"),
    _page("settings-images", "Site images", "Upload the logo, icons and link-preview image.", "settings",
          "admin_images", "logo", "icon", "favicon", "preview", "upload", danger="rewrites"),
    _page("jobs", "Jobs and queue", "Background work waiting and running, and the cache warm-up.", "jobs",
          "admin_queue", "queue", "workers", "cache", "warm", "redis"),
    _page("stats", "Stats", "Counts, names, lockouts and sign-in failures.", "data", "admin_stats",
          "statistics", "lockouts", "login failures", "duplicate names"),
    _page("facility-types", "Facility types", "The facility types the generator may place.", "data",
          "admin_facility_types", "facilities", "types", danger="rewrites"),
    _page("account", "Account", "Your sign-in, password, two-factor codes and API keys.", "account", "account",
          "password", "two factor", "2fa", "api key", "profile"),
)

EXEMPT = {
    "generate_status": "JSON for the Current job panel",
    "generate_job": "A job's own page, reached from History",
    "generate_job_stream": "Server-sent events for a job",
    "generate_job_log": "A job's raw log",
    "generate_sectors": "JSON for the sector picker",
    "generate_shape_report": "JSON for the line under the shape dials",
    "generate_command": "JSON for the palette's typed commands",
    "generate_system_download": "A file download",
    "admin_queue_tree": "A queue entry's detail, reached from the queue",
    "admin_queue_confirm": "A confirmation step, reached from the queue",
    "admin_queue_action": "A POST",
    "admin_backup": "One backup's detail, reached from Settings backups",
    "admin_backup_restore": "A POST",
    "admin_settings_export": "A file download",
    "admin_settings_import": "A POST",
    "admin_settings_import_save": "A POST",
    "admin_search": "The search page itself",
    "admin_palette": "JSON for the palette",
    "reset_generation_stats": "A POST",
    "change_naming_key": "A POST",
    "lift_lockout": "A POST",
    "account_two_factor": "A POST",
}
"""dict: Admin routes that are not pages of their own, with why. The registry test reads this."""

def tasks():
    """The Generate page's tasks as entries of kind `task` (`web/generate_flow.TASKS`)."""
    from planetgen.web import generate_flow

    return generate_flow.task_entries()


def setting_entries(role="admin", environ=None):
    """One `setting` entry per option the role may change, built from the settings model."""
    from planetgen.util import settings_view

    entries = []
    for option in settings_view.describe(role, environ=environ)["options"]:
        if not option["allowed"]:
            continue
        entries.append(Entry(
            key=f"setting:{option['path']}", kind="setting", label=option["title"], summary=option["summary"],
            destination="settings", endpoint="admin_settings", keywords=(option["path"].replace(".", " "),),
            query={"q": option["title"]}, anchor=f"opt-{option['path'].replace('.', '-')}",
            danger="rewrites" if option["restart"] else "safe",
        ))
    return entries


def entries(role="admin", environ=None):
    """Every entry the role can open: pages, tasks, then settings."""
    return [*PAGES, *tasks(), *setting_entries(role, environ=environ)]


def url(entry, url_for):
    """The entry's address; `url_for` is the caller's (`helpers.page_url` takes the endpoint without `web.`)."""
    address = url_for(entry.endpoint, **entry.params)
    if entry.query:
        address += "?" + urlencode(entry.query)
    return address + (f"#{entry.anchor}" if entry.anchor else "")


def search(query, found, limit=12):
    """
    The entries matching a typed query, best first.

    Args:
        query (str): What was typed; every word must appear somewhere in an entry's text.
        found (list): The `Entry` list to search (`entries()`).
        limit (int): At most this many.

    A word at the start of the label ranks above one inside it, which ranks above one in the summary or keywords.
    An empty query returns the pages.
    """
    words = query.lower().split()
    if not words:
        return [e for e in found if e.kind == "page"][:limit]
    scored = []
    for position, entry in enumerate(found):
        text = entry.text()
        if not all(word in text for word in words):
            continue
        label = entry.label.lower()
        score = sum(0 if label.startswith(w) else 1 if w in label else 2 for w in words)
        scored.append((score, ("page", "task", "setting").index(entry.kind), position, entry))
    scored.sort(key=lambda item: item[:3])
    return [item[3] for item in scored[:limit]]


def shell(section, page_url, admin=None):
    """
    What the shell draws around an admin page, or `None` for a page that is not under it or an admin who must first
    change the seeded credentials (the rail's links would only bounce them back).

    Args:
        section (str): The page's `render_page(section=...)`.
        page_url (callable): `helpers.page_url`.
        admin (dict, optional): The signed-in admin; without one there is no shell.

    Returns:
        dict: `destinations` (`[{"key", "label", "url", "current"}]`, the rail), `search_url`, `palette_url` and `command_url` (typed commands).
    """
    current = SECTION_DESTINATION.get(section)
    if current is None or not admin or admin.get("must_change_credentials"):
        return None
    return {
        "destinations": [
            {"key": key, "label": label, "url": page_url(endpoint), "current": key == current}
            for key, label, endpoint in DESTINATIONS
        ],
        "search_url": page_url("admin_search"),
        "palette_url": page_url("admin_palette"),
        "command_url": page_url("generate_command"),
    }
