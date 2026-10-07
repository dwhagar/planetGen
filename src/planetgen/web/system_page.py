# planetgen/web/system_page.py

"""
`/admin/generate/system`: a one-off star system from the browser, with
every option `generate.py system` takes, shown as Markdown (or wikitext)
on the page and never saved to the database.

The page runs the real command-line generator, `generate.py system
--output FILE`, in a throwaway directory and reads the page back, so the
browser and the terminal can never disagree on what an option means. A
separate process also keeps the flavor-text overrides (which set
`program_constants` globals) out of the web server. A system takes about
a second, so the page waits for it rather than starting a background job.

The result's code box sits inside a second form: its Download button
posts the text back to `/admin/generate/system/download`, which returns
it as a file, and Copy uses `static/copycode.js` like the system page.

Admins only, like the rest of `/admin/generate` (`generate_page.py`).
"""

import json
import os
import re
import subprocess
import tempfile

from flask import Response, request

from . import bp, jobs
from .generate_page import FormError, _admin_or_403, _admin_or_redirect, _no_store, _number
from .helpers import crumb, render_page, trusted_html

from planetgen.web.lib.mdconvert import markdown_to_html
from planetgen.generation.limits import MAX_NUM_ORBITS

TRISTATE_FIELDS = (
    # (generate.py option name, label) -- the same ten, in the same order,
    # as generate.py's TRISTATE_OPTIONS (a test checks they still match).
    ("habitable_world", "Habitable world"),
    ("asteroid_belt", "Asteroid belt"),
    ("comets", "Star-bound comet"),
    ("large_star", "Large, massive star"),
    ("moons", "Moons on every planet"),
    ("max_planets", "Maximum number of orbits"),
    ("intelligent_life", "Intelligent life"),
    ("binary_system", "Binary star"),
    ("wide_binary", "Wide (S-type) binary"),
    ("planets", "At least one planet or belt"),
)

TRISTATE_CHOICES = (("", "Generator decides"), ("yes", "Yes"), ("no", "No"))

FORMATS = (("markdown", "Markdown"), ("wikitext", "Wikitext"))

SYSTEM_FILE_EXAMPLE = """{
  "num_orbits": 5,
  "slots": [
    {"type": "planet", "planet_class": "M", "moons": 1},
    {"type": "asteroid_belt"},
    null,
    {"type": "planet", "planet_class": "J", "moons": 4},
    null
  ]
}"""

TIMEOUT_S = 120
"""int: How long the page waits for `generate.py` before giving up."""

MAX_LOG_BYTES = 256 * 1024
"""int: The most generator output (and debug log) the page shows."""


def _tristate(form, name):
    value = form.get(name) or ""
    if value not in ("", "yes", "no"):
        raise FormError("Choose Yes, No or Generator decides.")
    return {"": None, "yes": True, "no": False}[value]


def _text(form, name):
    return (form.get(name) or "").strip()


def system_request(form):
    """
    The `generate.py system` options a form asks for, checked the way
    `generate.py`'s `validate_system_args` checks them.

    Returns:
        dict: `argv` (the options after `system`, without `--output` or
        the logging flags), `format`, `debug` and `system_file` (the
        parsed JSON, or `None`).

    Raises:
        FormError: A bad value, or a combination the generator refuses.
    """
    states = {name: _tristate(form, name) for name, _label in TRISTATE_FIELDS}
    name = _text(form, "name")
    star_type = _text(form, "star_type")
    age = _text(form, "age")
    if age not in ("", "young", "old"):
        raise FormError("Age must be young, old or left to the generator.")
    num_orbits = _number(form, "num_orbits", "Orbital slots", int, minimum=0, maximum=MAX_NUM_ORBITS)
    chance_system = _number(form, "flavor_chance_system", "System flavor chance", float, minimum=0.0, maximum=1.0)
    chance_planet = _number(form, "flavor_chance_planet", "Planet flavor chance", float, minimum=0.0, maximum=1.0)
    max_flavor = bool(form.get("max_planet_flavor"))
    fmt = form.get("format") or "markdown"
    if fmt not in dict(FORMATS):
        raise FormError("Choose Markdown or Wikitext.")

    if states["planets"] is False and (states["moons"] or states["max_planets"] or states["habitable_world"]):
        raise FormError("No planets can't be combined with moons, maximum orbits or a habitable world.")
    if star_type and states["large_star"]:
        raise FormError("A star type can't be combined with forcing a large star.")
    if states["intelligent_life"] is not None and states["habitable_world"] is False:
        raise FormError("Intelligent life (yes or no) can't be combined with no habitable world.")
    if num_orbits is not None and states["planets"] is False:
        raise FormError("Orbital slots can't be combined with no planets.")
    if states["habitable_world"] and states["asteroid_belt"] and states["large_star"] is False:
        raise FormError("A habitable world and an asteroid belt together need a large star.")

    system_file = None
    raw_file = _text(form, "system_file")
    if raw_file:
        try:
            system_file = json.loads(raw_file)
        except ValueError as exc:
            raise FormError(f"The system file isn't valid JSON: {exc}.") from None
        if not isinstance(system_file, dict):
            raise FormError("The system file must be a JSON object ({ ... }).")

    argv = []
    for option, _label in TRISTATE_FIELDS:
        if states[option] is not None:
            argv.append(("+" if states[option] else "-") + option)
    if fmt == "markdown":
        argv.append("--markdown")
    # `--flag=value`, so a name starting with "-" is never read as an option.
    if name:
        argv.append(f"--name={name}")
    if star_type:
        argv.append(f"--star-type={star_type}")
    if age:
        argv.append(f"--age={age}")
    if num_orbits is not None:
        argv.append(f"--num-orbits={num_orbits}")
    if chance_system is not None:
        argv.append(f"--flavor-chance-system={chance_system}")
    if chance_planet is not None:
        argv.append(f"--flavor-chance-planet={chance_planet}")
    if max_flavor:
        argv.append("--max-planet-flavor")
    return {"argv": argv, "format": fmt, "debug": bool(form.get("debug")), "system_file": system_file}


def _read(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            return f.read()
    except FileNotFoundError:
        return ""


def _tail(text, limit=MAX_LOG_BYTES):
    return text[-limit:] if len(text) > limit else text


def run_generator(spec):
    """
    Runs `generate.py system --output` for a `system_request` and reads
    the page back. Nothing touches the database.

    Returns:
        dict: `ok`, `text` (the page), `output` (anything the generator
        printed) and `debug_log` (the same output, when `--debug` was
        asked for).
    """
    with tempfile.TemporaryDirectory(prefix="planetgen-system-") as tmp:
        page_path = os.path.join(tmp, "system.txt")
        # With --debug the generator narrates every choice to stdout, which
        # the page shows as the debug log; otherwise only errors.
        argv = [jobs.python_executable(), jobs.GENERATE_SCRIPT, "system", "--output", page_path,
                "--debug" if spec["debug"] else "--quiet"]
        if spec["system_file"] is not None:
            file_path = os.path.join(tmp, "system.json")
            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(spec["system_file"], f)
            argv.append("--system-file=" + file_path)
        argv += spec["argv"]
        try:
            proc = subprocess.run(argv, cwd=jobs.REPO_DIR, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, timeout=TIMEOUT_S)
            output = proc.stdout.decode("utf-8", errors="replace")
            ok = proc.returncode == 0
        except subprocess.TimeoutExpired as exc:
            output = (exc.stdout or b"").decode("utf-8", errors="replace")
            output += f"\nThe generator took longer than {TIMEOUT_S} seconds and was stopped."
            ok = False
        except OSError as exc:
            output, ok = f"The generator could not be started: {exc}", False
        text = _read(page_path)
    output = _tail(output.strip())
    return {"ok": ok and bool(text), "text": text, "output": output, "debug_log": output if spec["debug"] else ""}


def _title_of(text, fmt):
    """The system's name, from the page's first heading."""
    first = text.lstrip().split("\n", 1)[0] if text else ""
    match = re.match(r"^#\s+(.+)$", first) if fmt == "markdown" else re.match(r"^=\s*(.+?)\s*=$", first)
    return match.group(1).strip() if match else "System"


def download_name(title, fmt):
    """`"Ilkone.md"` / `"Ilkone.wiki"`, keeping only safe characters."""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "-", title).strip("-.") or "system"
    return f"{stem[:80]}.{'md' if fmt == 'markdown' else 'wiki'}"


def _page(admin, form=None, result=None, error=None, status=200):
    return _no_store(render_page(
        "generate_system.html",
        title="One-off system",
        section="admin_generate",
        breadcrumbs=[crumb("Admin", "admin"), crumb("Generate", "generate"), crumb("One-off system")],
        description="Generate one star system without saving it to the database.",
        admin=admin,
        form=form or {},
        tristate_fields=TRISTATE_FIELDS,
        tristate_choices=TRISTATE_CHOICES,
        formats=FORMATS,
        system_file_example=SYSTEM_FILE_EXAMPLE,
        max_num_orbits=MAX_NUM_ORBITS,
        result=result,
        error=error,
        status=status,
    ))


@bp.route("/admin/generate/system", methods=["GET", "POST"])
def generate_system():
    """The one-off system form (GET), and a generated system under it
    (POST, the same form filled in again so it can be re-rolled)."""
    if request.method in ("GET", "HEAD"):  # HEAD is a GET without the body, never the POST branch
        admin, response = _admin_or_redirect()
        return response or _page(admin)

    admin = _admin_or_403()
    try:
        spec = system_request(request.form)
    except FormError as exc:
        return _page(admin, form=request.form, error=str(exc), status=400)
    run = run_generator(spec)
    if not run["ok"]:
        return _page(admin, form=request.form, error="The generator failed; its output is below.",
                     result={"failed": True, **run}, status=500)
    title = _title_of(run["text"], spec["format"])
    result = {
        **run,
        "failed": False,
        "format": spec["format"],
        "format_label": dict(FORMATS)[spec["format"]],
        "title": title,
        "filename": download_name(title, spec["format"]),
        "rows": min(run["text"].count("\n") + 3, 30),
        "preview": trusted_html(markdown_to_html(run["text"])) if spec["format"] == "markdown" else None,
    }
    return _page(admin, form=request.form, result=result)


@bp.route("/admin/generate/system/download", methods=["POST"])
def generate_system_download():
    """The code box's text back as a file (the page's Download button)."""
    _admin_or_403()
    fmt = request.form.get("format") if request.form.get("format") in dict(FORMATS) else "markdown"
    text = (request.form.get("text") or "").replace("\r\n", "\n")
    filename = download_name(request.form.get("title") or "system", fmt)
    mimetype = "text/markdown" if fmt == "markdown" else "text/plain"
    response = Response(text, mimetype=mimetype)
    response.headers["Content-Disposition"] = f'attachment; filename="{filename}"'
    return _no_store(response)
