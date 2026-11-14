# planetgen/web/system_page.py

"""
`/admin/generate/system`: a one-off star system from the browser, with
every option `planetgen system` takes, shown as the page a saved system gets
and never saved to the database.

The page runs the real command-line generator, `planetgen system
--output FILE --format=document`, in a throwaway directory and reads the
system's Document back (`planetgen.describe`), so the
browser and the terminal can never disagree on what an option means. A
separate process also keeps the flavor-text overrides (which set
`program_constants` globals) out of the web server. A system takes about
a second, so the page waits for it rather than starting a background job.

The Document is shown as HTML (`describe.html`), and a Download button for
each text format posts that format's text to `/admin/generate/system/download`,
which returns it as a file.

Admins only, like the rest of `/admin/generate` (`generate_page.py`).
"""

import json
import os
import re
import tempfile

from flask import Response, request

from . import bp, jobs
from .generate_page import FormError, _admin_or_403, _admin_or_redirect, _no_store, _number
from .helpers import crumb, trusted_html

from planetgen.queue import api_jobs
from planetgen import describe
from planetgen.describe import DOCUMENT_FORMAT, html as describe_html
from planetgen.generation.limits import MAX_NUM_ORBITS
from planetgen.generation.run_system import ODDS_OPTIONS, forced_oddity_ids, max_stars_problem, stars_problem
from planetgen.physics import oddities
from planetgen.generation.star import STAR_TYPE_PATTERN
from planetgen.db.store import SYSTEM_NAME_MAX_LENGTH
from planetgen.util import log as planetgen_log
from planetgen.tuning import MAX_STARS_CEILING, MAX_STARS_DEFAULT

TRISTATE_FIELDS = (
    # (planetgen option name, label) -- the same ten, in the same order,
    # as planetgen's TRISTATE_OPTIONS (a test checks they still match).
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

ODDITY_IDS = tuple(oddities.BY_ID)
"""tuple: The oddity ids the Force oddities box accepts (UX.299)."""

STARLESS_RESULTS = (("", "A star system"), ("rogue-planet", "A rogue planet (no star)"), ("comet", "An interstellar comet (no star)"))
"""tuple: (UX.298) What the page can generate: a star system, or one of the starless results `planetgen phenomenon --type` makes,
described but not saved (`--output`). Only the name and the debug log apply to a starless result."""

TRISTATE_CHOICES = (("", "Generator decides"), ("yes", "Yes"), ("no", "No"))

ODDS_CHOICES = (("", "Generator decides"), ("yes", "Raise the odds"), ("no", "No"))
"""tuple: The choices for the options that only raise the odds (`ODDS_OPTIONS`, GEN.1130): Yes would promise what the generator does not force."""

FORMATS = (("markdown", "Markdown", "md", "text/markdown"),
           ("text", "Plain text", "txt", "text/plain"))
"""tuple: The text formats offered for download: `(format, label, file extension, mimetype)`."""

_FORMAT_EXTENSIONS = {fmt: extension for fmt, _label, extension, _mimetype in FORMATS}
_FORMAT_MIMETYPES = {fmt: mimetype for fmt, _label, _extension, mimetype in FORMATS}

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
"""int: How long the page waits for `planetgen` before giving up."""

MAX_LOG_BYTES = 256 * 1024
"""int: The most generator output (and debug log) the page shows."""


def _tristate(form, name):
    value = form.get(name) or ""
    if value not in ("", "yes", "no"):
        raise FormError("Choose Yes, No or Generator decides.")
    return {"": None, "yes": True, "no": False}[value]


def _text(form, name):
    return (form.get(name) or "").strip()


def _starless_request(form, result_kind):
    """The `planetgen phenomenon` options for a rogue planet or an interstellar comet (UX.298): its name and debug log only,
    so any other option that is filled in is refused instead of ignored."""
    others = [name for name, _label in TRISTATE_FIELDS if _tristate(form, name) is not None]
    others += [key for key in ("star_type", "age", "num_orbits", "stars", "max_stars", "flavor_chance_system", "flavor_chance_planet",
                               "system_file", "weirdness_dial", "force_oddities") if _text(form, key)]
    if form.get("max_planet_flavor"):
        others.append("max_planet_flavor")
    if others:
        raise FormError("A rogue planet or comet has no star, so only its name applies; clear: " + ", ".join(others).replace("_", " ") + ".")
    name = _text(form, "name")
    if len(name) > SYSTEM_NAME_MAX_LENGTH:
        raise FormError(f"The name is {len(name)} characters; a name can be at most {SYSTEM_NAME_MAX_LENGTH}.")
    argv = ["--type", result_kind, f"--format={DOCUMENT_FORMAT}"]
    if name:
        argv.append(f"--name={name}")
    if form.get("disregard_missing_data"):
        argv.append("--disregard-missing-data")
    return {"argv": argv, "debug": bool(form.get("debug")), "system_file": None, "command": "phenomenon"}


def system_request(form):
    """
    The `planetgen system` options a form asks for, checked the way
    `planetgen`'s `validate_system_args` checks them.

    Returns:
        dict: `argv` (the options after `system`, without `--output` or
        the logging flags), `debug` and `system_file` (the
        parsed JSON, or `None`).

    Raises:
        FormError: A bad value, or a combination the generator refuses.
    """
    result_kind = _text(form, "result_kind")
    if result_kind not in {value for value, _label in STARLESS_RESULTS}:
        raise FormError("The result must be a star system, a rogue planet or an interstellar comet.")
    if result_kind:
        return _starless_request(form, result_kind)
    states = {name: _tristate(form, name) for name, _label in TRISTATE_FIELDS}
    name = _text(form, "name")
    star_type = _text(form, "star_type")
    age = _text(form, "age")
    if age not in ("", "young", "old"):
        raise FormError("Age must be young, old, or left to the generator.")
    num_orbits = _number(form, "num_orbits", "Orbital slots", int, minimum=0, maximum=MAX_NUM_ORBITS)
    chance_system = _number(form, "flavor_chance_system", "System flavor chance", float, minimum=0.0, maximum=1.0)
    chance_planet = _number(form, "flavor_chance_planet", "Planet flavor chance", float, minimum=0.0, maximum=1.0)
    max_flavor = bool(form.get("max_planet_flavor"))
    stars = _number(form, "stars", "Number of stars", int, minimum=1, maximum=MAX_STARS_CEILING)
    max_stars = _number(form, "max_stars", "Most stars", int, minimum=1, maximum=MAX_STARS_CEILING)
    weirdness = _number(form, "weirdness_dial", "Weirdness dial", float, minimum=0.0)
    oddity_ids = _text(form, "force_oddities")
    if oddity_ids:
        try:
            oddity_ids = ",".join(forced_oddity_ids(oddity_ids))
        except ValueError as exc:
            raise FormError(f"Force oddities: {exc}.") from None

    if len(name) > SYSTEM_NAME_MAX_LENGTH:
        raise FormError(f"The name is {len(name)} characters; a system name can be at most {SYSTEM_NAME_MAX_LENGTH}.")
    if star_type and not STAR_TYPE_PATTERN.fullmatch(star_type.upper()):
        raise FormError("The star type must be a spectral type: a class letter (OBAFGKM), a subclass digit (0-9) "
                        "and a Yerkes class (0, IA+, IA, IAB, IB, II, III, IV, V, VI, VII or D), e.g. G2V.")
    if states["planets"] is False and (states["moons"] or states["max_planets"] or states["habitable_world"]
                                       or states["asteroid_belt"]):
        raise FormError("No planets can't be combined with moons, maximum orbits, an asteroid belt, or a habitable world.")
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

    # A form value wins over the system file's, as on the command line.
    resolved = lambda key, value: (system_file or {}).get(key) if value is None else value  # noqa: E731
    problem = max_stars_problem(resolved("max_stars", max_stars)) or stars_problem(
        resolved("stars", stars), resolved("max_stars", max_stars),
        resolved("binary_system", states["binary_system"]), resolved("wide_binary", states["wide_binary"]))
    if problem:
        raise FormError(problem[0].upper() + problem[1:])

    argv = []
    for option, _label in TRISTATE_FIELDS:
        if states[option] is not None:
            argv.append(("+" if states[option] else "-") + option)
    argv.append(f"--format={DOCUMENT_FORMAT}")
    # `--flag=value`, so a name starting with "-" is never read as an option.
    if name:
        argv.append(f"--name={name}")
    if star_type:
        argv.append(f"--star-type={star_type}")
    if age:
        argv.append(f"--age={age}")
    if num_orbits is not None:
        argv.append(f"--num-orbits={num_orbits}")
    if stars is not None:
        argv.append(f"--stars={stars}")
    if max_stars is not None:
        argv.append(f"--max-stars={max_stars}")
    if weirdness is not None:
        argv.append(f"--weirdness-dial={weirdness}")
    if oddity_ids:
        argv.append(f"--force-oddities={oddity_ids}")
    if chance_system is not None:
        argv.append(f"--flavor-chance-system={chance_system}")
    if chance_planet is not None:
        argv.append(f"--flavor-chance-planet={chance_planet}")
    if max_flavor:
        argv.append("--max-planet-flavor")
    if form.get("disregard_missing_data"):
        argv.append("--disregard-missing-data")
    debug = bool(form.get("debug"))
    areas = _text(form, "debug_areas")
    if areas:
        try:
            planetgen_log.parse_areas(areas)
        except ValueError as exc:
            raise FormError(f"Debug areas: {exc}.") from None
        if debug:
            argv.append(f"--debug-areas={areas}")
    return {"argv": argv, "debug": debug, "system_file": system_file}


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
    Runs `planetgen system --output` for a `system_request` and reads
    the Document back. Nothing touches the database.

    Returns:
        dict: `ok`, `text` (the Document as JSON), `output` (anything the generator
        printed) and `debug_log` (the same output, when `--debug` was
        asked for).
    """
    with tempfile.TemporaryDirectory(prefix="planetgen-system-") as tmp:
        page_path = os.path.join(tmp, "system.json")
        # With --debug the generator narrates every choice to stdout, which
        # the page shows as the debug log; otherwise only errors.
        argv = [jobs.python_executable(), *jobs.GENERATE_COMMAND, spec.get("command", "system"), "--output", page_path,
                "--debug" if spec["debug"] else "--quiet"]
        if spec["system_file"] is not None:
            file_path = os.path.join(tmp, "system.json")
            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(spec["system_file"], f)
            argv.append("--system-file=" + file_path)
        argv += spec["argv"]
        done = api_jobs.command_and_wait(argv, dict(os.environ), jobs.REPO_DIR, TIMEOUT_S, merge_stderr=True)
        output = done["stdout"]
        ok = done["returncode"] == 0
        if done["timed_out"]:
            output += f"\nThe generator took longer than {TIMEOUT_S} seconds and was stopped."
        elif done["error"]:
            output, ok = f"The generator could not be started: {done['error']}", False
        text = _read(page_path)
    output = _tail(output.strip())
    return {"ok": ok and bool(text), "text": text, "output": output, "debug_log": output if spec["debug"] else ""}


def download_name(title, fmt):
    """`"Ilkone.md"` / `"Ilkone.txt"`, keeping only safe characters."""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "-", title).strip("-.") or "system"
    return f"{stem[:80]}.{_FORMAT_EXTENSIONS[fmt]}"


def _downloads(document, title):
    """One entry per text format: `format`, `label`, `filename` and the Document written in it as `text`."""
    return [{"format": fmt, "label": label, "filename": download_name(title, fmt), "text": describe.render(document, fmt)}
            for fmt, label, _extension, _mimetype in FORMATS]


def _page(admin, form=None, result=None, error=None, status=200):
    from .generate_flow import render_screen

    return render_screen(
        "generate_system.html", admin, "One-off system", [crumb("One-off system")], form=form, error=error,
        status=status, description="Generate one star system without saving it to the database.",
        result_choices=[{"value": value, "label": label} for value, label in STARLESS_RESULTS], tristate_fields=TRISTATE_FIELDS, tristate_choices=TRISTATE_CHOICES, odds_choices=ODDS_CHOICES,
        odds_options=ODDS_OPTIONS, max_stars_ceiling=MAX_STARS_CEILING, max_stars_default=MAX_STARS_DEFAULT, system_file_example=SYSTEM_FILE_EXAMPLE,
        max_num_orbits=MAX_NUM_ORBITS, oddity_ids=", ".join(ODDITY_IDS), debug_areas=planetgen_log.AREAS, result=result)


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
    document = json.loads(run["text"])
    title = document["title"] or "System"
    result = {
        **run,
        "failed": False,
        "title": title,
        "page": trusted_html(describe_html.render(document)),
        "downloads": _downloads(document, title),
    }
    return _page(admin, form=request.form, result=result)


@bp.route("/admin/generate/system/download", methods=["POST"])
def generate_system_download():
    """The code box's text back as a file (the page's Download button)."""
    _admin_or_403()
    fmt = request.form.get("format") if request.form.get("format") in _FORMAT_EXTENSIONS else "markdown"
    text = (request.form.get("text") or "").replace("\r\n", "\n")
    filename = download_name(request.form.get("title") or "system", fmt)
    response = Response(text, mimetype=_FORMAT_MIMETYPES[fmt])
    response.headers["Content-Disposition"] = f'attachment; filename="{filename}"'
    return _no_store(response)
