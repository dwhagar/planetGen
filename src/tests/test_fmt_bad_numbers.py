# tests/test_fmt_bad_numbers.py

"""
TEST.53: the page formatters (`planetgen/web/lib/fmt.py`, `planetgen/web/lib/tabledisplay.py`)
fed bad numbers -- NaN, the infinities, negatives, zero, `None`, 1e300 and
an int past float range -- never raise, and never print "nan"/"inf": a
value that isn't a number shows a dash (`&ndash;` from the distance
formatters, which feed HTML; "–" from the rest, like their `None`).

The number formatters are found with `inspect` (every `format_*` callable
either module exposes), so a new one is swept without touching this file;
every other public function must be listed in `NOT_NUMBER_FORMATTERS`
with its own check below, so a new one can't slip past unnoticed.

Plus: empty tables (no neighbors, no rows, an empty location, an empty
pager) render as empty markup rather than raising.
"""

import inspect
import math
import re

import pytest

from planetgen.web.lib import fmt  # noqa: E402
from planetgen.web.lib import pagination  # noqa: E402
from planetgen.web.lib import tabledisplay  # noqa: E402

BAD_NUMBERS = [
    float("nan"), float("inf"), float("-inf"), -1, -1.5, -5e10, -1e300, 0, 0.0, -0.0, None, 1e300, 10 ** 400,
]
NON_NUMBERS = [None, float("nan"), float("inf"), float("-inf"), 10 ** 400]
"""Values that are no finite number: each must format as a dash."""

DASHES = ("&ndash;", "–")
_GARBAGE = re.compile(r"\b(?:nan|inf|infinity|none)\b", re.IGNORECASE)

NOT_NUMBER_FORMATTERS = {
    (fmt, "dash_unless_finite"),  # the guard itself
    (fmt, "esc"),  # escapes any value's str(): "nan" is the value's own text
    (fmt, "static_url"),
    (fmt, "quote"),  # urllib's, imported
    (fmt, "linkify_location"),
    (fmt, "nearest_neighbors_location"),
    (fmt, "nearest_systems_html"),
    (fmt, "utc_time_html"),
    (fmt, "inside_text"),
    (fmt, "runaway_text"),
    (tabledisplay, "to_plain_text"),
    (tabledisplay, "dash_unless_finite"),  # fmt's, imported
}
"""set: Public functions that aren't one-number formatters; each has its
own test below."""

_UTILS_HELPERS = {"format_body_radius_km", "format_relative_to_sol", "format_distance_km", "format_period_years",
                  "to_scientific_notation"}
"""set: `planetgen.util.format` functions `tabledisplay` imports for its own
use (they take a `SystemConfig`, or are wrapped by `tabledisplay`'s own
formatters), not part of its interface."""


def _public_functions(module):
    return {name: fn for name, fn in inspect.getmembers(module, callable)
            if not name.startswith("_") and inspect.isfunction(fn)}


def _number_formatters():
    found = []
    for module in (fmt, tabledisplay):
        for name, fn in sorted(_public_functions(module).items()):
            if not name.startswith("format_"):
                continue
            if module is tabledisplay and name in _UTILS_HELPERS:
                continue
            found.append((module, name, fn))
    return found


FORMATTERS = _number_formatters()
_POSITIONAL = (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD)


def _required_positional(fn):
    return sum(p.default is p.empty and p.kind in _POSITIONAL for p in inspect.signature(fn).parameters.values())


ONE_ARG = [(m, n, f) for m, n, f in FORMATTERS if _required_positional(f) == 1]
_IDS = [f"{m.__name__}.{n}" for m, n, _f in ONE_ARG]


def _id(value):
    return "10**400" if isinstance(value, int) and value > 10 ** 300 else repr(value)


def _is_number(value):
    """Whether `value` is a finite number a float can hold."""
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError, OverflowError):
        return False


def _clean(text):
    assert isinstance(text, str)
    assert not _GARBAGE.search(text), text
    assert len(text) < 200, text  # no 300-digit number spelled out


def test_the_sweep_finds_the_formatters():
    names = {(m.__name__.rpartition(".")[2], n) for m, n, _f in FORMATTERS}
    assert {("fmt", "format_number"), ("fmt", "format_distance_km"), ("fmt", "format_density"),
            ("fmt", "format_temperature_k"), ("tabledisplay", "format_star_mass"),
            ("tabledisplay", "format_body_distance")} <= names
    assert {n for m, n, _f in FORMATTERS} - {n for m, n, _f in ONE_ARG} == {"format_density"}


@pytest.mark.parametrize("module", [fmt, tabledisplay], ids=lambda m: m.__name__)
def test_every_public_function_is_covered(module):
    """A new public function is either a `format_*` (swept below) or
    listed in `NOT_NUMBER_FORMATTERS` with a test of its own."""
    swept = {n for m, n, _f in FORMATTERS if m is module}
    for name in _public_functions(module):
        covered = name in swept or (module, name) in NOT_NUMBER_FORMATTERS
        covered = covered or (module is tabledisplay and name in _UTILS_HELPERS)
        assert covered, f"{module.__name__}.{name} is not swept by this file"


@pytest.mark.parametrize("module,name,formatter", ONE_ARG, ids=_IDS)
@pytest.mark.parametrize("value", BAD_NUMBERS, ids=_id)
def test_formatter_never_raises_or_prints_garbage(module, name, formatter, value):
    _clean(formatter(value))


@pytest.mark.parametrize("module,name,formatter", ONE_ARG, ids=_IDS)
@pytest.mark.parametrize("value", NON_NUMBERS, ids=_id)
def test_no_number_is_a_dash(module, name, formatter, value):
    assert formatter(value) in DASHES


@pytest.mark.parametrize("module,name,formatter", ONE_ARG, ids=_IDS)
def test_zero_and_ordinary_values_still_format(module, name, formatter):
    for value in (0, 0.0, 1, 1234.5, 1e30):
        out = formatter(value)
        assert out not in DASHES and re.search(r"\d", out), (value, out)


@pytest.mark.parametrize("value", [-1.0, -5e10, -1e300])
@pytest.mark.parametrize("formatter", [tabledisplay.format_star_mass, tabledisplay.format_star_luminosity])
def test_a_negative_mass_or_luminosity_is_a_dash(formatter, value):
    assert formatter(value) == "–"


def test_distance_overflowing_to_infinite_km_is_a_dash():
    # 1e300 AU/ly/pc is finite, but infinite in km.
    for formatter in (fmt.format_distance_au, fmt.format_distance_ly, fmt.format_distance_pc):
        assert formatter(1e300) == "&ndash;"
    assert fmt.format_distance_km(1e300).endswith("ly)")


def test_format_number_keeps_ints_and_specs():
    assert fmt.format_number(1234) == "1,234"
    assert fmt.format_number(12_345) == "12,345"
    assert fmt.format_number(1_234_567) == "1.23 × 10⁶"
    assert fmt.format_number(1.5, ".1f") == "1.5"
    assert fmt.format_number(float("nan"), ".1f") == "–"


@pytest.mark.parametrize("edge", BAD_NUMBERS, ids=_id)
@pytest.mark.parametrize("count", BAD_NUMBERS, ids=_id)
def test_format_density_never_raises_or_prints_garbage(edge, count):
    out = fmt.format_density(edge, count)
    _clean(out)
    if not (_is_number(edge) and edge > 0 and _is_number(count)):
        assert out == "n/a"


def test_format_density_still_reads_for_a_real_sector():
    assert fmt.format_density(10.0, 3).startswith("0.00300 systems/ly&sup3;")


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf"), 1e300, -1e300, 10 ** 400, None, ""])
def test_utc_time_html_of_no_time_is_empty(value):
    assert fmt.utc_time_html(value) == ""


@pytest.mark.parametrize("value", [0, 0.0, -1, 1_790_000_000])
def test_utc_time_html_of_a_time(value):
    out = fmt.utc_time_html(value)
    assert out.startswith('<time datetime="') and out.endswith(" UTC</time>")


@pytest.mark.parametrize("value", BAD_NUMBERS, ids=_id)
def test_esc_and_static_url_take_any_value(value):
    assert isinstance(fmt.esc(value), str)
    assert fmt.esc(None) == ""
    assert fmt.static_url("style.css").startswith("static/style.css?v=")


@pytest.mark.parametrize("distance", [float("nan"), float("inf"), float("-inf"), None, 10 ** 400])
def test_nearest_systems_html_drops_a_bad_distance(distance):
    out = fmt.nearest_systems_html([{"id": 1, "name": "A", "distance_ly": distance}], lambda i: f"/system/{i}")
    assert out == '<a href="/system/1">A</a>'
    good = fmt.nearest_systems_html([{"id": 1, "name": "A", "distance_ly": 4.25}], lambda i: f"/system/{i}")
    assert good == '<a href="/system/1">A</a> (4.2 ly)'


def test_inside_and_runaway_text_with_bad_speeds():
    assert fmt.inside_text({}) is None and fmt.inside_text({"inside": None}) is None
    assert fmt.runaway_text({}) is None
    assert fmt.runaway_text({"runaway_class": "runaway", "runaway_speed_kms": None}) == "Runaway star"
    for speed in (float("nan"), float("inf"), 10 ** 400):
        assert fmt.runaway_text({"runaway_class": "runaway", "runaway_speed_kms": speed}) == "Runaway star, –"
    _clean(fmt.runaway_text({"runaway_class": "hypervelocity", "runaway_speed_kms": 1e300}))


def test_to_plain_text_of_dashes_and_empty():
    assert tabledisplay.to_plain_text("") == ""
    assert tabledisplay.to_plain_text("–") == "–"
    assert tabledisplay.to_plain_text(tabledisplay.format_star_radius(1e6)) == "1.00 × 10⁶ km"


# --- Empty tables -------------------------------------------------------------

def test_empty_tables_render():
    url = lambda i: f"/system/{i}"  # noqa: E731
    assert fmt.nearest_systems_html([], url) == ""
    assert fmt.nearest_neighbors_location(None, [], url) == " -- nearest: "
    assert fmt.linkify_location("", {}, url) == "" and fmt.linkify_location(None, {}, url) == ""
    assert pagination.render_pagination("/sectors", {}, "sectors_page", 1, 0) == ""
    assert pagination.page_slice([], 7) == ([], 1)
    assert pagination.page_count(0) == 1
