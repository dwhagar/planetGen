# tests/test_fuzz_text_rendering.py

"""
Property-based / brute-force tests for the text/HTML rendering layer:
`planetgen/describe/html.py`, `planetgen/web/lib/fmt.py`, `planetgen/web/lib/pagination.py`,
`planetgen/web/lib/tabledisplay.py`, and
`planetgen/db/render.py`.

Main invariants, checked structurally with `html.parser` rather than by
string search alone:
  - hostile text never produces markup: every tag in the output is one the
    renderer itself emits, with only its own attributes (no `<script`, no
    `href`/`on*` attribute, no `javascript:` URL);
  - pagination math is consistent for any total/page/page size;
  - number formatting never crashes on NaN/inf/huge/negative input.

Same `sys.path` setup as `test_pagination.py`/`test_page_shell.py`. See
`fuzz_support.py` for the `ci`/`deep` profiles.
"""

import html
import math
import re
from html.parser import HTMLParser
from unittest import mock

import pytest

from planetgen.describe.html import blocks_html
from planetgen.describe.model import split_by_body
from planetgen.util.format import format_number
from hypothesis import assume, example, given, settings
from hypothesis import strategies as st

from planetgen.web.lib import fmt  # noqa: E402
from planetgen.describe import html as describe_html  # noqa: E402
from planetgen.web.lib import pagination  # noqa: E402
from planetgen.web.lib import tabledisplay  # noqa: E402
from planetgen.db import render as systemRender  # noqa: E402
from planetgen.generation.config import SystemConfig  # noqa: E402
from planetgen.generation.system import StarSystem  # noqa: E402
from tests.fuzz_support import any_float, finite, hostile_text, non_finite, scaled  # noqa: E402

_MARKDOWN_FRAGMENTS = [
    "# ", "## ", "###### ", "####### ", "| a | b |", "|---|---|", "|:--|--:|", "| x |", "|", "||",
    "<sup>2</sup>", "<sup>", "</sup>", "<sup><script>x</script></sup>", "&lt;sup&gt;", "&amp;",
    "<img src=x onerror=alert(1)>", '"><svg onload=alert(1)>', "[x](javascript:alert(1))",
    "<a href='javascript:alert(1)'>x</a>", "\n", "\n\n", "\r\n\r\n", " ", "\t",
]
markdown_text = st.one_of(
    hostile_text,
    st.lists(st.one_of(st.sampled_from(_MARKDOWN_FRAGMENTS), hostile_text), max_size=12).map("".join),
)


class _TagCollector(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tags = []      # (tag, attrs) for every start/startend tag
        self.text = []

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))

    def handle_startendtag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))

    def handle_data(self, data):
        self.text.append(data)


def _parse(markup):
    collector = _TagCollector()
    collector.feed(markup)
    collector.close()
    return collector


# Attributes whose value is a caller-supplied page target (`action`, the
# GET-mode `href` built from it) -- always a constant script name/URL path
# in real callers, so only their escaping is checked, not their scheme.
_TRUSTED_TARGET_ATTRS = {"action", "href"}
# Attributes a browser follows as a URL -- the only place a `javascript:`/
# `data:` value is dangerous (as a hidden field's value it's inert text).
_URL_ATTRS = {"href", "src", "action", "formaction", "xlink:href"}


def _assert_no_active_content(markup, allowed_tags, allowed_attrs, trusted_targets=False):
    lowered = markup.lower()
    assert "<script" not in lowered
    parsed = _parse(markup)
    for tag, attrs in parsed.tags:
        assert tag in allowed_tags, (tag, markup)
        for name, value in attrs.items():
            assert name in allowed_attrs, (tag, name, markup)
            assert not name.startswith("on"), (tag, name)
            if (value is not None and name in _URL_ATTRS
                    and not (trusted_targets and name in _TRUSTED_TARGET_ATTRS)):
                assert not value.strip().lower().startswith(("javascript:", "data:", "vbscript:")), (name, value)
    return parsed


# ---------------------------------------------------------------------------
# The Document's HTML
# ---------------------------------------------------------------------------

_BLOCK_TAGS = {"p", "table", "caption", "thead", "tbody", "tr", "th", "td", "div", "span", "strong", "a", "sup"}
_BLOCK_ATTRS = {"class", "tabindex", "scope", "href", "data-q", "data-si", "data-nc", "data-colour"}


def _hostile_document(text):
    return {"title": text, "sections": [{"kind": "star", "title": text, "ref": None, "children": [], "blocks": [
        {"type": "facts", "title": text, "rows": [{"key": "a", "label": text, "value": [text]}]},
        {"type": "prose", "role": "summary", "text": [[text]]}]}]}


_DOC_TAGS = {"section", "h1", "h2", "h3", "h4", "h5", "h6", "p", "div", "table", "caption", "thead", "tbody", "tr", "th",
             "td", "span", "strong", "em", "a", "sup", "ul", "li", "br"}


class _DocAttrs(set):
    """The attributes the Document's HTML uses: a few fixed ones and any `data-*`."""

    def __contains__(self, name):
        return super().__contains__(name) or name.startswith("data-")


_DOC_ATTRS = _DocAttrs({"class", "tabindex", "scope", "href", "title"})


@given(text=markdown_text)
def test_the_document_html_contains_only_the_renderers_own_markup(text):
    _assert_no_active_content(describe_html.render(_hostile_document(text)), _DOC_TAGS, _DOC_ATTRS)


@given(text=markdown_text)
def test_the_document_html_never_loses_visible_text(text):
    visible = "".join(_parse(describe_html.render(_hostile_document(text))).text)
    squash = lambda value: re.sub(r"\s", "", value)  # noqa: E731
    assert squash(text) in squash(visible) or not text.strip()


# ---------------------------------------------------------------------------
# fmt
# ---------------------------------------------------------------------------

@given(value=st.one_of(hostile_text, st.none(), finite, any_float, st.integers()))
def test_esc_round_trips_and_leaves_no_markup_characters(value):
    out = fmt.esc(value)
    assert not set(out) & set("<>\"'")
    assert html.unescape(out) == ("" if value is None else str(value))


@given(prefix=hostile_text,
       entries=st.lists(st.tuples(hostile_text, st.floats(0, 1e6)), max_size=4),
       link_some=st.booleans())
def test_linkify_nearest_escapes_everything_and_links_only_known_names(prefix, entries, link_some):
    assume(fmt._LOCATION_NEIGHBOR_MARKER not in prefix)
    neighbors = ", ".join(f"{name} ({dist:.1f} ly)" for name, dist in entries)
    location = f"{prefix} -- nearest: {neighbors}" if entries else prefix
    name_to_id = {name: i + 1 for i, (name, _d) in enumerate(entries)} if link_some else {}
    out = fmt.linkify_nearest(location, name_to_id, lambda i: f"/system/{i}")
    _assert_no_active_content(out, {"a"}, {"href"})
    if not entries:
        assert out == ""
    assert len([t for t, _a in _assert_no_active_content(out, {"a"}, {"href"}).tags if t == "a"]) <= fmt.NEAREST_SHOWN


@given(location=st.one_of(st.none(), hostile_text),
       neighbors=st.lists(st.fixed_dictionaries({"id": st.integers(1, 10**9), "name": hostile_text,
                                                  "distance_ly": st.floats(0, 1e9)}), max_size=4))
def test_nearest_systems_html_escapes_every_name(location, neighbors):
    out = fmt.nearest_systems_html(neighbors, lambda i: f"/system/{i}")
    parsed = _assert_no_active_content(out, {"a"}, {"href"})
    assert len([t for t, _a in parsed.tags if t == "a"]) == min(len(neighbors), fmt.NEAREST_SHOWN)


@given(distance=st.one_of(st.none(), any_float))
def test_format_distance_ly_never_raises(distance):
    out = fmt.format_distance_ly(distance)
    units = (" km", " AU", " ly", " ly)", " AU)", " km)")
    # No value (None, NaN, an infinity, or one that overflows on the way to
    # km) is a dash, never "nan km"/"inf km" (TEST.53).
    assert out == "&ndash;" if distance is None else (out == "&ndash;" or out.endswith(units))
    assert not re.search(r"\b(?:nan|inf)\b", out, re.IGNORECASE)


@given(edge=st.floats(min_value=1e-3, max_value=1e6), count=st.integers(0, 10**7))
def test_format_density_for_realistic_sectors(edge, count):
    out = fmt.format_density(edge, count)
    assert out.startswith(f"{format_number(count / edge / edge / edge, '.5f')} systems/ly&sup3;")


@given(edge=st.one_of(st.just(0), st.just(0.0), st.just(None), st.just(-0.0), st.floats(max_value=-1e-3),
                      non_finite),
       count=st.integers(-10, 10**6))
def test_format_density_degenerate_edges_do_not_raise(edge, count):
    assert isinstance(fmt.format_density(edge, count), str)


@given(edge=any_float, count=st.integers(-10, 10**9))
@example(edge=1e-110, count=5)   # edge**3 underflowed to 0.0 -> ZeroDivisionError
@example(edge=1e300, count=5)    # edge**3 overflowed -> OverflowError
@example(edge=-5.643803094122362e+102, count=0)
def test_format_density_never_raises(edge, count):
    out = fmt.format_density(edge, count)
    assert out == "n/a" or out.endswith(("systems/ly&sup3;", "of local average)"))
    assert "inf" not in out


@given(name=st.text(alphabet=st.characters(blacklist_categories=("Cs",)), max_size=40))
def test_static_url_quotes_the_version(name):
    with mock.patch.object(fmt, "STATIC_VERSION", '1.0"<x> &y'):
        out = fmt.static_url(name)
    assert out.startswith(f"static/{name}?v=")
    assert not set(out.split("?v=", 1)[1]) & set("<>\"& ")


# ---------------------------------------------------------------------------
# pagination
# ---------------------------------------------------------------------------

totals = st.integers(0, 10**12)
page_sizes = st.integers(1, 10**4)
pages = st.integers(-10**12, 10**12)


@given(total=totals, size=page_sizes)
def test_page_count_is_ceiling_and_at_least_one(total, size):
    count = pagination.page_count(total, size)
    assert count == max(1, math.ceil(total / size)) if total < 2**50 else count >= 1
    assert (count - 1) * size < max(total, 1) <= count * size


@given(page_=pages, total=totals, size=page_sizes)
def test_clamp_page_and_offset_are_in_range(page_, total, size):
    clamped = pagination.clamp_page(page_, total, size)
    assert 1 <= clamped <= pagination.page_count(total, size)
    if 1 <= page_ <= pagination.page_count(total, size):
        assert clamped == page_
    offset = pagination.page_offset(clamped, size)
    assert offset >= 0 and offset % size == 0
    assert offset < max(total, 1)


@given(items=st.lists(st.integers(), max_size=300), page_=pages, size=st.integers(1, 60))
def test_page_slice_matches_list_slicing(items, page_, size):
    rows, clamped = pagination.page_slice(items, page_, size)
    start = (clamped - 1) * size
    assert rows == items[start:start + size]
    assert len(rows) <= size
    assert bool(rows) == bool(items)


@given(items=st.lists(st.integers(), max_size=300), size=st.integers(1, 60))
def test_every_page_together_reconstructs_the_list(items, size):
    out = []
    for p in range(1, pagination.page_count(len(items), size) + 1):
        rows, clamped = pagination.page_slice(items, p, size)
        assert clamped == p
        out.extend(rows)
    assert out == items


@given(total=st.integers(0, 10**6), page_=pages, size=page_sizes)
def test_fetch_page_never_asks_for_a_negative_offset_after_clamping(total, page_, size):
    calls = []

    def fetch(limit, offset):
        calls.append((limit, offset))
        return {"items": [], "total": total}

    envelope, clamped = pagination.fetch_page(fetch, page_, size)
    assert envelope["total"] == total
    assert 1 <= clamped <= pagination.page_count(total, size)
    assert calls[-1] == (size, (clamped - 1) * size)
    assert calls[-1][1] >= 0
    assert len(calls) == (1 if clamped == page_ else 2)


@given(raw=st.one_of(st.none(), hostile_text, st.integers(), st.binary(max_size=10)))
def test_parse_page_never_raises_and_is_positive(raw):
    value = pagination.parse_page(raw)
    assert isinstance(value, int) and value >= 1


@given(page_=st.integers(1, 10**6), last=st.integers(1, 10**6))
def test_page_numbers_window(page_, last):
    assume(page_ <= last)
    numbers = pagination._page_numbers(page_, last)
    shown = [n for n in numbers if n is not None]
    assert shown == sorted(set(shown))
    assert shown[0] == 1 and shown[-1] == last and page_ in shown
    assert all(1 <= n <= last for n in shown)
    for a, b in zip(numbers, numbers[1:]):
        assert not (a is None and b is None)
    # Gaps only where more than one page is hidden.
    for i, n in enumerate(numbers):
        if n is None:
            assert numbers[i + 1] - numbers[i - 1] > 2
    assert len(numbers) <= 2 * pagination._WINDOW + 5


@settings(max_examples=scaled(60))
@given(total=st.integers(0, 10**9), size=st.integers(1, 500), data=st.data(),
       params=st.dictionaries(hostile_text, hostile_text, max_size=4),
       anchor=st.one_of(st.none(), hostile_text), label=hostile_text,
       action=hostile_text)
def test_render_pagination_is_safe_and_consistent(total, size, data, params, anchor, label, action):
    last = pagination.page_count(total, size)
    page_ = data.draw(st.integers(1, last))
    out = pagination.render_pagination(action, params, "p", page_, total, page_size=size, anchor=anchor,
                                       label=label)
    if total <= size:
        assert out == ""
        return
    parsed = _assert_no_active_content(
        out, {"nav", "span", "a", "form", "input", "button"},
        {"class", "aria-label", "aria-current", "aria-hidden", "aria-disabled", "href", "method", "action",
         "type", "name", "value"}, trusted_targets=True)
    summary = re.search(r"Showing ([\d,]+)&ndash;([\d,]+) of ([\d,]+)", out)
    first, last_row, shown_total = (int(g.replace(",", "")) for g in summary.groups())
    assert shown_total == total and 1 <= first <= last_row <= total
    assert last_row - first + 1 <= size
    current = [t for t, a in parsed.tags if a.get("aria-current") == "page"]
    assert len(current) == 1
    for _t, attrs in parsed.tags:
        if "href" in attrs:
            assert attrs["href"].startswith(action)


# ---------------------------------------------------------------------------
# tabledisplay
# ---------------------------------------------------------------------------

@given(exponent=st.integers(-400, 400), coeff=st.floats(1, 9.99), extra=hostile_text)
def test_to_plain_text_replaces_every_sup_exponent(exponent, coeff, extra):
    assume("<sup>" not in extra)
    formatted = f"{coeff:.2f} × 10<sup>{exponent}</sup> {extra}"
    out = tabledisplay.to_plain_text(formatted)
    assert "<sup>" not in out
    assert out.endswith(extra)
    assert tabledisplay.to_plain_text(out) == out
    digits = str(exponent).translate(tabledisplay._SUPERSCRIPT_DIGITS)
    assert f"10{digits} " in out


@given(text=hostile_text)
def test_to_plain_text_is_a_no_op_without_sup(text):
    assume("<sup>" not in text)
    assert tabledisplay.to_plain_text(text) == text


nonzero_finite = st.floats(allow_nan=False, allow_infinity=False, min_value=-1e300, max_value=1e300).filter(
    lambda x: x == 0 or abs(x) > 1e-300)


@given(value=nonzero_finite)
def test_star_formatters_never_raise_on_finite_values(value):
    for formatter in (tabledisplay.format_mass, tabledisplay.format_luminosity,
                      tabledisplay.format_star_radius):
        out = formatter(value)
        assert isinstance(out, str) and out
        assert "<sup>" not in tabledisplay.to_plain_text(out)


@given(value=nonzero_finite, is_moon=st.booleans())
def test_format_body_distance_never_raises_on_finite_values(value, is_moon):
    out = tabledisplay.format_body_distance(value, is_moon)
    assert out.endswith((" km", " AU", " ly", " ly)", " AU)", " km)"))


@given(years=st.floats(min_value=0, max_value=1e9))
def test_format_period_never_raises(years):
    assert isinstance(tabledisplay.format_period(years), str)


@given(value=any_float)
@example(value=math.nan)
@example(value=math.inf)
@example(value=-math.inf)
@example(value=5e-324)   # 10**-324 == 0.0 made the old coefficient a division by zero
def test_formatters_survive_non_finite_and_subnormal_values(value):
    for formatter in (tabledisplay.format_mass, tabledisplay.format_luminosity,
                      tabledisplay.format_star_radius, lambda v: tabledisplay.format_body_distance(v, True)):
        assert isinstance(formatter(value), str)


# ---------------------------------------------------------------------------
# systemRender
# ---------------------------------------------------------------------------

@given(fmt_=hostile_text)
def test_render_star_system_rejects_unknown_formats_before_describing(fmt_):
    assume(fmt_ not in systemRender.FORMATS)
    with pytest.raises(ValueError):
        systemRender.render_star_system(object(), fmt_)


@given(fmt_=hostile_text)
def test_render_system_text_rejects_unknown_formats_before_touching_the_database(fmt_):
    assume(fmt_ not in systemRender.FORMATS)
    with pytest.raises(ValueError):
        systemRender.render_system_text(None, 1, fmt_)


_HOSTILE_NAMES = ["<script>alert(1)</script>", '"><img src=x onerror=alert(1)>', "# Heading\n\n| a |\n|--|",
                  "javascript:alert(1)", "<sup onmouseover=alert(1)>2</sup>", "𝔘𝔫𝔦𝔠𝔬𝔡𝔢 ‮", "\x00\t\r\n", ""]


@pytest.mark.parametrize("name", _HOSTILE_NAMES)
def test_generated_system_with_hostile_name_renders_to_inert_html(name):
    config = SystemConfig()
    config.NAME = name
    system = StarSystem(system_config=config)
    for fmt_ in systemRender.FORMATS:
        text = systemRender.render_star_system(system, fmt_)
        assert isinstance(text, str) and text
    from planetgen import describe
    _assert_no_active_content(describe_html.render(describe.describe(system)), _DOC_TAGS, _DOC_ATTRS)


def test_the_split_system_description_with_hostile_names_is_inert(mysql_config):
    from planetgen.db import store
    for name in _HOSTILE_NAMES[:4]:
        config = SystemConfig()
        config.NAME = name
        config.MOONS = True
        config.COMETS = True
        config.ASTEROID_BELT = True
        config.BINARY_SYSTEM = False
        system_id = store.save_system(StarSystem(system_config=config), config, config=mysql_config)
        conn = store.get_connection(mysql_config)
        try:
            sections = split_by_body(systemRender.describe_system(conn, system_id))
        finally:
            conn.close()
        chunks = [sections["overview"]] + [t for key in ("stars", "planets", "moons", "belts", "comets")
                                           for t in sections[key].values()]
        assert all(chunks)
        for chunk in chunks:
            parsed = _assert_no_active_content(blocks_html(chunk, lambda ref: f"#{ref}"), _BLOCK_TAGS, _BLOCK_ATTRS)
            for _tag, attrs in parsed.tags:        # a swatch's colour is a hex colour, nothing else (GEN.420)
                assert re.fullmatch(r"#[0-9a-f]{6}", attrs["data-colour"]) if "data-colour" in attrs else True
