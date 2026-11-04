# tests/test_references_crawl.py

"""
DOC.87 (References R1): the site speaks in universe, and the out-of-universe words live only on the References pages.

Three things are held here.

* The addresses (`web/lib/reference_urls.py`) and the one link helper (`web.helpers.reference_link`): every kind of entry has one address,
  a typo in a kind or an id is an error, and no template or module writes a `/references` address by hand.
* The word list (`plainwords.OUT_OF_UNIVERSE`, rule 16 of the style guide): what it flags and what it leaves alone, and the pages it is
  exempt on (`plainwords.exempt`).
* The page crawl: every in-universe page this suite can render without a database (the Help pages, the class pages, the sign-in and
  preferences pages), every sentence the Document of ten seeded systems writes and every property label, read for the words. Nothing is
  excused: there is no list of exceptions, so text with a word fails the test the moment it is written (DOC.94 locked the crawl when the
  last rewrite, DOC.93, landed).
* The lock around it (DOC.94): the visible text of every in-universe template is read too, so a page that needs data to render is covered,
  and every in-universe page route must be crawled or named with its reason, so a new page cannot skip the lock.

Every link into `/references` found by the crawl must resolve, and every page under `/references` must carry the Out of universe banner
(`id="out-of-universe-banner"`; the section is DOC.88, `tests/test_references_pages.py` holds its pages).
"""

import re
from html.parser import HTMLParser
from pathlib import Path

import pytest

from planetgen.describe import sheet
from planetgen.util import plainwords
from planetgen.web import helpers
from planetgen.web.lib import classref, life_help, reference_urls

from tests.test_class_pages import app, client  # noqa: F401
from tests.test_plain_words import SEEDS, _document_text

ROOT = Path(__file__).resolve().parents[1]
GUIDE = ROOT.parent / "docs" / "design" / "interface" / "text-style.md"

BANNER = 'id="out-of-universe-banner"'
"""str: What a References page carries so a reader always knows they left the story (DOC.88 builds it)."""

ADVICE = ("\nSay it in universe and, where a reader may want the why, link to the References entry with reference_link(kind, id) "
          "(docs/design/interface/text-style.md, rule 16). There is no list of exceptions.")


# --- The addresses and the link helper -------------------------------------------------------------------------------------------

@pytest.mark.parametrize("kind, entry, path", [
    ("overview", None, "/references"),
    ("papers", None, "/references/papers"),
    ("papers", "R0123", "/references/papers#paper-r0123"),
    ("data", None, "/references/data"),
    ("data", "R4821", "/references/data#paper-r4821"),
    ("method", "atmospheres", "/references/method/atmospheres"),
    ("method", "atmospheres/surface_pressure", "/references/method/atmospheres#surface-pressure"),
    ("stand-in", None, "/references/stand-ins"),
    ("stand-in", "5.1:GR1", "/references/stand-ins/5-1-gr1"),
    ("fitted", None, "/references/fitted"),
    ("fitted", "AM10", "/references/fitted#am10"),
    ("speculative", None, "/references/speculative"),
    ("speculative", "3a", "/references/speculative#class-3a"),
    ("gaps", "life", "/references/gaps#life"),
    ("scorecard", None, "/references/scorecard"),
    ("about", None, "/references/about"),
])
def test_each_kind_of_entry_has_one_address(kind, entry, path):
    assert reference_urls.reference_path(kind, entry) == path
    assert reference_urls.is_reference_path(path.split("#")[0])


@pytest.mark.parametrize("kind, entry", [
    ("nope", None), ("method", None), ("method", "wizardry"), ("method", "wizardry/surface_pressure"), ("about", "x"), ("overview", "x"),
    ("papers", "!!"), ("papers", "10.1086/151682"), ("papers", "R12"), ("data", "Isochrones"), ("stand-in", ""), ("stand-in", "GR1"),
])
def test_a_typo_in_a_kind_or_an_id_is_an_error_not_a_dead_link(kind, entry):
    with pytest.raises(ValueError):
        reference_urls.reference_path(kind, entry)


def test_every_kind_is_documented_and_every_area_has_a_method_address():
    assert set(reference_urls.KINDS) == {"overview", "papers", "data", "method", "stand-in", "fitted", "speculative", "gaps", "scorecard", "about"}
    for area in reference_urls.AREAS:
        assert reference_urls.reference_path("method", area) == f"/references/method/{area}"
    module = Path(reference_urls.__file__).read_text(encoding="utf-8")
    for kind in reference_urls.KINDS:
        assert f"* `{kind}`:" in module, f"{kind} is not described in the module's docstring"


def test_the_anchor_is_lower_case_with_one_hyphen_for_each_run_of_other_characters():
    assert reference_urls.anchor("10.1086/151682") == "10-1086-151682"
    assert reference_urls.anchor("  Surface__Pressure (bar) ") == "surface-pressure-bar"
    assert reference_urls.anchor(42) == "42"
    with pytest.raises(ValueError):
        reference_urls.anchor("---")


def test_the_word_rule_is_exempt_where_the_addresses_begin():
    assert reference_urls.ROOT in plainwords.EXEMPT_PATHS


def test_only_the_references_addresses_are_references_paths():
    assert not reference_urls.is_reference_path("/referencesx")
    assert not reference_urls.is_reference_path("/help/references")
    assert reference_urls.is_reference_path("/references/stand-ins/5-1-gr1")


def test_the_link_helper_is_a_template_global_and_escapes_its_label(app):
    from flask import render_template_string
    with app.test_request_context("/"):
        html = render_template_string("{{ reference_link('stand-in', '5.1:GR1') }} {{ reference_link('about', label='<b>Credits</b>') }}")
    assert html == ('<a class="reference-link" href="/references/stand-ins/5-1-gr1">Where this comes from</a> '
                    '<a class="reference-link" href="/references/about">&lt;b&gt;Credits&lt;/b&gt;</a>')
    with pytest.raises(ValueError):
        helpers.reference_link("method", "wizardry")


def test_no_template_or_module_writes_a_references_address_by_hand():
    """Pages never hard-code a References address: they call `reference_link`, so the layout lives in `reference_urls.py` alone. The routes
    that serve the pages (DOC.88) and the module docstrings that describe the layout are the only other places the path appears."""
    allowed = {"reference_urls.py"}
    offenders = []
    for base in (ROOT / "planetgen" / "web", ROOT / "planetgen" / "describe"):
        for path in sorted(base.rglob("*")):
            if path.suffix not in (".py", ".html") or path.name in allowed or "__pycache__" in path.parts:
                continue
            for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                if re.search(r"""["'(]/references\b""", line) and "route(" not in line and not line.lstrip().startswith("#"):
                    offenders.append(f"{path.relative_to(ROOT)}:{number}: {line.strip()[:100]}")
    assert not offenders, "write the link with reference_link(kind, id):\n  " + "\n  ".join(offenders)


# --- The word list ---------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("text, words", [
    ("A stand-in was used.", ["stand-in"]),
    ("The values are fitted to a result; the model is a tweak of the last.", ["fitted", "model", "tweak"]),
    ("Speculative chemistry, tier T5, tweak L3.", ["speculative", "tier", "tweak"]),
    ("An educated guess from the random seed of the generator.", ["educated guess", "random seed", "generator"]),
    ("A seeded run keeps its seed value and its seed number.", ["seeded", "seed value", "seed number"]),
    ("This version was generated and simulated.", ["version", "generated", "simulated"]),
    ("The Model atmosphere and the modeled and modelled models.", ["model", "modeled", "modelled", "models"]),
    ("Sea level and ground level are in universe.", []),
    ("The seed-based and tier-less words stay, as do the modelling-free ones and a frontier.", []),
    ("A seed vault, the seeds of a plant, the seedling and a tierce are in universe.", []),
    ("The star is 5 pc away and the planet is about 1 AU out.", []),
])
def test_the_checker_flags_the_out_of_universe_words_and_only_them(text, words):
    found = [detail.split("'")[1] for rule, detail in plainwords.out_of_universe(text)]
    assert found == words
    assert all(rule == "out-of-universe" for rule, _detail in plainwords.out_of_universe(text))


@pytest.mark.parametrize("path, free", [
    ("/references", True), ("/references/papers", True), ("/references/stand-ins/5-1-gr1", True), ("/admin", True), ("/admin/settings", True),
    ("/referencesx", False), ("/help/life", False), ("/help/references", False), ("/system/ABC", False), ("/", False),
])
def test_the_references_and_admin_pages_may_use_the_words(path, free):
    assert plainwords.exempt(path) is free


def test_the_guide_names_the_list_the_exemption_and_the_link_helper():
    text = GUIDE.read_text(encoding="utf-8")
    for needle in ("`OUT_OF_UNIVERSE`", "`exempt(path)`", "`reference_link(kind, id)`", "tests/test_references_crawl.py", "id=\"out-of-universe-banner\"", "DOC.87"):
        assert needle in text, needle


# --- The crawl -------------------------------------------------------------------------------------------------------------------

class _Visible(HTMLParser):
    """The text a reader sees on a page: its text nodes outside scripts and styles, and the `title`, `alt`, `aria-label` and `placeholder` of its tags."""

    def __init__(self):
        super().__init__()
        self.parts, self.hidden, self.links = [], 0, []

    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style"):
            self.hidden += 1
        for name, value in attrs:
            if name in ("title", "alt", "aria-label", "placeholder") and value:
                self.parts.append(value)
            if name == "href" and value and reference_urls.is_reference_path(value.split("#")[0].split("?")[0]):
                self.links.append(value)

    def handle_endtag(self, tag):
        if tag in ("script", "style") and self.hidden:
            self.hidden -= 1

    def handle_data(self, data):
        if not self.hidden:
            self.parts.append(data)


def _read(html):
    """`(visible text, [addresses into /references])` of a page."""
    reader = _Visible()
    reader.feed(html)
    return " ".join(" ".join(reader.parts).split()), reader.links


def _paths():
    """The in-universe pages this suite can render without a database."""
    paths = ["/help", "/help/habitability", "/help/life", "/classes", "/references", "/login", "/preferences", "/polities", "/species"]
    paths += [f"/help/life/{slug}" for slug, _title, _summary in life_help.TOPICS if slug != "kinds"]
    paths += [f"/help/life/classes/{number}" for number in range(1, 10)]
    if classref.preload():
        for slug, kind in classref.catalog().items():
            paths.append(f"/classes/{slug}")
            paths += [f"/classes/{slug}/{code}" for code in kind["classes"]]
    return paths


@pytest.fixture(scope="module")
def crawl(app):
    """`[(path, endpoint, visible text, links into /references)]` for every page of `_paths`, read once."""
    client = app.test_client()
    adapter = app.url_map.bind("localhost")
    pages = []
    for path in _paths():
        response = client.get(path)
        assert response.status_code == 200, path
        text, links = _read(response.get_data(as_text=True))
        pages.append((path, adapter.match(path)[0], text, links))
    return pages


def _found(crawl):
    """`{(word, where): [example]}` for the out-of-universe words on the pages outside References and admin, in the Document of ten seeded
    systems and in the sheet's labels and source lines."""
    found = {}
    for path, endpoint, text, _links in crawl:
        if plainwords.exempt(path):
            continue
        for _rule, detail in plainwords.out_of_universe(text):
            found.setdefault((detail.split("'")[1], endpoint), []).append(path)
    for seed in SEEDS:
        for _owner, _kind, text in _document_text(seed):
            for _rule, detail in plainwords.out_of_universe(text):
                found.setdefault((detail.split("'")[1], "document"), []).append(text[:100])
    labels = {spec.label for kind in ("star", "planet", "moon", "belt", "comet") for spec in sheet.specs_for(kind)}
    for label in sorted(labels):
        for _rule, detail in plainwords.out_of_universe(label):
            found.setdefault((detail.split("'")[1], "sheet label"), []).append(label)
    return found


def test_the_crawl_reads_a_real_corpus(crawl):
    paths = {path for path, _endpoint, _text, _links in crawl}
    assert len(paths) > 15 and "/help/life/classes/1" in paths
    assert sum(len(text.split()) for _path, _endpoint, text, _links in crawl) > 5000


def test_no_page_outside_references_uses_the_out_of_universe_words(crawl):
    found = [f"{word!r} on {where}: {examples[0]}" for (word, where), examples in sorted(_found(crawl).items())]
    assert not found, "out-of-universe words on an in-universe page:\n  " + "\n  ".join(found) + ADVICE


def test_new_text_with_a_word_that_is_not_listed_fails():
    assert plainwords.out_of_universe("A brand new sentence that mentions its tier.") and not plainwords.out_of_universe("A brand new sentence.")


# --- The lock (DOC.94) ----------------------------------------------------------------------------------------------------------
#
# The crawl above reads what this suite can render without a database. The lock closes the rest: no in-universe template may carry a
# word of the list in the text a reader sees, even a page that needs data to render, and no new page can arrive without being either
# crawled or named here. Both run in the gate for changed pages (`scripts/test.py changed`): the scan opens every template and the
# crawl renders the Help and class pages, so a changed template, page builder or `plainwords.py` selects them.

TEMPLATES = ROOT / "planetgen" / "web" / "templates"

EXEMPT_TEMPLATE_PREFIXES = ("admin", "generate", "references")
"""tuple: The templates (by file name) that serve the admin pages (`/admin/...`, `admin_*.html`, `generate*.html`) and the References
pages (`references*.html`): the pages where `plainwords.exempt` allows the out-of-universe words."""

NEEDS_DATA = {
    "web.account": "behind the sign-in",
    "web.galaxy": "needs the galaxy",
    "web.index": "needs the database",
    "web.lists": "needs the database",
    "web.logout": "behind the sign-in",
    "web.maintenance_page": "behind the sign-in",
    "web.nav": "needs the database",
    "web.nav_chart": "needs the database",
    "web.nearby": "needs the database",
    "web.object_list": "needs the database",
    "web.object_page": "needs the database",
    "web.phenomena": "needs the database",
    "web.phenomenon": "needs the database",
    "web.polity_page": "needs the database",
    "web.search": "needs the database",
    "web.sector": "needs the database",
    "web.sector_at": "needs the database",
    "web.sector_on_galaxy_map": "a redirect to the Galaxy Map",
    "web.sectors": "needs the database",
    "web.species_page": "needs the database",
    "web.system": "needs the database",
    "web.system_on_galaxy_map": "a redirect to the Galaxy Map",
    "web.systems": "needs the database",
    "web.uncharted_phenomenon": "needs the database",
    "web.uncharted_star": "needs the database",
}
"""dict: `{endpoint: why the crawl cannot render it here}` for the in-universe pages that need data or a sign-in. The words in these pages'
templates are still read, by `_template_text`; the sentences the Document writes for a real system are read by tests/test_plain_words.py.
A new page is crawled (add its path to `_paths`) or listed here with the reason."""


def _is_exempt_template(path):
    return path.name.startswith(EXEMPT_TEMPLATE_PREFIXES)


def _templates():
    """`[(name relative to the template folder, path)]` for every template, in order."""
    return [(str(path.relative_to(TEMPLATES)), path) for path in sorted(TEMPLATES.rglob("*.html"))]


def _template_text(path):
    """The text a reader sees in a template, without the Jinja tags and comments (what a tag prints is read by the crawl)."""
    source = path.read_text(encoding="utf-8")
    for pattern in (r"\{#.*?#\}", r"\{%.*?%\}", r"\{\{.*?\}\}"):
        source = re.sub(pattern, " ", source, flags=re.S)
    return _read(source)[0]


def test_the_scan_reads_every_template_and_the_exempt_ones_are_the_admin_and_references_pages():
    names = [name for name, _path in _templates()]
    assert len(names) > 40 and "base.html" in names and "partials/contents.html" in names
    exempt = {name for name, path in _templates() if _is_exempt_template(path)}
    assert "references.html" in exempt and "admin.html" in exempt and "generate_task.html" in exempt
    assert not exempt & {"base.html", "system.html", "sector.html", "help_index.html", "species.html"}


def test_no_template_outside_references_and_admin_has_an_out_of_universe_word():
    found = []
    for name, path in _templates():
        if _is_exempt_template(path):
            continue
        for _rule, detail in plainwords.out_of_universe(_template_text(path)):
            found.append(f"{detail.split(chr(39))[1]!r} in templates/{name}")
    assert not found, "out-of-universe words in the text of an in-universe template:\n  " + "\n  ".join(found) + ADVICE


def test_every_in_universe_page_is_crawled_or_listed_with_its_reason(app, crawl):
    crawled = {endpoint for _path, endpoint, _text, _links in crawl}
    pages = {}
    for rule in app.url_map.iter_rules():
        view = app.view_functions[rule.endpoint]
        if rule.endpoint.startswith("web.") and "GET" in rule.methods and not getattr(view, "json_only", False) and not plainwords.exempt(rule.rule):
            pages.setdefault(rule.endpoint, rule.rule)
    new = sorted(f"{endpoint} ({rule})" for endpoint, rule in pages.items() if endpoint not in crawled and endpoint not in NEEDS_DATA)
    assert not new, f"pages the crawl does not read and NEEDS_DATA does not name (add the path to _paths, or the endpoint with its reason): {new}"
    gone = sorted(endpoint for endpoint in NEEDS_DATA if endpoint not in pages)
    assert not gone, f"NEEDS_DATA names pages that are gone: {gone}"
    both = sorted(endpoint for endpoint in NEEDS_DATA if endpoint in crawled)
    assert not both, f"crawled now, delete from NEEDS_DATA: {both}"


def _references_routes(app):
    return sorted({rule.rule for rule in app.url_map.iter_rules() if reference_urls.is_reference_path(rule.rule) and "GET" in rule.methods
                   and "<" not in rule.rule})


def test_every_link_into_references_on_the_crawled_pages_resolves(app, crawl):
    client = app.test_client()
    links = sorted({link for _path, _endpoint, _text, found in crawl for link in found})
    broken = [link for link in links if client.get(link.split("#")[0]).status_code != 200]
    assert not broken, f"links into /references that do not resolve: {broken}"


def test_every_references_page_carries_the_out_of_universe_banner(app):
    routes = _references_routes(app)
    assert "/references" in routes, "the start page of the section is gone"
    client = app.test_client()
    missing = [route for route in routes if BANNER not in client.get(route).get_data(as_text=True)]
    assert not missing, f"References pages without the Out of universe banner ({BANNER}): {missing}"
