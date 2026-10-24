"""
The architecture map (DOC.67, `scripts/build_architecture_map.py`, `docs/generated/architecture-map.html`): the page
lists every module of the package, resolves calls the way its notes say, traces its flows from named entry functions,
and opens offline and drills down from the project map to a function's source.
"""

import importlib.util
import os

import pytest

pytest.importorskip("playwright.sync_api")

from tests.test_web_a11y import browser  # noqa: E402,F401 -- fixture

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_spec = importlib.util.spec_from_file_location("build_architecture_map", os.path.join(ROOT, "scripts", "build_architecture_map.py"))
amap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(amap)
PAGE = amap.PAGE


@pytest.fixture(scope="module")
def collected():
    return amap.collect()


@pytest.fixture(scope="module")
def data(collected):
    return amap.build_data(collected)


def _page():
    with open(PAGE, encoding="utf-8") as handle:
        return handle.read()


@pytest.mark.globaldocs
def test_the_page_lists_every_module_of_the_package(collected):
    problems = amap.module_problems(_page(), collected)
    assert not problems, "\n".join(problems) + "\nthe TODO thread runs python scripts/build_global_docs.py"


@pytest.mark.globaldocs
def test_a_new_test_module_needs_no_regeneration(collected):
    assert not amap.module_problems(_page().replace('"name":"tests.test_architecture_map"', '"name":"tests.other"'), collected)


def test_a_module_missing_from_the_page_is_reported(collected):
    light = amap.page_light(_page())
    victim = next(m["name"] for m in light["modules"] if m["name"].startswith("planetgen.queue."))
    text = _page().replace(f'"name":"{victim}"', '"name":"planetgen.gone_module"')
    problems = amap.module_problems(text, collected)
    assert f"missing from the page: {victim}" in problems
    assert "on the page but not in the code: planetgen.gone_module" in problems
    assert amap.module_problems("<html></html>", collected)


def test_every_area_is_named_and_every_package_module_is_in_one(data):
    light, _heavy = data
    keys = {area["key"] for area in light["areas"]}
    assert {"web", "api", "queue", "generation", "galaxy", "population", "physics", "describe", "db", "admin",
            "backup", "cli", "names", "util", "tuning", "tests"} <= keys
    placed = sorted(i for area in light["areas"] for i in area["modules"])
    assert placed == list(range(len(light["modules"])))
    assert all(area["purpose"] for area in light["areas"]), "an area has no purpose line"


def test_every_function_of_the_package_is_reachable_from_its_module(data):
    light, heavy = data
    modules = {i for i, m in enumerate(light["modules"]) if m["area"] != "tests"}
    assert {f[0] for f in heavy["functions"]} <= modules
    for module in modules:
        source = heavy["sources"].get(str(module))
        assert source is not None, light["modules"][module]["name"]


def test_calls_are_resolved_the_way_the_notes_say(data):
    light, heavy = data
    names = [m["name"] for m in light["modules"]]
    functions = heavy["functions"]
    by_key = {(names[f[0]], f[1]): i for i, f in enumerate(functions)}
    submit = functions[by_key[("planetgen.queue.api_jobs", "submit")]]
    resolved = {(names[functions[t][0]], functions[t][1]) for t, _kind in submit[6]}
    assert ("planetgen.queue.redisqueue", "redis_url") in resolved        # `redisqueue.redis_url()`: a module alias
    assert ("planetgen.queue.api_jobs", "new_job_id") in resolved         # a function of the module itself
    assert ("planetgen.queue.workers", "ensure_workers") in resolved      # `workers.ensure_workers()`
    assert ".enqueue" in submit[7]                                        # a method on an unknown value stays text
    assert not any(t in ("len", "int") for t in submit[7])                # builtins are left out


def test_a_method_call_on_an_unknown_value_is_only_resolved_when_the_name_is_unique(data):
    _light, heavy = data
    probable = [c for f in heavy["functions"] for c in f[6] if c[1] == 1]
    assert probable, "no call resolved by name at all"
    unique = {}
    for f in heavy["functions"]:
        if f[2] in ("function", "method"):
            unique.setdefault(f[1].rpartition(".")[2], []).append(f)
    for f in heavy["functions"]:
        for target, kind in f[6]:
            if kind == 1:
                assert len(unique[heavy["functions"][target][1].rpartition(".")[2]]) == 1


def test_the_flows_trace_the_five_paths_of_the_notes(data):
    _light, heavy = data
    keys = {flow["key"] for flow in heavy["flows"]}
    assert {"generate", "scatter", "page", "api-job", "tile", "backup"} <= keys
    for flow in heavy["flows"]:
        assert flow["entries"], flow["key"]
        assert any(heavy["functions"][e][6] for e in flow["entries"]), f"flow {flow['key']} traces no call"


def test_a_flow_naming_a_function_that_is_gone_fails_the_build(monkeypatch, collected):
    monkeypatch.setattr(amap, "FLOWS", [{"key": "x", "title": "x", "text": "x", "entries": ["planetgen.nowhere:gone"]}])
    with pytest.raises(SystemExit, match="no function"):
        amap.build_data(collected)


@pytest.mark.globaldocs
def test_the_check_command_reports_a_stale_module_list(tmp_path, monkeypatch, capsys):
    stale = tmp_path / "architecture-map.html"
    stale.write_text(_page().replace('"name":"planetgen.tuning"', '"name":"planetgen.tuning_old"'), encoding="utf-8")
    monkeypatch.setattr(amap, "PAGE", str(stale))
    assert amap.main(["--check"]) == 1
    assert "out of date" in capsys.readouterr().out
    monkeypatch.setattr(amap, "PAGE", PAGE)
    assert amap.main(["--check"]) == 0


def _heading(page, text):
    """Waits for the view to be the one named `text` (a click changes the hash, and the page renders after it)."""
    page.wait_for_function("t => (document.querySelector('h2') || {textContent: ''}).textContent.trim().startsWith(t)", arg=text)


def test_the_page_opens_offline_and_drills_down_to_a_functions_source(browser):
    page = browser.new_page()
    outside, errors = [], []
    page.on("request", lambda r: outside.append(r.url) if not r.url.startswith(("file:", "data:", "blob:")) else None)
    page.on("pageerror", lambda e: errors.append(str(e)))
    try:
        page.goto("file://" + PAGE)
        page.wait_for_selector("body[data-ready='1']", timeout=60000)
        assert page.locator("rect.box").count() >= 15
        page.click("g[data-area='queue'] rect", force=True)
        _heading(page, "queue")
        page.click("td > a[href='#/mod/planetgen.queue.units']")
        _heading(page, "queue.units")
        page.click("tbody tr a:text-is('run_unit')")
        _heading(page, "run_unit")
        assert "def run_unit" in page.inner_text("pre.src")
        assert page.locator("ul.edges li").count() > 0
        page.goto("file://" + PAGE + "#/flow")
        _heading(page, "Flows")
        page.wait_for_selector("section")
        assert page.locator("section").count() >= 5
        page.fill("#search", "ensure_workers")
        assert page.locator("#results a").count() >= 1
    finally:
        page.close()
    assert not outside and not errors
