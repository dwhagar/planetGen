# tests/test_web_browser_object_tree.py

"""
UX.172: `static/objecttree.js` in a real browser, on a small served tree with no database: a closed row loads its
children when opened and again, a page at a time, from "Show more"; a failed load says so and tries again; the arrow
keys move through the visible rows as the ARIA tree pattern says. Skipped without Playwright and Chromium.
"""

import os

import flask
import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from planetgen.web import STATIC_DIR  # noqa: E402
from planetgen.web.lib import objecttree  # noqa: E402
from tests.served import served  # noqa: E402

CONTEXT = objecttree.TreeContext(item_url=lambda ref: "/object/" + ref, tree_url=lambda ref: "/tree/" + ref)


def _planet(number):
    return {"ref": f"planet:{number}", "kind": "planet", "id": str(number), "name": f"Planet {number}",
            "child_count": 1 if number == 1 else 0, "counts": {"moon": 1} if number == 1 else {}, "summary": {}}


def _app():
    app = flask.Flask(__name__, static_folder=STATIC_DIR, static_url_path="/static")
    app.config["fail"] = 0

    @app.route("/")
    def index():
        star = {"ref": "star:0", "kind": "star", "id": "0", "name": "Sun", "child_count": 3, "counts": {"planet": 3},
                "summary": {}}
        return ('<!doctype html><html lang="en"><head><meta charset="utf-8"><title>t</title>'
                '<link rel="stylesheet" href="/static/style.css"></head><body><main id="main"><h1>Bodies</h1>'
                + objecttree.tree_html([star], CONTEXT, "Bodies")
                + '</main><script type="module" src="/static/objecttree.js"></script></body></html>')

    @app.route("/tree/<ref>")
    def tree(ref):
        if app.config["fail"]:
            app.config["fail"] -= 1
            return flask.Response("no", status=500)
        if ref == "star:0":
            after = flask.request.args.get("after")
            items = [_planet(1), _planet(2)] if not after else [_planet(3)]
            token = None if after else "more"
            return flask.jsonify({"html": objecttree.nodes_html(items, CONTEXT), "next": token})
        return flask.jsonify({"html": objecttree.nodes_html([{"ref": "moon:9", "kind": "moon", "id": "9", "name": "Luna",
                                                              "child_count": 0, "counts": {}, "summary": {}}], CONTEXT),
                              "next": None})

    return app


@pytest.fixture(scope="module")
def server():
    app = _app()
    with served(app) as url:
        yield url, app


@pytest.fixture(scope="module")
def browser():
    try:
        playwright = sync_api.sync_playwright().start()
        chromium = playwright.chromium.launch()
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"No Chromium for Playwright ({str(exc).splitlines()[0]})")
    try:
        yield chromium
    finally:
        chromium.close()
        playwright.stop()


@pytest.fixture()
def page(browser, server):
    url, app = server
    app.config["fail"] = 0
    page = browser.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.goto(url + "/")
    page.wait_for_selector("[role=tree]")
    yield page
    page.close()
    assert not errors, errors


def _names(page):
    return page.eval_on_selector_all("a.tree-name", "els => els.map(e => e.textContent)")


def test_opening_a_row_loads_its_children_then_show_more_loads_the_rest(page):
    page.click("button.tree-toggle >> nth=0")
    page.wait_for_selector("a.tree-name:text('Planet 2')")
    assert _names(page) == ["Sun", "Planet 1", "Planet 2"]
    page.click("[data-tree-more]")
    page.wait_for_selector("a.tree-name:text('Planet 3')")
    assert _names(page) == ["Sun", "Planet 1", "Planet 2", "Planet 3"]
    assert page.locator("[data-tree-more]").count() == 0
    assert page.locator("li.tree-pending").count() == 1  # Planet 1's own, still unopened


def test_a_failed_load_says_so_and_the_row_tries_again_when_reopened(page, server):
    _url, app = server
    app.config["fail"] = 1
    toggle = page.locator("button.tree-toggle").first
    toggle.click()
    page.wait_for_selector(".tree-error")
    toggle.click()
    toggle.click()
    page.wait_for_selector("a.tree-name:text('Planet 2')")
    assert page.locator(".tree-error").count() == 0


def test_the_arrow_keys_open_close_and_move_between_rows(page):
    page.focus("a.tree-name:text('Sun')")
    page.keyboard.press("ArrowRight")
    page.wait_for_selector("a.tree-name:text('Planet 2')")
    page.keyboard.press("ArrowRight")
    assert page.evaluate("document.activeElement.textContent") == "Planet 1"
    page.keyboard.press("ArrowRight")
    page.wait_for_selector("a.tree-name:text('Luna')")
    page.keyboard.press("ArrowDown")
    assert page.evaluate("document.activeElement.textContent") == "Luna"
    page.keyboard.press("ArrowLeft")
    assert page.evaluate("document.activeElement.textContent") == "Planet 1"
    page.keyboard.press("End")  # the last row of the level is its "Show more" button
    assert page.evaluate("document.activeElement.textContent") == "Show more"
    page.keyboard.press("Home")
    assert page.evaluate("document.activeElement.textContent") == "Sun"
    assert page.locator("a.tree-name[tabindex='0']").count() == 1
    assert page.get_attribute("li.tree-node", "aria-expanded") == "true"
    # A closed row's children are out of the keyboard's reach until it is opened again.
    page.focus("a.tree-name:text('Sun')")
    page.keyboard.press("ArrowLeft")
    assert page.is_hidden("a.tree-name:text('Planet 2')")
    page.keyboard.press("ArrowDown")
    assert page.evaluate("document.activeElement.textContent") == "Sun"


def test_an_open_tree_has_no_serious_accessibility_violation(page):
    """The tree with an opened row, its links and its toggles, checked with axe-core (WCAG 2.0/2.1 A and AA): a link inside
    a `<summary>` or a button is a nested control, and a link in the summary line must not rely on colour alone."""
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vendor", "axe-core", "axe.min.js"), encoding="utf-8") as axe:
        page.add_script_tag(content=axe.read())
    page.click("button.tree-toggle >> nth=0")
    page.wait_for_selector("a.tree-name:text('Planet 2')")
    page.click("button.tree-toggle >> nth=1")
    page.wait_for_selector("a.tree-name:text('Luna')")
    found = page.evaluate("""async () => {
        const result = await axe.run(document, {runOnly: {type: "tag", values: ["wcag2a", "wcag2aa", "wcag21a", "wcag21aa"]}});
        return result.violations.filter(v => v.impact === "serious" || v.impact === "critical")
            .map(v => v.id + ": " + v.nodes.map(n => n.target.join(" ")).join(" | "));
    }""")
    assert not found, found
