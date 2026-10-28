# tests/test_web_browser_bodyview.py

"""UX.171: `static/bodyview.js` in a real browser: a look record in a canvas's `data-look` is painted (lit sphere, clear corners)
and the hidden figure around it is shown. Skipped without Playwright and Chromium."""

import flask
import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.web import STATIC_DIR  # noqa: E402
from tests.served import served  # noqa: E402
from tests.test_bodyview import PLANET, STAR  # noqa: E402
from tests.test_web_browser_object_tree import browser  # noqa: E402,F401 -- fixture
from planetgen.web.lib import bodyview  # noqa: E402


def _app():
    app = flask.Flask(__name__, static_folder=STATIC_DIR, static_url_path="/static")

    @app.route("/<kind>")
    def page(kind):
        # MAP.328: "turning" is a planet with a stored 2-hour day, "backwards" the same with its axis tilted 180 degrees, "still"
        # one with no stored period, "star" a star with a day of 600 hours.
        row = STAR if kind in ("star", "starspin") else PLANET
        if kind == "starspin":
            row = dict(STAR, rotation_period_hours=600.0)
        if kind in ("turning", "backwards"):
            row = dict(PLANET, rotation_period_hours=2.0, axial_tilt_deg=0.0 if kind == "turning" else 180.0)
        shown = bodyview.view("star" if kind in ("star", "starspin") else "planet", row, [STAR])
        return ('<!doctype html><title>t</title><figure hidden><canvas width="320" height="320" data-look="'
                + shown["look"].replace('"', "&quot;")
                + '"></canvas><button type="button" data-spin aria-pressed="true" hidden>Spin</button></figure>'
                + '<script type="module" src="/static/bodyview.js"></script>')

    return app


@pytest.fixture(scope="module")
def url():
    with served(_app()) as base:
        yield base


@pytest.mark.parametrize("kind", ["star", "planet"])
def test_the_canvas_is_painted_and_its_figure_shown(browser, url, kind):
    page = browser.new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    try:
        page.goto(f"{url}/{kind}")
        page.wait_for_selector("canvas[data-drawn]")
        assert page.is_visible("figure")
        stats = page.evaluate("""() => {
            const c = document.querySelector('canvas'), d = c.getContext('2d').getImageData(0, 0, 320, 320).data;
            const at = (x, y) => d[(y * 320 + x) * 4 + 3];
            let opaque = 0;
            for (let i = 3; i < d.length; i += 4) if (d[i] > 250) opaque++;
            return {corner: at(2, 2), middle: at(160, 160), opaque};
        }""")
        assert stats["corner"] == 0 and stats["middle"] == 255 and stats["opaque"] > 20000, stats
    finally:
        page.close()
    assert not errors, errors


MIDDLE_ROW = """() => {
    const c = document.querySelector('canvas'), d = c.getContext('2d').getImageData(0, 160, 320, 1).data;
    const row = [];
    for (let x = 0; x < 320; x++) row.push(d[x * 4] + d[x * 4 + 1] + d[x * 4 + 2]);
    return {row, angle: Number(c.dataset.angle)};
}"""


def _sample(page):
    return page.evaluate(MIDDLE_ROW)


def _shift(before, after):
    """How many pixels the equator's pattern moved to the right between two samples (negative: to the left): the shift that
    makes the first row best match the second, over the sphere's middle."""
    best = None
    for s in range(-70, 71):
        error = sum(abs(before[x] - after[x + s]) for x in range(110, 210))
        if best is None or error < best[0]:
            best = (error, s)
    return best[1]


def _opened(browser, url, kind, **context):
    context.setdefault("viewport", {"width": 800, "height": 700})
    page = browser.new_context(**context).new_page()
    errors = []
    page.on("pageerror", lambda exc: errors.append(str(exc)))
    page.goto(f"{url}/{kind}")
    page.wait_for_selector("canvas[data-drawn]")
    return page, errors


def test_a_body_with_a_stored_day_turns_at_the_shared_rate_and_a_tilted_one_turns_backwards(browser, url):
    """MAP.328: the rendered example turns about its own axis: its angle follows spin.js's rate for its stored period (a 2-hour
    day is the fastest turn, one in 20 s), the surface moves to the right on a body with its axis up and to the left on one
    tilted 180 degrees."""
    shifts = {}
    for kind in ("turning", "backwards"):
        page, errors = _opened(browser, url, kind)
        try:
            first = _sample(page)
            page.wait_for_timeout(500)
            second = _sample(page)
            while second["angle"] - first["angle"] < 0.05 and not errors:
                page.wait_for_timeout(100)
                second = _sample(page)
            rate = 2 * 3.141592653589793 / 20
            elapsed_rate = (second["angle"] - first["angle"])
            assert elapsed_rate > 0, (first["angle"], second["angle"])
            assert second["angle"] - first["angle"] <= rate * 1.5 + 0.1, (first["angle"], second["angle"])
            shifts[kind] = _shift(first["row"], second["row"])
        finally:
            page.context.close()
        assert not errors, errors
    assert shifts["turning"] > 0 > shifts["backwards"], shifts


def test_a_star_with_a_stored_day_turns_and_a_body_with_none_stays_still(browser, url):
    page, errors = _opened(browser, url, "starspin")
    try:
        first = _sample(page)["angle"]
        page.wait_for_timeout(1200)
        assert _sample(page)["angle"] > first
        assert page.is_visible("[data-spin]")
    finally:
        page.context.close()
    still, errors2 = _opened(browser, url, "planet")
    try:
        first = _sample(still)
        still.wait_for_timeout(1000)
        assert _sample(still)["angle"] == first["angle"] and _sample(still)["row"] == first["row"]
        assert not still.is_visible("[data-spin]"), "a body with no stored period has no Spin switch"
    finally:
        still.context.close()
    assert not errors and not errors2, (errors, errors2)


def test_the_spin_switch_and_reduced_motion_hold_the_picture_still(browser, url):
    page, errors = _opened(browser, url, "turning")
    try:
        assert page.get_attribute("[data-spin]", "aria-pressed") == "true"
        page.click("[data-spin]")
        assert page.get_attribute("[data-spin]", "aria-pressed") == "false"
        page.wait_for_timeout(300)
        held = _sample(page)
        page.wait_for_timeout(1000)
        assert _sample(page)["angle"] == held["angle"] and _sample(page)["row"] == held["row"]
        page.click("[data-spin]")
        page.wait_for_timeout(1000)
        assert _sample(page)["angle"] > held["angle"]
    finally:
        page.context.close()
    quiet, errors2 = _opened(browser, url, "turning", reduced_motion="reduce")
    try:
        first = _sample(quiet)
        quiet.wait_for_timeout(1000)
        assert _sample(quiet)["angle"] == first["angle"]
        assert quiet.is_disabled("[data-spin]") and quiet.get_attribute("[data-spin]", "aria-pressed") == "false"
    finally:
        quiet.context.close()
    assert not errors and not errors2, (errors, errors2)
