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
        row = STAR if kind == "star" else PLANET
        shown = bodyview.view(kind, row, [STAR])
        return ('<!doctype html><title>t</title><figure hidden><canvas width="320" height="320" data-look="'
                + shown["look"].replace('"', "&quot;")
                + '"></canvas></figure><script type="module" src="/static/bodyview.js"></script>')

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
