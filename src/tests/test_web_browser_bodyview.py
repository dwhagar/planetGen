# tests/test_web_browser_bodyview.py

"""UX.171: `static/bodyview.js` in a real browser: a look record in a canvas's `data-look` is painted (lit sphere, clear corners)
and the hidden figure around it is shown. Skipped without Playwright and Chromium."""

import json

import flask
import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.web import STATIC_DIR  # noqa: E402
from tests.served import served  # noqa: E402
from tests.test_bodyview import PLANET, STAR  # noqa: E402
from tests.test_web_browser_object_tree import browser  # noqa: E402,F401 -- fixture
from planetgen.web.lib import bodyview  # noqa: E402


AIRS = {
    "noair": {},
    "thinair": {"atm": True, "glow": "#88aaff", "glowStrength": 0.02, "glowWidth": 0.03, "twilight": 0.005},
    "thickair": {"atm": True, "glow": "#88aaff", "glowStrength": 0.8, "glowWidth": 0.2, "twilight": 0.2},
}
"""The air of a still world in `_app`: none has no glow fields, a thin one almost none and a thick one the most (MAP.277)."""


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
        # UX.244: "locked" turns once an orbit (a 0.1-year orbit of 876.6 hours), "slow" has a 46.8-day rotation, "tilted" is Uranus-like.
        if kind == "locked":
            row = dict(PLANET, rotation_period_hours=0.1 * 365.25 * 24.0, period_years=0.1, axial_tilt_deg=3.0)
        if kind == "slow":
            row = dict(PLANET, rotation_period_hours=46.8 * 24.0, period_years=1.0, axial_tilt_deg=0.0)
        if kind in ("tilted", "upright"):
            row = dict(PLANET, rotation_period_hours=17.2, period_years=84.0, axial_tilt_deg=97.8 if kind == "tilted" else 0.0)
        shown = bodyview.view("star" if kind in ("star", "starspin") else "planet", row, [STAR])
        look = shown["look"]
        if kind in AIRS:     # MAP.277: the same still world with each kind of air (the glow fields the look record gets from `bodylook.air_glow`)
            record = json.loads(look)
            for key in ("glow", "glowStrength", "glowWidth", "twilight", "atm"):
                record.pop(key, None)
            record.update(AIRS[kind])
            look = json.dumps(record)
        return ('<!doctype html><title>t</title><figure hidden><canvas width="320" height="320" data-look="'
                + look.replace('"', "&quot;")
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
    const c = document.querySelector('canvas'), context = c.getContext('2d');
    const rows = [100, 115, 130, 145, 160, 175, 190, 205].map((y) => {
        const d = context.getImageData(0, y, 320, 1).data, row = [];
        for (let x = 0; x < 320; x++) row.push(d[x * 4] + d[x * 4 + 1] + d[x * 4 + 2]);
        return row;
    });
    return {row: rows[4], rows, angle: Number(c.dataset.angle)};
}"""


def _sample(page):
    return page.evaluate(MIDDLE_ROW)


def _shift(before, after):
    """How many pixels the surface's pattern moved to the right between two samples (negative: to the left), on the lit half of the
    sphere: the shift that makes the row of the first sample with the most pattern best match the second. The rows are compared by their
    first differences, so the slope of the shading toward the day and night line does not count as a pattern."""
    def slope(row):
        return [row[x + 1] - row[x] for x in range(len(row) - 1)]

    index = max(range(len(before["rows"])), key=lambda n: sum(abs(v) for v in slope(before["rows"][n])[75:135]))
    first, second = slope(before["rows"][index]), slope(after["rows"][index])
    best = None
    for s in range(-25, 26):
        error = sum(abs(first[x] - second[x + s]) for x in range(75, 135))
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


TURN_MIN = 0.08
TURN_MAX = 0.17
"""float: The turn, in radians, that two samples of a turning body must be apart: enough to see which way the surface moved, and not so much that it moved further than `_shift` looks (about 25 px, 0.18 rad on the 320 px sphere), where a pattern that repeats reads as moving the other way."""


def _turned(page, errors, attempts=20):
    """Two samples of the canvas whose angles differ by `TURN_MIN` to `TURN_MAX` radians (TEST.518). The interval is chosen by the angle the
    page reports, not by a sleep: a fixed wait turned the body 0.17 radians on a quiet machine, which is as far as `_shift` can see, and
    further on a busy one, so the shift came out the wrong way. A pair that turned too far (a stall between two polls) is dropped and the
    wait starts again from the newest sample."""
    first = second = _sample(page)
    turned = 0.0
    for _ in range(attempts):
        while not errors:
            page.wait_for_timeout(20)
            second = _sample(page)
            turned = second["angle"] - first["angle"]
            if turned >= TURN_MIN:
                break
        if errors or turned <= TURN_MAX:
            return first, second
        first = second
    raise AssertionError(f"no two samples {TURN_MIN} to {TURN_MAX} radians apart in {attempts} tries")


def test_a_body_with_a_stored_day_turns_at_the_shared_rate_and_a_tilted_one_turns_backwards(browser, url):
    """MAP.328: the rendered example turns about its own axis: its angle follows spin.js's rate for its stored period (a 2-hour
    day is the fastest turn, one in 20 s), the surface moves to the right on a body with its axis up and to the left on one
    tilted 180 degrees."""
    shifts = {}
    for kind in ("turning", "backwards"):
        page, errors = _opened(browser, url, kind)
        try:
            first, second = _turned(page, errors)
            rate = 2 * 3.141592653589793 / 20
            elapsed_rate = (second["angle"] - first["angle"])
            assert elapsed_rate > 0, (first["angle"], second["angle"])
            assert second["angle"] - first["angle"] <= rate * 1.5 + 0.1, (first["angle"], second["angle"])
            shifts[kind] = _shift(first, second)
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


def _halves(page, y=160):
    """The brightness of the left and the right part of the disc, along its middle row `y` (a locked world's disc is fitted into the picture
    right of the arrow, so its middle row and the parts are those of its disc)."""
    return page.evaluate("""(y) => {
        const c = document.querySelector('canvas'), d = c.getContext('2d').getImageData(0, y, 320, 1).data;
        const disc = c.dataset.disc ? JSON.parse(c.dataset.disc) : {x: 160, radius: 154};
        const from = (a, b) => { let sum = 0; for (let x = Math.round(disc.x + a * disc.radius); x < Math.round(disc.x + b * disc.radius); x++) sum += d[x * 4] + d[x * 4 + 1] + d[x * 4 + 2]; return sum; };
        return {left: from(-0.85, -0.1), right: from(0.1, 0.85)};
    }""", y)


def _brightness(page, xs, y=160):
    """The brightness (red plus green plus blue) of the pixel at each x along row `y`."""
    return page.evaluate("""([xs, y]) => {
        const d = document.querySelector('canvas').getContext('2d').getImageData(0, y, 320, 1).data;
        return xs.map((x) => d[x * 4] + d[x * 4 + 1] + d[x * 4 + 2]);
    }""", [xs, y])


@pytest.mark.parametrize("kind", ["upright", "tilted", "slow"])
def test_a_world_shows_a_curved_day_and_night_line_with_more_than_half_lit_and_keeps_turning(browser, url, kind):
    """UX.250 (was UX.244): the star is to the left and a little in front, so the day and night line is a curve bowing into the night side with
    more than half the disc lit: along the middle row the disc is still bright past its centre (the line was the centre) and dark by its right
    limb, and the surface turns past it."""
    page, errors = _opened(browser, url, kind)
    try:
        halves = _halves(page)
        assert halves["left"] > 1.5 * halves["right"], halves
        left_limb, past_centre, right_limb = _brightness(page, [114, 206, 275])
        assert past_centre > 1.8 * right_limb and left_limb > past_centre, (left_limb, past_centre, right_limb)
        first = _sample(page)["angle"]
        page.wait_for_timeout(1500)
        assert _sample(page)["angle"] > first, kind
        assert not page.evaluate("document.querySelector('canvas').dataset.locked")
    finally:
        page.context.close()
    assert not errors, errors


def test_a_slow_rotator_turns_visibly_not_at_the_floor_of_the_scale(browser, url):
    """UX.244: a 46.8-day rotation turns about once in four minutes, not at the one-in-ten-minutes floor the old scale gave every period over 1000 hours."""
    page, errors = _opened(browser, url, "slow")
    try:
        first = _sample(page)["angle"]
        page.wait_for_timeout(2000)
        turned = _sample(page)["angle"] - first
        assert turned > 2 * 3.141592653589793 / 600 * 1.5, turned
    finally:
        page.context.close()
    assert not errors, errors


def test_a_locked_world_holds_its_day_side_toward_the_star_and_says_so(browser, url):
    """UX.244, UX.245: a tidally locked world is drawn still with its star-facing point at the lit limb, has no Spin switch, says it is locked in its caption and on the picture."""
    page, errors = _opened(browser, url, "locked")
    try:
        first = _sample(page)
        halves = _halves(page, y=146)
        assert halves["left"] > 1.5 * halves["right"], halves
        page.wait_for_timeout(1200)
        again = _sample(page)
        assert again["angle"] == first["angle"] and again["row"] == first["row"]
        assert page.evaluate("document.querySelector('canvas').dataset.locked") == "1"
        assert not page.is_visible("[data-spin]")
        # UX.245: the picture itself carries the label (the strip under the disc is opaque where the label is) and the arrow on the left.
        marks = page.evaluate("""() => {
            const d = document.querySelector('canvas').getContext('2d').getImageData(0, 0, 320, 320).data;
            const at = (x, y) => d[(y * 320 + x) * 4 + 3];
            const disc = JSON.parse(document.querySelector('canvas').dataset.disc);
            return {label: at(160, 304), arrow: at(10, Math.round(disc.y)), free: at(2, 2), disc};
        }""")
        assert marks["label"] > 100 and marks["arrow"] > 100 and marks["free"] == 0, marks
        # UX.250: the arrow and its word are left of the disc and the label is under it, so neither covers it.
        disc = marks["disc"]
        assert disc["x"] - disc["radius"] >= 50 and disc["y"] + disc["radius"] <= 320 - 26, disc
    finally:
        page.context.close()
    assert not errors, errors
    row = dict(PLANET, rotation_period_hours=0.1 * 365.25 * 24.0, period_years=0.1)
    assert "Tidally locked" in bodyview.view("planet", row, [STAR])["note"]
    assert "turning much faster" in bodyview.view("planet", dict(row, rotation_period_hours=17.0), [STAR])["note"]


def _rim_and_dusk(page):
    """`(alpha just outside the limb on the lit side, brightness a little past the day and night line)` of the picture."""
    return page.evaluate("""() => {
        const c = document.querySelector('canvas'), disc = JSON.parse(c.dataset.disc), d = c.getContext('2d').getImageData(0, 0, c.width, c.height).data;
        const px = (nx, ny) => { const x = Math.round(disc.x + nx * disc.radius), y = Math.round(disc.y + ny * disc.radius); return (y * c.width + x) * 4; };
        const outside = d[px(-1.03, 0) + 3];
        const t = px(0.6, 0);                        // right of the centre: past the line, in the night side (the star is to the left and in front)
        return [outside, d[t] + d[t + 1] + d[t + 2]];
    }""")


def test_the_glow_and_the_dusk_follow_the_air_of_the_world(browser, url):
    """MAP.277: with no air there is no glow past the limb and the day and night line is hard; a thin air has hardly any glow; a thick air glows and
    carries the light a little past the line, so the night side just past it is brighter."""
    seen = {}
    for kind in AIRS:
        page, errors = _opened(browser, url, kind)
        try:
            seen[kind] = _rim_and_dusk(page)
        finally:
            page.context.close()
        assert not errors, errors
    assert seen["noair"][0] == 0, seen
    assert seen["thinair"][0] < 10 and seen["thickair"][0] > 40, seen
    assert seen["thickair"][1] > 1.5 * seen["noair"][1], seen
