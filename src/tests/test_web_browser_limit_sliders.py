"""
UX.116: the Generate page's limit sliders work in a browser. The slider moves its number box and a typed value moves the
slider; the core limits follow the galaxy-wide ones until one is set; a segment's Reset returns it to the defaults.
Skipped without Playwright, Chromium or the MySQL test server, like `test_web_a11y.py`.
"""

import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.api.authz import SESSION_COOKIE_NAME  # noqa: E402

from tests.test_web_a11y import (  # noqa: E402,F401 -- fixtures
    admin_token,
    base_url,
    browser,
    site_app,
    site_db,
)

PLAN = "plan-mass-limit"
"""The id prefix of the Plan form's scatter limits."""


@pytest.fixture()
def page(browser, base_url, admin_token):
    context = browser.new_context(viewport={"width": 1280, "height": 1400}, reduced_motion="reduce")
    context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
    page = context.new_page()
    page.goto(f"{base_url}/admin/generate")
    # The page's modules run in order after the load event under load: wait for the last of them to finish, or the
    # folds' own script closes the sections opened below.
    page.wait_for_function("() => document.querySelector('[data-ceiling-cost]')?.textContent.trim() !== ''")
    page.evaluate("() => document.querySelectorAll('details.fold').forEach((d) => { d.open = true; })")
    yield page
    context.close()


def _box(page, name):
    return page.locator(f"#{PLAN}-{name}-box" if name else f"#{PLAN}-box")


def _slide(slider, position):
    """Move a range input as a drag would (Playwright's `fill` does not take this one)."""
    slider.evaluate("(el, value) => { el.value = value; el.dispatchEvent(new Event('input', {bubbles: true})); }",
                    str(position))


def _value(locator):
    return float(locator.input_value())


def test_moving_a_slider_sets_its_box(page):
    floor = page.locator(f"#{PLAN}-floor")
    _slide(floor, 0)                                  # the first stop of the presets
    assert _value(_box(page, "floor")) == 2500.0
    _slide(floor, floor.get_attribute("max"))
    assert _value(_box(page, "floor")) == 4_000_000.0


def test_a_typed_value_is_kept_and_moves_the_slider_to_the_nearest_stop(page):
    box = _box(page, "floor")
    box.fill("9100")
    box.dispatch_event("change")
    assert _value(box) == 9100.0                                 # not snapped: any value inside the range is valid
    stops = page.locator(f"#{PLAN}-floor").evaluate("(el) => el.closest('[data-limit-slider]').dataset.stops")
    index = int(page.locator(f"#{PLAN}-floor").input_value())
    assert float(stops.strip("[]").split(",")[index]) == 9000.0


def test_a_typed_value_outside_the_range_is_clamped(page):
    box = _box(page, "floor")
    box.fill("10")
    box.dispatch_event("change")
    assert _value(box) == 2500.0


def test_a_limit_of_presets_snaps_to_one(page):
    box = _box(page, "")
    box.fill("13")
    box.dispatch_event("change")
    assert _value(box) in (12.0, 14.0)
    box.fill("14")
    box.dispatch_event("change")
    assert _value(box) == 14.0


def test_the_core_limits_follow_the_galaxy_wide_ones_until_one_is_set(page):
    core_floor, core_mass = _box(page, "core-floor"), _box(page, "core-mass")
    follow = page.locator(f"#{PLAN}-core input[name=core_follow]")
    assert _value(core_floor) == 9900.0 and _value(core_mass) == 15.4 and follow.is_checked()
    _box(page, "floor").fill("20000")
    _box(page, "floor").dispatch_event("change")
    assert _value(core_floor) == 22000.0
    _box(page, "").fill("10")
    _box(page, "").dispatch_event("change")
    assert _value(core_mass) == 11.0
    core_floor.fill("30000")                                     # set by hand: it stops following, the box un-ticks
    assert not follow.is_checked()
    _box(page, "floor").fill("5000")
    _box(page, "floor").dispatch_event("change")
    assert _value(core_floor) == 30000.0 and _value(core_mass) == 11.0


def test_ticking_follow_again_returns_the_core_to_the_galaxy_wide_limits(page):
    core_floor = _box(page, "core-floor")
    core_floor.fill("30000")
    follow = page.locator(f"#{PLAN}-core input[name=core_follow]")
    assert not follow.is_checked()
    follow.check()
    assert _value(core_floor) == 9900.0


def test_a_core_floor_under_the_galaxy_wide_one_is_warned_of(page):
    warning = page.locator(f"#{PLAN}-core [data-core-warning]")
    assert warning.is_hidden()
    _box(page, "core-floor").fill("5000")
    assert warning.is_visible()
    _box(page, "core-floor").fill("12000")
    assert warning.is_hidden()


def test_reset_returns_a_segment_to_its_defaults(page):
    halo, ceiling = _box(page, "halo-floor"), _box(page, "dense-ceiling")
    halo.fill("80")
    halo.dispatch_event("change")
    assert _value(halo) == 80.0
    page.locator(f"#{PLAN}-halo [data-limit-reset]").click()
    assert _value(halo) == 500.0
    assert "Reset" in page.locator(f"#{PLAN}-halo [data-limit-status]").inner_text()
    ceiling.fill("20000")
    ceiling.dispatch_event("change")
    page.locator(f"#{PLAN}-dense [data-limit-reset]").click()
    assert _value(ceiling) == 5000.0 and _value(_box(page, "dense-cap")) == 100.0


def test_the_ceiling_says_what_a_full_cell_costs(page):
    note = page.locator(f"#{PLAN}-dense [data-ceiling-cost]")
    assert "5,000 systems" in note.inner_text()
    _box(page, "dense-ceiling").fill("500000")
    assert "500,000 systems" in note.inner_text() and "hours" in note.inner_text()
