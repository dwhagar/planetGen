"""
UX.78: the reader's unit system (Automatic, Metric only, Customary) is applied
to every measurement the system page and a phenomenon page show, from the
`<span class="qty" data-q data-si>` the server writes, and the "Units" button
in the settings menu cycles through the three and remembers the choice.

Skipped without Playwright, Chromium or the MySQL test server, like
`test_web_browser_layout.py`.
"""

import re

import pytest

pytest.importorskip("playwright.sync_api")

from tests.test_web_a11y import (  # noqa: E402,F401 -- fixtures
    admin_token,
    base_url,
    browser,
    page_targets,
    sample_job,
    sample_job_tree,
    sample_params,
    site_app,
    site_db,
)

QUANTITIES = """() => [...document.querySelectorAll(".qty[data-q][data-si]")]
  .map((el) => ({ q: el.dataset.q, si: el.dataset.si, text: el.textContent }))"""

LARGE_UNIT = re.compile(r"\b(AU|ly|pc|kpc|Mpc|Gpc|mpc|cpc)\b")
METRIC_ONLY_UNIT = re.compile(r"\b(ft|mi|lb|mph|psi|°F|G)\b|\bmiles?\b")


def _open(browser, base_url, path, saved=None):
    context = browser.new_context(viewport={"width": 1100, "height": 900}, reduced_motion="reduce")
    if saved:
        context.add_init_script(f"window.localStorage.setItem('planetgen.units', '{saved}')")
    page = context.new_page()
    page.goto(base_url + path)
    page.wait_for_load_state("networkidle")
    return context, page


def _system_path(page_targets):
    path, _needs_admin = page_targets["web.system"]
    assert path, "the sample database has a system page"
    return path


def test_the_server_text_is_the_automatic_text(browser, base_url, page_targets):
    context, page = _open(browser, base_url, _system_path(page_targets))
    try:
        spans = page.evaluate(QUANTITIES)
        assert spans, "the system page shows measurements in .qty spans"
        assert page.evaluate("document.documentElement.dataset.units") == "auto"
        # Re-rendering in Automatic is the identity: the browser ladders agree with the server's text.
        page.evaluate("window.planetgenRenderQuantities(document)")
        assert page.evaluate(QUANTITIES) == spans
    finally:
        context.close()


def test_a_saved_preference_is_applied_on_load(browser, base_url, page_targets):
    path = _system_path(page_targets)
    context, page = _open(browser, base_url, path, "metric")
    try:
        assert page.evaluate("document.documentElement.dataset.units") == "metric"
        for span in page.evaluate(QUANTITIES):
            assert not LARGE_UNIT.search(span["text"]), span
            assert not METRIC_ONLY_UNIT.search(span["text"]), span
    finally:
        context.close()
    context, page = _open(browser, base_url, path, "customary")
    try:
        spans = page.evaluate(QUANTITIES)
        assert not any(re.search(r"\bkm\b|\bkg\b|\bPa\b", span["text"]) for span in spans
                       if span["q"] in ("length", "mass", "pressure")), spans
    finally:
        context.close()


def test_the_units_button_cycles_and_remembers(browser, base_url, page_targets):
    path = _system_path(page_targets)
    context, page = _open(browser, base_url, path)
    try:
        # The button lives in the settings menu: open it, press it.
        page.click("sl-dropdown.site-gear sl-icon-button")
        button = page.locator("[data-units-toggle]")
        assert button.inner_text() == "Units: Automatic"
        seen = []
        for expected in ("metric", "customary", "auto"):
            button.click()
            assert page.evaluate("document.documentElement.dataset.units") == expected
            seen.append(page.evaluate("window.localStorage.getItem('planetgen.units')"))
        assert seen == ["metric", "customary", None]
        button.click()  # metric again, then reload
        page.reload()
        page.wait_for_load_state("networkidle")
        assert page.evaluate("document.documentElement.dataset.units") == "metric"
    finally:
        context.close()
