# tests/test_web_generate_flow_browser.py

"""
The Generate screens (ADM.107, ADM.108) in a real browser: the Start screen, every task's page and the Review step each
pass the accessibility and phone-fit checks in both themes, a task's Review shows the command line and Back returns the
form with the draft filled in. Needs Playwright with Chromium and the MySQL test server, like test_web_a11y.py; skipped
without them.

Run with: pytest tests/test_web_generate_flow_browser.py
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_a11y import (  # noqa: E402,F401 -- the browser and axe fixtures, and the checks
    SCHEMES, VIEWPORTS, _axe_blocking, _overflow, axe_source, browser,
)
from tests.test_web_settings_browser import _page, config, server  # noqa: E402,F401

TASK_KEYS = ("add-sectors", "new-galaxy", "scatters", "check", "wipe")


@pytest.mark.parametrize("scheme", SCHEMES)
@pytest.mark.parametrize("viewport", sorted(VIEWPORTS))
def test_every_screen_is_accessible_and_fits(browser, axe_source, server, config, viewport, scheme):
    context, page, problems, base_url = _page(browser, server, viewport, scheme)
    try:
        paths = ["/admin/generate", "/admin/generate/history", *(f"/admin/generate/task/{key}" for key in TASK_KEYS)]
        for path in paths:
            assert page.goto(base_url + path, wait_until="load").status == 200, path
            page.evaluate("() => document.querySelectorAll('details').forEach((d) => { d.open = true; })")
            problems += [f"{path}: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]
        page.goto(base_url + "/admin/generate/task/check", wait_until="load")
        page.click("button:has-text('Review')")
        page.wait_for_selector("#review")
        problems += [f"review: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]
    finally:
        context.close()
    assert not problems, f"at {viewport}/{scheme}:\n  " + "\n  ".join(problems)


def test_review_shows_the_command_and_back_keeps_the_draft(browser, server, config):
    context, page, problems, base_url = _page(browser, server)
    try:
        page.goto(base_url + "/admin/generate/task/scatters", wait_until="load")
        page.check("input[name='action'][value='bright_band']")
        page.fill("input[name='down_to']", "250")
        page.click("button:has-text('Review')")
        page.wait_for_selector("#review")
        assert "250" in page.inner_text("#review pre.command-line >> nth=-1") or "250" in page.inner_text("#review")
        page.click("#review button:has-text('Back')")
        page.wait_for_selector("#task-scatters")
        assert page.input_value("input[name='down_to']") == "250"
        assert page.is_checked("input[name='action'][value='bright_band']")
        page.fill("input[name='down_to']", "")
        page.click("button:has-text('Review')")
        page.wait_for_selector("[role=alert]")
        assert "task-scatters" in page.inner_html("main")
    finally:
        context.close()
    # the empty box is refused on purpose (a 400 the browser logs as a failed resource)
    problems = [p for p in problems if "400 (BAD REQUEST)" not in p]
    assert not problems, "\n".join(problems)
