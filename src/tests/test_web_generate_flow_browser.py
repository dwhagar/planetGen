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
from planetgen.web import generate_page  # noqa: E402
from tests.test_web_settings_browser import _page, config, server  # noqa: E402,F401

ESTIMATE = {"sectors": 9, "systems": 120, "stars": 156, "bytes": 8_400_000, "seconds": 75.0, "workers": 2,
            "measured": True, "refused": False, "refusal": None,
            "disk": {"path": "/var/lib/mysql/", "total_bytes": 10 ** 12, "free_bytes": 4 * 10 ** 11},
            "summary": "About 8.4 MB and 1 m 15 s for 9 sectors.", "what": "shell 40"}
"""dict: A made-up size and time, so a galaxy job's Review step does not run the real estimate."""

TASK_KEYS = ("add-sectors", "new-galaxy", "scatters", "check", "wipe")


@pytest.mark.parametrize("scheme", SCHEMES)
@pytest.mark.parametrize("viewport", sorted(VIEWPORTS))
def test_every_screen_is_accessible_and_fits(browser, axe_source, server, config, viewport, scheme, monkeypatch):
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: dict(ESTIMATE))
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


def test_add_sectors_shows_one_panel_at_a_time_and_each_choice_reaches_its_review(browser, server, config, monkeypatch):
    """ADM.109: the two choices and the four place tabs show only the panel chosen; each tab's boxes build the
    command the old form built."""
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: dict(ESTIMATE))
    context, page, problems, base_url = _page(browser, server)
    try:
        page.goto(base_url + "/admin/generate/task/add-sectors", wait_until="load")
        visible = "[data-place-panel]:not([hidden])"
        assert page.locator(visible).count() == 1
        assert page.is_hidden("[data-what-panel='grid']")
        for tab, field, value, expected in (
            ("position", "center_x_pc", "10", "--radius-pc 20.0"),
            ("address", "center_shell", "40", "--shell 40"),
        ):
            page.click(f"label.place-tab:has(input[name='place'][value='{tab}'])")
            assert page.locator(visible).count() == 1
            assert page.get_attribute(visible, "data-place-panel") == tab
            page.fill(f"[data-place-panel='{tab}'] input[name='{field}']", value)
            page.fill("input[name='center_radius_pc']", "20")
            if tab == "position":
                page.fill("input[name='center_y_pc']", "0")
            else:
                page.fill("input[name='center_wedge']", "7")
            page.click("button:has-text('Review')")
            page.wait_for_selector("#review, [role=alert]")
            assert page.locator("[role=alert]").count() == 0, page.inner_text("[role=alert]")
            assert expected in page.text_content("#review")
            page.click("#review button:has-text('Back')")
            page.wait_for_selector("#task-add-sectors")
            assert page.is_checked(f"input[name='place'][value='{tab}']")
        page.check("input[name='what'][value='grid']")
        assert page.is_visible("[data-what-panel='grid']") and page.is_hidden("[data-what-panel='place']")
        page.click("label.place-tab:has(input[name='mode'][value='span'])")
        page.fill("input[name='span_shells']", "10:20")
        page.click("button:has-text('Review')")
        page.wait_for_selector("#review")
        assert "--shells=10:20" in page.text_content("#review")
    finally:
        context.close()
    problems = [p for p in problems if "400 (BAD REQUEST)" not in p]
    assert not problems, "\n".join(problems)


def test_a_typed_command_in_the_palette_opens_the_review_and_runs_nothing(browser, server, config, monkeypatch):
    """ADM.121: Ctrl+K, a typed line, Enter: the Review step for it, filled in, with Run still to press."""
    monkeypatch.setattr(generate_page, "run_estimate", lambda argv, env: dict(ESTIMATE))
    context, page, problems, base_url = _page(browser, server)
    try:
        page.goto(base_url + "/admin/generate", wait_until="load")
        page.keyboard.press("Control+k")
        page.fill("dialog.admin-palette input[type=search]", "add sectors shell 40 limit 200")
        page.wait_for_selector("dialog.admin-palette li[role=option] >> text=Review: Add sectors")
        page.keyboard.press("Enter")
        page.wait_for_selector("#review")
        assert "--shell 40" in page.text_content("#review") and "--limit 200" in page.text_content("#review")
        assert page.is_visible("#review button:has-text('Run')")
    finally:
        context.close()
    assert not problems, "\n".join(problems)


def test_the_palette_finds_a_page_and_the_keys_open_help_and_close_it(browser, server, config):
    """ADM.120, ADM.122: Ctrl+K and a word lands on that page; `?` lists the keys and Esc closes them."""
    context, page, problems, base_url = _page(browser, server)
    try:
        page.goto(base_url + "/admin/generate", wait_until="load")
        page.keyboard.press("Control+k")
        page.fill("dialog.admin-palette input[type=search]", "settings")
        page.wait_for_selector("dialog.admin-palette li[role=option]")
        page.keyboard.press("Enter")
        page.wait_for_url("**/admin/settings*")
        page.keyboard.press("?")
        page.wait_for_selector("dialog.admin-keys[open]")
        page.keyboard.press("Escape")
        page.wait_for_selector("dialog.admin-keys[open]", state="detached")
    finally:
        context.close()
    assert not problems, "\n".join(problems)


def test_a_review_can_be_saved_as_a_recipe_run_renamed_and_deleted(browser, axe_source, server, config):
    """ADM.125: Save as recipe on the Review step keeps the settings under a name; the Start screen lists it (and passes
    the accessibility and phone-fit checks with one listed); opening it shows the same Review step, which runs nothing;
    it can be renamed and deleted."""
    problems = []
    for viewport, scheme in (("desktop", "light"), ("phone", "dark")):
        context, page, found, base_url = _page(browser, server, viewport, scheme)
        try:
            page.goto(base_url + "/admin/generate/task/scatters", wait_until="load")
            page.check("input[name='action'][value='bright_band']")
            page.fill("input[name='down_to']", "250")
            page.click("button:has-text('Review')")
            page.wait_for_selector("#review")
            page.fill("#recipe-name", "Dim stars")
            page.click("#review button:has-text('Save as recipe')")
            page.wait_for_selector("#recipes")
            assert "Saved the recipe" in page.inner_text("#recipes")
            assert page.locator("#recipes a", has_text="Dim stars").count() == 1
            found += [f"start: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]

            page.click("#recipes a:has-text('Dim stars')")
            page.wait_for_selector("#review")
            assert "250" in page.inner_text("#review")
            assert page.locator("#review button:has-text('Save as recipe')").count() == 1

            page.fill("#recipe-name", "Dim stars")
            page.click("#review button:has-text('Save as recipe')")
            page.wait_for_selector("#recipe-error")
            assert "already have a recipe" in page.inner_text("#recipe-error")
            found += [f"review with error: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]

            page.goto(base_url + "/admin/generate", wait_until="load")
            page.fill("#recipes input[name='recipe_name']", "Dimmer stars")
            page.click("#recipes button:has-text('Rename')")
            page.wait_for_selector("#recipes")
            assert page.locator("#recipes a", has_text="Dimmer stars").count() == 1
            page.click("#recipes button:has-text('Delete')")
            page.wait_for_selector("#recipes")
            assert "Deleted the recipe" in page.inner_text("#recipes")
            assert page.locator("#recipes a").count() == 0
        finally:
            context.close()
        problems += [f"{viewport}/{scheme}: {p}" for p in found if "400 (BAD REQUEST)" not in p]
    assert not problems, "\n".join(problems)


def _sized(argv, env):
    """An estimate whose sector count is the form's `--limit`, so a change of the box is a change of the line."""
    limit = int(argv[argv.index("--limit") + 1]) if "--limit" in argv else 9
    return {**ESTIMATE, "sectors": limit, "systems": limit * 10}


def test_the_size_line_follows_the_form_and_the_palette_command(browser, axe_source, server, config, monkeypatch):
    """ADM.126: a moment after a box changes, the line under the form shows the sectors, systems and stars the run would
    make, or why it cannot be sized yet; a typed command in the palette shows its size too. Nothing is started."""
    monkeypatch.setattr(generate_page, "run_estimate", _sized)
    problems = []
    for viewport, scheme in (("desktop", "light"), ("phone", "dark")):
        context, page, found, base_url = _page(browser, server, viewport, scheme)
        try:
            page.goto(base_url + "/admin/generate/task/add-sectors", wait_until="load")
            page.check("input[name='what'][value='grid']")
            page.fill("input[name='shell_index']", "40")
            page.fill("input[name='shell_limit']", "200")
            page.wait_for_selector("#size-report:has-text('About 200 sectors')", timeout=1500)
            assert "about 2,000 star systems" in page.inner_text("#size-report")
            found += [f"size line: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]
            page.fill("input[name='shell_limit']", "300")
            page.wait_for_selector("#size-report:has-text('About 300 sectors')", timeout=1500)
            page.fill("input[name='shell_index']", "")
            page.wait_for_selector("#size-report:has-text('Not sized yet: ')", timeout=1500)
            assert page.locator("[role=alert]").count() == 0
            page.goto(base_url + "/admin/generate", wait_until="load")
            page.keyboard.press("Control+k")
            page.fill("dialog.admin-palette input[type=search]", "add sectors shell 40 limit 250")
            page.wait_for_selector("dialog.admin-palette li[role=option] >> text=About 250 sectors", timeout=3000)
        finally:
            context.close()
        problems += [f"{viewport}/{scheme}: {p}" for p in found if "400 (BAD REQUEST)" not in p]
    assert not problems, "\n".join(problems)
