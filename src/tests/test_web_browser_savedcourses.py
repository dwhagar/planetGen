"""
NAV.4, NAV.18: saved courses on the NAV page, in a real browser. The page is the real Flask NAV page (with the faked API of
`test_web_sector_nav`), served with the site's own static files. Save course asks for a name and keeps the course in this browser; the
Courses panel lists it with its Direct and Route forms, renames it, deletes it, and says what changed when the course is planned again.
Skipped without Playwright or Chromium.
"""

import json
import mimetypes
import os

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.test_web_sector_nav import _adjusted, app, client, fake  # noqa: E402,F401 -- fixtures and the faked course

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_STATIC = os.path.join(_SRC_DIR, "html", "static")
_ORIGIN = "http://planetgen.test"
_COURSE = "/nav?from=system:1001&to=system:1002"


@pytest.fixture(scope="module")
def browser():
    try:
        playwright = sync_api.sync_playwright().start()
    except Exception as exc:  # noqa: BLE001 -- any failure here means "no browser"
        pytest.skip(f"Playwright could not start ({exc})")
    try:
        chromium = playwright.chromium.launch()
    except Exception as exc:  # noqa: BLE001
        playwright.stop()
        pytest.skip(f"No Chromium for Playwright ({str(exc).splitlines()[0]})")
    try:
        yield chromium
    finally:
        chromium.close()
        playwright.stop()


@pytest.fixture
def page(browser, client, fake):
    fake.nav_galaxy_scope = True
    fake.nav_adjusted = _adjusted()
    context = browser.new_context()
    tab = context.new_page()
    tab.dialogs = []

    def serve(route):
        url = route.request.url[len(_ORIGIN):]
        path = url.split("?")[0]
        if path.startswith("/static/"):
            local = os.path.join(_STATIC, path[len("/static/"):])
            if os.path.isfile(local):
                with open(local, "rb") as handle:
                    route.fulfill(status=200, body=handle.read(), content_type=mimetypes.guess_type(local)[0] or "application/octet-stream")
                return
            route.fulfill(status=404, body="")
            return
        answer = client.get(url)
        route.fulfill(status=answer.status_code, body=answer.get_data(), content_type=answer.content_type)

    context.route(_ORIGIN + "/**", serve)
    yield tab
    context.close()


def _answer_prompt(tab, text):
    """Answer the next Name this course prompt with `text` (None cancels); alerts are kept in `tab.dialogs`."""
    def handle(dialog):
        tab.dialogs.append((dialog.type, dialog.message))
        if dialog.type == "prompt":
            dialog.accept(text) if text is not None else dialog.dismiss()
        else:
            dialog.accept()
    tab.once("dialog", handle)


def _stored(tab):
    return tab.evaluate("""() => {
        const key = Object.keys(localStorage).find(k => k.startsWith('planetgen.courses.'));
        return key ? JSON.parse(localStorage.getItem(key)) : null;
    }""")


def test_save_course_names_it_and_keeps_it_in_this_browser(page):
    page.goto(_ORIGIN + _COURSE, wait_until="load")
    button = page.locator("[data-course-save-button]")
    button.wait_for(state="visible")
    assert button.inner_text() == "Save course"
    assert "No saved courses yet" in page.locator("[data-courses-panel]").inner_text()
    _answer_prompt(page, "Home run")
    button.click()
    page.wait_for_function("document.querySelector('[data-course-save-button]').textContent === 'Update saved course'")
    stored = _stored(page)
    assert [entry["name"] for entry in stored] == ["Home run"]
    assert stored[0]["from"] == "system:1001" and stored[0]["to"] == "system:1002"
    assert "Saved as “Home run”." in page.locator("[data-course-status]").inner_text()
    assert "Nothing has changed since you saved it." in page.locator("[data-course-changes]").inner_text()
    panel = page.locator("[data-courses-panel]")
    assert "Home run" in panel.inner_text() and "Start 1001" in panel.inner_text()
    link = panel.locator("a.bookmarks-open")
    assert link.get_attribute("href") == "/nav?from=system:1001&to=system:1002"
    assert panel.get_by_role("link", name="Show Home run on the Galaxy Map").get_attribute("href") == "/galaxy?course=system:1001,system:1002"


def test_canceling_the_name_saves_nothing(page):
    page.goto(_ORIGIN + _COURSE, wait_until="load")
    page.locator("[data-course-save-button]").wait_for(state="visible")
    _answer_prompt(page, None)
    page.locator("[data-course-save-button]").click()
    page.wait_for_timeout(200)
    assert _stored(page) in (None, [])
    assert page.locator("[data-course-save-button]").inner_text() == "Save course"


def test_the_courses_list_shows_either_form_renames_and_deletes(page):
    page.goto(_ORIGIN + _COURSE, wait_until="load")
    page.locator("[data-course-save-button]").wait_for(state="visible")
    _answer_prompt(page, "Home run")
    page.locator("[data-course-save-button]").click()
    panel = page.locator("[data-courses-panel]")
    panel.locator("li.bookmarks-item").wait_for()
    direct = panel.locator(".courses-form-text").inner_text()
    assert "light-years" in direct or "ly" in direct
    assert panel.get_by_role("button", name="Show the direct form of Home run").get_attribute("aria-pressed") == "true"
    panel.get_by_role("button", name="Show the route form of Home run").click()
    route = panel.locator(".courses-form-text").inner_text()
    assert route != direct and "between the ends" in route
    assert panel.get_by_role("button", name="Show the route form of Home run").get_attribute("aria-pressed") == "true"
    assert panel.get_by_role("button", name="Show the direct form of Home run").get_attribute("aria-pressed") == "false"

    panel.get_by_role("button", name="Rename saved course Home run").click()
    box = panel.locator("form.bookmarks-rename input")
    box.fill("Long way home")
    box.press("Enter")
    panel.locator("a.bookmarks-open", has_text="Long way home").wait_for()
    assert [entry["name"] for entry in _stored(page)] == ["Long way home"]

    panel.get_by_role("button", name="Rename saved course Long way home").click()
    panel.locator("form.bookmarks-rename input").press("Escape")
    assert [entry["name"] for entry in _stored(page)] == ["Long way home"]

    panel.get_by_role("button", name="Delete saved course Long way home").click()
    page.wait_for_function("document.querySelector('[data-courses-panel]').textContent.includes('No saved courses yet')")
    assert _stored(page) == []
    assert page.locator("[data-course-save-button]").inner_text() == "Save course"


def test_a_saved_course_says_what_changed_when_it_is_planned_again(page, fake):
    page.goto(_ORIGIN + _COURSE, wait_until="load")
    page.locator("[data-course-save-button]").wait_for(state="visible")
    _answer_prompt(page, "Home run")
    page.locator("[data-course-save-button]").click()
    page.wait_for_function("document.querySelector('[data-course-save-button]').textContent === 'Update saved course'")
    # The route changed since it was saved: the browser keeps the old record, so alter it in storage to look older.
    page.evaluate("""() => {
        const key = Object.keys(localStorage).find(k => k.startsWith('planetgen.courses.'));
        const list = JSON.parse(localStorage.getItem(key));
        list[0].direct.distance_ly = list[0].direct.distance_ly + 1;
        localStorage.setItem(key, JSON.stringify(list));
    }""")
    page.reload(wait_until="load")
    page.locator("[data-course-changes]:not([hidden])").wait_for()
    assert "Nothing has changed" not in page.locator("[data-course-changes]").inner_text()
    page.locator("[data-course-save-button]").click()
    page.wait_for_function("document.querySelector('[data-course-changes]').textContent.includes('Nothing has changed')")
    assert [entry["name"] for entry in _stored(page)] == ["Home run"]
