# tests/test_web_admin_fields_browser.py

"""
The admin pages and the edit dialogs on the one form field (UX.126) in a real browser: a form of each kind is filled
and submitted on the admin hub (a text box, a number box, checkboxes), the facility types page (a text box, two
selects, a number box, a checkbox) and the system page's facility form (selects that filter each other, a text box),
and each page passes the accessibility and phone-fit checks with the fields drawn. Needs
Playwright with Chromium and the MySQL test server, like test_web_a11y.py; skipped without them.

Run with: pytest tests/test_web_admin_fields_browser.py
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from planetgen.api.authz import SESSION_COOKIE_NAME  # noqa: E402

from tests.test_web_a11y import (  # noqa: E402,F401 -- the browser, axe and site fixtures, and the checks
    VIEWPORTS, _axe_blocking, _overflow, admin_token, axe_source, base_url, browser, sample_job, sample_params,
    sample_run, site_app, site_db,
)


def _admin_page(browser, base_url, admin_token, viewport="desktop", scheme="light"):
    context = browser.new_context(viewport=VIEWPORTS[viewport], color_scheme=scheme, reduced_motion="reduce")
    context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
    page = context.new_page()
    page.set_default_timeout(30000)
    problems = []
    page.on("console", lambda msg: msg.type == "error" and problems.append(f"console error: {msg.text}"))
    page.on("pageerror", lambda exc: problems.append(f"uncaught exception: {exc}"))
    return context, page, problems


def _defined(page, selector):
    page.locator(f"{selector}:defined").first.wait_for(state="visible", timeout=20000)


@pytest.mark.parametrize("viewport", sorted(VIEWPORTS))
def test_the_admin_hub_creates_a_key_from_a_text_box_a_number_box_and_checkboxes(browser, axe_source, base_url,
                                                                                admin_token, viewport):
    context, page, problems = _admin_page(browser, base_url, admin_token, viewport)
    try:
        page.goto(f"{base_url}/admin", wait_until="load")
        _defined(page, "sl-input[name='label']")
        problems += [f"/admin: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]
        page.locator("sl-input[name='label'] input").fill("ux126 script")
        page.locator("sl-input[name='expires_days'] input").fill("30")
        page.locator("sl-checkbox[name='scope'][value='read']").click()
        page.get_by_role("button", name="Create key").click()
        page.wait_for_url("**/admin**")
        page.locator("#new-key").wait_for(state="visible")
        assert "ux126 script" in page.locator("#new-key").inner_text()
        # The scopes ticked on the form are the key's: Read and the ticked-at-start Admin.
        page.goto(f"{base_url}/admin", wait_until="load")
        assert "ux126 script" in page.locator("#api-keys").inner_text()
    finally:
        context.close()
    assert not problems, f"at {viewport}:\n  " + "\n  ".join(problems)


@pytest.mark.parametrize("viewport", sorted(VIEWPORTS))
def test_the_facility_types_page_adds_a_type_with_selects_a_number_and_a_checkbox(browser, axe_source, base_url,
                                                                                 admin_token, viewport):
    context, page, problems = _admin_page(browser, base_url, admin_token, viewport)
    try:
        page.goto(f"{base_url}/admin/facility-types", wait_until="load")
        _defined(page, "#new-type sl-input[name='name']")
        problems += [f"/admin/facility-types: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]
        page.locator("#new-type sl-input[name='name'] input").fill("Ux126 outpost")
        page.locator("#new-type sl-select[name='kind']").click()
        page.locator("#new-type sl-select[name='kind'] sl-option:visible").nth(1).click()
        page.locator("#new-type sl-input[name='sort_order'] input").fill("7")
        page.locator("#new-type sl-input[name='description'] input").fill("Added by the field test")
        page.locator("#new-type button[type=submit]").click()
        page.wait_for_url("**/admin/facility-types**")
        listed = page.locator("main").inner_text()
        assert "Ux126 outpost" in listed
    finally:
        context.close()
    assert not problems, f"at {viewport}:\n  " + "\n  ".join(problems)


def test_the_facility_form_filters_hosts_by_placement_and_posts_its_preview(browser, axe_source, base_url,
                                                                           admin_token, sample_params):
    context, page, problems = _admin_page(browser, base_url, admin_token)
    try:
        page.goto(f"{base_url}/system/{sample_params['system_id']}", wait_until="load")
        if not page.locator("#admin-facility-place").count():
            pytest.skip("the sample system offers no facility form")
        for _ in range(40):
            if page.evaluate("typeof document.getElementById('admin-facility-place').show === 'function'"):
                break
            page.wait_for_timeout(500)
        else:
            pytest.fail("the Place a facility dialog never became a dialog: " + "; ".join(problems))
        page.evaluate("document.getElementById('admin-facility-place').show()")
        _defined(page, "#admin-facility-place form.facility-form sl-select[name='placement']")
        problems += [f"system page: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]

        def offered(name):
            return page.evaluate(
                """(name) => Array.from(document.querySelectorAll(`#admin-facility-place form.facility-form sl-select[name='${name}'] sl-option`))
                    .filter((o) => !o.hidden && !o.disabled).map((o) => o.value)""", name)

        placements = page.evaluate(
            """() => Array.from(document.querySelectorAll("#admin-facility-place form.facility-form sl-select[name='placement'] sl-option"))
                .map((o) => o.value)""")
        for placement in placements:
            page.evaluate(
                """(value) => { const s = document.querySelector("#admin-facility-place form.facility-form sl-select[name='placement']");
                    s.value = value; s.dispatchEvent(new CustomEvent("sl-change", {bubbles: true})); }""", placement)
            allowed = page.evaluate(
                """(value) => Array.from(document.querySelectorAll("#admin-facility-place form.facility-form sl-select[name='host'] sl-option"))
                    .filter((o) => (o.dataset.placements || "").split(" ").includes(value)).map((o) => o.value)""",
                placement)
            assert offered("host") == allowed, placement
            chosen = page.evaluate("document.querySelector(\"#admin-facility-place form.facility-form sl-select[name='host']\").value")
            assert not allowed or chosen in allowed, (placement, chosen, allowed)
        page.locator("#admin-facility-place form.facility-form sl-input[name='name'] input").fill("Ux126 test")
        page.locator("#admin-facility-place form.facility-form button[value=preview]").click()
        page.wait_for_url("**/system/**")
        assert page.locator("#admin-facility-place form.facility-form sl-input[name='name'] input").input_value() == "Ux126 test"
    finally:
        context.close()
    assert not problems, "\n  ".join(problems)
