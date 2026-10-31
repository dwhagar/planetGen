"""
UX.146: the top menu bar in a real browser.

Lists is under the gear, the search is one icon-only link to the Search page (UX.160), Nearby is under Nav,
which opens Plot Course then Nearby Search, and Help (UX.180), left of the search icon, opens Classes then Habitability. Loaded in headless Chromium at five window sizes in both themes: the header
shows no Lists or Nearby link, every control lies inside the window with nothing wider than it, the search icon has
an accessible name and no text, the gear lists Lists, and the Nav menu (or, where the sections are folded into the
Menu, the Menu's Nav group) opens its two entries in order and works from the keyboard. Menus open under their button
and inside the window (UX.145).
Skipped without Playwright, Chromium or the MySQL test server, like `test_web_a11y.py`.
"""

import pytest

pytest.importorskip("playwright.sync_api")

from tests.test_web_a11y import (  # noqa: E402,F401 -- fixtures
    base_url,
    browser,
    site_app,
    site_db,
)

SIZES = ((1280, 720), (1024, 768), (820, 1024), (600, 900), (390, 844))
"""tuple: The window sizes checked, as (width, height): a laptop, a small laptop, a tablet, a wide phone and a phone."""

SLACK_PX = 1.0
"""float: How far a box may lie outside the window and still count as inside it (rounding)."""

GAP_PX = 16.0
"""float: How far a menu may be from its button."""

HEADER_STATE = """() => {
    const header = document.querySelector("header.site-header");
    const inner = header.querySelector(".site-header-inner");
    const visible = (el) => { const r = el.getBoundingClientRect(); return r.width > 0 && r.height > 0; };
    const links = [...header.querySelectorAll("a")].filter(visible).map((a) => a.textContent.trim());
    return {
        overflow: inner.scrollWidth - inner.clientWidth,
        pageOverflow: document.documentElement.scrollWidth - document.documentElement.clientWidth,
        links,
        wide: visible(header.querySelector(".site-wide")),
    };
}"""


def _inside(box, width, height):
    return [side for side, over in (("left", -box["x"]), ("top", -box["y"]), ("right", box["x"] + box["width"] - width),
                                    ("bottom", box["y"] + box["height"] - height)) if over > SLACK_PX]


def _attached(panel, trigger):
    below = abs(panel["y"] - (trigger["y"] + trigger["height"])) <= GAP_PX
    above = abs((panel["y"] + panel["height"]) - trigger["y"]) <= GAP_PX
    return below or above


def _load(page, url):
    """Opens `url` and waits until the header's dropdowns are Shoelace elements (a key pressed before is lost)."""
    page.goto(url, wait_until="load")
    page.locator("header sl-dropdown.site-gear:defined").wait_for()
    page.locator("header sl-icon-button.site-search-link:defined").wait_for()


def _close(page, dropdown):
    """Escape, and wait until the dropdown `dropdown` (a selector) has finished closing."""
    page.keyboard.press("Escape")
    page.wait_for_function("(sel) => !document.querySelector(sel).open", arg=dropdown)
    page.wait_for_timeout(100)


def _open(page, trigger, panel):
    box = trigger.bounding_box()
    trigger.click()
    panel.wait_for(state="visible")
    return box


@pytest.mark.parametrize("scheme", ["light", "dark"])
@pytest.mark.parametrize("size", SIZES, ids=lambda size: f"{size[0]}x{size[1]}")
def test_the_top_bar_menus(browser, base_url, size, scheme):
    width, height = size
    context = browser.new_context(viewport={"width": width, "height": height}, color_scheme=scheme,
                                  reduced_motion="reduce")
    problems = []
    try:
        page = context.new_page()
        _load(page, base_url + "/")
        state = page.evaluate(HEADER_STATE)

        # The bar: no Lists or Nearby link, and nothing wider than the window.
        for label in ("Lists", "Nearby", "Nearby Search"):
            if label in state["links"]:
                problems.append(f"the bar shows a {label} link: {state['links']}")
        if state["overflow"] > 0 or state["pageOverflow"] > 0:
            problems.append(f"the bar is wider than the window: {state}")

        # The search icon (UX.160): icon only, with a name, and a plain link to the Search page; no box drops down.
        search = page.locator("header sl-icon-button.site-search-link")
        if search.inner_text().strip() or search.get_attribute("label") != "Search":
            problems.append("the search icon is not an icon with the name Search")
        if page.locator("header").get_by_role("link", name="Search", exact=True).count() != 1:
            problems.append("no single link named Search in the bar")
        if page.locator("header input[type=search]").count():
            problems.append("the bar holds a search box")
        box = search.bounding_box()
        if _inside(box, width, height):
            problems.append(f"the search icon {box} lies outside the window")
        search.click()
        page.wait_for_url("**/search")
        _load(page, base_url + "/")

        # The gear: Lists is its first entry.
        gear = page.locator("header sl-dropdown.site-gear > [slot=trigger]")
        box = _open(page, gear, page.locator("header .site-gear-panel a", has_text="Lists"))
        lists = page.locator("header .site-gear-panel a", has_text="Lists").bounding_box()
        panel = page.locator("header sl-dropdown.site-gear .dropdown__panel").bounding_box()
        if _inside(panel, width, height) or _inside(lists, width, height) or not _attached(panel, box):
            problems.append(f"the gear panel {panel} or its Lists link {lists} is not at its button {box} inside the window")
        first = page.locator("header .site-gear-panel nav a").first.inner_text().strip()
        if first != "Lists":
            problems.append(f"the gear's first entry is {first!r}, not Lists")
        _close(page, "header sl-dropdown.site-gear")

        # Help (UX.180): Classes and Habitability, beside the search button (a book icon on the narrowest phones).
        help_button = page.locator("header sl-dropdown.site-help-menu > [slot=trigger]")
        items = page.locator("header sl-dropdown.site-help-menu sl-menu-item")
        if help_button.get_attribute("aria-haspopup") != "menu" or help_button.inner_text().strip() != "Help":
            if width >= 352:
                problems.append("the Help button does not read Help")
        help_box = help_button.bounding_box()
        search_box = search.bounding_box()
        if not help_box or help_box["x"] + help_box["width"] > search_box["x"] + SLACK_PX:
            problems.append(f"the Help button {help_box} is not left of the search button {search_box}")
        box = _open(page, help_button, items.first)
        labels = [item.inner_text().strip() for item in items.all()]
        if labels != ["Classes", "Habitability", "Properties", "Life"]:
            problems.append(f"the Help menu holds {labels}")
        panel = page.locator("header sl-dropdown.site-help-menu .dropdown__panel").bounding_box()
        if _inside(panel, width, height) or not _attached(panel, box):
            problems.append(f"the Help panel {panel} is not at its button {box} inside the window")
        items.nth(1).click()
        page.wait_for_url("**/help/habitability")
        _load(page, base_url + "/")

        # Nav.
        if state["wide"]:
            trigger = page.locator("header .site-nav-trigger")
            items = page.locator("header sl-dropdown.site-nav-menu sl-menu-item")
            box = _open(page, trigger, items.first)
            labels = [item.inner_text().strip() for item in items.all()]
            if labels != ["Plot Course", "Nearby Search"]:
                problems.append(f"the Nav menu holds {labels}")
            panel = page.locator("header sl-dropdown.site-nav-menu .dropdown__panel").bounding_box()
            if _inside(panel, width, height) or not _attached(panel, box):
                problems.append(f"the Nav panel {panel} is not at its button {box} inside the window")
            _close(page, "header sl-dropdown.site-nav-menu")
            if not page.evaluate("document.activeElement && document.activeElement.classList.contains('site-nav-trigger')"):
                problems.append("focus did not return to the Nav button after Escape")

            # From the keyboard alone: Enter opens the menu, the arrow keys move through its two entries, Enter goes there.
            trigger.focus()
            page.keyboard.press("Enter")
            items.first.wait_for(state="visible")
            page.keyboard.press("ArrowDown")
            if page.evaluate("document.activeElement.textContent.trim()") != "Plot Course":
                problems.append("the first arrow key does not reach Plot Course")
            page.keyboard.press("Enter")
            page.wait_for_url("**/nav")
            _load(page, base_url + "/")
            trigger = page.locator("header .site-nav-trigger")
            trigger.focus()
            page.keyboard.press("ArrowDown")
            items.first.wait_for(state="visible")
            page.keyboard.press("ArrowDown")
            if page.evaluate("document.activeElement.textContent.trim()") != "Nearby Search":
                problems.append("the second arrow key does not reach Nearby Search")
            page.keyboard.press("Enter")
            page.wait_for_url("**/nearby")
        else:
            menu = page.locator("header sl-dropdown.site-menu > [slot=trigger]")
            group = page.locator("header .site-menu-panel .site-nav-group")
            box = _open(page, menu, group)
            if group.locator(".site-nav-label").inner_text().strip() != "Nav":
                problems.append("the Menu has no Nav group")
            labels = [a.inner_text().strip() for a in group.locator("a").all()]
            if labels != ["Plot Course", "Nearby Search"]:
                problems.append(f"the Menu's Nav group holds {labels}")
            for link in group.locator("a").all():
                link_box = link.bounding_box()
                if _inside(link_box, width, height):
                    problems.append(f"{link.inner_text()!r} lies outside the window: {link_box}")
            panel = page.locator("header sl-dropdown.site-menu .dropdown__panel").bounding_box()
            if _inside(panel, width, height) or not _attached(panel, box):
                problems.append(f"the Menu panel {panel} is not at its button {box} inside the window")
            # Tab reaches the second entry right after the first.
            group.locator("a").first.focus()
            page.keyboard.press("Tab")
            if page.evaluate("document.activeElement.textContent.trim()") != "Nearby Search":
                problems.append("Tab from Plot Course does not reach Nearby Search")
            page.keyboard.press("Enter")
            page.wait_for_url("**/nearby")
    finally:
        context.close()
    assert not problems, f"{width}x{height} {scheme}:\n  " + "\n  ".join(problems)


def test_the_bar_fits_at_every_width(browser, base_url):
    """The sections fold into the Menu before they would overflow: nothing in the bar is wider than the window from a
    phone to a wide screen, in steps."""
    context = browser.new_context(viewport={"width": 1400, "height": 800}, reduced_motion="reduce")
    problems = []
    try:
        page = context.new_page()
        _load(page, base_url + "/")
        for width in range(320, 1401, 20):
            page.set_viewport_size({"width": width, "height": 800})
            page.wait_for_timeout(30)
            state = page.evaluate(HEADER_STATE)
            if state["overflow"] > 0 or state["pageOverflow"] > 0:
                problems.append(f"{width} px: {state}")
            names = page.evaluate("""() => [...document.querySelectorAll("header .site-name, header .site-sections a, header .site-nav-trigger, header sl-dropdown > [slot=trigger]")]
                .filter((el) => el.getBoundingClientRect().width > 0)
                .map((el) => { const r = el.getBoundingClientRect(); return [r.left, r.right]; })""")
            if any(left < -SLACK_PX or right > width + SLACK_PX for left, right in names):
                problems.append(f"{width} px: a control lies outside the window: {names}")
            tops = page.evaluate("""() => [...document.querySelectorAll("header .site-name, header .site-nav-trigger, header sl-dropdown.site-gear")]
                .filter((el) => el.getBoundingClientRect().width > 0).map((el) => Math.round(el.getBoundingClientRect().top))""")
            if len(set(tops)) > 1 and max(tops) - min(tops) > 12:
                problems.append(f"{width} px: the bar wraps onto a second row: {tops}")
    finally:
        context.close()
    assert not problems, "\n  ".join(problems)
