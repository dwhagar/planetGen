"""
TEST.56: no two controls overlap, and none runs off the screen.

Every page `test_web_a11y.py` checks (every GET page of the `web`
blueprint, against the same small generated database) is loaded in
headless Chromium at 390, 600, 820 and 1280 px wide, in the light and
the dark color scheme, and the bounding boxes of every visible control
(buttons, links styled as buttons, form fields, menus' summaries, map
buttons) are compared:

- no two intersect (one inside the other, such as an input in its
  label, doesn't count; touching edges don't either);
- none reaches past the left or right edge of the viewport, or above
  the top of the page, unless it sits in a box that scrolls sideways
  (a wide table), where scrolling reaches it.

At 390 px the header's menus are opened and checked too, as the
accessibility test does.

Each page is checked as it opens and again with every folded section
(`<details>`) open. Skipped without Playwright, Chromium or the MySQL
test server, like `test_web_a11y.py`.
"""

import pytest

pytest.importorskip("playwright.sync_api")

from planetgen.api.authz import SESSION_COOKIE_NAME  # noqa: E402

from tests.test_web_a11y import (  # noqa: E402,F401 -- fixtures
    PAGE_ENDPOINTS,
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

WIDTHS = (390, 600, 820, 1280)
HEIGHT = 900
SCHEMES = ("light", "dark")

# The maps' menus (Menu, Bookmarks, the history panel) are popovers that
# drop over what is below them, like the header's (TEST.106).
POPOVER_MENUS = "details.galaxy-menu, details.bookmarks-menu, details.galaxy-steps"

CONTROLS = ", ".join([
    "button", "a.btn", "input:not([type=hidden])", "select", "textarea", "summary",
    "sl-button", "sl-icon-button", "[role=button]", ".starmap-btn",
])

# Collects every visible control's box and finds the problems, in the
# page. A control is skipped when it is not rendered (no box, hidden,
# fully transparent), is a visually hidden screen-reader helper (1 px or
# clipped away) or is the skip link (off screen until focused).
FIND_PROBLEMS = """
(selector) => {
  const view = document.documentElement.clientWidth;
  const label = (el) => {
    let s = el.tagName.toLowerCase();
    if (el.id) s += "#" + el.id;
    if (typeof el.className === "string" && el.className.trim()) s += "." + el.className.trim().split(/\\s+/).join(".");
    for (const a of ["data-action", "name", "type", "href"]) if (el.getAttribute(a)) { s += `[${a}="${el.getAttribute(a)}"]`; break; }
    const text = (el.textContent || el.value || el.getAttribute("aria-label") || "").trim().replace(/\\s+/g, " ");
    return s + (text ? ` "${text.slice(0, 30)}"` : "");
  };
  const shown = (el) => {
    // Also false inside a closed <details> (its content is skipped, not
    // display: none).
    if (el.checkVisibility && !el.checkVisibility({ opacityProperty: true, visibilityProperty: true, contentVisibilityAuto: true })) return false;
    for (let up = el; up; up = up.parentElement) {
      const cs = getComputedStyle(up);
      if (cs.display === "none" || cs.visibility === "hidden" || cs.opacity === "0") return false;
      if (cs.clip === "rect(0px, 0px, 0px, 0px)" || cs.clipPath === "inset(50%)") return false;
    }
    return true;
  };
  // The part of the box any clipping ancestor leaves visible, and whether
  // one of those ancestors scrolls sideways.
  const visibleBox = (el) => {
    const r = el.getBoundingClientRect();
    let box = { left: r.left, top: r.top, right: r.right, bottom: r.bottom };
    let scrolls = false;
    for (let up = el.parentElement; up && up !== document.body; up = up.parentElement) {
      const cs = getComputedStyle(up);
      if (cs.overflowX !== "visible" || cs.overflowY !== "visible") {
        const u = up.getBoundingClientRect();
        if (cs.overflowX !== "visible") {
          if (cs.overflowX === "auto" || cs.overflowX === "scroll") scrolls = true;
          box.left = Math.max(box.left, u.left);
          box.right = Math.min(box.right, u.right);
        }
        if (cs.overflowY !== "visible") {
          box.top = Math.max(box.top, u.top);
          box.bottom = Math.min(box.bottom, u.bottom);
        }
      }
    }
    return { box, scrolls };
  };
  const items = [];
  for (const el of document.querySelectorAll(selector)) {
    // Map markers inside an <svg> overlap by design (a belt rings its
    // star); the map's own buttons are HTML.
    if (el.closest(".skip-link") || el.closest("svg") || !shown(el)) continue;
    const r = el.getBoundingClientRect();
    if (r.width <= 1 || r.height <= 1) continue;
    const { box, scrolls } = visibleBox(el);
    if (box.right - box.left <= 1 || box.bottom - box.top <= 1) continue;
    items.push({ el, r, box, scrolls });
  }
  const problems = [];
  const EPS = 0.5;
  for (const it of items) {
    if (it.scrolls) continue;
    const r = it.r;
    if (r.left < -EPS || r.right > view + EPS || r.top + window.scrollY < -EPS) {
      problems.push(`off screen: ${label(it.el)} at ${Math.round(r.left)}..${Math.round(r.right)} of ${view}`);
    }
  }
  for (let i = 0; i < items.length; i++) {
    for (let j = i + 1; j < items.length; j++) {
      const a = items[i], b = items[j];
      if (a.el.contains(b.el) || b.el.contains(a.el)) continue;
      const w = Math.min(a.box.right, b.box.right) - Math.max(a.box.left, b.box.left);
      const h = Math.min(a.box.bottom, b.box.bottom) - Math.max(a.box.top, b.box.top);
      if (w > EPS && h > EPS) {
        problems.push(`overlap ${Math.round(w)}x${Math.round(h)} px: ${label(a.el)} and ${label(b.el)}`);
      }
    }
  }
  return problems;
}
"""


SETTLE_SAMPLES = 4
SETTLE_STEP_MS = 60
SETTLE_LIMIT = 40


def _problems(page):
    """The layout problems of the page once it has settled (TEST.103): a map
    that is still placing its controls, or a menu mid-transition, shows a
    layout that isn't the one a visitor sees, and how long that takes depends
    on the machine. Sampled until the answer is empty or the same
    `SETTLE_SAMPLES` times in a row, so only a lasting overlap is reported."""
    steady, last = 0, None
    for _ in range(SETTLE_LIMIT):
        found = page.evaluate(FIND_PROBLEMS, CONTROLS)
        if not found:
            return found
        steady = steady + 1 if found == last else 1
        if steady >= SETTLE_SAMPLES:
            return found
        last = found
        page.wait_for_timeout(SETTLE_STEP_MS)
    return found


def test_the_check_finds_overlaps_and_controls_off_screen(browser):
    """The check itself, on a page made to fail it (and one that
    passes)."""
    page = browser.new_page(viewport={"width": 400, "height": 300})
    try:
        page.set_content("""
            <style>body { margin: 0 } .a { position: absolute; left: 10px; top: 10px; width: 80px; height: 30px }</style>
            <button class="a" id="one">One</button>
            <button class="a" id="two" style="left: 60px">Two</button>
            <label style="position: absolute; top: 100px; left: 10px">Name <input id="inside"></label>
            <a class="btn" style="position: absolute; top: 200px; left: 350px; width: 100px">Wide</a>
            <span style="position: absolute; top: 250px; left: 0; display: block; width: 400px; overflow-x: auto">
              <button style="margin-left: 500px">In a scroller</button></span>
            <details style="position: absolute; top: 150px"><summary>More</summary><button>Hidden</button></details>
            <button class="sr-only" style="position: absolute; width: 1px; height: 1px; overflow: hidden">sr</button>
            <svg style="position: absolute; top: 10px; left: 10px" width="50" height="50"><g role="button"><circle r="20"/></g></svg>
        """)
        problems = _problems(page)
    finally:
        page.close()
    assert len(problems) == 2, problems
    assert problems[0].startswith('off screen: a.btn "Wide" at 350..450 of 400')
    assert problems[1].startswith('overlap 30x30 px: button#one.a "One" and button#two.a "Two"')


@pytest.mark.parametrize("width", WIDTHS)
@pytest.mark.parametrize("endpoint", PAGE_ENDPOINTS)
def test_controls_do_not_overlap(browser, base_url, page_targets, admin_token, endpoint, width):
    path, needs_admin = page_targets[endpoint]
    if path is None:
        pytest.skip(needs_admin)

    problems = []
    for scheme in SCHEMES:
        context = browser.new_context(viewport={"width": width, "height": HEIGHT}, color_scheme=scheme,
                                      device_scale_factor=1, reduced_motion="reduce")
        if needs_admin:
            context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
        page = context.new_page()
        try:
            response = page.goto(base_url + path, wait_until="load")
            assert response.status == 200, f"{path} answered {response.status}"
            # The maps draw their own buttons (breadcrumb, slab slider)
            # once their data is in.
            page.wait_for_timeout(400)
            problems += [f"{scheme}: {p}" for p in _problems(page)]
            # Then with every folded section open (the header's menus
            # and the maps' popover menus drop down over the page by
            # design, so the header's open one at a time below and the
            # maps' stay closed: TEST.106).
            opened = page.evaluate("""(POPOVERS) => {
                let n = 0;
                for (const d of document.querySelectorAll("details:not([open])")) {
                    if (d.closest("header") || d.matches(POPOVERS)) continue;
                    d.open = true;
                    n++;
                }
                return n;
            }""", POPOVER_MENUS)
            if opened:
                page.wait_for_timeout(50)
                problems += [f"{scheme}, sections open: {p}" for p in _problems(page)]
            if width == 390:
                menus = ["sl-dropdown.site-gear"]
                if page.locator("sl-dropdown.site-menu > [slot=trigger]").is_visible():
                    menus.insert(0, "sl-dropdown.site-menu")
                # A header menu drops down over the page's first rows (the
                # breadcrumb, a page's action bar), by design (NAV.14: its
                # "..." can sit under it), so the page's own controls are
                # hidden while the menu's are checked (TEST.106).
                page.evaluate("""() => document.querySelectorAll("main").forEach((n) => { n.style.visibility = "hidden"; })""")
                for menu in menus:
                    page.locator(f"{menu} > [slot=trigger]").click()
                    page.wait_for_timeout(50)
                    problems += [f"{scheme}, {menu} open: {p}" for p in _problems(page)]
                    page.locator(f"{menu} > [slot=trigger]").click()
        finally:
            context.close()

    assert not problems, f"{path} at {width} px:\n  " + "\n  ".join(sorted(set(problems)))
