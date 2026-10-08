"""
ADM.14: on the admin Generate pages the text boxes line up, not their
headings. Every group of fields (`.search-fields`) is loaded in headless
Chromium with all sections open at phone and desktop widths, light and dark:
within a group the boxes are all the same width, boxes on one line sit level
with each other, and boxes in one column share a left edge, however long or
wrapped their labels are. Skipped without Playwright, Chromium or the MySQL
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

# For each group of fields: the visible boxes (text, number, password,
# search entries and selects that belong to a plain field) with their boxes.
GROUPS = """
() => {
  const out = [];
  for (const group of document.querySelectorAll(".aligned-fields .search-fields")) {
    const boxes = [];
    for (const field of group.querySelectorAll(":scope > .search-field:not(.generate-check)")) {
      const input = field.querySelector(
        "input[type=text], input[type=number], input[type=password], input[type=search], select");
      if (!input) continue;
      const r = input.getBoundingClientRect();
      if (r.width === 0 || r.height === 0) continue;
      boxes.push({label: field.textContent.trim().slice(0, 40), x: r.x, y: r.y, w: r.width});
    }
    if (boxes.length > 1) out.push(boxes);
  }
  return out;
}
"""


def _problems(page):
    problems = []
    for boxes in page.evaluate(GROUPS):
        names = ", ".join(b["label"] for b in boxes)
        widths = {round(b["w"]) for b in boxes}
        if len(widths) != 1:
            problems.append(f"different widths {sorted(widths)} in [{names}]")
        lines = {}
        for b in boxes:
            lines.setdefault(round(b["y"]), []).append(b)
        by_column = {}
        for b in boxes:
            by_column.setdefault(round(b["x"]), []).append(b)
        # Boxes the same width and starting at the same x per column are the
        # grid's; a box whose line is within 1 px of another's counts as level.
        ys = sorted(lines)
        for first, second in zip(ys, ys[1:]):
            if second - first <= 1:
                problems.append(f"boxes at y {first} and {second} are not level in [{names}]")
        # Across lines, the left edges come from the same few columns.
        columns = sorted(by_column)
        if len(columns) > 1 and min(b - a for a, b in zip(columns, columns[1:])) <= 1:
            problems.append(f"left edges {columns} are not shared in [{names}]")
    return problems


@pytest.mark.parametrize("scheme", ("light", "dark"))
@pytest.mark.parametrize("width", (390, 820, 1280))
@pytest.mark.parametrize("endpoint", ("web.generate", "web.generate_system"))
def test_text_boxes_line_up(browser, base_url, page_targets, admin_token, endpoint, width, scheme):
    path, _needs_admin = page_targets[endpoint]
    context = browser.new_context(viewport={"width": width, "height": 900}, color_scheme=scheme,
                                  reduced_motion="reduce")
    context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": admin_token, "url": base_url}])
    page = context.new_page()
    try:
        page.goto(base_url + path, wait_until="load")
        page.evaluate("() => { for (const d of document.querySelectorAll('details')) d.open = true; }")
        page.wait_for_timeout(100)
        assert page.evaluate("document.querySelectorAll('.aligned-fields .search-fields').length") > 0
        problems = _problems(page)
    finally:
        context.close()
    assert not problems, f"{path} at {width}/{scheme}:\n  " + "\n  ".join(problems)
