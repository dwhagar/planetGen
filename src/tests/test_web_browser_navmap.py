"""
NAV.41: the NAV page's course map in a real browser. Its stop names, the
bearing label and the scale label are at least the page's body text size
at phone and desktop widths, the map spans the page's width, and the
labels shown on a crowded route don't overlap.

The panel comes straight from `planetgen.web.maps.navmap.render_nav_map_panel` with the
site's own `style.css`, so no database is needed. Skipped without
Playwright or Chromium.
"""

import os

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from planetgen.web.maps.navmap import render_nav_map_panel  # noqa: E402

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

_STYLE = os.path.join(_SRC_DIR, "html", "static", "style.css")


def _link(name, **params):
    return f"/{name}?" + "&".join(f"{key}={value}" for key, value in sorted(params.items()))


def _crowded_route():
    """Twelve stops, several bunched together, with long names."""
    positions = [(0.0, 0.0), (0.4, 0.1), (0.8, 0.3), (1.0, 0.2), (3.0, 1.0), (3.2, 1.1),
                 (3.3, 0.9), (6.0, 2.0), (7.5, 1.0), (7.6, 1.2), (9.0, 3.0), (12.0, 4.0)]
    waypoints = []
    for index, (x, y) in enumerate(positions):
        role = "origin" if index == 0 else "destination" if index == len(positions) - 1 else "hop"
        waypoints.append({"id": index + 1, "name": f"Kepler Station {index + 1}", "position": (x, y, 0.0),
                          "role": role})
    return waypoints


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


MEASURE = """
() => {
  const body = parseFloat(getComputedStyle(document.body).fontSize);
  const viewport = document.querySelector(".navmap-viewport").getBoundingClientRect();
  const labels = [...document.querySelectorAll(".navmap-label")]
    .filter((el) => getComputedStyle(el).display !== "none")
    .map((el) => {
      const r = el.getBoundingClientRect();
      return { text: el.textContent, size: parseFloat(getComputedStyle(el).fontSize),
               left: r.left, right: r.right, top: r.top, bottom: r.bottom };
    });
  const panel = document.querySelector(".panel");
  const pad = getComputedStyle(panel);
  const panelInner = panel.clientWidth - parseFloat(pad.paddingLeft) - parseFloat(pad.paddingRight);
  return { body, viewportWidth: viewport.width, panelInner, labels };
}
"""


@pytest.mark.parametrize("width", [390, 1280])
def test_course_map_labels_are_readable_and_do_not_overlap(browser, width):
    with open(_STYLE, encoding="utf-8") as handle:
        css = handle.read()
    panel = render_nav_map_panel(_link, _crowded_route(), has_route=True)
    context = browser.new_context(viewport={"width": width, "height": 900}, device_scale_factor=1)
    try:
        page = context.new_page()
        page.set_content(f"<!doctype html><html><head><style>{css}</style></head>"
                         f"<body><main>{panel}</main></body></html>")
        result = page.evaluate(MEASURE)
    finally:
        context.close()

    assert result["viewportWidth"] == pytest.approx(result["panelInner"], abs=1), "the map should span the panel"
    labels = result["labels"]
    texts = [label["text"] for label in labels]
    assert "Kepler Station 1" in texts and "Kepler Station 12" in texts
    assert "Bearing 000" in texts
    for label in labels:
        assert label["size"] >= result["body"], f"{label['text']} is {label['size']}px, body is {result['body']}px"
    stops = [label for label in labels if label["text"].startswith("Kepler")]
    for i, a in enumerate(stops):
        for b in stops[i + 1:]:
            overlap = a["left"] < b["right"] - 0.5 and b["left"] < a["right"] - 0.5 \
                and a["top"] < b["bottom"] - 0.5 and b["top"] < a["bottom"] - 0.5
            assert not overlap, f"{a['text']} overlaps {b['text']} at {width}px"
