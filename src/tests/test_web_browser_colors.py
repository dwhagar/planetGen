"""UX.76: buttons and icons read in both color schemes, in every state.

Playwright fixes the OS scheme and `data-theme` picks the site's own, so the
four combinations (OS dark or light, with the site following or overriding) are
all checked. Hover and focus are forced through the DevTools protocol (Chromium
won't force `:visited`, so that state is checked against the stylesheet below).
"""
import pytest

from tests.map_site_support import SECTOR_ID, map_site  # noqa: F401
from tests.test_web_browser_fixture_maps import gl_browser, page  # noqa: F401

# One of each look of button, as the site's own CSS draws them.
FIXTURE = """
() => {
  const box = document.createElement('div');
  box.id = 'ux76';
  box.style.cssText = 'position:fixed;left:0;top:0;z-index:9999;background:var(--bg);padding:8px';
  box.innerHTML = `
    <button class="btn" id="b-primary">Primary</button>
    <a class="btn" id="b-link" href="/sector/1">Primary link</a>
    <a class="btn btn-secondary" id="b-secondary-link" href="/sector/1">Secondary link</a>
    <a class="btn btn-secondary btn-active" id="b-secondary-active" href="/sector/1">Active link</a>
    <button class="btn btn-secondary" id="b-secondary">Secondary</button>
    <button class="btn btn-secondary btn-bookmark" id="b-bookmark" aria-pressed="true">★</button>
    <button class="btn btn-danger" id="b-danger">Danger</button>
    <button class="starmap-btn starmap-toggle" id="b-toggle" aria-pressed="true">Toggle</button>
    <button class="starmap-btn" id="b-starmap">Starmap</button>
    <button class="icon-btn" id="b-icon"><svg class="icon" id="i-icon" width="24" height="24"><use href="${ICONS}#show-on-map"></use></svg></button>`;
  document.body.append(box);
}
""".replace("${ICONS}", "' + document.querySelector('link[href*=\"style.css\"]').href.replace(/style\\.css.*/, 'icons.svg') + '")

CONTRAST = """
(id) => {
  const el = document.getElementById(id);
  const parse = (c) => { const m = c.match(/rgba?\\(([^)]+)\\)/) || c.match(/color\\(srgb ([^)]+)\\)/);
    const p = m[1].split(/[ ,/]+/).map(parseFloat);
    return c.startsWith('color(') ? {r: p[0] * 255, g: p[1] * 255, b: p[2] * 255, a: p.length > 3 ? p[3] : 1}
                                  : {r: p[0], g: p[1], b: p[2], a: p.length > 3 ? p[3] : 1}; };
  const lum = (c) => { const f = (v) => { v /= 255; return v <= 0.03928 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); };
    return 0.2126 * f(c.r) + 0.7152 * f(c.g) + 0.0722 * f(c.b); };
  const over = (top, under) => ({r: top.r * top.a + under.r * (1 - top.a), g: top.g * top.a + under.g * (1 - top.a),
                                 b: top.b * top.a + under.b * (1 - top.a), a: 1});
  const layers = [];
  for (let e = el; e; e = e.parentElement) {
    const c = parse(getComputedStyle(e).backgroundColor);
    if (c.a > 0) { layers.push(c); if (c.a >= 0.99) break; }
  }
  let bg = {r: 255, g: 255, b: 255, a: 1};
  while (layers.length) bg = over(layers.pop(), bg);
  const fg = parse(getComputedStyle(el).color);
  const text = fg.a < 1 ? over(fg, bg) : fg;
  return (Math.max(lum(text), lum(bg)) + 0.05) / (Math.min(lum(text), lum(bg)) + 0.05);
}
"""

TEXT_IDS = ["b-primary", "b-link", "b-secondary-link", "b-secondary-active", "b-secondary", "b-bookmark",
            "b-danger", "b-toggle", "b-starmap"]


@pytest.mark.parametrize("os_scheme,theme", [("dark", None), ("light", None), ("light", "dark"), ("dark", "light")])
def test_buttons_and_icons_read_in_every_state(page, map_site, os_scheme, theme):
    page.emulate_media(color_scheme=os_scheme)
    page.goto(f"{map_site}/sector/{SECTOR_ID}", wait_until="load")
    if theme:
        page.evaluate("t => document.documentElement.setAttribute('data-theme', t)", theme)
    page.evaluate(FIXTURE)
    session = page.context.new_cdp_session(page)
    session.send("DOM.enable")
    session.send("CSS.enable")
    document = session.send("DOM.getDocument", {"depth": 0})
    failures = []
    for state in (None, "hover", "focus"):
        for ident in TEXT_IDS + ["i-icon"]:
            node = session.send("DOM.querySelector", {"nodeId": document["root"]["nodeId"], "selector": f"#{ident}"})
            session.send("CSS.forcePseudoState", {"nodeId": node["nodeId"],
                                                  "forcedPseudoClasses": [state] if state else []})
            ratio = page.evaluate(CONTRAST, ident)
            need = 3.0 if ident == "i-icon" else 4.5
            if ratio < need:
                failures.append(f"{ident} {state or 'rest'}: {ratio:.2f} < {need}")
            session.send("CSS.forcePseudoState", {"nodeId": node["nodeId"], "forcedPseudoClasses": []})
    assert not failures, f"{os_scheme}/{theme}: " + "; ".join(failures)


def test_visited_secondary_links_keep_their_color():
    from pathlib import Path
    css = (Path(__file__).resolve().parents[1] / "html" / "static" / "style.css").read_text()
    # `a.btn:visited` sets the filled button's dark text; the outlined look must outrank it.
    assert "a.btn.btn-secondary:visited" in css
    assert "a.btn.btn-secondary.btn-active:visited" in css
