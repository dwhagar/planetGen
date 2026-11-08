"""MAP.397: every control the Map View menu shows is something the visitor can use.

The menu lists the map's layers (`html/static/maplayers.js` `renderMenu`): a button for each, a second "Forward" button for the points of
interest, a row of three that sets a whole group, and headings. A browser test asks `check_menu_controls` about the menu of each map; it
fails on text that is not a heading or a control (a word with nothing behind it), a control that is disabled, one with no name, a "Forward"
button that is not on the same line as, and just after, its layer's button, and a control that does nothing when pressed."""

LIST_JS = """(container) => {
  const root = document.querySelector(container);
  const shown = (e) => e.getClientRects().length > 0;
  const out = [];
  const control = 'button, input, select, textarea, a[href], summary, label';
  root.querySelectorAll('*').forEach((e) => {
    if (!shown(e)) return;
    const own = [...e.childNodes].filter((n) => n.nodeType === 3 && n.textContent.trim()).map((n) => n.textContent.trim()).join(' ');
    const r = e.getBoundingClientRect();
    const inControl = !!e.closest(control) && root.contains(e.closest(control));
    out.push({tag: e.tagName.toLowerCase(), type: e.type || '', cls: e.className || '', text: own, name: (e.getAttribute('aria-label') || (e.labels && e.labels[0] && e.labels[0].textContent) || e.textContent || '').trim(),
              disabled: !!e.disabled, layer: e.dataset.layer || '', forward: e.dataset.forward || '', bulk: e.dataset.bulk || '',
              pressed: e.getAttribute('aria-pressed'), y: Math.round(r.y), x: Math.round(r.x), w: Math.round(r.width), inControl: inControl,
              parentTag: e.parentElement ? e.parentElement.tagName.toLowerCase() : '', parentCls: e.parentElement ? e.parentElement.className : ''});
  });
  return out;
}"""

HEADINGS = ("galaxy-kinds-heading", "view-menu-group-name", "gravity-hint", "view-menu-count", "gravity-slab")
"""tuple: The classes of the text the menu may show that is not a control: a group heading (the older `galaxy-kinds-heading` and the View menu's
`view-menu-group-name`, UX.277) and the one-line note the Gravity settings carry (TEST.696)."""

CONTROLS = ("button", "input", "select", "summary", "label", "a")


def check_menu_controls(page, container, press=True):
    """Checks every visible element of the menu in `container` (a selector); with `press` also presses each Forward button (twice) and each
    layer button (twice) and checks that it changed. Returns the count of Forward buttons and of controls seen."""
    items = page.evaluate(LIST_JS, container)
    assert items, f"{container}: the menu shows nothing"
    problems = []
    layer_at = {}
    for item in items:
        where = f"{item['tag']}[{item['layer'] or item['forward'] or item['bulk'] or item['text'] or item['name']}]"
        if item["text"] and not item["inControl"] and not any(heading in str(item["cls"]) for heading in HEADINGS):
            problems.append(f"{where}: text that is not a heading or a control: {item['text']!r} (class {item['cls']!r}, parent {item['parentTag']} {item['parentCls']!r})")
        if item["tag"] in CONTROLS and item["tag"] != "label":
            if item["disabled"]:
                problems.append(f"{where}: disabled")
            if not item["name"]:
                problems.append(f"{where}: no name ({item['type']!r}, class {item['cls']!r}, in {item['parentTag']} {item['parentCls']!r})")
            if item["tag"] == "button" and item["w"] < 8:
                problems.append(f"{where}: no width")
        if item["tag"] == "button" and item["layer"]:
            layer_at[item["layer"]] = item
    forwards = [item for item in items if item["tag"] == "button" and item["forward"]]
    for item in forwards:
        mine = layer_at.get(item["forward"])
        if mine is None:
            problems.append(f"Forward[{item['forward']}]: no button for its layer")
        elif mine["y"] != item["y"] or not mine["x"] < item["x"]:
            problems.append(f"Forward[{item['forward']}]: not on its layer's line, after it ({mine['x']},{mine['y']} / {item['x']},{item['y']})")
    assert not problems, f"{container}:\n  " + "\n  ".join(problems)
    if press:
        for item in forwards:
            button = page.locator(f"{container} button[data-forward=\"{item['forward']}\"]")
            layer = page.locator(f"{container} button[data-layer=\"{item['forward']}\"]")
            before = (button.get_attribute("aria-pressed"), layer.get_attribute("aria-pressed"))
            button.click()
            first = (button.get_attribute("aria-pressed"), layer.get_attribute("aria-pressed"))
            button.click()
            after = (button.get_attribute("aria-pressed"), layer.get_attribute("aria-pressed"))
            assert first != before or before == ("true", "true"), f"Forward[{item['forward']}] did nothing: {before} -> {first}"
            assert after != first, f"Forward[{item['forward']}] did nothing the second time: {first} -> {after}"
            assert "true" in (first[0],), f"Forward[{item['forward']}] is not pulled forward after the press: {first}"
    return len(forwards), len(items)
