# tests/test_poi_icons_browser.py

"""
MAP.320: what a browser draws from the point-of-interest symbols (`static/poi/*.svg`, `static/style.css` tokens), on the
fixture site of `tests/habitability_site.py`:

1. each of the seven symbols is in the key of the habitability page with its accessible name, drawn (solid pixels in its
   colour), in the light and the dark theme;
2. the inhabited world is drawn in the habitability colour of each tier, the uninhabitable one in brown (never red), the
   asteroid field in grey, and brown and grey stand off every tier in lightness as the page really draws them;
3. the seven silhouettes differ from one another, so no symbol depends on its colour to be told from the rest;
4. the Galaxy Map's own key lists the seven, and the habitable-world marker is drawn with its symbol over the map and in
   its legend line.
"""

import io
import re

import numpy as np
import pytest
from PIL import Image

from planetgen.physics.habitability_world import TIER_NAMES
from planetgen.web.lib import poi_icons

from tests.map_site_support import WORLD_FIELD, map_site  # noqa: F401
from tests.test_web_browser_fixture_maps import (  # noqa: F401 -- fixtures and the map drivers the lanes' own tests share
    _markers_until, _open_menu, _open_window, _screen_of, gl_browser, page)

pytestmark = pytest.mark.browser

THEMES = ("light", "dark")
ICON_IDS = tuple(icon.id for icon in poi_icons.ICONS)
SCORED = {icon.id for icon in poi_icons.ICONS if icon.scored}
MIN_LIGHTNESS_GAP = 12.0
"""float: As `tests/test_poi_icons.py`: how far brown and grey sit from every habitability colour in CIE L*."""
SHAPE_PX = 64
"""int: How big the symbols are drawn to compare their silhouettes."""
MAX_SHAPE_OVERLAP = 0.9
"""float: The most two symbols' silhouettes may overlap (intersection over union) and still count as different shapes."""


def _parse(color):
    return tuple(float(value) for value in re.findall(r"[\d.]+", color)[:3])


def _lightness(rgb):
    def linear(value):
        value /= 255.0
        return value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4
    red, green, blue = (linear(value) for value in rgb)
    y = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return 116 * y ** (1 / 3) - 16 if y > 0.008856 else 903.3 * y


def _hue_and_saturation(rgb):
    import colorsys
    hue, _light, saturation = colorsys.rgb_to_hls(*(value / 255.0 for value in rgb))
    return hue * 360, saturation


def _set_theme(page, theme):
    page.emulate_media(color_scheme=theme)
    page.evaluate("t => document.documentElement.setAttribute('data-theme', t)", theme)


def _open_key(page, base, theme):
    page.goto(f"{base}/help/habitability", wait_until="load")
    page.wait_for_selector("#poi-key svg[data-poi]", state="attached")
    _set_theme(page, theme)
    page.locator("#map-symbols").scroll_into_view_if_needed()


def _symbols(page, icon_id):
    return page.locator(f'#poi-key svg[data-poi="{icon_id}"]')


def _shot(locator):
    return np.asarray(Image.open(io.BytesIO(locator.screenshot())).convert("RGB")).astype(float)


def _background(page):
    return _parse(page.evaluate("getComputedStyle(document.querySelector('#map-symbols')).backgroundColor "
                                "=== 'rgba(0, 0, 0, 0)' ? getComputedStyle(document.body).backgroundColor "
                                ": getComputedStyle(document.querySelector('#map-symbols')).backgroundColor"))


def _token_colour(page, token):
    """The colour a custom property of the stylesheet resolves to in the page's theme now, as `(r, g, b)`."""
    return _parse(page.evaluate(
        "t => { const e = document.createElement('i'); e.style.color = `var(${t})`; document.body.appendChild(e);"
        " const c = getComputedStyle(e).color; e.remove(); return c; }", token))


def _colour_of(symbol):
    return _parse(symbol.evaluate("s => getComputedStyle(s).color"))


def _solid(image, rgb, tolerance=14):
    """How many pixels of `image` are `rgb` to within `tolerance` in each channel."""
    return int((np.abs(image - np.array(rgb)).max(axis=2) <= tolerance).sum())


# --- 1. Every symbol is there, named and drawn ----------------------------------------------------------------------

@pytest.mark.parametrize("theme", THEMES)
@pytest.mark.parametrize("icon_id", ICON_IDS)
def test_each_symbol_is_in_the_key_with_its_name_and_is_drawn_in_its_colour(page, map_site, icon_id, theme):
    _open_key(page, map_site, theme)
    symbols = _symbols(page, icon_id)
    expected = len(TIER_NAMES) if icon_id in SCORED else 1
    assert symbols.count() == expected
    label = poi_icons.BY_ID[icon_id].label
    for index in range(expected):
        symbol = symbols.nth(index)
        assert symbol.get_attribute("role") == "img"
        name = symbol.get_attribute("aria-label")
        assert name == (f"{label}: {TIER_NAMES[index]}" if icon_id in SCORED else label), name
        colour = _colour_of(symbol)
        assert _solid(_shot(symbol), colour) >= 30, (icon_id, theme, index, colour)     # the file's shape was drawn, filled with it
    assert page.locator(f"#poi-key-{icon_id}").inner_text().strip().startswith(label)


# --- 2. The colours ---------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("theme", THEMES)
def test_the_inhabited_world_wears_the_habitability_colour_and_the_uninhabitable_one_is_brown_never_red(page, map_site, theme):
    _open_key(page, map_site, theme)
    tier_colours = []
    for index, name in enumerate(TIER_NAMES):
        tier_colours.append(_token_colour(page, f"--tier-{name.lower()}"))
        for icon_id in sorted(SCORED):
            assert _colour_of(_symbols(page, icon_id).nth(index)) == tier_colours[index], (icon_id, name, theme)
    brown = _colour_of(_symbols(page, "uninhabitable-world"))
    grey = _colour_of(_symbols(page, "asteroid-field"))
    hue, saturation = _hue_and_saturation(brown)
    assert 20 <= hue <= 50 and saturation >= 0.3, ("brown", brown, hue, saturation)       # never red
    assert _hue_and_saturation(grey)[1] <= 0.2, ("grey", grey)
    for swatch, colour in (("brown", brown), ("grey", grey)):
        for name, tier in zip(TIER_NAMES, tier_colours):
            assert abs(_lightness(colour) - _lightness(tier)) >= MIN_LIGHTNESS_GAP, (theme, swatch, name)
    assert _colour_of(_symbols(page, "uninhabitable-world")) != _colour_of(_symbols(page, "inhabited-world").nth(3))


# --- 3. The shapes ---------------------------------------------------------------------------------------------------

def _silhouette(page, icon_id):
    """The symbol drawn `SHAPE_PX` across in the key as a boolean mask of the pixels that are not the page behind it."""
    symbol = _symbols(page, icon_id).first
    symbol.evaluate("(s, n) => { s.style.width = n + 'px'; s.style.height = n + 'px'; s.style.stroke = 'none'; }", SHAPE_PX)
    image = _shot(symbol)
    return np.abs(image - np.array(_background(page))).max(axis=2) > 40


def test_the_seven_silhouettes_differ_so_colour_is_never_the_only_sign(page, map_site):
    _open_key(page, map_site, "light")
    masks = {icon_id: _silhouette(page, icon_id) for icon_id in ICON_IDS}
    assert all(mask.sum() > 150 for mask in masks.values()), {icon_id: int(mask.sum()) for icon_id, mask in masks.items()}
    for first_index, first in enumerate(ICON_IDS):
        for second in ICON_IDS[first_index + 1:]:
            if masks[first].shape != masks[second].shape:
                continue
            union = (masks[first] | masks[second]).sum()
            overlap = (masks[first] & masks[second]).sum() / union
            assert overlap < MAX_SHAPE_OVERLAP, (first, second, round(float(overlap), 3))


# --- 4. The Galaxy Map ------------------------------------------------------------------------------------------------

def test_the_galaxy_map_keys_all_seven_and_draws_the_habitable_world_marker_with_its_symbol(page, map_site):
    _open_window(page, map_site, 27, 27, 0, 0, g=3)
    page.wait_for_timeout(1000)
    key = page.locator("#galaxy-map-symbols")
    assert key.count() == 1
    assert {svg.get_attribute("data-poi") for svg in key.locator("svg[data-poi]").all()} == set(ICON_IDS)
    _open_menu(page)
    page.locator('#galaxymap3d-kinds button[data-layer="marker-habitable"]').click()
    found = _markers_until(page, lambda entries: any(e["kind"] == "marker-habitable" for e in entries or []))
    assert {e["name"] for e in found if e["kind"] == "marker-habitable"} == {system["name"] for system in WORLD_FIELD}
    if page.locator("#galaxymap3d-menu").evaluate("m => m.open"):
        page.locator("#galaxymap3d-menu > summary").click()
    page.wait_for_timeout(300)

    # The legend line carries the symbol, in the colour of its token.
    legend = page.locator('#galaxymap3d-layer-legend [data-layer="marker-habitable"]')
    assert legend.count() == 1 and "Habitable worlds" in legend.inner_text()
    symbol = legend.locator('svg[data-poi="habitable-world"]')
    assert symbol.count() == 1
    green = _token_colour(page, "--tier-green")
    assert _colour_of(symbol) == green
    assert _solid(_shot(symbol), green) >= 30

    # And the marker over the map is the symbol, in that colour: a patch round the system holds its pixels.
    page.evaluate("document.querySelector('#galaxymap3d-canvas').scrollIntoView({block: 'center'})")
    page.wait_for_timeout(400)
    eden = next(system for system in WORLD_FIELD if system["name"] == "Fixture Eden")
    spot = _screen_of(page, (eden["x"], eden["y"], eden["z"]))
    assert spot
    box = page.locator("#galaxymap3d-canvas").bounding_box()
    deadline_shots = []
    for _ in range(20):
        image = _shot(page.locator("#galaxymap3d-canvas"))
        x, y = int(spot[0] - box["x"]), int(spot[1] - box["y"])
        patch = image[max(y - 14, 0):y + 15, max(x - 14, 0):x + 15]
        deadline_shots.append(_solid(patch, green, tolerance=45))
        if deadline_shots[-1] >= 25:
            break
        page.wait_for_timeout(150)
    assert deadline_shots[-1] >= 25, deadline_shots      # the symbol's fill (the texture is drawn once its file is read)
