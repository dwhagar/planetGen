# tests/test_map_region_labels_browser.py

"""
MAP.314: the names of charted regions on the Galaxy Map by zoom, on the fixture site of `tests/map_site_support.py`, whose two
charted regions (Fixture Beez Region and Fixture Zarh Region, the 9-blocks 9.0.0.0 and 9.1.0.0) lie in one family (the 81-block
81.0.0.0, Fixture Reach):

1. far out, where the view is cut into 81- and 243-blocks, the family is named, once, and no region by its own name;
2. closer, in 27- and 9-blocks, each charted region is named and the family no longer is;
3. closer still, in blocks of 3 sectors and below, no region is named (the info card of MAP.315 says it);
4. where nothing is charted no region is named, and the Region names switch of the Menu takes the names away;
5. a click on a name flies the camera to its block.

The distances are the camera's, from the middle of the fixture's regions, at which the node cut of `galaxynodecut.js` holds those
sizes on the fixture's sector edge.
"""

import pytest

from tests.map_site_support import REGION_FAMILY, REGION_NAMES, charted_regions, map_site  # noqa: F401
from tests.test_web_browser_fixture_maps import (  # noqa: F401 -- fixtures and the map drivers the lanes' own tests share
    LABELS, _label_spot, _open_galaxy, _settle, _window_size, gl_browser, page)

pytestmark = pytest.mark.browser

BIG_NODES = 700.0
"""float: A camera distance (pc) at which the nodes over the fixture's regions are 81-blocks (those of 243 from 3,000)."""
BIGGEST_NODES = 3000.0
MIDDLE_NODES = 150.0
"""float: A camera distance at which they are 27- and 9-blocks."""
FINE_NODES = 15.0
"""float: A camera distance at which they are blocks of 3 sectors and below."""
SWITCH = "label-region"
REGION_LABELS = "(() => document.querySelector('#galaxymap3d-canvas').__mapLabels.shown().filter((l) => l.kind === 'region').map((l) => l.text).sort())"
NODES_IN = ("() => document.querySelector('#galaxymap3d-canvas').galaxyNodeKeys()"
            ".filter((n) => Number(n.key.split('.')[0]) >= 9).every((n) => n.loaded)")


def _centre():
    """The middle of the fixture's regions, in pc."""
    regions = list(charted_regions().values())
    return [sum(region["centre"][axis] for region in regions) / len(regions) for axis in range(3)]


def _look(page, base, centre, distance):
    """Opens the Galaxy Map looking at `centre` (pc) from `distance` pc, in the orientation the map opens in, and waits for the
    nodes of 9 sectors and up to be in."""
    _open_galaxy(page, base)
    cam = next(part for part in page.evaluate("location.search").lstrip("?").split("&") if part.startswith("cam="))
    quaternion = cam[len("cam="):].split(",")[4:]
    _open_galaxy(page, base, "?cam=" + ",".join(str(value) for value in centre + [distance] + quaternion))
    page.wait_for_function(NODES_IN, timeout=60000)
    page.wait_for_timeout(600)


def _names(page):
    return page.evaluate(f"{REGION_LABELS}()")


def _wait_for_names(page, count=1):
    page.wait_for_function(f"() => {REGION_LABELS}().length >= {count}", timeout=30000)


def test_far_out_the_family_is_named_once_and_no_region_by_its_own_name(page, map_site):
    for distance in (BIG_NODES, BIGGEST_NODES):
        _look(page, map_site, _centre(), distance)
        _wait_for_names(page)
        assert _names(page) == [REGION_FAMILY], (distance, _names(page))


def test_closer_each_charted_region_is_named_and_the_family_is_not(page, map_site):
    _look(page, map_site, _centre(), MIDDLE_NODES)
    _wait_for_names(page, 2)
    assert _names(page) == sorted(REGION_NAMES.values()), _names(page)


def test_in_blocks_of_three_sectors_and_below_no_region_is_named(page, map_site):
    _look(page, map_site, _centre(), FINE_NODES)
    assert _names(page) == []


def test_where_nothing_is_charted_no_region_is_named(page, map_site):
    nowhere = [40000.0, 0.0, 0.0]     # far from every charted region: the views there are cut from blocks with none
    for distance in (BIG_NODES, MIDDLE_NODES):
        _look(page, map_site, nowhere, distance)
        assert _names(page) == [], (distance, _names(page))


def test_the_region_names_switch_takes_the_names_away_and_the_address_keeps_it(page, map_site):
    _look(page, map_site, _centre(), MIDDLE_NODES)
    _wait_for_names(page, 2)
    page.click("#galaxymap3d-menu > summary")
    toggle = page.locator('#galaxymap3d-kinds [data-label="region"]')
    assert toggle.inner_text() == "Region names"
    toggle.click()
    page.wait_for_function(f"() => {REGION_LABELS}().length === 0", timeout=10000)
    assert SWITCH in page.evaluate("new URLSearchParams(location.search).get('hide') || ''")
    toggle.click()
    _wait_for_names(page, 2)


def test_a_click_on_the_family_flies_to_its_block_and_on_a_region_to_the_region(page, map_site):
    _look(page, map_site, _centre(), BIG_NODES)
    _wait_for_names(page)
    family = next(label for label in page.evaluate(LABELS) if label["kind"] == "region")
    page.mouse.click(*_label_spot(page, family))
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=30000)
    _settle(page)
    assert _window_size(page) == 81
    _look(page, map_site, _centre(), MIDDLE_NODES)
    _wait_for_names(page, 2)
    region = next(label for label in page.evaluate(LABELS) if label["kind"] == "region" and label["text"] == REGION_NAMES[(0, 0, 0)])
    page.mouse.click(*_label_spot(page, region))
    page.wait_for_function("() => !document.querySelector('#galaxymap3d-canvas').galaxyMoving()", timeout=30000)
    _settle(page)
    assert _window_size(page) == 9
