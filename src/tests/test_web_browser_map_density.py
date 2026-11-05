# tests/test_web_browser_map_density.py

"""
Star density across the Galaxy Map's tiles (MAP.366): Boss saw square chunks of the star map missing while a tile beside
them at about the same distance kept its stars. The cut gives every tile the same number of stars (up to 96 bright and 96
generated) however big it is on screen, so a small tile beside a big one drew as many stars in a tenth of the area.
`static/starthin.js` keeps the density on screen of the tiles that reached their cap near that of the bigger tiles in the
view. These tests dive into the large synthetic galaxy (`large_site`, 34,000 seeded stars, every tile with its real budget)
and compare, at the same camera, the density of the capped tiles with the rule on and off.

Run with: pytest tests/test_web_browser_map_density.py
"""

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.map_site_support import LARGE_CENTER_PC, large_site  # noqa: E402,F401
from tests.test_web_browser_fixture_maps import GALAXY_CANVAS, _open, _settle, gl_browser, page  # noqa: E402,F401
from tests.test_web_browser_map_perf import INSTRUMENT, QUIET_MS, STEP, _wheel  # noqa: E402,F401

pytestmark = pytest.mark.timeout(900)

CAP = 90
"""int: Stars of one list from which a tile counts as capped (the server's cap is 96)."""
PX_CAP = 1200.0
"""float: Tiles over this many pixels count as this big (static/starthin.js PX_CAP)."""
STEPS = 8
"""int: Wheel steps of the dive; the density is read after each."""

NODE_STARS = """() => document.querySelector('#galaxymap3d-canvas').galaxyNodeStars()"""
THINNING = """(on) => document.querySelector('#galaxymap3d-canvas').galaxyStarThinning(on)"""


def _spread(nodes):
    """`(p90 / p10 of the stars shown per square pixel among the capped tiles, their number, the stars shown)`; the spread is
    `None` with fewer than five capped tiles."""
    capped = [n for n in nodes if n["level"] > 1 and (n.get("brightListed", 0) >= CAP or n.get("generatedListed", 0) >= CAP)]
    shown = sum(n["shown"] for n in nodes)
    if len(capped) < 5:
        return None, len(capped), shown
    density = sorted(n["shown"] / min(n["px"], PX_CAP) ** 2 for n in capped)
    return density[int(0.9 * len(density))] / max(density[int(0.1 * len(density))], 1e-12), len(capped), shown


def _read(page):
    page.wait_for_timeout(300)
    return page.evaluate(NODE_STARS)


def test_capped_tiles_at_one_view_differ_less_in_density_than_their_size_alone_would_make_them(page, large_site):
    _open(page, f"{large_site}/galaxy", GALAXY_CANVAS)
    page.evaluate(INSTRUMENT)
    _settle(page)
    page.evaluate("([x, y, z]) => document.querySelector('#galaxymap3d-canvas').galaxyCenterOn(x, y, z)", list(LARGE_CENTER_PC))
    rows = []
    for i in range(STEPS):
        page.evaluate(STEP, [QUIET_MS])
        page.evaluate(THINNING, False)
        plain = _spread(_read(page))
        page.evaluate(THINNING, True)
        thinned = _spread(_read(page))
        rows.append((i, plain, thinned))
        _wheel(page, -400)
    measured = [(i, p, t) for i, p, t in rows if p[0] is not None and t[0] is not None]
    assert len(measured) >= 4, rows
    problems = []
    for i, plain, thinned in measured:
        if not thinned[0] < plain[0]:
            problems.append(f"step {i}: the density spread did not fall ({plain[0]:.1f} to {thinned[0]:.1f})")
        if thinned[0] > 12.0:
            problems.append(f"step {i}: the capped tiles still differ {thinned[0]:.1f} times in density")
        if thinned[2] < 0.25 * plain[2]:
            problems.append(f"step {i}: the rule left {thinned[2]:.0f} of {plain[2]:.0f} stars shown")
    assert not problems, "\n".join(problems)
    assert max(p[0] for _i, p, _t in measured) >= 10.0, "no step showed the chunks the rule is for: the fixture no longer makes them"


def test_the_stars_shown_change_smoothly_as_the_camera_moves(page, large_site):
    _open(page, f"{large_site}/galaxy", GALAXY_CANVAS)
    page.evaluate(INSTRUMENT)
    _settle(page)
    page.evaluate("([x, y, z]) => document.querySelector('#galaxymap3d-canvas').galaxyCenterOn(x, y, z)", list(LARGE_CENTER_PC))
    page.evaluate(STEP, [QUIET_MS])
    for _ in range(3):
        _wheel(page, -400)
        page.evaluate(STEP, [QUIET_MS])
    previous = sum(n["shown"] for n in _read(page))
    for i in range(6):
        _wheel(page, -120)
        page.evaluate(STEP, [QUIET_MS])
        now = sum(n["shown"] for n in _read(page))
        assert abs(now - previous) <= 0.35 * previous, f"small step {i}: {previous:.0f} stars shown became {now:.0f}"
        previous = now
