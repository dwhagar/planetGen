"""
The page scripts' own unit tests (TEST.57, TEST.58): every
`tests/js/*.test.mjs` file, run with node's built-in test runner
(`node --test`), skipped where node isn't installed.

Those files load `html/static/*.js` against `tests/js/fakedom.mjs`, a
small stand-in for the browser's DOM (no browser, no WebGL), and cover
the Galaxy Map's stage view and controls, the phenomenon diagram's zoom,
the Generate page's job panel and the facility form. What they need from
the Python side (the galaxy's density shape, sector designations, the
Galaxy Map panel's buttons) is computed here and handed over as JSON in
`PLANETGEN_JS_FIXTURES`, so the two languages are checked against each
other rather than against copied numbers.

Run one file by hand from the repository root with
`PLANETGEN_JS_FIXTURES="$(python -m tests.test_js_unit)" node --test
src/tests/js/<name>.test.mjs` (from `src/`).
"""

import glob
import json
import math
import os
import re
import shutil
import subprocess

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))

NODE = shutil.which("node")
JS_TESTS = sorted(glob.glob(os.path.join(HERE, "js", "*.test.mjs")))

EDGE_PC = 4.0
GALAXY_RADIUS_PC = 15000.0


def fixtures():
    """The values the JavaScript tests compare against, as a dict."""
    from planetgen.web.maps.galaxymap3d import render_galaxy_map3d_panel
    from planetgen.galaxy.density import build_galaxy_shape, shape_with_terms
    from planetgen.galaxy.spherical import band_range, provisional_sector_designation, wedge_count
    from planetgen.galaxy.skeleton import expected_system_count_at_density_1
    from planetgen.physics.units import pc_to_ly

    shape = build_galaxy_shape(2800.0, 350.0, 200.0, 1.0, 2, math.radians(15.0), 0.4)
    shape = {**shape_with_terms(shape),
             "sector_min_density": 1.0 / expected_system_count_at_density_1(pc_to_ly(EDGE_PC))}

    designations = []
    # Shells up to the designation's 12 bits (past 32 bits of packing, which is why the browser needs BigInt); bands within the
    # designation's 12-bit field (|band| < 2040 here).
    for shell in (0, 1, 2, 7, 100, 3749, 4095):
        low, high = band_range(shell)
        for band in sorted({max(low, -2040), max(low, -1), 0, min(high, 1), min(high, 2040)}):
            if not low <= band <= high:
                continue
            count = wedge_count(shell, band)
            for wedge in sorted({0, count // 2, count - 1}):
                designations.append({"shell": shell, "band": band, "wedge": wedge,
                                     "designation": provisional_sector_designation(shell, band, wedge)})

    def actions(html):
        # The row and its Menu (MAP.55), up to the side panel's end.
        controls = html[html.index('id="galaxymap3d-controls"'):html.index('<script type="application/json"')]
        return re.findall(r'data-action="([^"]+)"', controls)

    view = {"stamp": "0123456789abcdef", "tiles": {}, "edge_pc": EDGE_PC, "has_shape": True}
    return {
        "starLight": star_light(),
        "sectorMap": sector_map(),
        "edgePc": EDGE_PC,
        "galaxyRadius": GALAXY_RADIUS_PC,
        "shape": shape,
        "designations": designations,
        "orbitPositions": orbit_positions(),
        "tilePack": tile_pack(),
        "tileCut": tile_cut(),
        "galaxyControls": {
            "public": actions(render_galaxy_map3d_panel("db", None, EDGE_PC, view)),
            "admin": actions(render_galaxy_map3d_panel("db", None, EDGE_PC, view, generate={"url": "/x"})),
        },
    }


def orbit_positions():
    """MAP.70's `positions_at` over a hand-built scene at several times, for
    static/orbitpositions.js to match."""
    from planetgen.physics.body_positions import positions_at
    from tests.test_body_positions import sample_scene

    scene = sample_scene()
    return {"scene": scene, "cases": [{"years": y, "positions": positions_at(scene, y)}
                                      for y in (0.0, 0.05, 1.8, 37.5, -2.2, 400.0, 3000.0, -2.7e4, 1.0e5)],
            "turns": wrapped_turns(), "stumpff": stumpff_values(),
            "steps": universal_steps()}


def wrapped_turns():
    """GEN.374's `kepler.wrapped_turn_deg` for steps from a day to a billion days, for `wrappedTurnDeg` to match."""
    from planetgen.physics.kepler import wrapped_turn_deg

    return [{"elapsed": elapsed, "period": period, "degrees": wrapped_turn_deg(elapsed, period)}
            for period in (0.319 / 365.25, 1.0, 11.86, 1.3e-3, 0.0)
            for elapsed in (0.0, 1.0 / 365.25, 2.5, -7.25, 3000.0, 2.7e6, 1.0e9 / 365.25 * 1000.0, -4.0e6)]


def stumpff_values():
    """GEN.140's `kepler.stumpff_c2_c3` from deep in the series to well past the limit, for the compiled-for-the-browser
    copy to match."""
    from planetgen.physics.kepler import stumpff_c2_c3

    return [{"z": z, "c": list(stumpff_c2_c3(z))}
            for z in (0.0, 1e-9, 1e-3, -1e-3, 0.3, -0.3, 0.999, -0.999, 1.0, -1.0, 5.0, -5.0, 40.0, -40.0)]


def universal_steps():
    """GEN.289's `kepler.universal_step` over a circular, an eccentric, a thin and an open state for steps of every
    size and sign, for `universalStep` in static/orbitpositions.js to match."""
    from planetgen.physics.kepler import universal_step

    mu = 1.327e11
    states = [((1.5e8, 0.0, 0.0), (0.0, mu ** 0.5 / 1.5e8 ** 0.5, 0.0)),
              ((1.5e8, 2.0e7, -1.0e6), (-5.0, 29.0, 1.0)),
              ((1.5e8, 0.0, 0.0), (0.0, 12.0, 0.0)),
              ((1.5e8, 0.0, 0.0), (0.0, 45.0, 0.0))]
    cases = []
    for position, velocity in states:
        for dt in (0.0, 1.0e3, 1.0e7, 3.15576e7, -1.0e6, 3.15576e10):
            moved = universal_step(position, velocity, mu, dt)
            cases.append({"position": list(position), "velocity": list(velocity), "mu": mu, "dt": dt,
                          "position_after": list(moved[0]), "velocity_after": list(moved[1])})
    return cases


def tile_pack():
    """MAP.159's packed `/galaxy/tiles` answer (`tilepack.pack_response`) of three
    tiles -- stars with and without a temperature, radius, class and
    system, an id past 32 bits, generated stars with printed IDs, clouds,
    points, an incomplete tile, an empty one and one with mixed ids -- as base64, with the tiles
    `tilepack.unpack_response` reads from it, for static/tilepack.js to
    match."""
    import base64
    from planetgen.web.lib import tilepack

    def star(star_id, x, y, z, luminosity, **rest):
        return {"id": star_id, "x": x, "y": y, "z": z, "luminosity_sol": luminosity, **rest}

    full = {
        "clouds": [{"type": "nebula", "id": 4, "x": 1.0, "y": 2.0, "z": 3.0, "radius_pc": 5.5}],
        "points": [{"type": "black_hole", "id": "0D01F3", "x": -1.5, "y": 0.25, "z": 9.0}],
        "stars": [
            star(7, 1234.567, -8000.125, 40.0, 1.2e5, temperature_k=21500, radius_sol=14.0, star_type="B", system_id=3),
            star(2 ** 40 + 3, 1230.0, -7999.5, 41.5, 0.0004, star_type="M"),
            star(9, 1240.25, -7990.0, 39.0, 3.0, temperature_k=5778, radius_sol=1.0, star_type="G", tier=2),
        ],
        "generated": [star("FE81000A2B-0000005-000", 1.0, 2.0, 3.0, 0.85, temperature_k=5200, star_type="K", system_id=12),
                      star("FE81000A2B-0000006-000", 1.5, 2.0, 3.0, 0.5, star_type="M")],
    }
    parts = [("3/1/2/3", tilepack.pack_tile(full)),
             ("3/1/2/4", tilepack.pack_tile({"stars": [star(1, 5.0, 5.0, 5.0, 1.0)], "incomplete": ["stars"]})),
             ("3/1/2/5", tilepack.pack_tile({})),
             ("3/1/2/6", tilepack.pack_tile({"stars": [star(3, 1.0, 1.0, 1.0, 1.0), star("odd", 2.0, 2.0, 2.0, 1.0)]}))]
    body = tilepack.pack_response({"stamp": "0123456789abcdef", "generation": "g1", "edge_pc": 4.0}, parts)
    header, tiles = tilepack.unpack_response(body)
    return {"body": base64.b64encode(body).decode("ascii"), "header": header, "tiles": tiles,
            "sizes": {key: len(raw) for key, raw in parts}}


def tile_cut():
    """MAP.152: the distance cut and the tiles to fetch for a few cameras, which
    `static/galaxytilecut.js` must reproduce key for key."""
    from planetgen.galaxy.viewport import (
        DISTANCE_CUT_EDGES, DISTANCE_CUT_FADE_BAND, TILE_MAX_LEVEL, TILE_ROOT_EDGE_PC,
        distance_cut_tiles, fetch_cut_tiles,
    )

    bounds = ((-12600.0, -12600.0, -960.0), (12600.0, 12600.0, 960.0))
    cases = []
    for camera, look in (
        ((34000.0, -20000.0, 28000.0), (0.0, 0.0, 0.0)),
        ((8000.0, 0.0, 0.0), (8000.0, 0.0, -100.0)),
        ((8012.0, 21.0, 9.0), (8000.0, 4.0, 4.0)),
        ((0.0, 0.0, 0.0), (1.0, 1.0, 0.0)),
        ((-5000.0, 3000.0, 700.0), (-5100.0, 3000.0, 700.0)),
    ):
        cases.append({
            "camera": camera, "look": look, "bounds": bounds,
            "cut": distance_cut_tiles(camera, look, bounds),
            "fetch": fetch_cut_tiles(camera, look, bounds),
        })
    return {"edges": DISTANCE_CUT_EDGES, "band": DISTANCE_CUT_FADE_BAND, "root": TILE_ROOT_EDGE_PC,
            "finestLevel": TILE_MAX_LEVEL - 1, "cases": cases}


def star_light():
    """MAP.87's brightness boost (`starmap.star_light_boost`) and a
    10 px, 0.5-strength halo boosted by it (`_boost_light`), from 1e-6 to
    1e7 L_sun, for static/starlight.js to match."""
    from planetgen.web.maps.starmap import _boost_light, star_light_boost

    rows = []
    for exponent in range(-60, 71, 3):
        luminosity = 10 ** (exponent / 10)
        boost = star_light_boost(luminosity)
        size, glow = _boost_light(10.0, 0.5, boost)
        rows.append({"luminositySol": luminosity, "boost": boost, "sizePx": size, "glow": glow})
    return rows


def sector_map():
    """The Sector Map's scene for two systems, a nebula, a rogue planet and a
    neighbour (test_starmap.py's own fixtures): the JSON
    `/sector/<id>/scene` answers."""
    from planetgen.web.maps.starmap import map_scene_data
    from tests.test_starmap import _link, _make_system, _phenomenon

    near = _make_system(100.0, 50.0, -30.0)
    far = dict(_make_system(-200.0, -120.0, 80.0), id=2, name="Far System")
    data = map_scene_data(
        _link, 1000.0, (5, 0, 17), (500.0, 200.0, -100.0), [near, far],
        phenomena=[_phenomenon(), dict(_phenomenon(type_="rogue_planet", descriptor=None, radius_ly=0.0,
                                                   offset=(-3.0, 1.0, 2.0)), id=2, name="Lonely")],
        neighbors=[{
            "shell_index": 6, "band_index": -2, "wedge_index": 9, "direction_pc": [1.0, 0.0, 0.0],
            "exists": True, "sector_id": 42, "sector_name": "Next Door", "designation": "600000063",
        }],
    )
    return {"data": data}


@pytest.fixture(scope="module")
def fixture_json():
    return json.dumps(fixtures())


@pytest.mark.skipif(NODE is None, reason="node is not installed")
@pytest.mark.parametrize("path", JS_TESTS, ids=[os.path.basename(p) for p in JS_TESTS])
def test_js_file(path, fixture_json):
    env = dict(os.environ, PLANETGEN_JS_FIXTURES=fixture_json)
    result = subprocess.run([NODE, "--test", "--test-reporter=spec", path], capture_output=True, text=True,
                            env=env, timeout=300)
    if result.returncode != 0:
        pytest.fail(f"node --test {os.path.basename(path)} failed:\n{result.stdout[-20000:]}\n{result.stderr[-5000:]}")


def test_every_js_test_file_is_collected():
    """A misnamed file (not *.test.mjs) would silently never run."""
    names = [n for n in os.listdir(os.path.join(HERE, "js")) if n.endswith((".mjs", ".js"))]
    stray = [n for n in names if not n.endswith(".test.mjs") and n != "fakedom.mjs" and not n.startswith("_")]
    assert stray == [], f"tests/js files that would never run: {stray}"


if __name__ == "__main__":
    print(json.dumps(fixtures()))
