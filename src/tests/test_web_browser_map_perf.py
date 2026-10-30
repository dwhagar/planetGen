# tests/test_web_browser_map_perf.py

"""
The Galaxy Map's node path under load (MAP.249): a dive from the whole galaxy to a sector, a flight across a level and the climb back out over a large
synthetic galaxy (`large_site`: tens of thousands of stars, every node carrying its real budget of stars), measuring at
each wheel step how long the view takes to settle (the flight over, the nodes it asked for in, the star buffer
stable), what the step asked for (node requests and bytes), what it drew (stars written and moved in the buffer) and
how the page's frames ran (the longest frame, the longest task on the page's thread).

The figures that vary with the machine (milliseconds under software WebGL) are held to generous limits that catch a
regression of a different order, a dive that stalls or a step that rewrites the whole buffer; the figures that do not
(requests, bytes, stars written) are held close. Set `PLANETGEN_MAP_PERF_REPORT` to a file path to write every step's
figures there as JSON, and `PLANETGEN_MAP_PERF_PROFILE` to one to write the functions the run spent most time in (the numbers in docs/design/interface/galaxy-map-performance.md come from it). A step over its
limit fails with the step and the figure.

Run with: pytest tests/test_web_browser_map_perf.py
"""

import json
import math
import os

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from tests.map_site_support import LARGE_CENTRE_PC, large_site  # noqa: E402,F401
from tests.test_web_browser_fixture_maps import GALAXY_CANVAS, _open, _settle, gl_browser, page  # noqa: E402,F401

pytestmark = pytest.mark.timeout(900)

INSTRUMENT = """() => {
    if (window.__perf) return;
    const perf = window.__perf = {requests: 0, bytes: 0, inflight: 0, frames: [], longTasks: [], last: 0};
    const fetchNow = window.fetch.bind(window);
    window.fetch = function (input, init) {
        const url = typeof input === "string" ? input : (input.url || String(input));   // a Request has .url, a URL (an icon file, MAP.320) is its own string
        const counted = url.indexOf("/nodes") >= 0;
        if (counted) { perf.requests++; perf.inflight++; }
        return fetchNow(input, init).then((response) => {
            if (!counted) return response;
            return response.clone().arrayBuffer().then((buffer) => { perf.bytes += buffer.byteLength; perf.inflight--; return response; });
        }, (error) => { if (counted) perf.inflight--; throw error; });
    };
    try {
        new PerformanceObserver((list) => { for (const e of list.getEntries()) perf.longTasks.push(e.duration); })
            .observe({entryTypes: ["longtask"]});
    } catch (e) { /* no long task entries in this browser */ }
    const frame = (now) => { if (perf.last) perf.frames.push(now - perf.last); perf.last = now; requestAnimationFrame(frame); };
    requestAnimationFrame(frame);
}"""
"""str: Run once per page (`page.evaluate(INSTRUMENT)`): counts the node requests and their bytes, collects the page's
frame intervals and its long tasks (over 50 ms on the main thread)."""

STEP = """([wait_ms]) => new Promise((resolve) => {
    const canvas = document.querySelector("#galaxymap3d-canvas");
    const perf = window.__perf;
    const start = performance.now();
    const before = {requests: perf.requests, bytes: perf.bytes, frames: perf.frames.length, tasks: perf.longTasks.length};
    const fill0 = canvas.galaxyStarFill();
    let syncs = fill0.syncs, quiet = performance.now();
    const poll = () => {
        const now = performance.now();
        const fill = canvas.galaxyStarFill();
        if (fill.syncs !== syncs || canvas.galaxyMoving() || perf.inflight) { syncs = fill.syncs; quiet = now; }
        if (now - quiet >= wait_ms || now - start > 60000) {
            const frames = perf.frames.slice(before.frames).sort((a, b) => a - b);
            const tasks = perf.longTasks.slice(before.tasks);
            const stage = canvas.galaxyStage();
            resolve({
                settleMs: Math.round(quiet - start), timedOut: now - start > 60000,
                requests: perf.requests - before.requests, bytes: perf.bytes - before.bytes,
                written: fill.written - fill0.written, moved: fill.moved - fill0.moved, syncs: fill.syncs - fill0.syncs,
                stars: fill.count, frames: frames.length,
                frameMedianMs: frames.length ? Math.round(frames[frames.length >> 1]) : 0,
                frameMaxMs: frames.length ? Math.round(frames[frames.length - 1]) : 0,
                longestTaskMs: tasks.length ? Math.round(Math.max.apply(null, tasks)) : 0,
                level: stage.win ? stage.win.centre.m : null, dist: canvas.galaxyCamera().dist,
            });
            return;
        }
        setTimeout(poll, 40);
    };
    poll();
})"""
"""str: Waits for the view to settle after what the test just did (no flight, no node request in the air and no change to
the star buffer for `wait_ms`) and returns the figures of the step since it began."""

QUIET_MS = 1000
"""int: How long the view must stay still for a step to count as settled: longer than the map's 300 ms fetch debounce plus the stored-node lookup and the request's start on a busy machine (at 400 ms a climb step could end before its fetch had begun and the next step took the blame, TEST.439), shorter than its 1.5 s prefetch idle."""
WHEEL_DELTA = 400
DIVE_STEPS_MAX = 30
"""The wheel steps of the dive in, at most (it stops at the closest zoom)."""
PAN_PC = 150.0
PAN_STOPS = 6
PAN_AT_PC = 4.0 * PAN_PC
"""The camera distance the dive pans at: the zoom the 3-block window used to be reached at (MAP.337: the wheel no longer steps through windows, the zoom is the camera's alone)."""
"""The flight across: `PAN_STOPS` jumps round a circle of `PAN_PC` about the dive's target, at the zoom the dive reached
the 3-block level at."""

LIMITS = {
    # Node requests and bytes of a step: a step asks for a handful of nodes, each at most NODE_STARS stars (the dive
    # measured 8 requests and 260 kB at most).
    "requests": 12,
    "bytes": 1_000_000,
    # Milliseconds, under software WebGL on a busy machine: generous, for an order-of-magnitude stall (the dive
    # measured 3.2 s to settle, 0.2 s for the longest task, 0.25 s for the longest frame).
    "settleMs": 15000,
    "longestTaskMs": 1500,
    "frameMaxMs": 2000,
    # Stars written or moved in the buffer by a step, as a multiple of the stars the buffer holds after it: a step
    # that rewrote the whole buffer at each of its syncs would write several times that (the climb measured 1.3, as
    # each level it climbs to replaces nearly all the stars on screen).
    "buffer_multiple": 2.0,
}
"""Per step; see the module docstring."""


def _wheel(page, delta):
    box = page.locator(GALAXY_CANVAS).bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2)
    page.mouse.wheel(0, delta)


def _step(page, name, act, steps):
    act()
    figures = page.evaluate(STEP, [QUIET_MS])
    figures["step"] = name
    steps.append(figures)
    return figures


def _profile(page, out_path, run):
    """Runs `run()` under the browser's sampling profiler and writes the 30 functions with the most self time, in
    milliseconds, to `out_path` (for finding what a slow step spends its time on)."""
    session = page.context.new_cdp_session(page)
    session.send("Profiler.enable")
    session.send("Profiler.setSamplingInterval", {"interval": 500})
    session.send("Profiler.start")
    try:
        run()
    finally:
        profile = session.send("Profiler.stop")["profile"]
    nodes = {node["id"]: node["callFrame"] for node in profile["nodes"]}
    own = {}
    for node_id, delta in zip(profile["samples"], profile["timeDeltas"]):
        frame = nodes[node_id]
        key = f"{frame['functionName'] or '(anonymous)'} {frame['url'].rsplit('/', 1)[-1]}:{frame['lineNumber'] + 1}"
        own[key] = own.get(key, 0) + delta / 1000.0
    top = sorted(own.items(), key=lambda item: -item[1])[:30]
    with open(out_path, "w") as out:
        json.dump([[name, round(ms, 1)] for name, ms in top], out, indent=1)


def _dive_and_climb(page, steps):
    page.evaluate("([x, y, z]) => document.querySelector('#galaxymap3d-canvas').galaxyCentreOn(x, y, z)", list(LARGE_CENTRE_PC))
    _step(page, "centre", lambda: None, steps)
    for i in range(DIVE_STEPS_MAX):
        figures = _step(page, f"in {i + 1}", lambda: _wheel(page, -WHEEL_DELTA), steps)
        if figures["dist"] <= PAN_AT_PC and not any(step["step"] == "pan 1" for step in steps):
            # Across the level the dive just reached: jump the centre round a circle.
            for k in range(PAN_STOPS):
                angle = 2.0 * math.pi * (k + 1) / PAN_STOPS
                point = [LARGE_CENTRE_PC[0] + PAN_PC * math.cos(angle), LARGE_CENTRE_PC[1] + PAN_PC * math.sin(angle), LARGE_CENTRE_PC[2]]
                _step(page, f"pan {k + 1}", lambda: page.evaluate(
                    "([x, y, z]) => document.querySelector('#galaxymap3d-canvas').galaxyCentreOn(x, y, z)", point), steps)
            _step(page, "pan back", lambda: page.evaluate(
                "([x, y, z]) => document.querySelector('#galaxymap3d-canvas').galaxyCentreOn(x, y, z)", list(LARGE_CENTRE_PC)), steps)
        if figures["dist"] <= 2.0 + 1e-6 and i > 2:
            break
    page.locator(GALAXY_CANVAS).focus()
    for i in range(12):
        before = page.evaluate("() => document.querySelectorAll('#galaxymap3d-steps [data-steps-panel] li').length")
        _step(page, f"up {i + 1}", lambda: page.keyboard.press("Escape"), steps)
        if page.evaluate("() => document.querySelectorAll('#galaxymap3d-steps [data-steps-panel] li').length") >= before:
            break


def test_the_node_path_stays_within_budget_on_a_dive_a_flight_across_and_a_climb_over_a_large_galaxy(page, large_site):
    _open(page, f"{large_site}/galaxy", GALAXY_CANVAS)
    page.evaluate(INSTRUMENT)
    _settle(page)
    steps = []
    profile = os.environ.get("PLANETGEN_MAP_PERF_PROFILE")
    if profile:
        _profile(page, profile, lambda: _dive_and_climb(page, steps))
    else:
        _dive_and_climb(page, steps)

    report = os.environ.get("PLANETGEN_MAP_PERF_REPORT")
    if report:
        with open(report, "w") as out:
            json.dump(steps, out, indent=1)
    problems = []
    for step in steps:
        if step["timedOut"]:
            problems.append(f"{step['step']}: never settled")
        for key in ("requests", "bytes", "settleMs", "longestTaskMs", "frameMaxMs"):
            if step[key] > LIMITS[key]:
                problems.append(f"{step['step']}: {key} {step[key]} over {LIMITS[key]}")
        for key in ("written", "moved"):
            if step[key] > LIMITS["buffer_multiple"] * max(step["stars"], 1):
                problems.append(f"{step['step']}: {key} {step[key]} stars against a buffer of {step['stars']}")
    assert min(step["dist"] for step in steps) <= PAN_AT_PC, "the dive never came in to the pan's zoom"
    assert any(step["step"] == "pan 1" for step in steps), "the flight across never ran"
    assert steps[-1]["step"].startswith("up"), steps[-1]
    assert not problems, "\n".join(problems)
