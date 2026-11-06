# planetgen/galaxy/avoid.py

"""
The course planner (NAV.26)
===========================

A straight course that enters an object's keep-out sphere (NAV.54, `galaxy.keepout`) is bent around it. This module is the geometry only:
numpy, no storage, no generation. It takes a start, an end and the spheres a straight course meets (NAV.25, `db.obstacles`) and returns the
`adjusted` object of `docs/design/orbits/course-avoidance.md` section 7, the short detour around them. It works in any one unit (parsecs for
the galaxy, light-years in a sector, km in a system), as long as every number is in it.

The method is section 5 of that document:

1. A sphere that holds the start or the end is dropped (a planet lies inside its star's Hill sphere). If one sphere holds both, the course
   is the straight line.
2. A straight line that misses every other sphere is the course.
3. Otherwise a wrap is inserted at the first leg that still crosses a sphere, from the closed form for one sphere (the shortest path
   around one sphere lies in the plane of the two ends and its center: a tangent line, an arc of a great circle, a tangent line).
4. A validated block-coordinate polish then recomputes each wrap from its two neighbors, keeping a change only when the legs and the
   arc stay clear of every sphere and the length does not rise.
5. A sphere smaller than `tuning.COURSE_BELOW_RESOLUTION` of the straight length is still stepped around, with two tangent points and no
   arc samples; it is listed apart and not drawn.

Overlapping spheres break the closed form (the arc of one may enter another). When the plan crosses a sphere, the connected set of
overlapping spheres is replaced by its smallest enclosing sphere and the course is planned again; if that sphere swallows an end, or the
plan is still wrong, the course is the straight line with the warning `overlap_fallback`. A course that needs more than
`tuning.COURSE_MAX_WRAPS` wraps is the straight line with the warning `wrap_limit`. The result is a short detour, not a proven shortest
path.

Rounding: a hit is tested at `R (1 - HIT_MARGIN)` and a waypoint sits at `R (1 + WAYPOINT_MARGIN)`, so a tangent leg never fails on
rounding. A sphere holds a point when the center is no farther than `R (1 + CONTAIN_MARGIN)`.
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.util import detmath

HIT_MARGIN = 1e-7
"""float: A leg hits a sphere when it comes nearer the center than `R * (1 - HIT_MARGIN)`."""

WAYPOINT_MARGIN = 1e-6
"""float: A waypoint on a sphere sits at `R * (1 + WAYPOINT_MARGIN)` from its center."""

CONTAIN_MARGIN = 1e-9
"""float: A sphere contains a point when the center is no farther than `R * (1 + CONTAIN_MARGIN)` (the closed test of section 2.1)."""

KEEPOUT_RULESET = 1
PLANNER_RULESET = 1
"""int: The versions written into `adjusted.ruleset`, so a saved course made under another rule can be told from this one."""

_MAX_INSERTS = 400
_POLISH_SWEEPS = 60
_ROUNDS = 40
_MERGE_ROUNDS = 4


# ---------------------------------------------------------------------------------------------------------------------------------
# One leg against spheres
# ---------------------------------------------------------------------------------------------------------------------------------

def _closest(p, q, centers):
    d = q - p
    a = detmath.dot(d, d)
    f = centers - p
    tu = detmath.dot(f, d) / a
    t = np.clip(tu, 0.0, 1.0)
    dist = np.linalg.norm(f - np.outer(t, d), axis=1)
    return dist, tu, math.sqrt(a)


def _empty(p, q):
    return detmath.dot(q - p, q - p) == 0.0


def seg_hits(p, q, centers, radii):
    """Which spheres the segment p-q enters (an open ball: grazing is not a hit), as a boolean array."""
    if _empty(p, q) or len(radii) == 0:
        return np.zeros(len(radii), dtype=bool)
    dist, _tu, _dl = _closest(p, q, centers)
    return dist < radii


def seg_first_hit(p, q, centers, radii):
    """The index of the sphere the segment p-q enters first, or -1."""
    if _empty(p, q) or len(radii) == 0:
        return -1
    dist, tu, dl = _closest(p, q, centers)
    hit = dist < radii
    if not hit.any():
        return -1
    half = np.sqrt(np.maximum(radii * radii - dist * dist, 0.0)) / dl
    t_in = np.where(hit, np.maximum(tu - half, 0.0), np.inf)
    return int(np.argmin(t_in))


def tangent_wrap(a, b, c, radius):
    """
    The shortest path from `a` to `b` around the one sphere (`c`, `radius`), in the plane of the three points. `None` when the segment
    clears the sphere (or an end lies inside it), else `(entry, exit, arc_angle)`: the two tangent points and the angle between them at `c`.
    """
    u, v = a - c, b - c
    da, db = float(np.linalg.norm(u)), float(np.linalg.norm(v))
    if da < radius * (1 - CONTAIN_MARGIN) or db < radius * (1 - CONTAIN_MARGIN):
        return None
    e1 = u / da
    w = v - detmath.dot(v, e1) * e1
    wn = float(np.linalg.norm(w))
    if wn < 1e-12 * db:                       # collinear: any perpendicular gives the same length
        k = int(np.argmin(np.abs(e1)))
        t = np.zeros(3)
        t[k] = 1.0
        w = t - detmath.dot(t, e1) * e1
        wn = float(np.linalg.norm(w))
    e2 = w / wn
    phi = math.atan2(detmath.dot(v, e2), detmath.dot(v, e1))
    alpha_a = math.acos(min(1.0, radius / da))
    alpha_b = math.acos(min(1.0, radius / db))
    if alpha_a + alpha_b >= phi:
        return None
    ang_e, ang_x = alpha_a, phi - alpha_b
    entry = c + radius * (math.cos(ang_e) * e1 + math.sin(ang_e) * e2)
    exit_ = c + radius * (math.cos(ang_x) * e1 + math.sin(ang_x) * e2)
    return entry, exit_, ang_x - ang_e


def _arc_angle(c, entry, exit_):
    u, v = entry - c, exit_ - c
    nu, nv = float(np.linalg.norm(u)), float(np.linalg.norm(v))
    if nu == 0.0 or nv == 0.0:
        return 0.0
    return math.acos(max(-1.0, min(1.0, detmath.dot(u, v) / (nu * nv))))


def arc_points(c, entry, exit_, step):
    """Points along the shorter great-circle arc from `entry` to `exit_` about `c`, at most `step` radians apart, both ends included."""
    angle = _arc_angle(c, entry, exit_)
    if angle < 1e-12:
        return [entry, exit_]
    n = max(1, int(math.ceil(angle / step)))
    u, v = entry - c, exit_ - c
    sin_total = math.sin(angle)
    points = []
    for k in range(n + 1):
        f = k / n
        if sin_total < 1e-9:
            points.append(c + (1.0 - f) * u + f * v)
        else:
            points.append(c + (math.sin((1.0 - f) * angle) * u + math.sin(f * angle) * v) / sin_total)
    return points


# ---------------------------------------------------------------------------------------------------------------------------------
# A chain of points and wraps
# ---------------------------------------------------------------------------------------------------------------------------------
# An item is ("pt", p) or ("wrap", sphere_index, entry, exit).

def _end_of(item):
    return item[1] if item[0] == "pt" else item[3]


def _start_of(item):
    return item[1] if item[0] == "pt" else item[2]


def _wrap_arc_length(item, centers, radii):
    if item[0] != "wrap":
        return 0.0
    i = item[1]
    return radii[i] * (1.0 + WAYPOINT_MARGIN) * _arc_angle(centers[i], item[2], item[3])


def chain_length(items, centers, radii):
    """The length of a chain: its legs plus each wrap's arc."""
    total = 0.0
    for k in range(len(items) - 1):
        total += float(np.linalg.norm(_start_of(items[k + 1]) - _end_of(items[k])))
    for item in items:
        total += _wrap_arc_length(item, centers, radii)
    return total


def _pieces_length(items, k, centers, radii):
    """The length of what item `k` touches: its two legs and its arc."""
    return (float(np.linalg.norm(_start_of(items[k]) - _end_of(items[k - 1])))
            + float(np.linalg.norm(_start_of(items[k + 1]) - _end_of(items[k])))
            + _wrap_arc_length(items[k], centers, radii))


def violations(items, centers, radii):
    """
    The pieces of a chain that enter a sphere, as a set of `(own, other)` pairs (`own` is the wrap's sphere for an arc, `-1` for a leg).
    An empty set means the chain is clear.
    """
    hit_radii = radii * (1.0 - HIT_MARGIN)
    bad = set()
    for k in range(len(items) - 1):
        for j in np.nonzero(seg_hits(_end_of(items[k]), _start_of(items[k + 1]), centers, hit_radii))[0]:
            bad.add((-1, int(j)))
    for item in items:
        if item[0] != "wrap":
            continue
        i = item[1]
        for p in arc_points(centers[i], item[2], item[3], tuning.COURSE_ARC_SAMPLE_RAD)[1:-1]:
            inside = np.linalg.norm(centers - p, axis=1) < hit_radii
            inside[i] = False
            for j in np.nonzero(inside)[0]:
                bad.add((i, int(j)))
    return bad


def _wrap_clear(a, item, b, centers, radii):
    """Whether a wrap, with the legs from `a` into it and from it to `b`, enters no sphere but its own (the legs touch it only at a tangent)."""
    hit_radii = radii * (1.0 - HIT_MARGIN)
    if seg_hits(a, item[2], centers, hit_radii).any() or seg_hits(item[3], b, centers, hit_radii).any():
        return False
    i = item[1]
    for p in arc_points(centers[i], item[2], item[3], tuning.COURSE_ARC_SAMPLE_RAD)[1:-1]:
        inside = np.linalg.norm(centers - p, axis=1) < hit_radii
        inside[i] = False
        if inside.any():
            return False
    return True


def repair(items, centers, radii):
    """Inserts a wrap at the first leg that still crosses a sphere until no leg does (or `_MAX_INSERTS` wraps are in)."""
    hit_radii = radii * (1.0 - HIT_MARGIN)
    inflated = radii * (1.0 + WAYPOINT_MARGIN)
    for _ in range(_MAX_INSERTS):
        for k in range(len(items) - 1):
            a, b = _end_of(items[k]), _start_of(items[k + 1])
            j = seg_first_hit(a, b, centers, hit_radii)
            if j < 0:
                continue
            wrap = tangent_wrap(a, b, centers[j], inflated[j])
            if wrap is None:
                continue
            items.insert(k + 1, ("wrap", j, wrap[0], wrap[1]))
            break
        else:
            return items
    return items


def _merge_adjacent(items, centers, radii):
    """Joins two wraps in a row on the same sphere into one when the joined wrap is clear and no longer."""
    inflated = radii * (1.0 + WAYPOINT_MARGIN)
    k = 1
    while k < len(items) - 2:
        first, second = items[k], items[k + 1]
        if first[0] == "wrap" and second[0] == "wrap" and first[1] == second[1]:
            i = first[1]
            a, b = _end_of(items[k - 1]), _start_of(items[k + 2])
            wrap = tangent_wrap(a, b, centers[i], inflated[i])
            if wrap is not None:
                cand = ("wrap", i, wrap[0], wrap[1])
                before = chain_length(items[k - 1:k + 3], centers, radii)
                after = chain_length([items[k - 1], cand, items[k + 2]], centers, radii)
                if _wrap_clear(a, cand, b, centers, radii) and after <= before + 1e-12 * max(before, 1.0):
                    items[k:k + 2] = [cand]
                    continue
        k += 1
    return items


def polish(items, centers, radii, sweeps=_POLISH_SWEEPS):
    """
    Block-coordinate descent: each wrap is recomputed exactly from its two neighbors and kept only when clear of every sphere and no
    longer; a wrap whose neighbors now see each other past its sphere is dropped when the direct leg is clear of every sphere.
    """
    hit_radii = radii * (1.0 - HIT_MARGIN)
    inflated = radii * (1.0 + WAYPOINT_MARGIN)
    for _ in range(sweeps):
        before = chain_length(items, centers, radii)
        k = 1
        while k < len(items) - 1:
            if items[k][0] != "wrap":
                k += 1
                continue
            i = items[k][1]
            a, b = _end_of(items[k - 1]), _start_of(items[k + 1])
            wrap = tangent_wrap(a, b, centers[i], inflated[i])
            if wrap is None:
                if seg_first_hit(a, b, centers, hit_radii) < 0:
                    del items[k]
                    continue
            else:
                cand = ("wrap", i, wrap[0], wrap[1])
                old = _pieces_length(items, k, centers, radii)
                trial = items[:k] + [cand] + items[k + 1:]
                if _pieces_length(trial, k, centers, radii) <= old + 1e-12 * max(old, 1.0) and _wrap_clear(a, cand, b, centers, radii):
                    items[k] = cand
            k += 1
        _merge_adjacent(items, centers, radii)
        after = chain_length(items, centers, radii)
        if before - after < 1e-12 * max(before, 1.0):
            break
    return items


def repair_polish(items, centers, radii):
    """Alternates the polish and the repair until the chain stops growing."""
    for _ in range(_ROUNDS):
        items = polish(items, centers, radii)
        grown = repair(list(items), centers, radii)
        if len(grown) == len(items):
            return items
        items = grown
    return polish(items, centers, radii)


def plan_chain(start, end, centers, radii):
    """The chain of points and wraps from `start` to `end` around the spheres, clear of them or not (the caller checks `violations`)."""
    items = [("pt", start), ("pt", end)]
    items = repair(items, centers, radii)
    if any(item[0] == "wrap" for item in items):
        items = repair_polish(items, centers, radii)
    return items


# ---------------------------------------------------------------------------------------------------------------------------------
# Overlapping spheres
# ---------------------------------------------------------------------------------------------------------------------------------

def enclosing_sphere(centers, radii):
    """
    The smallest sphere holding the given spheres: exact for one and two, and for more a refinement (Badoiu and Clarkson) that is then
    grown to hold every one. Returns `(center, radius)`.
    """
    centers = np.asarray(centers, dtype=float)
    radii = np.asarray(radii, dtype=float)
    if len(radii) == 1:
        return centers[0].copy(), float(radii[0])
    if len(radii) == 2:
        d = float(np.linalg.norm(centers[1] - centers[0]))
        big, small = (0, 1) if radii[0] >= radii[1] else (1, 0)
        if d + radii[small] <= radii[big]:
            return centers[big].copy(), float(radii[big])
        half = (d + radii[0] + radii[1]) / 2.0
        direction = (centers[1] - centers[0]) / d
        return centers[0] + direction * (half - radii[0]), half
    center = (centers.min(axis=0) + centers.max(axis=0)) / 2.0
    for k in range(1, 400):
        far = int(np.argmax(np.linalg.norm(centers - center, axis=1) + radii))
        gap = centers[far] - center
        norm = float(np.linalg.norm(gap))
        if norm > 0.0:
            center = center + gap / (k + 1)
    radius = float(np.max(np.linalg.norm(centers - center, axis=1) + radii))
    return center, radius


def _clusters(centers, radii):
    """The connected sets of overlapping spheres, as lists of indexes (singletons left out)."""
    n = len(radii)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        gap = np.linalg.norm(centers[i + 1:] - centers[i], axis=1)
        for off in np.nonzero(gap < radii[i] + radii[i + 1:])[0]:
            ra, rb = find(i), find(i + 1 + int(off))
            if ra != rb:
                parent[ra] = rb
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return [members for members in groups.values() if len(members) > 1]


def merge_overlaps(centers, radii, refs, wanted):
    """
    Replaces each cluster of overlapping spheres that holds an index in `wanted` by its enclosing sphere. Returns the new
    `(centers, radii, refs)`, where a merged sphere's ref is the list of its members' refs.
    """
    drop = set()
    merged = []
    for members in _clusters(centers, radii):
        if not wanted.intersection(members):
            continue
        drop.update(members)
        c, r = enclosing_sphere(centers[members], radii[members])
        flat = []
        for m in members:
            flat.extend(refs[m] if isinstance(refs[m], list) else [refs[m]])
        merged.append((c, r, flat))
    keep = [i for i in range(len(radii)) if i not in drop]
    new_centers = [centers[i] for i in keep] + [m[0] for m in merged]
    new_radii = [radii[i] for i in keep] + [m[1] for m in merged]
    new_refs = [refs[i] for i in keep] + [m[2] for m in merged]
    return np.array(new_centers).reshape(-1, 3), np.array(new_radii, dtype=float), new_refs


# ---------------------------------------------------------------------------------------------------------------------------------
# The course
# ---------------------------------------------------------------------------------------------------------------------------------

def _as_list(ref):
    return list(ref) if isinstance(ref, list) else [ref]


def _straight(start, end, frame, unit, warnings, ignored=(), passed_through=None):
    length = float(np.linalg.norm(end - start))
    return {
        "frame": frame, "unit": unit,
        "length": length, "straight_length": length, "extra": 0.0,
        "waypoints": [{"p": start.tolist(), "kind": "start", "ref": None}, {"p": end.tolist(), "kind": "end", "ref": None}],
        "legs": [{"type": "line", "from": 0, "to": 1}],
        "polyline": [start.tolist(), end.tolist()],
        "avoided": [], "passed_through": list(passed_through or []),
        "ignored_below_resolution": list(ignored),
        "ruleset": {"keepout": KEEPOUT_RULESET, "planner": PLANNER_RULESET},
        "warnings": list(warnings),
    }


def _describe(items, start, end, centers, radii, refs, frame, unit, small, passed_through):
    waypoints = [{"p": start.tolist(), "kind": "start", "ref": None}]
    legs = []
    polyline = [start.tolist()]
    avoided = []
    ignored = []
    arc_step = 2.0 * math.acos(1.0 - tuning.COURSE_POLYLINE_SAGITTA)
    length = 0.0
    previous = start
    for item in items[1:-1]:
        i = item[1]
        label = _as_list(refs[i])
        ref = label[0] if len(label) == 1 else label
        radius = float(radii[i]) * (1.0 + WAYPOINT_MARGIN)
        entry, exit_ = item[2], item[3]
        waypoints.append({"p": entry.tolist(), "kind": "tangent_in", "ref": ref})
        legs.append({"type": "line", "from": len(waypoints) - 2, "to": len(waypoints) - 1})
        waypoints.append({"p": exit_.tolist(), "kind": "tangent_out", "ref": ref})
        legs.append({"type": "arc", "from": len(waypoints) - 2, "to": len(waypoints) - 1,
                     "center": centers[i].tolist(), "radius": radius, "obstacle": ref})
        length += float(np.linalg.norm(entry - previous)) + radius * _arc_angle(centers[i], entry, exit_)
        previous = exit_
        polyline.append(entry.tolist())
        if i in small:
            ignored.extend(label)
        else:
            polyline.extend(p.tolist() for p in arc_points(centers[i], entry, exit_, arc_step)[1:-1])
            avoided.extend(label)
        polyline.append(exit_.tolist())
    waypoints.append({"p": end.tolist(), "kind": "end", "ref": None})
    legs.append({"type": "line", "from": len(waypoints) - 2, "to": len(waypoints) - 1})
    polyline.append(end.tolist())
    length += float(np.linalg.norm(end - previous))
    straight = float(np.linalg.norm(end - start))
    return {
        "frame": frame, "unit": unit,
        "length": length, "straight_length": straight, "extra": length - straight,
        "waypoints": waypoints, "legs": legs, "polyline": polyline,
        "avoided": avoided, "passed_through": list(passed_through or []),
        "ignored_below_resolution": ignored,
        "ruleset": {"keepout": KEEPOUT_RULESET, "planner": PLANNER_RULESET},
        "warnings": [],
    }


def plan_course(start, end, spheres, frame="galactic", unit="pc", passed_through=None):
    """
    The `adjusted` course (course-avoidance.md section 7) from `start` to `end`, a short detour around `spheres`.

    `spheres` is a list of dicts with `ref` (the object reference, such as `"star:101"`), `center` (three numbers) and `radius`, all in
    the same unit as `start` and `end`. `passed_through` is copied into the result as given (the notes on gas and dust crossed; they do
    not bend the course). Spheres of no size are left out.
    """
    start = np.asarray(start, dtype=float)
    end = np.asarray(end, dtype=float)
    spheres = [s for s in spheres if float(s["radius"]) > 0.0]
    if not spheres or _empty(start, end):
        return _straight(start, end, frame, unit, [], passed_through=passed_through)
    centers = np.array([s["center"] for s in spheres], dtype=float).reshape(-1, 3)
    radii = np.array([s["radius"] for s in spheres], dtype=float)
    refs = [s["ref"] for s in spheres]
    straight = float(np.linalg.norm(end - start))

    near_start = np.linalg.norm(centers - start, axis=1) <= radii * (1.0 + CONTAIN_MARGIN)
    near_end = np.linalg.norm(centers - end, axis=1) <= radii * (1.0 + CONTAIN_MARGIN)
    if (near_start & near_end).any():
        return _straight(start, end, frame, unit, [], passed_through=passed_through)
    keep = ~(near_start | near_end)
    centers, radii, refs = centers[keep], radii[keep], [r for r, k in zip(refs, keep) if k]
    if len(radii) == 0:
        return _straight(start, end, frame, unit, [], passed_through=passed_through)

    items = plan_chain(start, end, centers, radii)
    bad = violations(items, centers, radii)
    for _ in range(_MERGE_ROUNDS):
        if not bad:
            break
        wanted = set()
        for own, other in bad:
            wanted.add(other)
            if own >= 0:
                wanted.add(own)
        if not any(_clusters(centers, radii)):
            break
        centers, radii, refs = merge_overlaps(centers, radii, refs, wanted)
        near_start = np.linalg.norm(centers - start, axis=1) <= radii * (1.0 + CONTAIN_MARGIN)
        near_end = np.linalg.norm(centers - end, axis=1) <= radii * (1.0 + CONTAIN_MARGIN)
        if (near_start | near_end).any():
            return _straight(start, end, frame, unit, ["overlap_fallback"], passed_through=passed_through)
        items = plan_chain(start, end, centers, radii)
        bad = violations(items, centers, radii)
    if bad:
        return _straight(start, end, frame, unit, ["overlap_fallback"], passed_through=passed_through)
    if sum(1 for item in items if item[0] == "wrap") > tuning.COURSE_MAX_WRAPS:
        return _straight(start, end, frame, unit, ["wrap_limit"], passed_through=passed_through)

    small = {i for i in range(len(radii)) if radii[i] < tuning.COURSE_BELOW_RESOLUTION * straight}
    return _describe(items, start, end, centers, radii, refs, frame, unit, small, passed_through)


def spheres_from_obstacles(obstacles):
    """The `spheres` argument from the dicts `db.obstacles.obstacles_along_segment` returns (parsecs)."""
    return [{"ref": o["ref"], "center": o["center_pc"], "radius": o["radius_pc"]} for o in obstacles]
