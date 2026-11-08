# tests/audit/prevalence_oddities.py

"""
The oddity census (TEST.492, part of the prevalence suite): the share of the bodies that can carry each row of the oddity table (`physics/data/oddities.json`) that
do carry it, from the same seeded systems as the system census (`prevalence_systems`, whose record holds each system's `oddities` counts).

**What is counted.** Every planet and moon of every system, once: it is *gated* for a row when the row's kind is the body's and the row's gate holds on the finished body
(`oddities.candidates` over `oddities.quantities_of`), and *carried* when the row's id is in `body.oddities`. A body that carries a row whose gate its later
life no longer meets (a tide that moved it, a branch it took) counts as gated too, so a share is never over 1.

**The range.** A row is a rate on a gated population (design 4.1), so the share that carry it is the row's own probability (`oddities.probability`, the rate times the
dial) within the binomial band of the gated count (4 sigma) and a pad of 0.05; the upper bound is the whole story for a row that can be declined or whose gate the body
leaves (`W2`, `W19` and `W20` decline the bodies their effect cannot act on), whose lower bound is 0. A row with fewer than `MIN_GATED` gated bodies in the run is *reported*, not judged
(the sample cannot say), and a row expected `MUST_APPEAR` times or more that never appears is a collapse whatever the band says.
"""

import math

from planetgen.physics import oddities
from tests.audit import prevalence as census

REPORT_NAME = "oddity-census-latest.md"
SIGMAS = 4.0
PAD_SHARE = 0.05
MIN_GATED = 20
MUST_APPEAR = 5
"""int: A share expected this many times or more in the run must appear at least once, whatever its band."""

SHARE_CAN_FALL_SHORT = frozenset({"W2", "W19", "W20"})
"""frozenset: Rows whose effect can decline a body the gate passes (the Venus-type spin wants a planet the tide leaves turning backwards; the nightside collapse wants a locked
planet whose nightside freezes its gas out; the other climate branch wants a bistable window), so the share carried is under the rate and only the upper bound judges it."""


def record_of(system):
    """`{row id: [gated bodies, carried bodies]}` over the planets and moons of `system`, for the rows that gate at least one of them."""
    counts = {}
    for planet in system.all_planets:
        if planet.body_type == "a":
            continue
        for body in (planet, *planet.moons):
            gated = {row.id for row in oddities.candidates(oddities.kind_of(body), oddities.quantities_of(body))}
            carried = set(getattr(body, "oddities", ()) or ())
            for row_id in gated | carried:
                entry = counts.setdefault(row_id, [0, 0])
                entry[0] += 1
                entry[1] += 1 if row_id in carried else 0
    return counts


def summarize(records):
    """`{row id: {"gated": n, "carried": n, "share": carried / gated}}` over the records that finished (`oddities` is each record's `record_of`)."""
    total = {}
    for record in records:
        if record.get("abandoned") or record.get("error"):
            continue
        for row_id, (gated, carried) in record.get("oddities", {}).items():
            entry = total.setdefault(row_id, {"gated": 0, "carried": 0})
            entry["gated"] += gated
            entry["carried"] += carried
    for entry in total.values():
        entry["share"] = entry["carried"] / entry["gated"] if entry["gated"] else 0.0
    return total


def expected_range(row, gated):
    """`(low, high)`: the share of `gated` bodies that should carry `row`, from its probability."""
    p = oddities.probability(row)
    half = SIGMAS * math.sqrt(p * (1.0 - p) / max(gated, 1)) + PAD_SHARE
    low = 0.0 if row.id in SHARE_CAN_FALL_SHORT else max(0.0, p - half)
    return round(low, 4), round(min(1.0, p + half), 4)


def compare(summary, rows=None, known=None):
    """The flags of the run: a row's share out of its range, or a row that gates enough bodies and never appears. A row with too few gated bodies is not judged
    (`report_lines` lists it). `known` maps a row id to the item that explains its flag."""
    known = known or {}
    findings = []
    for row in (oddities.ROWS if rows is None else rows):
        mine = summary.get(row.id)
        if not mine or mine["gated"] < MIN_GATED:
            continue
        low, high = expected_range(row, mine["gated"])
        share = mine["share"]
        p = oddities.probability(row)
        if share < low or share > high:
            kind = "missing" if mine["carried"] == 0 and p * mine["gated"] >= MUST_APPEAR else "out-of-range"
            findings.append(census.Finding(kind, "oddities", row.id, share, low, high, f"{row.rate.source[:140]}", known.get(row.id, "")))
        elif mine["carried"] == 0 and p * mine["gated"] >= MUST_APPEAR and row.id not in SHARE_CAN_FALL_SHORT:
            findings.append(census.Finding("missing", "oddities", row.id, share, low, high, "expected to appear in the run", known.get(row.id, "")))
    return findings


def format_report(summary, findings, header="", rows=None):
    lines = ["# Oddity census (TEST.492), latest run", ""]
    if header:
        lines += [header, ""]
    unexplained = [f for f in findings if not f.known]
    explained = [f for f in findings if f.known]
    lines += [f"Flags: {len(unexplained)} new, {len(explained)} known anomalies.", ""]
    lines += [f"- {finding.line()}" for finding in unexplained]
    if explained:
        lines += ["", "Known anomalies:", ""]
        lines += [f"- {finding.line()}" for finding in explained]
    lines += ["", "## Share of the gated bodies that carry each row", ""]
    for row in (oddities.ROWS if rows is None else rows):
        mine = summary.get(row.id, {"gated": 0, "carried": 0, "share": 0.0})
        note = "" if mine["gated"] >= MIN_GATED else " (too few gated bodies to judge)"
        low, high = expected_range(row, mine["gated"])
        lines.append(f"- {row.id} {row.name}: {mine['carried']} of {mine['gated']} = {mine['share']:.3f}, expected {low:.3f} to {high:.3f}{note}")
    return "\n".join(lines) + "\n"
