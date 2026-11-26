# planetgen/web/lib/life_air_chart.py

"""
The "Life and air" chart of a planet's or moon's page (UX.237, docs/design/habitability/life-classes.md section 6.2): the partial
pressure of oxygen, carbon dioxide and methane against the age of the star, as an SVG line chart with the stages of the life
marked, and the same figures as a table. The data is `GET /api/planets/<id>/life-air` (`physics.life_air.air_history`). Imports no Flask.
"""

import math
from html import escape

from markupsafe import Markup

from planetgen.generation import evolution
from planetgen.physics import air_events, life_ages
from planetgen.web.lib import reference_urls

WIDTH, HEIGHT = 640, 280
LEFT, RIGHT, TOP, BOTTOM = 64, 14, 24, 40
FLOOR_KPA = 1.0e-4
"""float: The bottom of the pressure axis, kPa; a gas at zero is drawn on it."""

SERIES = {"O2": "life-air-o2", "CO2": "life-air-co2", "CH4": "life-air-ch4"}
"""dict: The CSS class of each gas's line (a colour and a dash pattern, so the lines differ without colour)."""
OTHER_STYLES = 6
"""int: How many styles (`life-air-g0` to `life-air-g5`) the other gases of a stored history (API.51) take in turn, each a colour with a dash pattern of its own (a drawing choice, not a physical constant)."""


def _styles(gases):
    """`[(gas, css class)]`: the three gases of the first chart keep their own lines, every other gas takes the next of `OTHER_STYLES`, in the order the history lists them (the biggest first)."""
    styled, other = [], 0
    for gas in gases:
        if gas["formula"] in SERIES:
            styled.append((gas, SERIES[gas["formula"]]))
        else:
            styled.append((gas, f"life-air-g{other % OTHER_STYLES}"))
            other += 1
    return styled


EXPLAIN = (
    "Each line running across the chart is one gas in the air (the key under the chart says which). Its height is the pressure that gas alone makes, in kilopascals (kPa); for scale, the whole air of the Earth today is about 101 kPa, and its oxygen makes about 21 kPa of that.",
    "Time runs from the birth of the star on the left to now on the right, in billions of years (Gyr). Each gridline up is ten times the pressure of the one below it, and a gas on the bottom line is not in the air at all.",
    "The numbered vertical lines are the stages the life reached, listed under the chart. The dashed vertical lines without a number are mass extinctions: a flood of lava or a large impact that killed much of the life, and often changed the air.",
)
"""tuple: The plain sentences the page gives with the chart (UX.322): what a line is, what the axes are, and what the vertical lines are."""

CAUSES = {air_events.FLOOD_BASALT: "Flood basalt eruptions", air_events.IMPACT: "A large impact"}
"""dict: The words for what ended a life age (`air_events.FLOOD_BASALT`, `air_events.IMPACT`)."""


def _kpa_text(value):
    if value <= 0.0:
        return "none"
    return f"{value:.3g} kPa" if value >= 0.01 else f"{value * 1000.0:.3g} Pa"


def _axis_ticks(top):
    """The ages, Gyr, of the x axis ticks: about five, at a round step."""
    step = 10 ** math.floor(math.log10(top / 4.0))
    for factor in (1, 2, 5, 10):
        if top / (step * factor) <= 6:
            step *= factor
            break
    count = int(top / step)
    return [round(index * step, 6) for index in range(count + 1)]


def _class_name(entry):
    return "Class " + str(entry["life_class"]) + (entry.get("variant") or "")


def _makers(history):
    """What each class of a world with several living classes made of the air (GEN.1213): `[{"gas", "text"}]` with "Class 1a 70 percent, Class 3a 30 percent" for
    each gas that more than one class changed, from `air_history`'s `by_class`; empty when one class lives there."""
    if len(history.get("classes") or []) < 2:
        return []
    makers = []
    for gas in history["gases"]:
        by_class = gas.get("by_class") or []
        if len(by_class) > 1:
            parts = [f"{_class_name(entry)} {round(100.0 * entry['share'])} percent" for entry in sorted(by_class, key=lambda entry: -entry["share"])]
            makers.append({"gas": gas["name"], "text": ", ".join(parts)})
    return makers


def chart(history, class_id, grade=0):
    """
    The view of a chart for `history` (`air_history`'s dict, or `None`): `{"svg", "table", "stages", "note", "makers", "title"}`, or `None` when there is
    nothing to draw. `svg` and `note` are HTML (everything in them escaped), `table` is `{"header": [...], "rows": [[...]]}` with one row for
    each span of ages in which the air differs, and `stages` is a list of `{"number", "name", "age"}`.
    """
    if not history or not history.get("gases"):
        return None
    top = float(history["age_gyr"])
    values = [max(value, FLOOR_KPA) for gas in history["gases"] for _age, value in gas["steps"]]
    low, high = math.floor(math.log10(min(values))), math.ceil(math.log10(max(values)))
    if high <= low:
        high = low + 1
    plot_w, plot_h = WIDTH - LEFT - RIGHT, HEIGHT - TOP - BOTTOM

    def x(age):
        return LEFT + plot_w * min(max(age, 0.0), top) / top

    def y(value):
        return TOP + plot_h * (1.0 - (math.log10(max(value, FLOOR_KPA)) - low) / (high - low))

    parts = [f'<svg class="life-air-chart" viewBox="0 0 {WIDTH} {HEIGHT}" role="img" aria-labelledby="life-air-title life-air-desc">',
             '<title id="life-air-title">Life and air through the ages</title>']
    gas_names = ", ".join(gas["name"].lower() for gas in history["gases"])
    parts.append(f'<desc id="life-air-desc">The partial pressure of {escape(gas_names)} from the birth of the star to now, on a logarithmic axis; '
                 'the same figures are in the table below.</desc>')
    for exponent in range(low, high + 1):
        row = y(10.0 ** exponent)
        label = "none" if 10.0 ** exponent <= FLOOR_KPA * 1.0001 else f"{10.0 ** exponent:g}"
        parts.append(f'<line class="life-air-grid" x1="{LEFT}" x2="{WIDTH - RIGHT}" y1="{row:.1f}" y2="{row:.1f}"/>')
        parts.append(f'<text class="life-air-label" x="{LEFT - 6}" y="{row + 4:.1f}" text-anchor="end">{escape(label)}</text>')
    for age in _axis_ticks(top):
        column = x(age)
        parts.append(f'<line class="life-air-grid" x1="{column:.1f}" x2="{column:.1f}" y1="{TOP}" y2="{HEIGHT - BOTTOM}"/>')
        parts.append(f'<text class="life-air-label" x="{column:.1f}" y="{HEIGHT - BOTTOM + 16}" text-anchor="middle">{age:g}</text>')
    parts.append(f'<text class="life-air-label" x="{LEFT + plot_w / 2:.1f}" y="{HEIGHT - 6}" text-anchor="middle">Time since the star formed (Gyr, billions of years)</text>')
    parts.append(f'<text class="life-air-label" transform="translate(14 {TOP + plot_h / 2:.1f}) rotate(-90)" text-anchor="middle">Pressure of the gas (kPa)</text>')
    stages = []
    for number, stage in enumerate(history["stages"], 1):
        column = x(stage["age_gyr"])
        name = evolution.stage_name(stage["key"], class_id, grade if stage["key"] == "multicellularity" else 0)
        stages.append({"number": number, "name": name, "age": f"{stage['age_gyr']:.2f}"})
        parts.append(f'<line class="life-air-stage" x1="{column:.1f}" x2="{column:.1f}" y1="{TOP}" y2="{HEIGHT - BOTTOM}"/>')
        parts.append(f'<text class="life-air-stage-number" x="{column:.1f}" y="{TOP - 7}" text-anchor="middle">{number}</text>')
    extinctions = []
    for event in history.get("events") or []:
        if event["kill_share"] < life_ages.AGE_KILL_THRESHOLD:
            continue      # UX.321: a shock that killed less than the share that ends an age of life is no mass extinction, and the page lists and marks only those
        column = x(event["age_gyr"])
        extinctions.append({"age": f"{event['age_gyr']:.2f}", "cause": CAUSES.get(event["cause"], str(event["cause"])), "killed": f"{round(100.0 * event['kill_share'])} percent"})
        parts.append(f'<line class="life-air-event" x1="{column:.1f}" x2="{column:.1f}" y1="{TOP}" y2="{HEIGHT - BOTTOM}"/>')
    styled = _styles(history["gases"])
    for gas, css in styled:
        points, previous = [], None
        for age, value in gas["steps"]:
            if previous is not None:
                points.append((x(age), y(previous)))
            points.append((x(age), y(value)))
            previous = value
        points.append((x(top), y(previous)))
        path = " ".join(f"{px:.1f},{py:.1f}" for px, py in points)
        parts.append(f'<polyline class="life-air-line {css}" points="{path}"/>')
    parts.append("</svg>")
    legend = "".join(f'<li><span class="life-air-key {css}"></span> {escape(gas["name"])}</li>' for gas, css in styled)
    spans = sorted({age for gas in history["gases"] for age, _ in gas["steps"]} | {top})
    rows = []
    for start, end in zip(spans, spans[1:]):
        row = [f"{start:.2f} to {end:.2f}"]
        for gas in history["gases"]:
            current = [value for age, value in gas["steps"] if age <= start][-1]
            row.append(_kpa_text(current))
        rows.append(row)
    makers = _makers(history)
    note = ""
    if history.get("speculative"):
        note = Markup('This class\'s air numbers are less sure than the rest. <a class="reference-link" href="{}">More information</a>').format(
            reference_urls.reference_path("speculative", class_id))
    ages = [{"number": age["number"], "span": f"{age['start_gyr']:.2f} to " + ("now" if age["end_gyr"] is None else f"{age['end_gyr']:.2f}"),
             "ended": "" if age["cause"] is None or age["end_gyr"] is None else f"the mass extinction at {age['end_gyr']:.2f} Gyr"}
            for age in history.get("life_ages") or []]
    marks = ""
    if stages:
        marks += '<li><span class="life-air-key life-air-stage-key"></span> A stage the life reached (numbered, listed below)</li>'
    if extinctions:
        marks += '<li><span class="life-air-key life-air-event-key"></span> A mass extinction</li>'
    return {"svg": "".join(parts), "legend": legend, "marks": marks, "explain": EXPLAIN, "stages": stages, "note": note, "makers": makers, "stored": history.get("source") == "stored",
            "extinctions": extinctions, "extinction_floor": f"{round(100.0 * life_ages.AGE_KILL_THRESHOLD)} percent", "life_ages": ages,
            "table": {"header": ["Age (Gyr)"] + [gas["name"] for gas in history["gases"]], "rows": rows}}
