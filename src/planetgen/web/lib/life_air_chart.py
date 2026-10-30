# planetgen/web/lib/life_air_chart.py

"""
The "Life and air" chart of a planet's or moon's page (UX.237, docs/design/habitability/life-classes.md section 6.2): the partial
pressure of oxygen, carbon dioxide and methane against the age of the star, as an SVG line chart with the stages of the life
marked, and the same figures as a table. The data is `GET /api/planets/<id>/life-air` (`physics.life_air.air_history`). Imports no Flask.
"""

import math
from html import escape

from planetgen.generation import evolution

WIDTH, HEIGHT = 640, 280
LEFT, RIGHT, TOP, BOTTOM = 64, 14, 24, 40
FLOOR_KPA = 1.0e-4
"""float: The bottom of the pressure axis, kPa; a gas at zero is drawn on it."""

SERIES = {"O2": "life-air-o2", "CO2": "life-air-co2", "CH4": "life-air-ch4"}
"""dict: The CSS class of each gas's line (a colour and a dash pattern, so the lines differ without colour)."""


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


def chart(history, class_id, grade=0):
    """
    The view of a chart for `history` (`air_history`'s dict, or `None`): `{"svg", "table", "stages", "note", "title"}`, or `None` when there is
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
    parts.append(f'<text class="life-air-label" x="{LEFT + plot_w / 2:.1f}" y="{HEIGHT - 6}" text-anchor="middle">Age of the star (Gyr)</text>')
    parts.append(f'<text class="life-air-label" transform="translate(14 {TOP + plot_h / 2:.1f}) rotate(-90)" text-anchor="middle">Partial pressure (kPa)</text>')
    stages = []
    for number, stage in enumerate(history["stages"], 1):
        column = x(stage["age_gyr"])
        name = evolution.stage_name(stage["key"], class_id, grade if stage["key"] == "multicellularity" else 0)
        stages.append({"number": number, "name": name, "age": f"{stage['age_gyr']:.2f}"})
        parts.append(f'<line class="life-air-stage" x1="{column:.1f}" x2="{column:.1f}" y1="{TOP}" y2="{HEIGHT - BOTTOM}"/>')
        parts.append(f'<text class="life-air-stage-number" x="{column:.1f}" y="{TOP - 7}" text-anchor="middle">{number}</text>')
    for gas in history["gases"]:
        points, previous = [], None
        for age, value in gas["steps"]:
            if previous is not None:
                points.append((x(age), y(previous)))
            points.append((x(age), y(value)))
            previous = value
        points.append((x(top), y(previous)))
        path = " ".join(f"{px:.1f},{py:.1f}" for px, py in points)
        parts.append(f'<polyline class="life-air-line {SERIES.get(gas["formula"], "life-air-other")}" points="{path}"/>')
    parts.append("</svg>")
    legend = "".join(f'<li><span class="life-air-key {SERIES.get(gas["formula"], "life-air-other")}"></span> {escape(gas["name"])}</li>'
                     for gas in history["gases"])
    spans = sorted({age for gas in history["gases"] for age, _ in gas["steps"]} | {top})
    rows = []
    for start, end in zip(spans, spans[1:]):
        row = [f"{start:.2f} to {end:.2f}"]
        for gas in history["gases"]:
            current = [value for age, value in gas["steps"] if age <= start][-1]
            row.append(_kpa_text(current))
        rows.append(row)
    note = ""
    if history.get("speculative"):
        note = escape(f"Speculative chemistry: this class's numbers are tier {history['speculative']['tier']}, tweak level "
                      f"{history['speculative']['tweak']} (stand-ins, not measured).")
    return {"svg": "".join(parts), "legend": legend, "stages": stages, "note": note,
            "table": {"header": ["Age (Gyr)"] + [gas["name"] for gas in history["gases"]], "rows": rows}}
