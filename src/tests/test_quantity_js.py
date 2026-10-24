# tests/test_quantity_js.py

"""
UX.128: the browser's Quantity (`html/static/quantity.js`) is the mirror of `planetgen.physics.quantity.Quantity`, and
`unitladder.js` shows a Quantity, or a value of a registry unit, the way Python's `format_measure` and its
`format_distance_km`-style names do (checked through Node, skipped without it). The browser sends no value through a
factor of its own: `formatDistanceKm(384400)` is `formatMeasure(held("km", 384400))`.
"""

import json
import math
import os
import shutil
import subprocess

import pytest

from planetgen.physics import unit_registry as registry
from planetgen.physics.quantity import Quantity
from planetgen.util import format as fmt
from planetgen.util import ladders

_STATIC = os.path.join(os.path.dirname(__file__), "..", "..", "src", "html", "static")

pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")

_VALUES = (0.0, 1.0, 1.5, 0.001, 3.7, 12345.678, 6.02e23, 1.2e-9, 273.15, 1e-30, 1e30, 384400.0, 0.378, 92.0)

_HELD = [  # (kind, unit) the site holds values in
    ("length", "km"), ("length", "au"), ("length", "ly"), ("length", "pc"), ("mass", "Msun"), ("mass", "Mearth"),
    ("duration", "yr"), ("duration", "d"), ("speed", "km_s"), ("temperature", "K"), ("temperature", "degC"),
    ("pressure", "atm"), ("acceleration", "g0"), ("power", "Lsun"), ("flux", "W_m2"), ("magnetic_field", "G"),
    ("number_density", "cm-3"),
]

_WRAPPERS = {  # the JavaScript name and the Python one it mirrors
    "formatDistanceKm": fmt.format_distance_km, "formatDistanceAu": fmt.format_distance_au,
    "formatDistanceLy": fmt.format_distance_ly, "formatDistancePc": fmt.format_distance_pc,
    "formatSpeedKms": fmt.format_speed_kms, "formatPeriodYears": fmt.format_period_years,
    "formatPressureAtm": fmt.format_pressure_atm, "formatGravityG": fmt.format_gravity_g,
    "formatMassSolar": fmt.format_mass_solar, "formatLuminositySol": fmt.format_luminosity_sol,
}


def _node(script, job):
    prelude = ("import fs from 'node:fs';"
               f"const q = await import('file://{os.path.abspath(os.path.join(_STATIC, 'quantity.js'))}');"
               f"const l = await import('file://{os.path.abspath(os.path.join(_STATIC, 'unitladder.js'))}');"
               "const job = JSON.parse(fs.readFileSync(0, 'utf8'));")
    run = subprocess.run(["node", "--input-type=module", "-e", prelude + script], input=json.dumps(job),
                         capture_output=True, text=True, check=True, timeout=120)
    return json.loads(run.stdout, parse_int=float)


def test_a_quantity_answers_every_unit_of_its_kind_as_the_python_one_does():
    cases = [[kind, value, unit, other.key] for kind, unit in _HELD for value in _VALUES
             for other in registry.units_of(kind)]
    answer = _node(
        "console.log(JSON.stringify(job.map(([k, v, u, t]) => { try { const x = new q.Quantity(k, v, u);"
        " return [x.as(t), x.base()]; } catch (e) { return null; } })));", cases)
    checked = 0
    for (kind, value, unit, target), got in zip(cases, answer):
        try:
            quantity = Quantity(kind, value, unit)
        except ValueError:
            assert got is None, (kind, value, unit)
            continue
        assert got == [quantity.as_(target), quantity.base()], (kind, value, unit, target)
        checked += 1
    assert checked > 1000


def test_a_quantity_refuses_what_the_python_one_refuses():
    cases = [("length", -1.0, "km"), ("length", math.inf, "km"), ("length", math.nan, "km"), ("speed", 3e8, "m_s"),
             ("temperature", -1.0, "K"), ("temperature", -300.0, "degC"), ("mass", -1.0, "kg"),
             ("length", 1.0, "kg"), ("nonsense", 1.0, "km"), ("length", 1.0, "nonsense")]
    answer = _node(
        "console.log(JSON.stringify(job.map(([k, v, u]) => { try { new q.Quantity(k, v === null ? Infinity : v, u);"
        " return 'accepted'; } catch (e) { return 'refused'; } })));",
        [[kind, None if not math.isfinite(value) else value, unit] for kind, value, unit in cases])
    for case, got in zip(cases, answer):
        try:
            Quantity(case[0], case[1], case[2])
            expected = "accepted"
        except ValueError:
            expected = "refused"
        assert got == expected, case


def test_set_changes_the_truth_and_every_form_follows():
    answer = _node(
        "const x = new q.Quantity('length', 1, 'km'); const before = x.as('m'); x.set(2, 'au');"
        " console.log(JSON.stringify([before, x.as('km'), x.unit, x.value]));", None)
    after = Quantity("length", 2, "au")
    assert answer == [1000.0, after.km, "au", 2.0]


def test_a_quantity_round_trips_through_json():
    answer = _node("const x = new q.Quantity('mass', 5.972e24, 'kg');"
                   " console.log(JSON.stringify(q.Quantity.fromJSON(JSON.parse(JSON.stringify(x))).toJSON()));", None)
    assert answer == {"kind": "mass", "value": 5.972e24, "unit": "kg"}


@pytest.mark.parametrize("kind, unit", _HELD, ids=[unit for _kind, unit in _HELD])
def test_format_measure_gives_the_text_python_gives(kind, unit):
    values = [value for value in _VALUES if not (kind == "temperature" and unit == "degC" and value < 0)]
    answer = _node("console.log(JSON.stringify(job.map(([k, v, u]) => { try { return l.formatMeasure(new q.Quantity(k, v, u)); }"
                   " catch (e) { return null; } })));", [[kind, value, unit] for value in values])
    for value, got in zip(values, answer):
        try:
            expected = fmt.format_measure(Quantity(kind, value, unit))
        except ValueError:
            expected = None
        assert got == expected, (kind, unit, value)


def test_the_format_names_give_the_text_python_gives_for_numbers_and_non_numbers():
    values = [None, *_VALUES, -5.0, math.inf, -math.inf, math.nan]
    encoded = [None if value is None else ("nan" if math.isnan(value) else "inf" if value == math.inf
                                          else "-inf" if value == -math.inf else value) for value in values]
    names = list(_WRAPPERS)
    answer = _node("const back = (v) => v === 'nan' ? NaN : v === 'inf' ? Infinity : v === '-inf' ? -Infinity : v;"
                   " console.log(JSON.stringify(job.names.map((n) => job.values.map((v) => l[n](back(v))))));",
                   {"names": names, "values": encoded})
    for name, row in zip(names, answer):
        for value, got in zip(values, row):
            assert got == _WRAPPERS[name](value), (name, value)


def test_a_ladder_is_found_for_every_registry_kind_python_shows():
    answer = _node("const d = await import('file://" + os.path.abspath(os.path.join(_STATIC, "unitladders.js")) + "');"
                   " console.log(JSON.stringify(d.LADDER_OF_KIND));", None)
    assert answer == ladders.LADDER_OF_KIND
