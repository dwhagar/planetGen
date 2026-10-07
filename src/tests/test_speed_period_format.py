# tests/test_speed_period_format.py

"""
The speed ladder (UX.13) and the time-period ladder (UX.14):
`stellarObjects.utils.format_speed_kms` shows a speed in one of km/h <
km/s < Mm/s < c, and `format_duration_seconds` / `format_period_years` a
period in the largest of µs < ms < s < minutes < hours < days < years < ky
< My < Gy it is at least 1 of. The browser copies (`html/static/speed.js`,
`html/static/period.js`) are checked against the same tables when Node is
installed.
"""

import json
import os
import re
import shutil
import subprocess

import pytest

from planetgen.physics import constants as pc
from stellarObjects.utils import (
    PERIOD_LADDER, SPEED_LADDER,
    format_duration_seconds, format_galactic_orbit, format_period_years, format_speed_kms, format_speed_ms,
)

_STATIC = os.path.join(os.path.dirname(__file__), "..", "html", "static")
_C = pc.SPEED_OF_LIGHT_KMS
_YEAR = pc.SECONDS_PER_YEAR

# (km/s, expected text): one value in each unit and each rung boundary.
SPEED_CASES = [
    (0.0, "0 km/h"),
    (0.001, "3.6 km/h"),
    (0.01, "36 km/h"),
    (0.5, "1,800 km/h"),
    (0.999, "3,596 km/h"),
    (1.0, "1 km/s"),
    (29.78, "29.8 km/s"),
    (205.7, "206 km/s"),
    (999.0, "999 km/s"),
    (1000.0, "1 Mm/s"),
    (4500.0, "4.5 Mm/s"),
    (0.0999 * _C, "29.9 Mm/s"),
    (0.1 * _C, "0.1 c"),
    (0.25 * _C, "0.25 c"),
    (_C, "1 c"),
    (3.4 * _C, "3.4 c"),
    (-29.78, "-29.8 km/s"),
]

# (seconds, expected text).
PERIOD_CASES = [
    (0.0, "0 s"),
    (3e-7, "0.3 µs"),
    (1e-6, "1 µs"),
    (9.47e-4, "947 µs"),
    (1.56e-3, "1.56 ms"),
    (0.5, "500 ms"),
    (1.0, "1 s"),
    (45.0, "45 s"),
    (60.0, "1 minute"),
    (90.0, "1.5 minutes"),
    (3600.0, "1 hour"),
    (12 * 3600.0, "12 hours"),
    (86400.0, "1 day"),
    (27.32 * 86400, "27.3 days"),
    (_YEAR, "1 year"),
    (1.88 * _YEAR, "1.88 years"),
    (999 * _YEAR, "999 years"),
    (1e3 * _YEAR, "1 ky"),
    (2.36e8 * _YEAR, "236 My"),
    (13.8e9 * _YEAR, "13.8 Gy"),
    (1.2e13 * _YEAR, "12,000 Gy"),
]


@pytest.mark.parametrize("kms,expected", SPEED_CASES)
def test_speed_ladder(kms, expected):
    assert format_speed_kms(kms) == expected


@pytest.mark.parametrize("seconds,expected", PERIOD_CASES)
def test_period_ladder(seconds, expected):
    assert format_duration_seconds(seconds) == expected


def test_speed_switches_to_c_at_a_tenth_of_light_speed():
    assert format_speed_kms(0.1 * _C * 0.999).endswith(" Mm/s")
    assert format_speed_kms(0.1 * _C).endswith(" c")


def test_unit_wrappers_and_none():
    assert format_speed_ms(29_780) == "29.8 km/s"
    assert format_speed_kms(None) == "–"
    assert format_speed_ms(None) == "–"
    assert format_period_years(1.0) == "1 year"
    assert format_period_years(1.88) == "1.88 years"
    assert format_period_years(0.5) == "183 days"
    assert format_period_years(None) == "–"
    assert format_duration_seconds(None) == "–"


def test_singular_only_when_exactly_one():
    assert format_period_years(1.001) == "1 year"
    assert format_period_years(1.01) == "1.01 years"
    assert format_duration_seconds(2 * 86400) == "2 days"


def test_galactic_orbit_text_uses_both_ladders():
    assert format_galactic_orbit(205.7, 0.23626) == "206 km/s (236 My per orbit)"


def test_ladders_are_sorted():
    assert [step for _, _, step in SPEED_LADDER] == sorted(step for _, _, step in SPEED_LADDER)
    assert [unit for *_, unit in PERIOD_LADDER] == sorted(unit for *_, unit in PERIOD_LADDER)


def test_year_is_julian():
    assert _YEAR == 365.25 * 86400


def _js_constants(name):
    with open(os.path.join(_STATIC, name), encoding="utf-8") as handle:
        return handle.read()


def test_browser_constants_match_python():
    speed_js = _js_constants("speed.js")
    assert re.search(r"SPEED_OF_LIGHT_KMS = ([0-9.]+);", speed_js).group(1) == repr(_C)
    for label, *_ in SPEED_LADDER:
        assert f'["{label}",' in speed_js
    period_js = _js_constants("period.js")
    assert "SECONDS_PER_YEAR = 365.25 * 24 * 3600;" in period_js
    for plural, singular, _ in PERIOD_LADDER:
        assert f'["{plural}", "{singular}",' in period_js


def _run_node(module, function, values):
    script = (
        "const m = await import(process.argv[1]);"
        "const values = JSON.parse(process.argv[2]);"
        f"console.log(JSON.stringify(values.map(v => m.{function}(v))));"
    )
    module_url = "file://" + os.path.abspath(os.path.join(_STATIC, module))
    out = subprocess.run(
        ["node", "--input-type=module", "-e", script, module_url, json.dumps(values)],
        check=True, capture_output=True, text=True,
    ).stdout
    return json.loads(out)


@pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")
def test_browser_speed_copy_matches_python():
    values = [kms for kms, _ in SPEED_CASES] + [None]
    assert _run_node("speed.js", "formatSpeedKms", values) == [format_speed_kms(v) for v in values]


@pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")
def test_browser_period_copy_matches_python():
    values = [s for s, _ in PERIOD_CASES] + [None]
    assert _run_node("period.js", "formatDurationSeconds", values) == [format_duration_seconds(v) for v in values]
    years = [0.5, 1.0, 1.88, 2.36e8]
    assert _run_node("period.js", "formatPeriodYears", years) == [format_period_years(v) for v in years]
