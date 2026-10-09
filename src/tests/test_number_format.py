# tests/test_number_format.py

"""
The site's one number formatter (UX.20): `planetgen.util.format.
format_number` shows a whole number with 7 or more digits, or a number
with decimals and 5 or more digits before the decimal point, in
scientific notation (UX.36; 3 significant figures, Unicode superscripts), and `html/static/numberformat.js` mirrors it, checked
against the same table when Node is installed.
"""

import json
import os
import shutil
import subprocess

import pytest

from planetgen.util.format import format_distance_km, format_number, format_period_years, scientific_text

_STATIC = os.path.join(os.path.dirname(__file__), "..", "html", "static")

# (value, Python spec, JS decimals, expected)
CASES = [
    (0, ",.0f", 0, "0"),
    (42, ",.0f", 0, "42"),
    (9_999, ",.0f", 0, "9,999"),
    (10_000, ",.0f", 0, "10,000"),
    (999_999, ",.0f", 0, "999,999"),  # UX.36: 6 whole digits stay plain
    (999_999.4, ",.0f", 0, "999,999"),
    (999_999.6, ",.0f", 0, "1.00 × 10⁶"),  # rounds to 7 digits
    (1_234_567, ",.0f", 0, "1.23 × 10⁶"),
    (-54_321, ",.0f", 0, "-54,321"),
    (-7_654_321, ",.0f", 0, "-7.65 × 10⁶"),
    (1234.567, ",.2f", 2, "1,234.57"),  # 4 whole digits with decimals: plain
    (12345.6, ",.1f", 1, "1.23 × 10⁴"),  # 5 with decimals: scientific
    (0.5, ",.2f", 2, "0.50"),
    # UX.79: ties round away from zero on the shortest decimal, in both copies.
    (9.995, ",.2f", 2, "10.00"),
    (1.005, ",.2f", 2, "1.01"),
    (2.675, ",.2f", 2, "2.68"),
    (0.285, ",.2f", 2, "0.29"),
    (-1.005, ",.2f", 2, "-1.01"),
    (0.5, ",.0f", 0, "1"),
    (1.5, ",.0f", 0, "2"),
    (2.5, ",.0f", 0, "3"),
    (-2.5, ",.0f", 0, "-3"),
    # UX.80: a negative that rounds to zero is "0", not "-0".
    (-0.4, ",.0f", 0, "0"),
    (-0.0, ",.0f", 0, "0"),
    (-0.004, ",.2f", 2, "0.00"),
    (-0.4, ",.1f", 1, "-0.4"),
]


@pytest.mark.parametrize("value,spec,_decimals,expected", CASES)
def test_format_number(value, spec, _decimals, expected):
    assert format_number(value, spec) == expected


def test_three_figures_rounds_half_up_and_has_no_negative_zero():
    from planetgen.util.format import _three_figures
    assert _three_figures(9.995) == "10"
    assert _three_figures(1.005) == "1.01"
    assert _three_figures(-0.0) == "0"
    assert _three_figures(0.0) == "0"
    assert _three_figures(2.675) == "2.68"


@pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")
def test_browser_three_figures_matches_python():
    from planetgen.util.format import _three_figures
    values = [9.995, 1.005, 2.675, -0.0, 0.0, 12345.678, 0.000499, 4.2, 495.5, -1.005]
    script = (
        "const m = await import(process.argv[1]);"
        "console.log(JSON.stringify(JSON.parse(process.argv[2]).map((v) => m.threeFigures(v))));"
    )
    module_url = "file://" + os.path.abspath(os.path.join(_STATIC, "numberformat.js"))
    out = subprocess.run(["node", "--input-type=module", "-e", script, module_url, json.dumps(values)],
                         check=True, capture_output=True, text=True).stdout
    assert json.loads(out) == [_three_figures(v) for v in values]


def test_scientific_text_negative_exponent():
    assert scientific_text(0.000123) == "1.23 × 10⁻⁴"


def test_distance_ladder_goes_scientific_past_six_digits():
    assert format_distance_km(9_999) == "9,999 km"
    assert format_distance_km(384_400) == "384,400 km"


def test_periods_go_scientific_past_six_digits_of_years():
    # Past 999 years the period ladder moves to ky, so a period only goes
    # scientific past 999,999 Gy.
    assert format_period_years(9_000) == "9 ky"
    assert format_period_years(1.2346e13) == "12,346 Gy"
    assert format_period_years(1.2346e15) == "1.23 × 10⁶ Gy"


@pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")
def test_browser_copy_matches_python():
    pairs = [[value, decimals] for value, _spec, decimals, _ in CASES]
    script = (
        "const m = await import(process.argv[1]);"
        "const pairs = JSON.parse(process.argv[2]);"
        "console.log(JSON.stringify(pairs.map(([v, d]) => m.formatNumber(v, d))));"
    )
    module_url = "file://" + os.path.abspath(os.path.join(_STATIC, "numberformat.js"))
    out = subprocess.run(
        ["node", "--input-type=module", "-e", script, module_url, json.dumps(pairs)],
        check=True, capture_output=True, text=True,
    ).stdout
    assert json.loads(out) == [expected for *_rest, expected in CASES]
