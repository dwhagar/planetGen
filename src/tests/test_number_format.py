# tests/test_number_format.py

"""
The site's one number formatter (UX.20): `stellarObjects.utils.
format_number` shows anything with 5 or more digits before the decimal
point in scientific notation (3 significant figures, Unicode
superscripts), and `html/static/numberformat.js` mirrors it, checked
against the same table when Node is installed.
"""

import json
import os
import shutil
import subprocess

import pytest

from stellarObjects.utils import format_distance_km, format_number, scientific_text, years_to_time_string

_STATIC = os.path.join(os.path.dirname(__file__), "..", "html", "static")

# (value, Python spec, JS decimals, expected)
CASES = [
    (0, ",.0f", 0, "0"),
    (42, ",.0f", 0, "42"),
    (9_999, ",.0f", 0, "9,999"),
    (9_999.4, ",.0f", 0, "9,999"),
    (9_999.6, ",.0f", 0, "1.00 × 10⁴"),  # rounds to 5 digits
    (10_000, ",.0f", 0, "1.00 × 10⁴"),
    (1_234_567, ",.0f", 0, "1.23 × 10⁶"),
    (-54_321, ",.0f", 0, "-5.43 × 10⁴"),
    (1234.567, ",.2f", 2, "1,234.57"),
    (12345.6, ",.1f", 1, "1.23 × 10⁴"),
    (0.5, ",.2f", 2, "0.50"),
]


@pytest.mark.parametrize("value,spec,_decimals,expected", CASES)
def test_format_number(value, spec, _decimals, expected):
    assert format_number(value, spec) == expected


def test_scientific_text_negative_exponent():
    assert scientific_text(0.000123) == "1.23 × 10⁻⁴"


def test_distance_ladder_goes_scientific_past_four_digits():
    assert format_distance_km(9_999) == "9,999 km"
    assert format_distance_km(384_400) == "3.84 × 10⁵ km"


def test_periods_go_scientific_past_four_digits_of_years():
    assert years_to_time_string(12_345.6) == "1.23 × 10⁴ years"
    assert years_to_time_string(9_000) == "9000 years"


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
