# tests/test_distance_format.py

"""
The distance ladder (UX.6): `stellarObjects.utils.
format_distance_m` picks the largest of km < AU < mpc < cpc < ly < pc <
kpc < Mpc < Gpc the value is at least 1 of, parsec values carry a ly/AU/km
parenthetical, and body radii are always km in scientific notation. The
browser copy (`html/static/distance.js`) is checked against the same table
when Node is installed.
"""

import json
import os
import shutil
import subprocess

import pytest

from planetgen.physics import constants as pc
from stellarObjects.config import SystemConfig
from stellarObjects.utils import (
    DISTANCE_PAREN_MIN_LY,
    distance_parenthetical, format_body_radius_km, format_distance_au, format_distance_km,
    format_distance_ly, format_distance_m, format_distance_pc,
)

_STATIC = os.path.join(os.path.dirname(__file__), "..", "html", "static")

# (meters, expected text): one value in each unit and each rung boundary.
CASES = [
    (0.0, "0 km"),
    (500.0, "0.5 km"),
    (12_000.0, "12 km"),
    (9_999e3, "9,999 km"),
    (384_400e3, "3.84 × 10⁵ km"),
    (pc.AU_M * 0.999, "1.49 × 10⁸ km"),
    (pc.AU_M, "1 AU"),
    (pc.AU_M * 5.2, "5.2 AU"),
    (pc.MILLIPARSEC_M * 0.999, "206 AU"),
    (pc.MILLIPARSEC_M, "1 mpc (206 AU)"),
    (pc.MILLIPARSEC_M * 2.4, "2.4 mpc (495 AU)"),
    (pc.CENTIPARSEC_M, "1 cpc (0.0326 ly)"),
    (pc.CENTIPARSEC_M * 15.3, "15.3 cpc (0.499 ly)"),
    (pc.LIGHTYEAR_M * 0.999, "30.6 cpc (0.999 ly)"),
    (pc.LIGHTYEAR_M, "1 ly"),
    (pc.LIGHTYEAR_M * 2.5, "2.5 ly"),
    (pc.PARSEC_M, "1 pc (3.26 ly)"),
    (pc.PARSEC_M * 4.2, "4.2 pc (13.7 ly)"),
    (pc.KILOPARSEC_M * 8, "8 kpc (2.61 × 10⁴ ly)"),
    (pc.MEGAPARSEC_M * 1.5, "1.5 Mpc (4.89 × 10⁶ ly)"),
    (pc.GIGAPARSEC_M * 2, "2 Gpc (6.52 × 10⁹ ly)"),
]


@pytest.mark.parametrize("meters,expected", CASES)
def test_ladder_picks_the_largest_unit_at_least_one(meters, expected):
    assert format_distance_m(meters) == expected


def test_parsec_parenthetical_switches_from_ly_to_au_at_a_hundredth_ly():
    just_over = pc.LIGHTYEAR_M * 0.0101
    just_under = pc.LIGHTYEAR_M * 0.0099
    assert format_distance_m(just_over) == "3.1 mpc (0.0101 ly)"
    assert format_distance_m(just_under) == "3.04 mpc (626 AU)"
    assert format_distance_m(pc.LIGHTYEAR_M * DISTANCE_PAREN_MIN_LY).endswith("(0.01 ly)")


def test_parenthetical_switches_from_au_to_km_at_a_hundredth_au():
    # No parsec value is this small, so the rule is checked on its own.
    assert distance_parenthetical(pc.AU_M * 0.01) == "0.01 AU"
    assert distance_parenthetical(pc.AU_M * 0.0099) == "1.48 × 10⁶ km"


def test_unit_wrappers_agree():
    assert format_distance_km(149_597_870.7) == "1 AU"
    assert format_distance_au(1.0) == "1 AU"
    assert format_distance_ly(1.0) == "1 ly"
    assert format_distance_pc(1.0) == "1 pc (3.26 ly)"
    assert format_distance_km(None) == "–"


def test_exact_constants():
    assert pc.AU_TO_KM == 149_597_870.7
    assert pc.LY_TO_M == 9_460_730_472_580_800
    assert pc.AU_PER_PARSEC == pytest.approx(648000 / 3.141592653589793, rel=1e-12)
    assert pc.LY_TO_AU == pytest.approx(63241.077, rel=1e-8)


def test_body_radius_is_always_scientific_km():
    config = SystemConfig()
    config.MARKDOWN = True
    assert format_body_radius_km(config, 6371) == "6.37 × 10<sup>3</sup> km"
    assert format_body_radius_km(config, 1737.4) == "1.74 × 10<sup>3</sup> km"
    assert format_body_radius_km(config, 696_000) == "6.96 × 10<sup>5</sup> km"


@pytest.mark.skipif(shutil.which("node") is None, reason="Node isn't installed")
def test_browser_copy_matches_python():
    meters = [m for m, _ in CASES] + [pc.LIGHTYEAR_M * 0.0101, pc.LIGHTYEAR_M * 0.0099]
    script = (
        "const m = await import(process.argv[1]);"
        "const values = JSON.parse(process.argv[2]);"
        "console.log(JSON.stringify(values.map(v => m.formatDistanceM(v))));"
    )
    module_url = "file://" + os.path.abspath(os.path.join(_STATIC, "distance.js"))
    out = subprocess.run(
        ["node", "--input-type=module", "-e", script, module_url, json.dumps(meters)],
        check=True, capture_output=True, text=True,
    ).stdout
    assert json.loads(out) == [format_distance_m(v) for v in meters]
