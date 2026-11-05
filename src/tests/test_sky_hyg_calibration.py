# tests/test_sky_hyg_calibration.py

"""
VIEW.30 (S12 of the system-sky plan, goal #4907): the generated sky held against the real one.

A developer-machine calibration. It is skipped when the HYG catalogue is not downloaded (`python3 -m planetgen.cli.datasets fetch hyg_stars`, which setup runs) and unless asked for with
`PLANETGEN_SKY_CALIBRATION`, because it plans a galaxy, scatters its bright stars and charts a sky. Run it when the sky tiers, the backfill, the star populations or the sky photometry change:

    PLANETGEN_SKY_CALIBRATION=small python3 -m pytest src/tests/test_sky_hyg_calibration.py -s        # about 25 minutes
    PLANETGEN_SKY_CALIBRATION=milky-way python3 -m pytest src/tests/test_sky_hyg_calibration.py -s    # about 2 hours 15 minutes; "all" runs both

The HYG catalogue is the oracle (`datasets/hyg_stars.py`: 8,920 stars at V <= 6.5, nothing in the generator reads it). Each case plans a galaxy, scatters its bright stars, charts the sky tiers of a cell
for 6.5 and counts the stars `sky_from_database` finds there, with no dust and with the default dust law. The reports are reports/system-sky-hyg-calibration-2026-11-05.md.

* `milky-way`: the galaxy of a Milky Way-sized disk (1e11 solar masses, a 2,600 parsec scale length, a thin disk 404 parsecs thick that is 69 percent of the stars at the Sun's place) at its reference
  radius of 3.15 scale lengths, the Sun's place, where the relative density is 0.91. The scatter places every star of 20 solar masses and 4,000,000 solar luminosities and up. Measured 2026-11-05: 5,733
  stars with no dust, 64 percent of the catalogue's 8,920 (70 percent of it at the place's density), 3,207 (36 percent) with the default dust. By absolute magnitude, no dust, against the catalogue: brighter than -6
  1,548 against 10 (dust takes the real ones), -6 to -3 456 against 318, -3 to 0 2,141 against 3,363, 0 to 2 1,083 against 3,616, 2 to 4 373 against 1,164, 4 to 6 119 against 219.
* `small`: the galaxy the default seed draws (a 480 parsec disk, 92 percent thick disk at the place where the density is 1), which is the cheap regression guard: 2,184 stars, 24.5 percent.
  Star 2's analysis (reports/gen1279-1280-sky-counts-cause-2026-11-05.md) puts most of that shortfall on the galaxy, not the generator, so its band is only there to catch a change.

The bands hold each case's dust-free count to a window round what was measured so that a change that moves it either way is seen and the band is then set again on purpose.
What is still short in the Milky Way-sized case is the middle of the range: a third of the real stars from absolute magnitude 0 to 4 (1,456 against 4,780).
"""

import argparse
import json
import math
import os

import numpy as np
import pytest

from planetgen import tuning
from planetgen.cli import generate as generate_cli
from planetgen.datasets import fetch, hyg_stars
from planetgen.db import store
from planetgen.galaxy import density, sky
from planetgen.generation import run_plan, sky_charting
from planetgen.queue import commands
from tests.test_galaxy_gen import _mysql_argv
from tests.test_work_units import control_database  # noqa: F401 -- the fixture

SEED_HEX = "5e7a1c3b9d2f4608a1b2c3d4e5f60718"
"""str: The galaxy seed of the calibration."""

LIMIT = 6.5
"""float: The limiting magnitude of the naked-eye sky."""

CASES = {
    "small": {"pins": (), "scatter": {}, "site": "density-one", "band": (0.15, 0.40)},
    "milky-way": {"pins": ("--mass", "1.0e11", "--disk-scale-length", "2600pc"), "scatter": {"--phenomenon-min-mass": "20", "--bright-star-min-luminosity": "4000000"},
                  "site": "reference", "band": (0.50, 0.80)},
}
"""dict: The calibrations: the plan's pins, the scatter's floors, where the observer stands and the share of the catalogue's count the dust-free sky must hold (measured 0.245 and 0.643)."""

WANTED = os.environ.get("PLANETGEN_SKY_CALIBRATION", "")

pytestmark = [
    pytest.mark.skipif(fetch.locate(hyg_stars.DATASET).missing, reason="the HYG catalogue is not downloaded (python3 -m planetgen.cli.datasets fetch hyg_stars)"),
    pytest.mark.fresh_database,
]


def _parsed(add_arguments, validate, mysql_config, *extra):
    parser = argparse.ArgumentParser(prefix_chars="-+")
    add_arguments(parser)
    args = parser.parse_args(list(extra) + _mysql_argv(mysql_config))
    validate(args, parser)
    return args


def _density_one_place(shape):
    """The point in the disk plane, between two arms, from the reference radius outward where the relative density falls to 1."""
    def at(radius):
        angle = density._interarm_angle(radius, shape)
        return (radius * math.cos(angle), radius * math.sin(angle), 0.0)

    low, high = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * shape.disk_scale_length_pc, 12.0 * shape.disk_scale_length_pc
    for _ in range(40):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if density.relative_density(at(middle), shape) > 1.0 else (low, middle)
    return at(0.5 * (low + high))


def _reference_place(shape):
    """The Sun's place in the model: the reference radius (3.15 scale lengths) in the disk plane, between two arms."""
    radius = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * shape.disk_scale_length_pc
    angle = density._interarm_angle(radius, shape)
    return (radius * math.cos(angle), radius * math.sin(angle), 0.0)


def _reference_place(shape):
    """The Sun's place in the model: the reference radius (3.15 scale lengths) in the disk plane, between two arms."""
    radius = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * shape.disk_scale_length_pc
    angle = density._interarm_angle(radius, shape)
    return (radius * math.cos(angle), radius * math.sin(angle), 0.0)


@pytest.mark.parametrize("name", sorted(CASES))
def test_the_generated_naked_eye_sky_holds_the_share_of_the_real_one_that_was_measured(name, control_database):  # noqa: F811
    if WANTED not in (name, "all") and not (WANTED == "1" and name == "small"):
        pytest.skip(f"the {name} sky calibration is run on request: PLANETGEN_SKY_CALIBRATION={name}")
    case = CASES[name]
    mysql_config = control_database
    real = hyg_stars.count_at_or_brighter(LIMIT)
    assert 8900 <= real <= 8940

    args = _parsed(generate_cli.add_plan_arguments, generate_cli.validate_plan_args, mysql_config, "--seed", SEED_HEX, *case["pins"])
    commands.run_alone("plan-skeleton", json.loads(json.dumps(run_plan.skeleton_params(args, mysql_config))))
    floors = case["scatter"]
    for kind, flag in (("--mass", "--phenomenon-min-mass"), ("--luminosity", "--bright-star-min-luminosity")):
        extra = [flag, floors[flag]] if flag in floors else []
        run_plan.scatter_bright_stars(_parsed(generate_cli.add_scatter_arguments, generate_cli.validate_scatter_args, mysql_config, kind, *extra))

    conn = store.get_connection(mysql_config)
    try:
        shape = store.get_galaxy_shape(conn).shape
        observer_pc = _density_one_place(shape) if case["site"] == "density-one" else _reference_place(shape)
        relative = density.relative_density(observer_pc, shape)
        assert 0.9 <= relative <= 1.1, relative
        sky_charting.chart(conn, observer_pc, LIMIT)
        observer = sky.Observer(observer_pc, 1.0e9)
        clear = sky.sky_from_database(conn, observer, 1.0e9, LIMIT, galaxy_key=SEED_HEX, dust={"av_per_kpc_at_sun": 0.0})
        dusty = sky.sky_from_database(conn, observer, 1.0e9, LIMIT, galaxy_key=SEED_HEX)
    finally:
        conn.close()

    share = len(clear.uid) / real
    absolute = clear.magnitude - 5.0 * np.log10(clear.distance_pc / 10.0)
    bins = [int(((absolute >= low) & (absolute < high)).sum()) for low, high in ((-12, -6), (-6, -3), (-3, 0), (0, 2), (2, 4), (4, 6), (6, 20))]
    print(f"\nsky calibration {name}: {len(clear.uid)} stars dust-free ({share:.3f} of the catalogue's {real}), {len(dusty.uid)} with the default dust;"
          f" relative density {relative:.2f}; by absolute magnitude (brighter than -6, -6 to -3, -3 to 0, 0 to 2, 2 to 4, 4 to 6, fainter) {bins}")
    low, high = case["band"]
    assert low <= share <= high, f"the dust-free sky holds {len(clear.uid)} stars at {LIMIT}, {share:.3f} of the catalogue's {real}; the band is {case['band']}"
    assert 0 < len(dusty.uid) < len(clear.uid)
    assert float(clear.magnitude.max()) <= LIMIT
