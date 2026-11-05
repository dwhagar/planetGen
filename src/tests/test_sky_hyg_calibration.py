# tests/test_sky_hyg_calibration.py

"""
VIEW.30 (S12 of the system-sky plan, goal #4907): the generated sky held against the real one.

A developer-machine calibration. It is skipped when the HYG catalogue is not downloaded (`python3 -m planetgen.cli.datasets fetch hyg_stars`, which setup runs) and unless
`PLANETGEN_SKY_CALIBRATION=1`, because it plans a galaxy, scatters its bright stars and charts a sky, about 25 to 40 minutes on a four-core machine. Run it when the sky tiers, the
backfill, the star populations or the sky photometry change:

    PLANETGEN_SKY_CALIBRATION=1 python3 -m pytest src/tests/test_sky_hyg_calibration.py -s

What it does: the HYG catalogue is the oracle (`datasets/hyg_stars.py`: 8,920 stars at V <= 6.5 and 25,791 at V <= 7.5, nothing in the generator reads it). The test plans the galaxy the seed
`SEED_HEX` makes with every lever on Auto, scatters its bright stars at the default thresholds, finds the place in the disc plane (between two arms, past the reference radius) where the model's relative
density is 1, the real local density, charts the sky tiers of that cell for 6.5 and counts the stars `sky_from_database` finds there.

The band (measured 2026-11-05 on this seed, reports/system-sky-hyg-calibration-2026-11-05.md): the dust-free sky holds 2,184 stars, 24 percent of the catalogue's 8,920; with the default dust law 586, 7 percent.
By absolute magnitude the model is right for M_V -6 to -3 (321 against 318) and for the faint dwarfs, but holds a fifth of the real stars from -3 to 4 (A and F stars, the giants) and 24 times
the real number brighter than -6 (237 against 10: with no dust the distant O and B stars all show). The test holds the dust-free count between 15 and 40 percent of the catalogue so that a change
that moves it either way is seen and the band is then set again on purpose; closing the gap is the follow-up, not a reason to widen the band.
"""

import argparse
import json
import math
import os

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
"""str: The galaxy seed of the calibration: the galaxy it makes with every lever on Auto has a 480 parsec disc, a density-1 place at 2,885 parsecs and 7,300 pre-placed stars."""

LIMIT = 6.5
"""float: The limiting magnitude of the naked-eye sky."""

BAND = (0.15, 0.40)
"""tuple: The share of the catalogue's count the dust-free generated sky must hold (measured 0.245; see the module docstring)."""

pytestmark = [
    pytest.mark.skipif(os.environ.get("PLANETGEN_SKY_CALIBRATION") != "1", reason="the sky calibration takes 25 to 40 minutes: PLANETGEN_SKY_CALIBRATION=1 runs it"),
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
    """The point in the disc plane, between two arms, from the reference radius outward where the relative density falls to 1."""
    def at(radius):
        angle = density._interarm_angle(radius, shape)
        return (radius * math.cos(angle), radius * math.sin(angle), 0.0)

    low, high = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * shape.disk_scale_length_pc, 12.0 * shape.disk_scale_length_pc
    for _ in range(40):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if density.relative_density(at(middle), shape) > 1.0 else (low, middle)
    return at(0.5 * (low + high))


def test_the_generated_naked_eye_sky_holds_the_share_of_the_real_one_that_was_measured(control_database):  # noqa: F811
    mysql_config = control_database
    real = hyg_stars.count_at_or_brighter(LIMIT)
    assert 8900 <= real <= 8940

    args = _parsed(generate_cli.add_plan_arguments, generate_cli.validate_plan_args, mysql_config, "--seed", SEED_HEX)
    commands.run_alone("plan-skeleton", json.loads(json.dumps(run_plan.skeleton_params(args, mysql_config))))
    for kind in ("--mass", "--luminosity"):
        run_plan.scatter_bright_stars(_parsed(generate_cli.add_scatter_arguments, generate_cli.validate_scatter_args, mysql_config, kind))

    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        observer_pc = _density_one_place(skeleton.shape)
        assert density.relative_density(observer_pc, skeleton.shape) == pytest.approx(1.0, abs=1e-6)
        sky_charting.chart(conn, observer_pc, LIMIT)
        observer = sky.Observer(observer_pc, 1.0e9)
        clear = sky.sky_from_database(conn, observer, 1.0e9, LIMIT, galaxy_key=SEED_HEX, dust={"av_per_kpc_at_sun": 0.0})
        dusty = sky.sky_from_database(conn, observer, 1.0e9, LIMIT, galaxy_key=SEED_HEX)
    finally:
        conn.close()

    share = len(clear.uid) / real
    print(f"\nsky calibration: {len(clear.uid)} stars dust-free ({share:.3f} of the catalogue's {real}), {len(dusty.uid)} with the default dust")
    assert BAND[0] <= share <= BAND[1], f"the dust-free sky holds {len(clear.uid)} stars at {LIMIT}, {share:.3f} of the catalogue's {real}; the band is {BAND}"
    assert 0 < len(dusty.uid) < len(clear.uid)
    assert float(clear.magnitude.max()) <= LIMIT
