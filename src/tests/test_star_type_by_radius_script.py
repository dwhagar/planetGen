# tests/test_star_type_by_radius_script.py

"""
GEN.542: `scripts/star_type_by_radius.py` runs to its tables with the shape the plan draws from the seed (it stopped with a TypeError when the plan's shape options defaulted to None).
"""

import pathlib
import subprocess
import sys

SCRIPT = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "star_type_by_radius.py"


def test_the_script_runs_with_a_few_stars_and_prints_every_table():
    run = subprocess.run([sys.executable, str(SCRIPT), "--stars", "20", "--seed", "3"], capture_output=True, text=True, timeout=300)
    assert run.returncode == 0, run.stderr[-800:]
    assert "## Class shares inside each population" in run.stdout
    assert run.stdout.count("| R (kpc) |") == 4                                     # the arm crest and midway, in the plane and 500 pc above it
    assert "spiral-arm crest, z = 500 pc" in run.stdout and "midway between arms, z = 0 pc" in run.stdout
