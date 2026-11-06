"""
PERF.232: the first fill of a galaxy works out a band share per luminosity threshold and population, and does it in less time
by three changes that must not move a number: the stretches of a track that clear a luminosity come from one pass over the
samples (`bright_piece_kernel`) instead of two dozen array operations, a star whose bright windows all lie outside its
population's age window measures exactly 0 without its windows being worked out (`star_population._outside_window`), and what
a share reads whatever the threshold (the living count, the cells with a star that can count) is worked out once
(`star_population._fraction_terms`). Each is held against the way it was worked out before, to the last bit.
"""

import math
import random

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy.age_law import AgeLaw
from planetgen.generation import star_population as sp
from planetgen.physics import alpha_tracks, bright_piece_kernel, post_main_sequence as pm

HORIZONTAL_BRANCH = (pm.HORIZONTAL_BRANCH,)


def _cases(seed, count):
    rng = random.Random(seed)
    for _ in range(count):
        mass = math.exp(rng.uniform(math.log(0.5), math.log(150.0)))
        fe_h = rng.choice([-2.0, -1.2, -0.4, 0.0, 0.07, 0.25])
        yield mass, fe_h, [10.0 ** rng.uniform(0.0, 6.5) for _ in range(5)]


def test_the_kernel_finds_the_stretches_the_array_version_finds_to_the_last_bit():
    found = 0
    for mass, fe_h, floors in _cases(5, 300):
        line = pm.timeline(mass, fe_h)
        for floor in floors:
            for skip in ((), HORIZONTAL_BRANCH):
                fast = line.bright_pieces(floor, skip=skip)
                assert repr(fast) == repr(line._bright_pieces_by_array(math.log10(floor), skip)), (mass, fe_h, floor, skip)
                found += bool(fast)
    assert found > 500       # most cases have a stretch, so the comparison is not of empty lists


def test_the_kernel_joins_stretches_that_touch_and_cuts_them_to_the_phase():
    frac = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    log_l = np.array([0.0, 2.0, 0.0, 2.0, 2.0, 0.0])
    counts, starts, ends = bright_piece_kernel.bright_ranges(frac, log_l, 1.0, np.array([0, 3], dtype=np.int64),
                                                              np.array([5, 5], dtype=np.int64), np.array([0.0, 3.5]),
                                                              np.array([5.0, 4.5]))
    # Phase 0 is bright from 0.5 to 1.5 and from 2.5 to 4.5; phase 1 (samples 3 to 5, cut to 3.5..4.5) is bright from 3.5 to 4.5.
    assert counts.tolist() == [2, 1]
    assert starts.tolist() == [0.5, 2.5, 3.5]
    assert ends.tolist() == [1.5, 4.5, 4.5]


def test_a_track_with_a_missing_sample_takes_the_array_path(monkeypatch):
    cached = pm.timeline(3.0, 0.0)
    line = pm.Timeline(cached.mass_sol, cached.fe_h, cached.main_sequence_gy, cached.windows)   # not the cached one: it keeps the holed light it is given (TEST.612)
    real = pm.isochrones.post_main_sequence_light

    def holed(mass, fe_h):
        frac, log_l = real(mass, fe_h)
        log_l = log_l.copy()
        log_l[2] = np.nan
        return frac, log_l

    monkeypatch.setattr(pm.isochrones, "post_main_sequence_light", holed)
    calls = []
    monkeypatch.setattr(bright_piece_kernel, "bright_ranges", lambda *args: calls.append(args))
    line.bright_pieces(100.0, skip=HORIZONTAL_BRANCH)
    assert not calls


def _laws():
    uniform = {name: AgeLaw.uniform(tuning.STELLAR_POPULATION_AGE_RANGES_GY[name]) for name in ("young", "old", "halo")}
    ramp = AgeLaw((0.0, 4.0, 10.0), (3.0, 1.0, 0.0))
    return list(uniform.values()) + [ramp]


def test_a_star_the_prune_skips_measures_exactly_what_the_full_sum_gives():
    skipped = measured = 0
    for mass, fe_h, floors in _cases(11, 250):
        for law in _laws():
            for floor in floors:
                full = sum(sp._overlap(start, end, law) * chance for start, end, chance in sp._bright_windows(mass, floor, fe_h))
                if sp._outside_window(mass, floor, law, fe_h):
                    assert full == 0.0, (mass, fe_h, floor, law)
                    skipped += 1
                else:
                    assert sp._bright_measure(mass, floor, law, fe_h) == full
                    measured += full > 0.0
    assert skipped > 1000 and measured > 100


def test_a_threshold_of_zero_asks_whether_the_windows_are_outside_whatever_the_threshold():
    for mass, fe_h, floors in _cases(3, 150):
        for law in _laws():
            if sp._outside_window(mass, 0.0, law, fe_h):
                assert all(sp._outside_window(mass, floor, law, fe_h) for floor in floors)


def test_the_fraction_is_the_sum_over_every_cell_and_star_as_it_was_worked_out_before():
    for population in ("old", "halo", "young"):
        law = sp._law(population, None)
        fe_bin = sp.reference_fe_h(population)
        for floor in (150.0, 4000.0):
            mass_low, mass_high = None, None
            bright = living = 0.0
            for a, b, alpha, norm in sp._mass_grid((), sp._grid_cells(fe_bin, True)):
                count = norm * sp.imf.power_law_integral(a, b, alpha)
                points = [sp._power_law_quantile(a, b, alpha, (j + 0.5) / sp._FRACTION_POINTS_PER_CELL)
                          for j in range(sp._FRACTION_POINTS_PER_CELL)]
                living += count * sum(sp._living_measure(m, law, fe_bin) for m in points) / len(points)
                bright += count * sum(sp._full_bright_measure(m, floor, law, fe_bin) for m in points) / len(points)
            expected = bright / living if living > 0 else 0.0
            assert sp._bright_fraction(floor, population, None, law, fe_bin) == expected


def test_the_birth_track_pieces_are_the_ones_the_prune_reads():
    # `_outside_window` reads the pieces of `alpha_tracks.track_segments`; a star with a second piece (the thick disk's) is checked too.
    assert any(len(alpha_tracks.track_segments(fe_h)) > 1 for fe_h in (-2.0, -1.2, -0.4, 0.0, 0.07, 0.25))


def test_the_missing_sample_test_leaves_the_cached_timeline_alone():   # TEST.612
    """`pm.timeline` is cached, so a test that gives its timeline a holed track must use a copy: the next test to ask for the same star got the hole.
    Each order of the two tests that met here is run on its own, in a fresh process."""
    import os
    import subprocess
    import sys
    here = os.path.dirname(os.path.abspath(__file__))
    holed = "tests/test_bright_measure_exactness.py::test_a_track_with_a_missing_sample_takes_the_array_path"
    exact = "tests/test_post_main_sequence.py::test_the_bright_pieces_are_exactly_the_ages_the_star_is_that_bright"
    for order in ((holed, exact), (exact, holed)):
        done = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", *order], cwd=os.path.dirname(here),
                              capture_output=True, text=True)
        assert done.returncode == 0, (order, done.stdout[-1500:])
