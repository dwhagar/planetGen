# tests/test_phenomenon_cuts.py

"""
GEN.377: the phenomenon scatter's mass cut by place. The neutron stars and black holes the scatter places are cut at the tier limit of
their band (core, baseline or halo: the limits GEN.333 gave the stars, from the same Generate page settings), the share below the cut is
left to the sector's fill at the same cut, and a galaxy whose tiers share one cut draws what one cut drew.
"""

import collections

import numpy as np
import pytest

from planetgen.galaxy import scatter_limits
from planetgen.galaxy.scatter_limits import BASELINE, CORE, HALO, ScatterLimits, TierLimit
from planetgen.generation import component_pass, phenomenon_reach as reach, phenomenon_scatter as scatter
from planetgen.generation.phenomenon_cuts import Cuts, as_cuts
from tests import test_gen_scatter_tiers as tiers
from tests.anchor_masses import anchor_mass
from tests.test_gen_component_pass import _cluster

SHAPE, EDGE = tiers.SHAPE, tiers.EDGE_PC
REACH_E_VALUE = 4000.0
PLACED_PER_STAR = {tiers.CORE_SHELL: 2.0, tiers.DISK_SHELL: 2.0, tiers.MIXED_SHELL: 400.0}
"""float: The galaxy's expected systems at density 1 the shell is drawn at: the toy galaxy's core and disk shells hold thousands of holes at 2, its mixed shell (mostly halo) a few hundred at 400."""
BASE, CORE_CUT, HALO_CUT = 14.0, 4.0, 8.0
COMPACT = ("black-hole", "neutron-star")


def _core_cuts():
    """The cuts the Generate page's settings give: a core mass limit of 4 under a baseline cut of 14; the halo takes the baseline's."""
    return Cuts.of(tiers._skeleton(core_mass=CORE_CUT), BASE)


def _three_cuts():
    """Cuts with a mass of their own in all three tiers (the halo's too), for the places where the baseline and the halo would otherwise agree."""
    limits = ScatterLimits(SHAPE, EDGE, (TierLimit(1.0, CORE_CUT), TierLimit(1.0, BASE), TierLimit(1.0, HALO_CUT)))
    return Cuts((CORE_CUT, BASE, HALO_CUT), limits)


def _compact(rows):
    return [row for row in rows if row[3] in COMPACT and row[4] != scatter.MAGNETAR_SUBTYPE]


def _shell(shell, cuts, seed=11):
    high = tiers._high(shell)
    return high, list(scatter.scatter_shell(SHAPE, shell, high, EDGE, PLACED_PER_STAR[shell], seed, min_mass_solar=cuts))


def test_the_cuts_are_the_generate_pages_three_limits():
    cuts = _core_cuts()
    assert cuts.tiers == (CORE_CUT, BASE, BASE) and cuts.baseline == BASE and cuts.lowest == CORE_CUT
    assert not cuts.is_uniform
    unset = Cuts.of(tiers._skeleton(), BASE)
    assert unset.tiers[CORE] == pytest.approx(BASE * 1.1) and unset.tiers[HALO] == BASE
    assert Cuts.uniform(BASE).is_uniform and as_cuts(BASE).tiers == (BASE, BASE, BASE)
    assert as_cuts(None).baseline > 0.0 and as_cuts(cuts) is cuts


def test_a_band_reads_the_cut_of_its_own_tier():
    cuts = _three_cuts()
    high = tiers._high(tiers.MIXED_SHELL)
    assert cuts.cut_at(tiers.CORE_SHELL, 0) == CORE_CUT
    assert cuts.cut_at(tiers.DISK_SHELL, 0) == BASE
    assert cuts.cut_at(tiers.MIXED_SHELL, -high) == HALO_CUT
    assert cuts.cut_at(tiers.MIXED_SHELL, 0) == BASE
    assert cuts.cut_at(None, None) == BASE, "a place with no band takes the baseline's"


def test_a_tier_keeps_the_share_its_cut_leaves_over_the_lowest_cuts():
    shares = scatter.tier_shares("black-hole", "stellar", _core_cuts())
    assert shares[CORE] == 1.0
    assert shares[BASELINE] == pytest.approx(scatter.share_above("black-hole", "stellar", BASE) / scatter.share_above("black-hole", "stellar", CORE_CUT))
    assert scatter.tier_shares("neutron-star", None, _core_cuts()) == [0.0, 0.0, 0.0], "no tier keeps a neutron star over the ordinary ones' 2.2 solar masses"


def test_a_core_band_places_the_compact_objects_down_to_its_own_cut_and_no_further():
    high, rows = _shell(tiers.CORE_SHELL, _core_cuts())
    holes = [row for row in rows if (row[3], row[4]) == ("black-hole", "stellar")]
    assert len(holes) > 30
    masses = np.array([anchor_mass(row, SHAPE) for row in holes])      # the cut is in the anchor law's masses (GEN.554)
    assert masses.min() >= CORE_CUT - 1e-6 and masses.min() < BASE, "the core places holes the baseline's cut would leave to the sectors"
    flat = [row for row in scatter.scatter_shell(SHAPE, tiers.CORE_SHELL, high, EDGE, PLACED_PER_STAR[tiers.CORE_SHELL], 11, min_mass_solar=BASE)
            if (row[3], row[4]) == ("black-hole", "stellar")]
    assert min(anchor_mass(row, SHAPE) for row in flat) >= BASE - 1e-6
    # the share of the law above each cut: more holes where the cut is lower
    assert len(holes) / len(flat) == pytest.approx(scatter.share_above("black-hole", "stellar", CORE_CUT) / scatter.share_above("black-hole", "stellar", BASE), rel=0.45)


def test_a_disk_band_keeps_the_baseline_cut_whatever_the_core_takes():
    _high, rows = _shell(tiers.DISK_SHELL, _core_cuts())
    holes = _compact(rows)
    assert holes and min(anchor_mass(row, SHAPE) for row in holes) >= BASE - 1e-6


def test_a_halo_band_and_a_baseline_band_each_cut_at_their_own_tier_limit():
    cuts = _three_cuts()
    shell = tiers.MIXED_SHELL
    high, rows = _shell(shell, cuts)
    tier_of = scatter_limits.tier_indices(SHAPE, shell, high, EDGE)
    by_tier = collections.defaultdict(list)
    for row in _compact(rows):
        by_tier[int(tier_of[row[1] + high])].append(anchor_mass(row, SHAPE))
    assert set(by_tier) == {BASELINE, HALO}
    assert min(by_tier[BASELINE]) >= BASE - 1e-6
    assert min(by_tier[HALO]) >= HALO_CUT - 1e-6 and min(by_tier[HALO]) < BASE, "the halo's lower cut places holes under the baseline's"


def test_cuts_that_agree_draw_what_one_cut_drew():
    for shell in (tiers.CORE_SHELL, tiers.DISK_SHELL, tiers.MIXED_SHELL):
        high = tiers._high(shell)
        number = list(scatter.scatter_shell(SHAPE, shell, high, EDGE, PLACED_PER_STAR[shell], 5, min_mass_solar=BASE))
        uniform = list(scatter.scatter_shell(SHAPE, shell, high, EDGE, PLACED_PER_STAR[shell], 5, min_mass_solar=Cuts.uniform(BASE)))
        equal = Cuts((BASE, BASE, BASE), ScatterLimits(SHAPE, EDGE, tuple(TierLimit(1.0, BASE) for _ in range(3))))
        assert number == uniform == list(scatter.scatter_shell(SHAPE, shell, high, EDGE, PLACED_PER_STAR[shell], 5, min_mass_solar=equal))


def test_the_expected_rows_of_a_shell_follow_the_cut_of_its_bands():
    args = (SHAPE, tiers.CORE_SHELL, tiers._high(tiers.CORE_SHELL), EDGE, PLACED_PER_STAR[tiers.CORE_SHELL])
    low, flat, tiered = (scatter.shell_expected(*args, cuts) for cuts in (CORE_CUT, BASE, _core_cuts()))
    assert tiered == pytest.approx(low) and tiered > flat, "a core shell expects what the core's lower cut places"


def test_beyond_the_outline_each_object_is_cut_at_the_tier_of_the_band_it_lands_in():
    cuts = _three_cuts()
    rows = list(reach.scatter_shell(SHAPE, 800, 60, EDGE, REACH_E_VALUE, 3, min_mass_solar=cuts, check=True))
    holes = _compact(rows)
    assert len(holes) > 20
    kinds = set()
    for row in holes:
        tier = scatter_limits.tier_index(SHAPE, row[0], row[1], EDGE)
        kinds.add(tier)
        assert row[12] >= cuts.tiers[tier] - 1e-9
    assert HALO in kinds
    assert any(row[12] < BASE for row in holes if scatter_limits.tier_index(SHAPE, row[0], row[1], EDGE) == HALO)
    flat = _compact(list(reach.scatter_shell(SHAPE, 800, 60, EDGE, REACH_E_VALUE, 3, min_mass_solar=BASE)))
    assert len(holes) > len(flat), "the halo's lower cut places more holes than the baseline's"


def test_a_clusters_phenomena_are_cut_at_the_tier_of_the_band_each_lands_in():
    cluster = _cluster(center=(20.0, 10.0, 5.0), mass=3.0e8, half=15.0)
    cuts = _core_cuts()
    rows = _compact(list(component_pass.scatter_phenomena([cluster], EDGE, 9, min_mass_solar=cuts, shape=SHAPE)))
    flat = _compact(list(component_pass.scatter_phenomena([cluster], EDGE, 9, min_mass_solar=BASE, shape=SHAPE)))
    assert rows and flat
    for row in rows:
        assert anchor_mass(row, SHAPE) >= cuts.tiers[scatter_limits.tier_index(SHAPE, row[0], row[1], EDGE)] - 1e-6
    assert min(anchor_mass(row, SHAPE) for row in rows) < BASE <= min(anchor_mass(row, SHAPE) for row in flat) + 1e-6


def test_the_stage_line_names_the_tiers_that_cut_elsewhere():
    assert scatter.cut_summary(BASE) == scatter.cut_summary(Cuts.uniform(BASE))
    line = scatter.cut_summary(_core_cuts())
    assert line.startswith(scatter.cut_summary(BASE)) and "the core's is 4 solar masses" in line
