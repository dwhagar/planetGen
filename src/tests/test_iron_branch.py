# tests/test_iron_branch.py

"""GEN.512: the iron branch of a rocky planet (`bulk_draw.IRON_BRANCH_TABLE`) and its tuning to the homogeneous samples of rocky planets."""

import dataclasses
import types

import pytest

from planetgen.generation import possibility_table as pt, tuning_harness as th
from planetgen.physics import bulk_draw as bd, iron_branch_tuning as it
from planetgen.util import draw

SAMPLE_TARGETS = ["t5-outlier-share", "t5-within-share"]


def _ratios():
    return bd.HostRatios(fe=0.9, si=0.9, al=0.08, ca=0.07, ni=0.05, c_over_o=0.5, y_envelope=0.27, z=0.014)


def _body(planet_class="M", is_moon=False):
    return types.SimpleNamespace(planet_class=planet_class, is_moon=is_moon)


@pytest.mark.parametrize("target_id", SAMPLE_TARGETS)
def test_the_shipped_branch_reproduces_the_published_sample(target_id):
    target = th.TARGETS.get(target_id)
    result = th.simulate_and_compare(target, it.sampler(bd.IRON_ENRICHED_WEIGHT, bd.IRON_SCATTER_DEX), n=3000)
    assert result.passed, (target_id, result.simulated, target.value, result.z)


def test_a_branch_that_is_far_too_common_is_refused_by_the_outlier_share():
    target = th.TARGETS.get("t5-outlier-share")
    result = th.simulate_and_compare(target, it.sampler(0.6, bd.IRON_SCATTER_DEX), n=3000)
    assert not result.passed and result.simulated > 0.25


def test_the_within_share_barely_moves_with_the_scatter_so_its_power_is_stated():
    """The scatter is not measured: 0 and 0.1 dex both pass the share within 1 sigma (the measurement error dominates)."""
    target = th.TARGETS.get("t5-within-share")
    for scatter in (0.0, 0.1):
        assert th.simulate_and_compare(target, it.sampler(bd.IRON_ENRICHED_WEIGHT, scatter), n=3000).passed


def test_the_table_has_two_tunable_rows_each_with_a_target_and_the_enriched_weight_is_the_constant():
    rows = {row.name: row for row in bd.IRON_BRANCH_TABLE.rows}
    assert set(rows) == {"star-matched", "iron-enriched"}
    assert all(row.tunable and row.state == pt.PARTLY_VERIFIED for row in rows.values())
    probabilities = dict((row.name, p) for row, p in bd.IRON_BRANCH_TABLE.probabilities(pt.Subject("M", {})))
    assert probabilities["iron-enriched"] == pytest.approx(bd.IRON_ENRICHED_WEIGHT)
    assert {t.row for t in th.TARGETS if t.fill_id == "GEN.512"} == {"bulk-iron-branch:star-matched", "bulk-iron-branch:iron-enriched"}


def test_a_moon_an_icy_small_or_carbon_body_has_no_iron_branch_and_draws_nothing():
    ratios = _ratios()
    draw.set_run_seed(1)
    for body in (_body("M", is_moon=True), _body("D"), _body("S"), _body("C")):
        assert bd.iron_branch_ratios(body, ratios, 0.08, 0.07, False, 0.0) is ratios
    assert "D" not in bd.IRON_BRANCH_CLASSES and "M" in bd.IRON_BRANCH_CLASSES


def test_star_matched_iron_is_the_stars_times_the_scatter_and_enriched_iron_lands_in_the_range():
    ratios = _ratios()
    draw.set_run_seed(5)
    same = bd.branch_ratios(ratios, "star-matched", 0.0, 0.08, 0.07, False, 0.0)
    assert same.fe == ratios.fe
    scattered = [bd.branch_ratios(ratios, "star-matched", 0.1, 0.08, 0.07, False, 0.0).fe / ratios.fe for _ in range(400)]
    logs = sorted(__import__("math").log10(x) for x in scattered)
    assert abs(logs[len(logs) // 2]) < 0.03 and 0.07 < (logs[340] - logs[60]) / 2 < 0.14, "a normal of 0.1 dex"
    low, high = bd.RANGES["enriched_cmf"]
    star_cmf = bd.dry_partition(ratios, 0.08, 0.07, False).metal
    for _ in range(100):
        enriched = bd.branch_ratios(ratios, "iron-enriched", 0.1, 0.08, 0.07, False, 0.0)
        cmf = bd.dry_partition(enriched, 0.08, 0.07, False).metal
        assert max(low, star_cmf) - 1e-6 <= cmf <= high + 1e-6
        assert enriched.fe > ratios.fe and dataclasses.replace(enriched, fe=ratios.fe) == ratios


def test_a_generated_rocky_planet_uses_the_branch(monkeypatch):
    calls = []
    original = bd.iron_branch_ratios

    def spy(*args, **kwargs):
        calls.append(args[0].planet_class)
        return original(*args, **kwargs)
    monkeypatch.setattr(bd, "iron_branch_ratios", spy)
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from tests.fuzz_support import deterministic_entropy
    for seed in range(6):
        with deterministic_entropy(seed):
            StarSystem(system_config=SystemConfig())
    assert calls, "no rocky planet reached the iron branch"
