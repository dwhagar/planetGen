# tests/test_estimate_sampling.py

"""
PERF.36: a run of more than `generationStats.SAMPLE_ABOVE` sectors is estimated from a sample, its sectors are
resolved as they are submitted and not held, and a web request is refused a ball of more than
`run_common.WEB_ENUMERATION_LIMIT` sectors.
"""

import math
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import run_common, run_galaxy
from planetgen.generation import stats as generationStats
from planetgen.generation.stats import GenerationStats, estimate
from planetgen.physics.units import ly_to_pc
from tests.test_api import admin_client, client, default_admin_client, first_admin_password  # noqa: F401
from tests.test_galaxy_gen import _all_sectors, _mysql_argv, _plan_wide_galaxy, _run_cli

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)


# --- The sample ----------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("count", [0, 1, 199, 200])
def test_a_run_no_larger_than_the_sample_is_sampled_whole(count):
    assert generationStats.sample_positions(count) == list(range(count))


@pytest.mark.parametrize("count", [201, 1000, 20_001, 7_455_000])
def test_the_sample_spreads_evenly_through_the_run(count):
    places = generationStats.sample_positions(count)
    assert len(places) == generationStats.SAMPLE_SIZE == len(set(places))
    assert places == sorted(places) and 0 <= places[0] and places[-1] < count
    stride = count / generationStats.SAMPLE_SIZE
    assert all(math.isclose(b - a, stride, abs_tol=1) for a, b in zip(places, places[1:]))
    assert places[0] < stride and count - places[-1] <= stride    # neither end of the run is left out


def test_an_estimate_from_a_sample_scales_up_and_says_so():
    stats = GenerationStats()
    count = 100_000
    expected = [1000.0 * math.exp(-n / 30_000) for n in range(count)]     # a smooth fall-off through the run
    full = estimate([(1.0, e) for e in expected], stats, "galaxy", workers=4)
    sample = [(1.0, expected[n]) for n in generationStats.sample_positions(count)]
    sampled = estimate(sample, stats, "galaxy", workers=4, total=count)
    assert sampled.sampled == generationStats.SAMPLE_SIZE and full.sampled == 0
    assert sampled.sectors == count == full.sectors
    for name in ("systems", "stars", "bytes", "seconds"):
        assert getattr(sampled, name) == pytest.approx(getattr(full, name), rel=0.01), name
    assert "worked out from 200 sampled sectors" in sampled.summary()
    assert sampled.as_dict()["sampled"] == generationStats.SAMPLE_SIZE


def test_an_estimate_without_a_total_is_unchanged():
    result = estimate([(1.0, 20.0), (1.0, 10.0)], GenerationStats(), "galaxy", workers=1)
    assert result.sectors == 2 and result.systems == 30.0 and result.sampled == 0
    assert "sampled" not in result.summary()


# --- The walk ------------------------------------------------------------------------------------------------------


class _Bounds:
    """Reaches shells below 50 only."""

    @staticmethod
    def contains(shell, band):
        return shell < 50


class _Density:
    bounds = _Bounds()

    @staticmethod
    def resolve(args, address, position_pc):
        return None if address[0] >= 50 else SimpleNamespace(address=address, position=position_pc)


def _pending(addresses, occupied=(), limit=None):
    return run_galaxy._Pending.of_addresses(SimpleNamespace(), lambda: iter(addresses), set(occupied), _Density(),
                                            EDGE_PC, limit=limit)


ADDRESSES = [(shell, 0, wedge) for shell in (1, 2, 60) for wedge in range(6)]


def test_a_walk_counts_what_exists_what_is_out_of_reach_and_what_is_left():
    pending = _pending(ADDRESSES, occupied=[(1, 0, 2), (2, 0, 5)])
    assert pending.survey() == 10
    assert (pending.count, pending.existed, pending.skipped) == (10, 2, 6)
    items = list(pending.items())
    assert [item[0] for item in items] == [a for a in ADDRESSES if a[0] < 50 and a not in {(1, 0, 2), (2, 0, 5)}]
    assert all(item[2].address == item[0] and item[3] == "" for item in items)


def test_a_limit_stops_the_walk_at_that_many_sectors():
    pending = _pending(ADDRESSES, occupied=[(1, 0, 0)], limit=4)
    assert pending.survey() == 4 and pending.existed == 1
    assert [item[0] for item in pending.items()] == [(1, 0, 1), (1, 0, 2), (1, 0, 3), (1, 0, 4)]


def test_the_sample_takes_sectors_spread_through_the_run(monkeypatch):
    monkeypatch.setattr(generationStats, "SAMPLE_SIZE", 3)
    pending = _pending([(shell, 0, 0) for shell in range(1, 31)])
    pending.survey()
    assert [arg.address[0] for arg in pending.sample_args()] == [1 + n for n in (5, 15, 25)]
    assert _pending([]).sample_args() == []


def test_a_small_run_is_a_list_estimated_whole_and_a_large_one_streams_and_is_sampled(monkeypatch):
    small = _pending([(shell, 0, 0) for shell in range(1, 11)])
    batch, estimate_args, total = run_galaxy._plan_pending(small)
    assert isinstance(batch, list) and len(batch) == len(estimate_args) == 10 and total is None
    monkeypatch.setattr(generationStats, "SAMPLE_ABOVE", 9)
    monkeypatch.setattr(generationStats, "SAMPLE_SIZE", 4)
    large = _pending([(shell, 0, 0) for shell in range(1, 11)])
    batch, estimate_args, total = run_galaxy._plan_pending(large)
    assert not isinstance(batch, list) and len(estimate_args) == 4 and total == 10
    assert [item[0] for item in batch] == [(shell, 0, 0) for shell in range(1, 11)]


def test_the_candidates_of_a_neighbourhood_walk_with_their_distance_in_the_suffix():
    candidates = [(1, 0, 0, 1.0, 2.0, 3.0, 4.5), (1, 0, 1, 1.0, 2.0, 3.0, 5.5), (60, 0, 0, 0.0, 0.0, 0.0, 6.5)]
    pending = run_galaxy._Pending.of_candidates(SimpleNamespace(), candidates, {(1, 0, 0)}, _Density(),
                                                ", {distance:.1f} pc away", skip={(9, 9, 9)})
    assert pending.survey() == 1 and (pending.existed, pending.skipped) == (1, 1)
    ((address, position, args, suffix),) = pending.items()
    assert (address, position, suffix) == ((1, 0, 1), (1.0, 2.0, 3.0), ", 5.5 pc away") and args.address == address


# --- In the run ----------------------------------------------------------------------------------------------------


FIRST_THREE = ["--shell", "1", "--band", "0", "--limit", "3", "--num-systems", "2"]


def test_a_streamed_shell_run_estimates_from_a_sample_and_generates_every_sector(mysql_config, monkeypatch, capsys):
    _plan_wide_galaxy(mysql_config)
    monkeypatch.setattr(generationStats, "SAMPLE_ABOVE", 1)
    monkeypatch.setattr(generationStats, "SAMPLE_SIZE", 2)
    _run_cli(FIRST_THREE + ["--estimate-only"] + _mysql_argv(mysql_config))
    line = [m for m in capsys.readouterr().out.splitlines() if m.startswith(run_common.ESTIMATE_PREFIX)]
    body = __import__("json").loads(line[0][len(run_common.ESTIMATE_PREFIX):])
    assert body["sectors"] == 3 and body["sampled"] == 2 and body["systems"] == 6
    assert _all_sectors(mysql_config) == []

    _run_cli(FIRST_THREE + _mysql_argv(mysql_config))
    assert [(r["shell_index"], r["band_index"], r["wedge_index"]) for r in _all_sectors(mysql_config)] \
        == [(1, 0, 0), (1, 0, 1), (1, 0, 2)]


# --- The web limit -------------------------------------------------------------------------------------------------


def test_a_web_neighbourhood_over_the_limit_is_refused_before_anything_is_listed(mysql_config, monkeypatch):
    _plan_wide_galaxy(mysql_config)
    _run_cli(FIRST_THREE + _mysql_argv(mysql_config))
    center = _all_sectors(mysql_config)[0]["id"]
    monkeypatch.setattr(run_common, "WEB_ENUMERATION_LIMIT", 5)
    monkeypatch.setattr(run_galaxy, "enumerate_sectors_within_radius",
                        lambda *a, **k: pytest.fail("the ball was enumerated"))
    with pytest.raises(run_galaxy.EnumerationTooLarge, match="more than the 5 a web request may list"):
        run_galaxy.generate_sector_neighborhood(center, radius_ly=40.0, config=mysql_config, estimate_only=True)


def test_the_neighbourhood_route_answers_413_for_a_ball_too_large_for_a_web_request(admin_client, monkeypatch):
    def too_large(*args, **kwargs):
        raise run_galaxy.EnumerationTooLarge("That ball holds 9,000,000 sectors")

    monkeypatch.setattr(run_galaxy, "generate_sector_neighborhood", too_large)
    sector = admin_client.post("/api/sectors", json={"name": "Hood", "edge_ly": 10.0}).get_json()["id"]
    response = admin_client.post(f"/api/sectors/{sector}/generate-neighborhood", json={"estimate_only": True})
    assert response.status_code == 413 and "9,000,000" in response.get_json()["error"]


def test_the_generate_pages_radial_estimate_stops_over_the_limit(mysql_config, monkeypatch, capsys):
    _plan_wide_galaxy(mysql_config)
    _run_cli(FIRST_THREE + _mysql_argv(mysql_config))
    center = _all_sectors(mysql_config)[0]["id"]
    monkeypatch.setattr(run_common, "WEB_ENUMERATION_LIMIT", 5)
    with pytest.raises(SystemExit):
        _run_cli(["--center-sector", str(center), "--radius-pc", "40", "--estimate-only"] + _mysql_argv(mysql_config))
    assert not [m for m in capsys.readouterr().out.splitlines() if m.startswith(run_common.ESTIMATE_PREFIX)]
