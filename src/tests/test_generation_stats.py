# tests/test_generation_stats.py

"""
`planetgen.generation.stats` (PERF.3, PERF.10): density buckets,
decaying averages, the stored speeds and sizes, the size and time
estimate, and the disk-space refusal.

Tests that take the `mysql_config` fixture (see `conftest.py`) need a
MySQL server and skip without one.
"""

import json
import math

import pytest

import generate
from planetgen.db import store
from planetgen.generation import stats as generationStats
from planetgen.generation.stats import (
    Bucket, DiskSpace, GenerationStats, bucket_bounds, bucket_index, check_disk, estimate, format_bytes,
    format_duration,
)
from tests.test_galaxy_gen import _all_sectors, _mysql_argv, _plan_wide_galaxy, _run_cli


# ---------------------------------------------------------------------------
# Buckets
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("density, index", [
    (None, 0), (float("nan"), 0), (0.0, 0), (0.005, 0), (0.01, 0), (0.03, 0), (0.0317, 1),
    (0.1, 2), (1.0, 4), (3.0, 4), (3.17, 5), (10.0, 6), (1e4, 12),
])
def test_densities_fall_in_half_decade_buckets(density, index):
    assert bucket_index(density) == index


def test_every_bucket_contains_its_own_bounds():
    for index in range(0, 14):
        low, high = bucket_bounds(index)
        assert bucket_index(low * 1.0001) == index
        assert bucket_index(high * 0.9999) == index
        assert math.isclose(high / low, 10 ** (1 / generationStats.BUCKETS_PER_DECADE))


def test_a_bucket_starts_at_its_first_sample_then_decays():
    bucket = Bucket("sector", 4)
    bucket.add(1.2, seconds=2.0, systems=10, stars=13)
    assert bucket.seconds_per_task == 2.0
    assert bucket.seconds_per_system == 0.2
    assert bucket.stars_per_system == 1.3
    bucket.add(1.5, seconds=4.0, systems=10, stars=10)
    decay = generationStats.DECAY
    assert math.isclose(bucket.seconds_per_task, 2.0 + decay * 2.0)
    assert math.isclose(bucket.seconds_per_system, 0.2 + decay * 0.2)
    assert math.isclose(bucket.stars_per_system, 1.3 - decay * 0.3)
    assert bucket.samples == 2 and bucket.max_density == 1.5


def test_an_empty_sector_counts_its_time_as_one_system():
    bucket = Bucket("sector", 0)
    bucket.add(0.02, seconds=0.5, systems=0, stars=0)
    assert bucket.seconds_per_system == 0.5
    assert bucket.stars_per_system == 1.0


# ---------------------------------------------------------------------------
# Estimates
# ---------------------------------------------------------------------------

def test_an_unmeasured_server_estimates_from_the_defaults():
    stats = GenerationStats()
    result = estimate([(1.0, 20.0), (1.0, 10.0)], stats, "galaxy", workers=1)
    assert result.sectors == 2 and result.systems == 30.0
    assert not result.measured
    assert math.isclose(result.seconds, 30 * generationStats.DEFAULT_SECONDS_PER_SYSTEM)
    assert result.bytes == math.ceil(30 * generationStats.DEFAULT_BYTES_PER_SYSTEM * 1.1)
    assert "not yet measured" in result.summary()


def test_measured_buckets_drive_the_time_and_workers_divide_it():
    stats = GenerationStats()
    stats.record("sector", 1.0, seconds=1.0, systems=10, stars=10)     # 0.1 s/system at density ~1
    stats.record("sector", 10.0, seconds=10.0, systems=20, stars=20)   # 0.5 s/system at density ~10
    result = estimate([(1.0, 10.0), (10.0, 20.0)], stats, "galaxy", workers=2)
    assert result.measured
    assert math.isclose(result.seconds, (10 * 0.1 + 20 * 0.5) / 2)
    # An unmeasured bucket borrows its nearest measured neighbor's speed.
    assert stats.seconds_per_system("sector", 30.0) == 0.5
    assert stats.seconds_per_system("sector", 0.05) == 0.1


def test_more_workers_than_sectors_dont_shorten_the_estimate():
    stats = GenerationStats()
    one = estimate([(1.0, 10.0)], stats, "galaxy", workers=1)
    many = estimate([(1.0, 10.0)], stats, "galaxy", workers=8)
    assert one.seconds == many.seconds


def test_a_measured_size_replaces_the_default():
    stats = GenerationStats()
    stats.sizes["galaxy"] = {"bytes_per_system": 1000.0, "systems": 1, "total_bytes": 1000}
    assert estimate([(1.0, 100.0)], stats, "galaxy").bytes == 110000
    assert estimate([(1.0, 100.0)], stats, "other").bytes > 110000


# ---------------------------------------------------------------------------
# Disk space
# ---------------------------------------------------------------------------

GB = 10 ** 9


def _sized(bytes_):
    result = generationStats.Estimate(sectors=1, systems=1.0)
    result.bytes = bytes_
    return result


def test_a_run_that_fits_is_not_refused():
    result = _sized(1 * GB)
    assert check_disk(result, DiskSpace("/var/lib/mysql", 100 * GB, 50 * GB)) is None
    assert result.as_dict()["refused"] is False


def test_more_than_a_quarter_of_the_disk_is_refused():
    result = _sized(26 * GB)
    refusal = check_disk(result, DiskSpace("/var/lib/mysql", 100 * GB, 90 * GB))
    assert refusal.startswith("Refused") and "quarter" in refusal
    assert result.as_dict()["refused"] is True


def test_leaving_less_than_5_gb_free_is_refused():
    result = _sized(2 * GB)
    refusal = check_disk(result, DiskSpace("/var/lib/mysql", 100 * GB, 6 * GB))
    assert refusal.startswith("Refused") and "5.0 GB" in refusal
    assert "1.0 GB more free space" in refusal


def test_an_unmeasurable_disk_refuses_nothing():
    result = _sized(10 ** 15)
    assert check_disk(result, None) is None


def test_a_remote_database_has_no_disk_to_check():
    class _Conn:
        def execute(self, *_args):
            raise AssertionError("not asked")

    assert generationStats.database_disk(_Conn(), "db.example.com") is None


@pytest.mark.parametrize("count, text", [
    (0, "0 bytes"), (999, "999 bytes"), (1500, "1.5 KB"), (12_000_000, "12 MB"), (2.5e9, "2.5 GB"), (3e12, "3.0 TB"),
])
def test_format_bytes(count, text):
    assert format_bytes(count) == text


@pytest.mark.parametrize("seconds, text", [
    (0, "0 s"), (12, "12 s"), (250, "4 m 10 s"), (7500, "2 h 05 m"), (90000, "1 day 1 h"), (3 * 86400, "3 days 0 h"),
])
def test_format_duration(seconds, text):
    assert format_duration(seconds) == text


# ---------------------------------------------------------------------------
# Stored in the control database
# ---------------------------------------------------------------------------

@pytest.fixture
def control_config(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    conn.close()
    return mysql_config


def test_buckets_and_sizes_survive_a_round_trip(control_config):
    stats = GenerationStats(control_config)
    assert stats.available
    stats.record("sector", 1.0, seconds=2.0, systems=10, stars=13)
    stats.record("sector", 12.0, seconds=5.0, systems=50, stars=60)
    stats.record("scatter", 0.5, seconds=3.0, systems=100, stars=100)
    stats.sizes["galaxy"] = {"bytes_per_system": 5000.0, "systems": 10, "total_bytes": 50000}
    stats._dirty_sizes.add("galaxy")
    stats.flush()

    again = GenerationStats(control_config)
    assert again.rows() == stats.rows()
    assert again.bytes_per_system("galaxy") == 5000.0
    again.record("sector", 1.0, seconds=4.0, systems=10, stars=13)
    again.flush()
    third = GenerationStats(control_config)
    bucket = third.buckets[("sector", bucket_index(1.0))]
    assert bucket.samples == 2
    assert math.isclose(bucket.seconds_per_task, 2.0 + generationStats.DECAY * 2.0)


def test_the_size_is_measured_from_the_galaxy_database(control_config):
    conn = store.get_connection(control_config)
    try:
        stats = GenerationStats(control_config)
        assert stats.measure_size(conn, control_config.database) is None   # no systems yet
        conn.execute("CREATE TEMPORARY TABLE star_systems (id INT)")   # hides the real, empty one
        conn.execute("INSERT INTO star_systems VALUES (1), (2)")
        size = stats.measure_size(conn, control_config.database)
    finally:
        conn.close()
    assert size["systems"] == 2 and size["bytes_per_system"] > 0
    assert size["bytes_per_system"] * 2 == size["total_bytes"]


def test_without_the_stats_tables_nothing_is_recorded(mysql_config):
    stats = GenerationStats(mysql_config)       # a galaxy database: no generation_stats table
    assert not stats.available
    stats.record("sector", 1.0, seconds=1.0, systems=1, stars=1)
    stats.flush()                                # no error
    assert stats.seconds_per_system("sector", 1.0) == 1.0


# ---------------------------------------------------------------------------
# generate.py: the estimate before a bulk run, the refusal, the question,
# and the speed recorded after it
# ---------------------------------------------------------------------------

RING_0 = ["--ring", "0", "--num-systems", "2"]   # ring 0 holds 3 slots


def _estimate_line(capsys):
    lines = [line for line in capsys.readouterr().out.splitlines() if line.startswith(generate.ESTIMATE_PREFIX)]
    assert len(lines) == 1
    return json.loads(lines[0][len(generate.ESTIMATE_PREFIX):])


def test_estimate_only_writes_nothing(mysql_config, capsys):
    _plan_wide_galaxy(mysql_config)
    _run_cli(RING_0 + ["--estimate-only"] + _mysql_argv(mysql_config))
    estimate_ = _estimate_line(capsys)
    assert estimate_["sectors"] == 3 and estimate_["systems"] == 6
    assert estimate_["what"] == "ring 0 layer 0" and estimate_["refused"] is False
    assert _all_sectors(mysql_config) == []


def test_a_run_the_disk_cant_hold_is_refused_before_anything_is_written_under_strict(mysql_config, monkeypatch,
                                                                                      capsys):
    _plan_wide_galaxy(mysql_config)
    monkeypatch.setattr(generationStats, "database_disk",
                        lambda conn, host: DiskSpace("/data", 100 * GB, 5 * GB))
    with pytest.raises(SystemExit) as exit_info:
        _run_cli(RING_0 + ["--yes", "--strict"] + _mysql_argv(mysql_config))
    assert exit_info.value.code == 1
    assert _all_sectors(mysql_config) == []
    assert "at least 5.0 GB must stay free" in capsys.readouterr().out


def test_a_run_the_disk_cant_hold_warns_and_goes_ahead(mysql_config, monkeypatch, capsys):
    # GEN.81: the console never says no; it warns and does what was asked.
    _plan_wide_galaxy(mysql_config)
    monkeypatch.setattr(generationStats, "database_disk",
                        lambda conn, host: DiskSpace("/data", 100 * GB, 5 * GB))
    _run_cli(RING_0 + ["--yes"] + _mysql_argv(mysql_config))
    assert len(_all_sectors(mysql_config)) == 3
    out = capsys.readouterr().out
    assert "WARNING:" in out and "at least 5.0 GB must stay free" in out


def test_estimate_only_reports_a_refusal(mysql_config, monkeypatch, capsys):
    _plan_wide_galaxy(mysql_config)
    monkeypatch.setattr(generationStats, "database_disk",
                        lambda conn, host: DiskSpace("/data", 100 * GB, 5 * GB))
    _run_cli(RING_0 + ["--estimate-only"] + _mysql_argv(mysql_config))
    assert _estimate_line(capsys)["refused"] is True


def test_a_terminal_is_asked_first_and_no_stops_the_run(mysql_config, monkeypatch):
    _plan_wide_galaxy(mysql_config)
    monkeypatch.setattr(generate, "_interactive", lambda: True)
    questions = []
    monkeypatch.setattr("builtins.input", lambda prompt: questions.append(prompt) or "n")
    with pytest.raises(SystemExit) as exit_info:
        _run_cli(RING_0 + _mysql_argv(mysql_config))
    assert exit_info.value.code == 0
    assert questions == ["Generate these 3 sectors? [y/N] "]
    assert _all_sectors(mysql_config) == []

    monkeypatch.setattr("builtins.input", lambda prompt: questions.append(prompt) or "y")
    _run_cli(RING_0 + _mysql_argv(mysql_config))
    assert len(_all_sectors(mysql_config)) == 3


def test_yes_skips_the_question(mysql_config, monkeypatch):
    _plan_wide_galaxy(mysql_config)
    monkeypatch.setattr(generate, "_interactive", lambda: True)
    monkeypatch.setattr("builtins.input", lambda prompt: pytest.fail("asked"))
    _run_cli(RING_0 + ["--yes"] + _mysql_argv(mysql_config))
    assert len(_all_sectors(mysql_config)) == 3


def test_a_run_records_its_speed_and_the_databases_size(control_config, monkeypatch):
    monkeypatch.setenv(generationStats.STATS_ENV_VAR, "1")
    monkeypatch.setenv("PLANETGEN_CONTROL_DATABASE", control_config.database)
    # Without the GEN.23 bright-star backfill each sector holds exactly
    # its --num-systems, so the per-task counts below are exact.
    monkeypatch.setattr(generate, "backfill_bright_stars", lambda *args, **kwargs: {"blocks": 0, "stars": 0})
    _plan_wide_galaxy(control_config)
    _run_cli(RING_0 + _mysql_argv(control_config))
    stats = GenerationStats(control_config)
    (bucket,) = stats.rows("sector")
    assert bucket["samples"] == 3 and bucket["seconds_per_task"] > 0
    assert bucket["systems_per_task"] == pytest.approx(2.0)
    assert stats.sizes[control_config.database]["systems"] == 6
    assert generate._RUN_STATS == {}


def test_the_web_neighborhood_can_be_estimated_and_refused(mysql_config, monkeypatch):
    _plan_wide_galaxy(mysql_config)
    _run_cli(RING_0 + _mysql_argv(mysql_config))
    center = _all_sectors(mysql_config)[0]["id"]
    result = generate.generate_sector_neighborhood(center, radius_ly=40.0, config=mysql_config, estimate_only=True)
    assert result["generated"] == 0 and result["estimate"]["sectors"] >= 1
    assert len(_all_sectors(mysql_config)) == 3

    monkeypatch.setattr(generationStats, "database_disk",
                        lambda conn, host: DiskSpace("/data", 100 * GB, 5 * GB))
    with pytest.raises(generate.GenerationRefused, match="must stay free"):
        generate.generate_sector_neighborhood(center, radius_ly=40.0, config=mysql_config)
    assert len(_all_sectors(mysql_config)) == 3


def test_absurd_densities_stay_finite():
    assert generationStats.bucket_index(float("inf")) == generationStats.MAX_BUCKET
    assert generationStats.bucket_index(1e308) == generationStats.MAX_BUCKET
    assert generationStats.bucket_index(float("nan")) == 0
    stats = generationStats.GenerationStats()
    result = generationStats.estimate([(float("inf"), float("inf"))], stats, "galaxy")
    assert result.systems == generationStats.MAX_SYSTEMS_PER_SECTOR
    assert math.isfinite(result.seconds) and result.bytes > 0
    assert generationStats.format_bytes(float("inf")) == "more than any disk holds"
    assert generationStats.format_duration(float("inf")) == "longer than anyone will wait"
