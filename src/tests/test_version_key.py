# tests/test_version_key.py

"""
DB.6: the 22-hex-digit version key, what made the galaxy (stored with
its seed), and the run history (`generation_runs`, one row per
`generate.py` run that changes the galaxy).
"""

import json
import secrets
import sys

import pytest

import generate
from stellarObjects import _db, galaxySeed, versionKey
from stellarObjects._version import __version__
from tests.test_galaxy_gen import _mysql_argv, _seed_skeleton

SEED = bytes.fromhex("00112233445566778899AABBCCDDEEFF")


# ---------------------------------------------------------------------------
# The key
# ---------------------------------------------------------------------------

def test_the_design_example():
    assert versionKey.version_key("7.127.352", (3, 12, 3), "Linux", "x86_64") == "0007007F000160030C0300"


def test_each_part_has_its_own_digits():
    key = versionKey.version_key("1.2.3", (3, 9, 20), "Windows", "AMD64")
    assert len(key) == versionKey.KEY_DIGITS
    assert (key[0:4], key[4:8], key[8:14]) == ("0001", "0002", "000003")
    assert (key[14:16], key[16:18], key[18:20]) == ("03", "09", "14")
    assert (key[20], key[21]) == ("1", "0")


def test_releases_whose_parts_sum_alike_get_different_keys():
    assert versionKey.version_key("7.127.352", (3, 12, 3), "Linux", "x86_64") != \
        versionKey.version_key("7.128.351", (3, 12, 3), "Linux", "x86_64")


@pytest.mark.parametrize("system, digit", [("Linux", 0), ("Windows", 1), ("Darwin", 3), ("FreeBSD", 2),
                                           ("OpenBSD", 2), ("SunOS", 2), ("Plan9", 0xF), ("", 0xF)])
def test_os_digits(system, digit):
    assert versionKey.os_digit(system) == digit


@pytest.mark.parametrize("machine, digit", [("x86_64", 0), ("AMD64", 0), ("arm64", 1), ("aarch64", 1),
                                            ("i686", 2), ("armv7l", 3), ("riscv64", 4), ("ppc64le", 0xF)])
def test_architecture_digits(machine, digit):
    assert versionKey.arch_digit(machine) == digit


@pytest.mark.parametrize("version", ["7.127", "7.127.352.1", "7.x.352", "65536.0.0", "1.0.16777216"])
def test_a_version_that_doesnt_fit_is_refused(version):
    with pytest.raises(ValueError):
        versionKey.version_key(version, (3, 12, 3), "Linux", "x86_64")


def test_the_running_code_has_a_key():
    current = versionKey.current()
    assert current["planetgen_version"] == __version__
    assert current["version_key"] == versionKey.version_key()
    major, revision, build = (int(part) for part in __version__.split("."))
    assert current["version_key"][:14] == f"{major:04X}{revision:04X}{build:06X}"
    assert current["version_key"][14:20] == "".join(f"{part:02X}" for part in sys.version_info[:3])


# ---------------------------------------------------------------------------
# What made the galaxy
# ---------------------------------------------------------------------------

def _maker(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        return _db.get_galaxy_maker(conn)
    finally:
        conn.close()


def test_a_new_galaxy_records_what_made_it_and_a_re_plan_keeps_it(mysql_config, monkeypatch):
    _seed_skeleton(mysql_config, galaxy_seed=SEED)
    made = _maker(mysql_config)
    assert made == _db.GalaxyMaker(**versionKey.current())

    monkeypatch.setattr(versionKey, "__version__", "99.0.0")
    _seed_skeleton(mysql_config)  # same seed kept
    assert _maker(mysql_config) == made
    _seed_skeleton(mysql_config, galaxy_seed=SEED)  # same seed given
    assert _maker(mysql_config) == made
    _seed_skeleton(mysql_config, galaxy_seed=b"\x01" * 16)  # a new galaxy
    assert _maker(mysql_config).planetgen_version == "99.0.0"


# ---------------------------------------------------------------------------
# The run history
# ---------------------------------------------------------------------------

def _runs(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        return conn.execute("SELECT * FROM generation_runs ORDER BY id").fetchall()
    finally:
        conn.close()


def _generate(mysql_config, monkeypatch, *argv):
    monkeypatch.setenv("PLANETGEN_CONTROL_DATABASE", mysql_config.database)
    monkeypatch.setattr(sys, "argv", ["generate.py", *argv] + _mysql_argv(mysql_config))
    generate.main()


def test_every_run_that_changes_the_galaxy_is_recorded(mysql_config, monkeypatch):
    monkeypatch.setattr(secrets, "randbits", lambda bits: 0xABC)
    _generate(mysql_config, monkeypatch, "plan", "--no-bright-stars", "--max-ring", "40",
              "--seed", galaxySeed.format_seed(SEED))
    _generate(mysql_config, monkeypatch, "galaxy", "--ring", "1", "--limit", "1", "--num-systems", "1",
              "--backfill-from", "none")
    runs = _runs(mysql_config)
    assert [run["command"] for run in runs] == ["plan", "galaxy"]
    plan, galaxy = runs
    assert json.loads(plan["arguments"]) == ["plan", "--no-bright-stars", "--max-ring", "40",
                                             "--seed", galaxySeed.format_seed(SEED)]
    assert "--mysql-password" not in galaxy["arguments"]
    # The plan ran before the galaxy had a seed; the galaxy run against it.
    assert bytes(galaxy["galaxy_seed"]) == SEED
    for run in runs:
        assert bytes(run["run_seed"]) == (0xABC).to_bytes(16, "big")
        assert run["version_key"] == versionKey.version_key()
        assert run["planetgen_version"] == __version__
        assert run["outcome"] == "ok"
        assert run["finished_at"] >= run["started_at"]


def test_a_failed_run_says_so_and_an_estimate_writes_nothing(mysql_config, monkeypatch):
    with pytest.raises(SystemExit):
        # No plan yet: the galaxy run refuses.
        _generate(mysql_config, monkeypatch, "galaxy", "--ring", "1", "--limit", "1", "--num-systems", "1")
    assert [(run["command"], run["outcome"]) for run in _runs(mysql_config)] == [("galaxy", "failed")]

    _seed_skeleton(mysql_config, layers=[(0, 5)])
    _generate(mysql_config, monkeypatch, "galaxy", "--ring", "1", "--limit", "1", "--num-systems", "1",
              "--estimate-only")
    assert len(_runs(mysql_config)) == 1


def test_the_history_never_stops_a_run(mysql_config, monkeypatch):
    def broken(*_args, **_kwargs):
        raise RuntimeError("no history table")

    monkeypatch.setattr(_db, "start_generation_run", broken)
    _seed_skeleton(mysql_config, layers=[(0, 5)])
    _generate(mysql_config, monkeypatch, "galaxy", "--ring", "1", "--limit", "1", "--num-systems", "1",
              "--backfill-from", "none")
    assert _runs(mysql_config) == []
