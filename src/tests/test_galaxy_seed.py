# tests/test_galaxy_seed.py

"""
GEN.39: the galaxy's 128-bit seed -- typed and shown as 32 hex digits,
stored in `galaxy_shape.galaxy_seed` bit for bit, kept by every later
plan, and the root of every unit's seed (SHA-256 of the seed and the
unit's `kind:address`). The same-sectors-at-any-worker-count checks are
in `test_parallel_generation.py`.
"""

import hashlib
import random
import sys

import pytest

import generate
from stellarObjects import _db
from planetgen.galaxy import seed as galaxySeed
from tests.test_galaxy_gen import _SKELETON_SHAPE, _mysql_argv, _seed_skeleton

SEED_HEX = "00112233445566778899AABBCCDDEEFF"
SEED = bytes.fromhex(SEED_HEX)


# ---------------------------------------------------------------------------
# The seed and the unit seeds
# ---------------------------------------------------------------------------

def test_a_seed_reads_and_prints_as_32_hex_digits():
    assert galaxySeed.parse_seed(SEED_HEX) == SEED
    assert galaxySeed.parse_seed(SEED_HEX.lower()) == SEED
    assert galaxySeed.format_seed(SEED) == SEED_HEX
    assert len(galaxySeed.new_seed()) == 16


@pytest.mark.parametrize("text", ["", SEED_HEX[:-1], SEED_HEX + "0", "0x" + SEED_HEX[2:], "G" + SEED_HEX[1:],
                                  SEED_HEX[:16] + " " + SEED_HEX[16:]])
def test_anything_but_32_hex_digits_is_refused(text):
    with pytest.raises(ValueError, match="32 hex digits"):
        galaxySeed.parse_seed(text)


def test_a_unit_seed_is_the_sha256_of_the_seed_and_its_kind_and_address():
    expected = int.from_bytes(hashlib.sha256(SEED + b"sector:12/3/0").digest(), "big")
    assert galaxySeed.unit_seed(SEED, "sector", (12, 3, 0)) == expected
    assert galaxySeed.unit_seed(SEED, "sector", "12/3/0") == expected
    assert galaxySeed.unit_seed(SEED, "sector", (12, 3, 1)) != expected
    assert galaxySeed.unit_seed(SEED, "bright-stars", (12, 3, 0)) != expected
    assert 0 <= galaxySeed.short_seed(SEED, "bright-stars", "scatter") < 2 ** 63
    with pytest.raises(ValueError):
        galaxySeed.unit_seed(SEED[:8], "sector", (0, 0, 0))


def test_every_bit_of_the_seed_reaches_a_unit_seed():
    base = galaxySeed.unit_seed(SEED, "sector", (1, 0, 0))
    for bit in range(128):
        flipped = bytearray(SEED)
        flipped[bit // 8] ^= 0x80 >> (bit % 8)
        assert galaxySeed.unit_seed(bytes(flipped), "sector", (1, 0, 0)) != base


def test_a_seeded_unit_draws_from_its_own_seed_and_puts_the_stream_back():
    random.seed(5)
    expected_after = random.random()
    random.seed(5)
    with galaxySeed.seeded(SEED, "sector", (1, 0, 0)):
        drawn = [random.random() for _ in range(3)]
    assert random.random() == expected_after
    random.seed(galaxySeed.unit_seed(SEED, "sector", (1, 0, 0)))
    assert drawn == [random.random() for _ in range(3)]


def test_with_no_galaxy_seed_a_unit_leaves_the_stream_alone():
    random.seed(9)
    expected = [random.random() for _ in range(2)]
    random.seed(9)
    with galaxySeed.seeded(None, "sector", (1, 0, 0)):
        first = random.random()
    assert [first, random.random()] == expected


def test_the_bright_star_seeds_come_from_the_galaxy_seed():
    skeleton = _db.GalaxySkeletonInfo(_SKELETON_SHAPE, 1.0, 1, 1.0, SEED)
    scatter = generate._bright_star_seed(skeleton, "scatter")
    assert scatter == galaxySeed.short_seed(SEED, "bright-stars", "scatter")
    assert generate._bright_star_seed(skeleton, "band/100-500") != scatter


# ---------------------------------------------------------------------------
# Storage and `plan --seed`
# ---------------------------------------------------------------------------

def _stored_seed(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        return _db.get_galaxy_seed(conn)
    finally:
        conn.close()


def _plan(mysql_config, monkeypatch, *extra):
    monkeypatch.setenv("PLANETGEN_CONTROL_DATABASE", mysql_config.database)
    monkeypatch.setattr(sys, "argv", ["generate.py", "plan", "--no-bright-stars", "--max-ring", "40", *extra]
                        + _mysql_argv(mysql_config))
    generate.main()


def test_the_seed_is_stored_bit_for_bit(mysql_config):
    for seed in (SEED, b"\xff" * 16, b"\x00" * 16, bytes(range(16))):
        _seed_skeleton(mysql_config, galaxy_seed=seed)
        assert _stored_seed(mysql_config) == seed
        conn = _db.get_connection(mysql_config)
        try:
            assert _db.get_galaxy_shape(conn).galaxy_seed == seed
        finally:
            conn.close()


def test_plan_stores_its_seed_and_logs_it(mysql_config, monkeypatch, capsys):
    _plan(mysql_config, monkeypatch, "--seed", SEED_HEX.lower())
    assert _stored_seed(mysql_config) == SEED
    assert f"Galaxy seed {SEED_HEX} (from --seed)" in capsys.readouterr().out


def test_a_plan_without_seed_draws_one_and_a_re_plan_keeps_it(mysql_config, monkeypatch, capsys):
    _plan(mysql_config, monkeypatch)
    drawn = _stored_seed(mysql_config)
    assert drawn is not None and len(drawn) == 16
    assert f"Galaxy seed {galaxySeed.format_seed(drawn)} (drawn at random)" in capsys.readouterr().out
    _plan(mysql_config, monkeypatch, "--arm-count", "3")
    assert _stored_seed(mysql_config) == drawn
    assert "(kept from the earlier plan)" in capsys.readouterr().out


def test_a_new_seed_replaces_the_old_only_before_any_sector(mysql_config, monkeypatch):
    _plan(mysql_config, monkeypatch, "--seed", SEED_HEX)
    other = "F" * 32
    _plan(mysql_config, monkeypatch, "--seed", other)
    assert _stored_seed(mysql_config) == bytes.fromhex(other)
    monkeypatch.setattr(sys, "argv", ["generate.py", "galaxy", "--ring", "1", "--limit", "1", "--num-systems", "1",
                                      "--backfill-from", "none"] + _mysql_argv(mysql_config))
    generate.main()
    with pytest.raises(SystemExit) as raised:
        _plan(mysql_config, monkeypatch, "--seed", SEED_HEX)
    assert raised.value.code == 1
    assert _stored_seed(mysql_config) == bytes.fromhex(other)
    # The same seed again is no change, so it plans.
    _plan(mysql_config, monkeypatch, "--seed", other)


def test_plan_refuses_a_malformed_seed(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["generate.py", "plan", "--seed", "1234"])
    with pytest.raises(SystemExit) as raised:
        generate.main()
    assert raised.value.code == 2
    assert "32 hex digits" in capsys.readouterr().err
