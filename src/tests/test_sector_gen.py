"""
generate.generate_sector_name regression tests (generate.py's `sector`
subcommand).

Covers the two-word invariant: `generate_sector_name` joins two independent
`generate_phoneme_salad_name` calls into one name, and that inner function
can itself split a long result into two words (`split_long_word`) -- left
unchecked, that silently produced 3-4 word sector names instead of 2.

Run with: pytest tests/test_sector_gen.py
"""
from generate import generate_sector_name

TRIALS = 500


def test_generate_sector_name_is_always_two_words():
    for _ in range(TRIALS):
        name = generate_sector_name()
        words = name.split(" ")
        assert len(words) == 2, f"expected exactly 2 words, got {len(words)}: {name!r}"
        assert all(word for word in words), f"unexpected empty word in {name!r}"


# ---------------------------------------------------------------------------
# generate_sector_phenomena -- sector-level exotic phenomena generation
# (science-based Poisson rates, see program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM).
# ---------------------------------------------------------------------------

from types import SimpleNamespace

import pytest

import generate as sectorGen
from stellarObjects import program_constants
from stellarObjects.compactRemnant import BlackHole
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


def _make_cheap_system(star_type="G2V"):
    """A fast-to-generate system (no planets/moons) for seeding a sector
    without the cost of full planet/life generation -- this file's own
    tests only care about phenomenon counts/types, not system content."""
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.PLANETS = False
    return StarSystem(system_config=cfg), cfg


def _seeded_sector(system_count):
    sector = SpaceSector("Phenomena Test Sector")
    for _ in range(system_count):
        system, cfg = _make_cheap_system()
        sector.add_system(system, system_config=cfg)
    return sector


def test_generate_sector_phenomena_passes_rate_times_system_count_as_the_poisson_mean(monkeypatch):
    system_count = 7
    sector = _seeded_sector(system_count)

    captured_means = []

    def fake_sample_poisson_count(mean):
        captured_means.append(mean)
        return 0  # no actual placement needed for this test

    monkeypatch.setattr(sectorGen, "_sample_poisson_count", fake_sample_poisson_count)

    args = SimpleNamespace(markdown=False)
    entries = sectorGen.generate_sector_phenomena(sector, args)

    assert entries == []
    expected_means = [rate * system_count for rate in program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM.values()]
    assert sorted(captured_means) == pytest.approx(sorted(expected_means))


def test_generate_sector_phenomena_builds_the_right_type_and_count(monkeypatch):
    # Force only "nebula" (a non-massive type -- no Hill-sphere placement
    # risk, so this is fully deterministic) to draw a count, every other
    # type draws 0.
    def fake_sample_poisson_count(mean):
        # nebula's own rate is the smallest of the seven (see
        # PHENOMENON_RATE_PER_STAR_SYSTEM's own docstring) -- comparing the
        # mean directly, rather than hardcoding its value here, keeps this
        # test correct even if the cited rate is retuned later.
        nebula_rate = program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM["nebula"]
        return 3 if mean == pytest.approx(nebula_rate * 4) else 0

    monkeypatch.setattr(sectorGen, "_sample_poisson_count", fake_sample_poisson_count)

    sector = _seeded_sector(4)
    args = SimpleNamespace(markdown=False)
    entries = sectorGen.generate_sector_phenomena(sector, args)

    assert len(entries) == 3
    assert all(e.phenomenon_type == "nebula" for e in entries)
    assert entries == sector.phenomena
    for entry in entries:
        assert entry.phenomenon.name  # a real, generated Nebula object


def test_generate_sector_phenomena_threads_galactic_center_dist_ly_to_compact_remnants(monkeypatch):
    # Force exactly one black hole and confirm its own Hill-sphere/orbit
    # calculation actually used the sector's real distance from the
    # galactic center, not the fallback constant.
    def fake_sample_poisson_count(mean):
        black_hole_rate = program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM["black-hole"]
        return 1 if mean == pytest.approx(black_hole_rate * 1) else 0

    monkeypatch.setattr(sectorGen, "_sample_poisson_count", fake_sample_poisson_count)

    sector = _seeded_sector(1)
    args = SimpleNamespace(markdown=False)
    galactic_center_dist_ly = 5000.0
    entries = sectorGen.generate_sector_phenomena(sector, args, galactic_center_dist_ly=galactic_center_dist_ly)

    assert len(entries) == 1
    black_hole = entries[0].phenomenon
    assert isinstance(black_hole, BlackHole)
    assert black_hole.galactic_center_dist_ly == galactic_center_dist_ly
    assert black_hole.system_perimeter == pytest.approx(
        black_hole.calculate_system_perimeter(galactic_center_dist_ly)
    )


def test_generate_sector_phenomena_skips_a_massive_draw_that_cannot_be_placed(monkeypatch):
    # If SpaceSector.add_phenomenon can't fit a massive phenomenon
    # (Hill-sphere placement failure), generate_sector_phenomena must
    # silently skip that one draw rather than crashing the whole sector.
    def fake_sample_poisson_count(mean):
        black_hole_rate = program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM["black-hole"]
        return 2 if mean == pytest.approx(black_hole_rate * 1) else 0

    monkeypatch.setattr(sectorGen, "_sample_poisson_count", fake_sample_poisson_count)

    def fake_add_phenomenon(self, phenomenon, phenomenon_type, position=None, min_separation_ly=None):
        raise ValueError("no room -- simulated placement failure")

    monkeypatch.setattr(SpaceSector, "add_phenomenon", fake_add_phenomenon)

    sector = _seeded_sector(1)
    args = SimpleNamespace(markdown=False)
    entries = sectorGen.generate_sector_phenomena(sector, args)

    assert entries == []  # both draws failed to place and were skipped, no crash


def test_generate_sector_phenomena_honors_markdown_flag(monkeypatch):
    # A mean of 30 makes P(count == 0) ~= e^-30 (effectively zero) while
    # staying well within _sample_poisson_count's own "small means" scope
    # (its docstring flags O(mean) draws as inefficient much above a few
    # dozen) -- every other type pinned to 0 so only nebula draws at all.
    inflated = {key: 0.0 for key in program_constants.PHENOMENON_RATE_PER_STAR_SYSTEM}
    inflated["nebula"] = 30.0
    monkeypatch.setattr(program_constants, "PHENOMENON_RATE_PER_STAR_SYSTEM", inflated)

    sector = _seeded_sector(1)
    args = SimpleNamespace(markdown=True)
    entries = sectorGen.generate_sector_phenomena(sector, args)

    assert entries
    assert all(e.phenomenon_type == "nebula" for e in entries)
    assert all(e.phenomenon.system_config.MARKDOWN is True for e in entries)


# ---------------------------------------------------------------------------
# generate_sector -- capacity handling (see generate_sector's own
# "Capacity" docstring note). A sector's cube can only physically hold so
# many systems before every point left is within some existing system's
# own Hill sphere; SpaceSector.add_system raises ValueError once placement
# fails. generate_sector must stop and return what fit instead of letting
# that exception propagate and discard every system already generated.
# ---------------------------------------------------------------------------

def _real_sector_args(num_systems, name, extra_argv=()):
    """A real `sector` subcommand namespace, built via generate.py's own
    parser/validators -- generate_sector reads far more of `args` (via
    build_sector_configs/build_system_config) than a hand-built
    SimpleNamespace could safely stand in for."""
    parser, command_parsers = sectorGen.build_parser()
    args = parser.parse_args(["sector", "--num-systems", str(num_systems), "--name", name, *extra_argv])
    command_parser = command_parsers["sector"]
    sectorGen.validate_shared_generation_args(args, command_parser)
    sectorGen.validate_sector_args(args, command_parser)
    return args


def test_generate_sector_stops_gracefully_once_a_system_cannot_be_placed(monkeypatch):
    # Mirrors test_generate_sector_phenomena_skips_a_massive_draw_that_cannot_be_placed
    # above, one level up: a placement failure partway through the system
    # loop must not discard the systems already placed (each one an
    # expensive full StarSystem generation) or crash the whole sector.
    original_add_system = SpaceSector.add_system
    call_count = {"n": 0}

    def flaky_add_system(self, *args, **kwargs):
        call_count["n"] += 1
        if call_count["n"] > 3:
            raise ValueError("no room -- simulated placement failure")
        return original_add_system(self, *args, **kwargs)

    monkeypatch.setattr(SpaceSector, "add_system", flaky_add_system)

    args = _real_sector_args(10, "CapacityTestSector", extra_argv=["-planets"])
    _sector_name, sector = sectorGen.generate_sector(args)

    assert len(sector.entries) == 3  # stopped right after the 3 successful placements
    assert call_count["n"] == 4  # the 4th (failing) attempt, then no more


def test_generate_sector_does_not_generate_systems_past_the_first_placement_failure(monkeypatch):
    # The whole point of stopping early (rather than catching the
    # ValueError but still looping through every remaining config) is to
    # avoid paying for a full StarSystem generation -- the expensive part
    # -- on a config that's never going to fit anyway.
    original_add_system = SpaceSector.add_system
    call_count = {"n": 0}

    def flaky_add_system(self, *args, **kwargs):
        call_count["n"] += 1
        if call_count["n"] > 2:
            raise ValueError("no room -- simulated placement failure")
        return original_add_system(self, *args, **kwargs)

    monkeypatch.setattr(SpaceSector, "add_system", flaky_add_system)

    generated_count = {"n": 0}
    original_star_system_init = StarSystem.__init__

    def counting_init(self, *args, **kwargs):
        generated_count["n"] += 1
        return original_star_system_init(self, *args, **kwargs)

    monkeypatch.setattr(StarSystem, "__init__", counting_init)

    args = _real_sector_args(10, "NoWastedGenerationSector", extra_argv=["-planets"])
    sectorGen.generate_sector(args)

    # 2 successful placements + the 1 that triggered the failure -- not
    # all 10 requested configs.
    assert generated_count["n"] == 3


def test_generate_sector_real_capacity_overflow_returns_a_partial_sector():
    # End-to-end (no monkeypatching): a sector this size physically cannot
    # hold 120 systems at real Hill-sphere spacing (see this project's own
    # local-stellar-density constants), which used to raise ValueError and
    # discard the whole sector. 120, not 60: a lucky draw occasionally fit
    # all 60 (sampled fills run about 20-59). -planets keeps this fast (skips
    # planet/moon generation) without touching the real placement logic
    # this test actually cares about.
    args = _real_sector_args(120, "RealCapacityOverflowSector", extra_argv=["-planets"])
    _sector_name, sector = sectorGen.generate_sector(args)

    assert 0 < len(sector.entries) < 120
