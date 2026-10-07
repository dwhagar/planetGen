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
# (research-based Poisson rates, see tuning.PHENOMENON_DENSITY_PC3).
# ---------------------------------------------------------------------------

from types import SimpleNamespace

import pytest

import generate as sectorGen
from planetgen import tuning
from planetgen.generation.phenomena.compact_remnant import BlackHole
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem


def _make_cheap_system(star_type="G2V"):
    """A fast-to-generate system (no planets/moons) for seeding a sector
    without the cost of full planet/life generation -- this file's own
    tests only care about phenomenon counts/types, not system content."""
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.PLANETS = False
    return StarSystem(system_config=cfg), cfg


def _seeded_sector(system_count, edge_ly=100.0):
    # Wider than the default sector: a forced phenomenon must always find
    # room clear of the seeded systems' Hill spheres, which in a default
    # 13 ly sector it failed to about once in 700 runs.
    sector = SpaceSector("Phenomena Test Sector", edge_ly=edge_ly)
    for _ in range(system_count):
        system, cfg = _make_cheap_system()
        sector.add_system(system, system_config=cfg)
    return sector


def _rate(kind):
    return tuning.phenomenon_rate_per_star(kind)


def test_generate_sector_phenomena_passes_rate_per_star_times_star_count_as_the_poisson_mean(monkeypatch):
    sector = _seeded_sector(7)
    star_count = sectorGen.sector_star_count(sector)
    assert star_count >= 7

    captured_means = []

    def fake_sample_poisson_count(mean):
        captured_means.append(mean)
        return 0  # no actual placement needed for this test

    monkeypatch.setattr(sectorGen, "_sample_poisson_count", fake_sample_poisson_count)

    args = SimpleNamespace(markdown=False)
    entries = sectorGen.generate_sector_phenomena(sector, args)

    assert entries == []
    expected_means = [_rate(kind) * star_count for kind, _type, _factory in sectorGen.SECTOR_PHENOMENON_KINDS]
    assert sorted(captured_means) == pytest.approx(sorted(expected_means))


def test_research_rates_per_star():
    """Boss's research densities at 0.14 stars per pc^3."""
    assert _rate("rogue-planet") == pytest.approx(6.5)
    assert _rate("brown-dwarf") == pytest.approx(0.03 / 0.14)
    assert _rate("neutron-star") == pytest.approx(0.005)
    assert _rate("black-hole") == pytest.approx(0.001)
    assert _rate("comet") == pytest.approx(0.05)
    assert _rate("asteroid-field") == 0.0
    assert _rate("runaway-star") == pytest.approx(0.015)


def test_rate_scale_dials_a_kind(monkeypatch):
    monkeypatch.setitem(tuning.PHENOMENON_RATE_SCALE, "rogue-planet", 0.1)
    assert _rate("rogue-planet") == pytest.approx(0.65)


def _only(kind, count):
    """A fake Poisson draw returning `count` for `kind` and 0 otherwise
    (matched by the mean, since each kind's rate differs)."""
    def fake(mean):
        return count if mean == pytest.approx(_rate(kind) * fake.stars) and mean > 0 else 0
    return fake


def test_generate_sector_phenomena_builds_the_right_type_and_count(monkeypatch):
    sector = _seeded_sector(4)
    fake = _only("planetary-nebula", 3)
    fake.stars = sectorGen.sector_star_count(sector)
    monkeypatch.setattr(sectorGen, "_sample_poisson_count", fake)

    args = SimpleNamespace(markdown=False)
    entries = sectorGen.generate_sector_phenomena(sector, args)

    assert len(entries) == 3
    assert all(e.phenomenon_type == "nebula" for e in entries)
    assert all(e.phenomenon.nebula_type == "planetary" for e in entries)
    assert entries == sector.phenomena
    # Each planetary nebula sits on its own new hot white dwarf system.
    assert len(sector.entries) == 4 + 3
    for entry in entries:
        host = next(e for e in sector.entries if e.position == entry.position)
        assert sectorGen._spectral_code(host.star_system.stars[0]) in tuning.PLANETARY_NEBULA_CENTRAL_STAR_TYPES


def test_molecular_clouds_are_dark_family_nebulae(monkeypatch):
    sector = _seeded_sector(2)
    fake = _only("molecular-cloud", 2)
    fake.stars = sectorGen.sector_star_count(sector)
    monkeypatch.setattr(sectorGen, "_sample_poisson_count", fake)

    entries = sectorGen.generate_sector_phenomena(sector, SimpleNamespace(markdown=False))
    assert [e.phenomenon.nebula_type for e in entries] == ["dark", "dark"]
    assert all(e.phenomenon.nebula_class in "MNPQ" for e in entries)


def test_o_stars_sit_in_h_ii_regions_and_cool_stars_light_nothing(monkeypatch):
    monkeypatch.setattr(sectorGen, "_sample_poisson_count", lambda mean: 0)
    sector = SpaceSector("Hosted Nebula Sector")
    for star_type in ("O5V", "G2V", "M3V"):
        system, cfg = _make_cheap_system(star_type)
        cfg.BINARY_SYSTEM = False
        sector.add_system(system, system_config=cfg)
    hot = sector.entries[0]
    assert sectorGen._spectral_code(hot.star_system.stars[0]) == "O5V"

    entries = sectorGen.generate_sector_phenomena(sector, SimpleNamespace(markdown=False))
    assert len(entries) == 1
    assert entries[0].position == hot.position
    assert entries[0].phenomenon.nebula_class in ("C", "D", "E")


def test_brown_dwarfs_are_rogue_planet_rows(monkeypatch):
    sector = _seeded_sector(2)
    fake = _only("brown-dwarf", 2)
    fake.stars = sectorGen.sector_star_count(sector)
    monkeypatch.setattr(sectorGen, "_sample_poisson_count", fake)

    entries = sectorGen.generate_sector_phenomena(sector, SimpleNamespace(markdown=False))
    assert [e.phenomenon_type for e in entries] == ["rogue-planet", "rogue-planet"]
    assert all(e.phenomenon.mass_bin == "brown-dwarf" for e in entries)


def test_generate_sector_phenomena_threads_galactic_center_dist_ly_to_compact_remnants(monkeypatch):
    # Force exactly one black hole and confirm its own Hill-sphere/orbit
    # calculation actually used the sector's real distance from the
    # galactic center, not the fallback constant.
    sector = _seeded_sector(1)
    fake = _only("black-hole", 1)
    fake.stars = sectorGen.sector_star_count(sector)
    monkeypatch.setattr(sectorGen, "_sample_poisson_count", fake)

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
    sector = _seeded_sector(1)
    fake = _only("black-hole", 2)
    fake.stars = sectorGen.sector_star_count(sector)
    monkeypatch.setattr(sectorGen, "_sample_poisson_count", fake)

    def fake_add_phenomenon(self, phenomenon, phenomenon_type, position=None, min_separation_ly=None):
        raise ValueError("no room -- simulated placement failure")

    monkeypatch.setattr(SpaceSector, "add_phenomenon", fake_add_phenomenon)

    args = SimpleNamespace(markdown=False)
    entries = sectorGen.generate_sector_phenomena(sector, args)

    assert entries == []  # both draws failed to place and were skipped, no crash


def test_generate_sector_phenomena_honors_markdown_flag(monkeypatch):
    sector = _seeded_sector(1)
    fake = _only("comet", 4)
    fake.stars = sectorGen.sector_star_count(sector)
    monkeypatch.setattr(sectorGen, "_sample_poisson_count", fake)

    args = SimpleNamespace(markdown=True)
    entries = sectorGen.generate_sector_phenomena(sector, args)

    assert len(entries) == 4
    assert all(e.phenomenon.system_config.MARKDOWN is True for e in entries)


def test_a_local_density_sector_gets_about_six_and_a_half_rogues_per_star():
    sector = _seeded_sector(9)
    stars = sectorGen.sector_star_count(sector)
    counts = []
    for _ in range(5):
        sector.phenomena.clear()
        entries = sectorGen.generate_sector_phenomena(sector, SimpleNamespace(markdown=False))
        counts.append(sum(1 for e in entries if e.phenomenon_type == "rogue-planet"))
    mean = sum(counts) / len(counts)
    expected = (_rate("rogue-planet") + _rate("brown-dwarf")) * stars
    assert 0.6 * expected < mean < 1.4 * expected


def test_flag_fast_stars(monkeypatch):
    sector = _seeded_sector(20)
    monkeypatch.setitem(tuning.PHENOMENON_RATE_SCALE, "runaway-star", 1 / 0.015)
    assert sectorGen.flag_fast_stars(sector, galactic_center_dist_ly=26000.0) == 20
    for entry in sector.entries:
        system = entry.star_system
        assert system.runaway_class in ("runaway", "hypervelocity")
        if system.runaway_class == "runaway":
            low, high = tuning.RUNAWAY_STAR_SPEED_RANGE_KMS
            assert low <= system.runaway_speed_kms <= high

    # Hypervelocity stars crowd the center (r^-2).
    sector = _seeded_sector(5)
    monkeypatch.setitem(tuning.PHENOMENON_RATE_SCALE, "runaway-star", 0.0)
    monkeypatch.setitem(tuning.PHENOMENON_RATE_SCALE, "hypervelocity-star", 1e12)
    assert sectorGen.flag_fast_stars(sector, galactic_center_dist_ly=1.0) == 5
    assert all(e.star_system.runaway_class == "hypervelocity" for e in sector.entries)
    assert all(500 <= e.star_system.runaway_speed_kms <= 1000 for e in sector.entries)


# ---------------------------------------------------------------------------
# generate_sector -- capacity handling (see generate_sector's own
# "Capacity" docstring note). A sector's cube can only physically hold so
# many systems before every point left is within some existing system's
# own Hill sphere; SpaceSector.add_system raises ValueError once placement
# fails. generate_sector must stop and return what fit instead of letting
# that exception propagate and discard every system already generated.
# ---------------------------------------------------------------------------

def _real_sector_args(num_systems, name, planets=None):
    """A real `sector` subcommand namespace, built via generate.py's own
    parser/validators -- generate_sector reads far more of `args` (via
    build_sector_configs/build_system_config) than a hand-built
    SimpleNamespace could safely stand in for. `planets=False` makes every
    system planet-less (cheap); `sector` takes no `-planets` option any
    more (GEN.51), so it is set on the namespace build_system_config reads."""
    parser, command_parsers = sectorGen.build_parser()
    args = parser.parse_args(["sector", "--num-systems", str(num_systems), "--name", name])
    command_parser = command_parsers["sector"]
    sectorGen.validate_shared_generation_args(args, command_parser)
    sectorGen.validate_sector_args(args, command_parser)
    args.planets = planets
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

    args = _real_sector_args(10, "CapacityTestSector", planets=False)
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

    args = _real_sector_args(10, "NoWastedGenerationSector", planets=False)
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
    args = _real_sector_args(120, "RealCapacityOverflowSector", planets=False)
    _sector_name, sector = sectorGen.generate_sector(args)

    assert 0 < len(sector.entries) < 120
