# tests/test_bughunt_sector_placement.py

"""
Tier 1 bug-hunt coverage: `spaceSector.SpaceSector.add_system` under
stress -- many real (large, `+large_star`-biased, wide-Hill-sphere)
systems crammed into a small sector, checking the module's own central
invariant (its docstring's "Minimum separation (Hill spheres)" section):
every successfully-placed pair must be at least `required_separation_ly`
apart, and a sector with no more room must raise the documented
`ValueError` (`SECTOR_MAX_PLACEMENT_ATTEMPTS` exhausted) rather than
silently placing an overlapping system or hanging (the placement loop is
already a bounded `for` loop, not `while True`, so a hang isn't expected
here -- this is confirming that bound actually holds under real stress,
not just trusting the source).
"""

import pytest

from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector, distance_between, required_separation_ly
from stellarObjects.systemData import StarSystem

from tests.bughunt_support import FUZZ_SEEDS, run_seeded

_SEEDS = FUZZ_SEEDS[:15]


def _make_large_system(star_type="B3V"):
    """A deliberately large-Hill-sphere system -- a hot, massive star
    forced large, to make overlap the *likely* outcome in a small sector
    rather than something a fuzz loop might never actually trigger."""
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    cfg.LARGE_STAR = True
    cfg.STAR_TYPE = star_type
    return StarSystem(system_config=cfg)


def test_placed_systems_never_violate_required_separation():
    """Fill a small sector with large-Hill-sphere systems until it either
    fills up (ValueError, the documented/expected outcome) or a fixed
    attempt cap is reached; whichever happens, every pair that DID get
    placed must satisfy required_separation_ly."""
    def check(seed):
        sector = SpaceSector(name=f"stress-{seed}", edge_ly=50.0)
        placed = 0
        try:
            for _ in range(12):
                system = _make_large_system()
                sector.add_system(system)
                placed += 1
        except ValueError:
            pass  # documented, expected failure mode once the sector fills

        entries = sector.entries
        for i in range(len(entries)):
            for j in range(i + 1, len(entries)):
                a, b = entries[i], entries[j]
                dist = distance_between(a.position, b.position)
                required = required_separation_ly(a.star_system, b.star_system)
                assert dist >= required - 1e-9, (
                    f"seed={seed}: entries {i}/{j} placed {dist} ly apart, "
                    f"required {required} ly -- Hill spheres overlap"
                )
        assert placed >= 1, f"seed={seed}: not even one system could be placed in a fresh 50 ly sector"

    run_seeded(check, seeds=_SEEDS)


def test_sector_full_raises_clean_value_error_not_hang():
    """A sector far too small for even two large systems must raise the
    documented ValueError quickly, not loop forever or silently accept an
    overlapping placement."""
    import time

    sector = SpaceSector(name="tiny", edge_ly=0.001)
    first = _make_large_system("O5V")
    sector.add_system(first)

    second = _make_large_system("O5V")
    start = time.monotonic()
    with pytest.raises(ValueError):
        sector.add_system(second)
    elapsed = time.monotonic() - start
    assert elapsed < 10.0, f"add_system exhaustion took {elapsed}s -- expected a fast, bounded failure"


def test_explicit_position_bypasses_hill_sphere_check_as_documented():
    """An explicit position is documented to skip the automatic
    separation check entirely -- confirms that's still true (a caller
    deliberately placing two systems close together, e.g. a hand-crafted
    scenario, must not be silently blocked)."""
    sector = SpaceSector(name="explicit", edge_ly=100.0)
    a = _make_large_system("O5V")
    b = _make_large_system("O5V")
    sector.add_system(a, position=(0.0, 0.0, 0.0))
    sector.add_system(b, position=(0.0001, 0.0, 0.0))  # deliberately overlapping
    assert len(sector.entries) == 2
