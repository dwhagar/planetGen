# tests/test_ejected_losses.py

"""GEN.1111 part 1: a system keeps the planet architecture each of its stars was built from, so the bodies it threw out can be read from it."""

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw


def _system(**overrides):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return StarSystem(cfg)


def test_a_system_keeps_the_architecture_of_its_star_and_reads_its_thrown_bodies_from_it():
    with draw.bound("ejected-losses-1"):
        system = _system()
    assert len(system.architectures) == 1
    star, architecture = next(iter(system.architectures.values()))
    assert star is system.star and architecture.planets is not None
    for owner, loss in system.ejected_losses():
        assert owner is system.star and loss.kind == "thrown" and loss.ejected > 0.0 and loss.perturber_au is not None


def test_a_system_with_no_planets_throws_nothing_out():
    with draw.bound("ejected-losses-2"):
        system = _system(PLANETS=False)
    assert system.ejected_losses() == []


def test_some_systems_throw_bodies_out_and_each_loss_names_where_and_by_whom():
    found = []
    for seed in range(16):
        with draw.bound(f"ejected-losses-pop-{seed}"):
            found.extend(_system().ejected_losses())
    assert found, "a giant-hosting G star in sixteen throws some body out"
    for _star, loss in found:
        assert loss.mass_kg > 0.0 and loss.distance_au > 0.0 and loss.perturber_mass_kg > 0.0 and loss.hill_radius_au > 0.0
        assert 0.0 < loss.ejected <= 1.0


def test_a_system_rebuilt_from_storage_has_no_plan_and_so_no_ejections():
    with draw.bound("ejected-losses-3"):
        system = _system()
    rebuilt = StarSystem.from_dict(system.to_dict())
    assert rebuilt.ejected_losses() == []
