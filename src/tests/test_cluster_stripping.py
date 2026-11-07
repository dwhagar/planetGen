# tests/test_cluster_stripping.py

"""
GEN.1288: a star born in a cluster has the planets the cluster pulls off their orbits taken out of its system and made again as rogues, with the note that its birth cluster unbound them.
`cluster_unbinding.history_of` and `unbinding_epoch_myr` give the fate of a planet from the cluster's density and speed through its life; `planet_architecture.strip` moves the planet from the
architecture's planets to its losses; `StarSystem` reads `CLUSTER_HISTORY`; `ejection_events` and `rogue_yield` make the rogue.
"""

import math

import pytest

from planetgen.generation import ejection_events, rogue_yield
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena import starless_host
from planetgen.generation.system import StarSystem
from planetgen.physics import cluster_unbinding as cu, constants, planet_architecture as pa, rogue_kick
from planetgen.util import draw

EARTH = constants.EARTH_MASS_TO_KG
SUN = constants.SOLAR_MASS_TO_KG
DENSE = ((100.0, 1.0e5, 20.0),)
"""A history that unbinds every planet beyond a few AU: 100 Myr at 1e5 solar masses per cubic parsec and 20 km/s."""


def test_a_clusters_history_is_pieces_from_the_end_of_the_disk_to_its_age_with_a_density_that_falls():
    history = cu.history_of(300.0, 3.0, 0.2)
    assert len(history) == cu.HISTORY_PIECES
    assert sum(piece[0] for piece in history) == pytest.approx(200.0 - cu.DISK_LIFETIME_MYR)
    densities = [piece[1] for piece in history]
    assert densities == sorted(densities, reverse=True) and densities[0] > 10.0 * densities[-1]
    assert history[-1][1] == pytest.approx(cu.half_mass_density(300.0, 3.0))       # the radius has stopped growing by then: the present one


def test_a_cluster_no_older_than_the_disk_has_no_history():
    assert cu.history_of(300.0, 3.0, 0.002) == ()


def test_the_virial_dispersion_of_a_plummer_cluster_is_the_textbook_one():
    assert cu.virial_dispersion_km_s(1.0e4, 3.0) == pytest.approx(math.sqrt(4.30091e-3 * 1.0e4 / (6.0 * 1.30 * 3.0)), rel=1e-3)


def test_the_draws_that_unbind_a_planet_are_the_share_the_history_unbinds_over_its_soft_pieces():
    history = cu.history_of(3000.0, 2.0, 0.1)
    a, mass = 1000.0, 1.0
    draws = [(k + 0.5) / 2000.0 for k in range(2000)]
    gone = sum(1 for u in draws if cu.unbinding_epoch_myr(a, mass, history, u) is not None) / len(draws)
    soft = tuple(piece for piece in history if cu.soft(a, mass, piece[2]))
    assert soft and gone == pytest.approx(cu.unbound_fraction(a, soft), abs=2e-3)


def test_the_time_a_planet_is_pulled_off_grows_with_its_draw_and_lies_inside_the_clusters_life():
    history = cu.history_of(3000.0, 2.0, 0.1)
    times = [cu.unbinding_epoch_myr(1000.0, 1.0, history, u) for u in (0.01, 0.1, 0.3, 0.5)]
    times = [t for t in times if t is not None]
    assert len(times) >= 2 and times == sorted(times) and all(cu.DISK_LIFETIME_MYR <= t <= 100.0 for t in times)


def test_a_tight_orbit_is_hard_against_the_cluster_and_is_never_pulled_off():
    history = cu.history_of(300.0, 3.0, 0.2)
    assert not cu.soft(1.0, 1.0, history[-1][2])
    assert all(cu.unbinding_epoch_myr(1.0, 1.0, history, u) is None for u in (0.1, 0.5, 0.99))


def test_a_draw_outside_zero_to_one_is_refused():
    with pytest.raises(cu.ClusterUnbindingError):
        cu.unbinding_epoch_myr(1000.0, 1.0, DENSE, 1.0)


def _seed(distance, mass_earth):
    return pa.Seed(distance, mass_earth * EARTH, 0.1 * mass_earth * EARTH, False, 1, False, distance)


def test_strip_moves_a_planet_from_the_planets_to_the_losses_and_the_chain_mass_still_adds_up():
    near, far = _seed(1.0, 1.0), _seed(40.0, 10.0)
    lost = pa.Lost(0.0, 0.0, 0, 0.0, 0, (), (), (), 2.7, SUN, ())
    architecture = pa.Architecture((near, far), 1.0, 12.0, 2, lost)
    stripped = pa.strip(architecture, [(far, 5.0e6)])
    assert stripped.planets == (near,) and stripped.lost.stripped == ((far, 5.0e6),)
    [loss] = [each for each in pa.losses(stripped) if each.kind == "stripped"]
    assert loss.count == 1 and loss.mass_kg == far.mass_kg and loss.distance_au == 40.0 and loss.time_yr == 5.0e6
    assert loss.ejected == 1.0 and loss.collided == 0.0 and loss.hill_radius_au > 0.0 and loss.perturber_au is None
    assert pa.strip(architecture, []) is architecture


def _system(history, seed):
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    cfg.CLUSTER_HISTORY = history
    with draw.bound(seed):
        return StarSystem(system_config=cfg)


def test_a_system_with_no_cluster_history_is_the_system_it_always_was_and_a_dense_one_loses_its_wide_planets():
    found = None
    for seed in range(40):
        plain, dense = _system(None, seed), _system(DENSE, seed)
        if not plain.ejected_losses() and any(kind == "stripped" for _star, kind in [(s, l.kind) for s, l in dense.ejected_losses()]):
            found = plain, dense
            break
    assert found, "no system of 40 seeds had a planet to strip"
    plain, dense = found
    assert len(dense.planets) < len(plain.planets)
    assert all(loss.kind == "stripped" for _star, loss in dense.ejected_losses())
    again = _system(None, seed)
    assert [p.distance for p in again.planets] == [p.distance for p in plain.planets]


def test_a_stripped_loss_makes_a_rogue_that_leaves_with_the_unbinding_speed_and_says_its_birth_cluster_unbound_it():
    donor = _system(DENSE, 3)
    loss = pa.Loss("stripped", 1, 100.0 * EARTH, 10.0 * EARTH, 40.0, 40.0, 5.0e6, None, None, None, 4.0, True, 1.0, 0.0, 0.0)
    [event] = ejection_events.events("uid", donor.star, [loss], draw.Stream("t"))
    assert event.cause == "cluster" and event.delay_gy == pytest.approx(5.0e-3)
    circular = rogue_kick.circular_speed_km_s(donor.star.mass, 40.0)
    assert math.sqrt(sum(c * c for c in event.kick_km_s)) == pytest.approx(cu.unbinding_kick_km_s(circular))
    assert event.moon_limit_km == pytest.approx(starless_host.stripped_moon_limit_km(loss)) and event.moon_limit_km > 0.0
    note = rogue_yield.note_of("Test", event)
    assert note.startswith("Unbound by its birth cluster from the Test system")
    assert rogue_yield.note_of("Test", event._replace(cause="thrown")).startswith("Thrown out of the Test system")
