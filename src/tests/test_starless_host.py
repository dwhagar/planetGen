# tests/test_starless_host.py

"""GEN.1111 slice 1: the host a rogue planet is built against and the planet built through the planet code with it."""

import copy

import pytest

from planetgen.generation.phenomena import starless_host
from planetgen.generation.planet import Planet
from planetgen.physics import constants
from planetgen.util import draw
from tests.test_gen_builder_internals import make_system

EARTH_KG = constants.EARTH_MASS_TO_KG


def _donor():
    system = make_system()
    return system, system.star


def test_the_host_gives_no_light_the_rogues_age_and_the_donors_everything_else():
    system, donor = _donor()
    host = starless_host.StarlessHost(donor, 3.5)
    assert host.luminosity == 0.0 and host.age == 3.5
    assert host.mass == donor.mass and host.snow_line_au == donor.snow_line_au
    assert host.habitable_zone == donor.habitable_zone
    assert donor.luminosity > 0.0, "the donor is untouched"


def test_a_host_can_be_copied_and_a_missing_attribute_is_an_attribute_error():
    _system, donor = _donor()
    host = starless_host.StarlessHost(donor, 1.0)
    assert copy.copy(host).age == 1.0
    with pytest.raises(AttributeError):
        host.no_such_attribute


_UNTIL = ("the air of a body with no star: the atmosphere table's shoreline row (d1) is a logistic in log10 of the insolation and of the host's luminosity, "
          "which a body with no light does not have (GEN.1111 lock wiring: a no-host rule for the table rows that read the star)")


@pytest.mark.xfail(strict=True, reason=_UNTIL)
def test_an_ejected_planet_is_a_planet_of_the_given_mass_with_no_light_on_it():
    system, donor = _donor()
    with draw.bound("starless-host-planet"):
        planet = starless_host.build_ejected_planet(system.system_config, donor, 5.0, 3.0 * EARTH_KG, moon_count=1)
    assert isinstance(planet, Planet) and isinstance(planet.star, starless_host.StarlessHost)
    assert planet.mass == 3.0 * EARTH_KG
    assert planet.energy_flux_w_m2 < 1e-6, "no starlight reaches it"
    assert len(planet.moons) == 1 and planet.moons[0].is_moon


@pytest.mark.xfail(strict=True, reason=_UNTIL)
def test_the_formation_distance_the_loss_gives_is_the_one_the_planet_keeps():
    system, donor = _donor()
    with draw.bound("starless-host-formation"):
        planet = starless_host.build_ejected_planet(system.system_config, donor, 4.0, 2.0 * EARTH_KG, formation_au=6.0, moon_count=0)
    assert planet.formation_distance_au == pytest.approx(6.0)


# --- the moons a thrown-out body keeps (RP4) -------------------------------------------------------------------------------------------

def _thrown(**over):
    from planetgen.physics.planet_architecture import Loss
    fields = dict(kind="thrown", count=1, mass_kg=300.0 * EARTH_KG, core_kg=10.0 * EARTH_KG, distance_au=5.0, formation_au=5.0, time_yr=1e7,
                  perturber_au=3.0, perturber_mass_kg=3000.0 * EARTH_KG, theta=8.0, hill_radius_au=0.35, ice_rich=True, ejected=1.0,
                  collided=0.0, bound=0.0)
    fields.update(over)
    return Loss(**fields)


def test_the_moon_limit_is_never_wider_than_the_zone_the_moons_had_before_the_ejection():
    from planetgen.physics import clearance_ledger
    loss = _thrown()
    window = loss.hill_radius_au * constants.AU_TO_KM * clearance_ledger.moon_hill_fraction(0.0)[0]
    for u in (0.0, 0.1, 0.5, 0.9, 1.0):
        limit = starless_host.ejected_moon_limit_km(loss, loss.mass_kg, 1.99e30, u)
        assert 0.0 < limit <= window * (1 + 1e-12)


def test_a_closer_encounter_leaves_a_narrower_moon_system():
    loss = _thrown()
    limits = [starless_host.ejected_moon_limit_km(loss, loss.mass_kg, 1.99e30, u) for u in (0.05, 0.3, 0.6, 0.95)]
    assert limits == sorted(limits) and limits[0] < limits[-1]


def test_the_limit_is_the_survival_distance_of_the_closest_approach_drawn():
    from planetgen.physics import rogue_moons
    loss = _thrown()
    mass = loss.mass_kg
    mutual = rogue_moons.mutual_hill_radius_m(5.0 * constants.AU_TO_M, 3.0 * constants.AU_TO_M, mass, loss.perturber_mass_kg, 1.99e30)
    r_min_m = rogue_moons.closest_approach_hill(0.2) * mutual
    expected_km = rogue_moons.survival_distance_m(r_min_m, mass, loss.perturber_mass_kg) / 1000.0
    window_km = loss.hill_radius_au * constants.AU_TO_KM * 0.4895
    got = starless_host.ejected_moon_limit_km(loss, mass, 1.99e30, 0.2)
    assert got == pytest.approx(min(expected_km, window_km), rel=1e-3)


def test_a_loss_with_no_scatterer_has_no_moon_limit():
    with pytest.raises(ValueError):
        starless_host.ejected_moon_limit_km(_thrown(kind="swept", perturber_au=None, perturber_mass_kg=None), 1e25, 1.99e30, 0.5)


# --- the rogue body of an ejection event (part 3) -------------------------------------------------------------------------------------

def _event(donor, **over):
    from planetgen.generation import ejection_events
    fields = dict(donor_uid="u", index=0, mass_kg=300.0 * EARTH_KG, core_kg=10.0 * EARTH_KG, formation_au=5.0, distance_au=5.0, perturber_au=3.0,
                  perturber_mass_kg=3000.0 * EARTH_KG, ice_rich=True, delay_gy=0.5, age_gy=donor.age, since_ejection_gy=donor.age - 0.5,
                  kick_km_s=(1.0, 0.0, 0.0), orbital_velocity_km_s=(0.0, 1.0, 0.0), moon_limit_km=2.0e6)
    fields.update(over)
    return ejection_events.EjectionEvent(**fields)


def test_the_host_of_an_event_is_the_donor_with_the_rogues_age():
    _system, donor = _donor()
    event = _event(donor, age_gy=4.0)
    host = starless_host.host_of(donor, event)
    assert host.age == 4.0 and host.luminosity == 0.0 and host.mass == donor.mass


def _moon_distances_km(donor, system, limit_km, seeds=6):
    found = []
    for seed in range(seeds):
        event = _event(donor, mass_kg=5.0 * EARTH_KG, core_kg=5.0 * EARTH_KG, moon_limit_km=limit_km)
        with draw.bound(f"rogue-body-moons-{seed}"):
            # a real star as host: the atmosphere table cannot read a body with no light yet
            planet = starless_host.build_rogue_body(system.system_config, donor, event, moon_count=2)
        found.extend(moon.distance * constants.AU_TO_KM for moon in planet.moons)   # a moon's distance is in AU
    return found


def test_the_moons_of_a_rogue_body_stay_inside_the_edge_the_ejection_left():
    from planetgen.physics import clearance_ledger
    system, donor = _donor()
    event = _event(donor, mass_kg=5.0 * EARTH_KG, core_kg=5.0 * EARTH_KG, moon_limit_km=2.0e5)
    with draw.bound("rogue-body-hill"):
        planet = starless_host.build_rogue_body(system.system_config, donor, event, moon_count=2)
    assert planet.hill_radius == pytest.approx(event.moon_limit_km / clearance_ledger.moon_hill_fraction(0.0)[0])
    distances = _moon_distances_km(donor, system, 2.0e5)
    assert distances, "a limit of 2e5 km leaves room for moons"
    assert max(distances) <= 2.0e5 * (1 + 1e-9)


def test_an_ejection_that_left_no_room_leaves_no_moons_and_a_wide_one_lets_them_be_wider():
    system, donor = _donor()
    assert _moon_distances_km(donor, system, 4.0e4) == [], "inside the planet's own clearance there is no room for a moon"
    narrow, wide = _moon_distances_km(donor, system, 2.0e5), _moon_distances_km(donor, system, 2.0e6)
    assert max(wide) > max(narrow)
