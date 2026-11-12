# tests/test_reparent.py

"""
GEN.583 (with TEST.244): the state-preserving reparent. A body keeps its absolute position and velocity when its parent changes; the total linear momentum and the orbital angular
momentum about the centre of mass do not change; the body goes to its new parent and back to where it was; a moon's velocity in its star's frame is the sum of its hops; a body that
would not be bound is raised before anything is changed; and a planet passed to a new parent keeps what it is and has the climate, the class and the spin of its new place.
"""

import math

import pytest

from planetgen.generation import reparent as gen_reparent
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, planets as planetPhysics, reparent as kin
from planetgen.util import draw

MSUN = constants.SOLAR_MASS_TO_KG
YEAR = constants.SECONDS_PER_YEAR


def circular(radius_au, primary_mass_kg, phase=0.0, tilt=0.0):
    """Position (AU) and velocity (AU/yr) of a body on a circle about a primary: in a plane tilted `tilt` about x."""
    speed = math.sqrt(kin.primary_mu_au3_per_year2(primary_mass_kg) / radius_au)
    x, y = radius_au * math.cos(phase), radius_au * math.sin(phase)
    vx, vy = -speed * math.sin(phase), speed * math.cos(phase)
    return (x, y * math.cos(tilt), y * math.sin(tilt)), (vx, vy * math.cos(tilt), vy * math.sin(tilt))


def close(a, b, rel=1e-12, scale=1.0):
    return all(abs(x - y) <= rel * scale for x, y in zip(a, b))


def test_the_new_relative_state_is_the_difference_and_the_absolute_state_is_kept():
    position, velocity = circular(1.0, MSUN, phase=0.7, tilt=0.2)
    parent_position, parent_velocity = (3.0, 0.5, -0.2), (0.0, -1.0, 0.1)
    state = kin.reparent_state(position, velocity, parent_position, parent_velocity, 0.8 * MSUN)
    assert close(kin.chain_sum([state.position, parent_position]), position)
    assert close(kin.chain_sum([state.velocity, parent_velocity]), velocity)
    assert state.mu == pytest.approx(kin.primary_mu_au3_per_year2(0.8 * MSUN))


def test_the_period_of_the_new_orbit_is_the_generators_own_kepler_law():
    position, velocity = circular(1.0, MSUN)
    state = kin.reparent_state(position, velocity, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), MSUN)
    assert state.bound and state.orbit["eccentricity"] == pytest.approx(0.0, abs=1e-9)
    assert state.orbit["period"] == pytest.approx(planetPhysics.calculate_orbital_period_years(state.orbit["semi_major_axis"], MSUN), rel=1e-12)


def test_going_to_a_new_parent_and_back_returns_the_state():
    position, velocity = circular(2.0, MSUN, phase=1.1, tilt=0.4)
    parent_position, parent_velocity = (5.0, -1.0, 0.3), (0.1, 0.9, -0.05)
    there = kin.reparent_state(position, velocity, parent_position, parent_velocity, 0.6 * MSUN)
    back = kin.reparent_state(there.position, there.velocity, tuple(-p for p in parent_position), tuple(-v for v in parent_velocity), MSUN)
    assert close(back.position, position, scale=2.0) and close(back.velocity, velocity, scale=6.0)


def test_a_moons_velocity_in_its_stars_frame_is_the_sum_of_its_hops():
    planet_position, planet_velocity = circular(1.0, MSUN, phase=0.3)
    moon_position, moon_velocity = circular(0.0026, 6.0e24, phase=2.0, tilt=0.1)
    assert kin.chain_sum([moon_velocity, planet_velocity]) == pytest.approx(tuple(a + b for a, b in zip(moon_velocity, planet_velocity)), rel=1e-15)
    assert kin.chain_sum([moon_position, planet_position]) == pytest.approx(tuple(a + b for a, b in zip(moon_position, planet_position)), rel=1e-15)


def test_a_body_that_would_leave_is_not_bound_and_is_reported():
    position, velocity = circular(1.0, MSUN)
    state = kin.reparent_state(position, velocity, (0.5, 0.0, 0.0), (0.0, -40.0, 0.0), 0.1 * MSUN)
    assert not state.bound and state.orbit["kind"] == "hyperbolic"


def test_the_barycentre_the_momentum_and_the_angular_momentum_do_not_change_with_the_parent_a_state_is_measured_from():
    masses = (MSUN, 0.5 * MSUN, 3.0e-6 * MSUN)
    positions = [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0), circular(1.0, MSUN, phase=0.4, tilt=0.1)[0]]
    velocities = [(0.0, 0.0, 0.0), (0.0, 2.0, 0.0), circular(1.0, MSUN, phase=0.4, tilt=0.1)[1]]
    # The planet now measured from the second star: the chain sum puts it back where it was.
    state = kin.reparent_state(positions[2], velocities[2], positions[1], velocities[1], 0.5 * MSUN)
    positions_again = [positions[0], positions[1], kin.chain_sum([state.position, positions[1]])]
    velocities_again = [velocities[0], velocities[1], kin.chain_sum([state.velocity, velocities[1]])]
    assert kin.barycentre(masses, positions_again) == pytest.approx(kin.barycentre(masses, positions), rel=1e-13)
    assert kin.momentum(masses, velocities_again) == pytest.approx(kin.momentum(masses, velocities), rel=1e-13)
    assert kin.angular_momentum(masses, positions_again, velocities_again) == pytest.approx(kin.angular_momentum(masses, positions, velocities), rel=1e-12)
    # and the angular momentum is that of the system's own motion: a boost of every velocity leaves it as it was
    boosted = [tuple(v + b for v, b in zip(velocity, (5.0, -3.0, 1.0))) for velocity in velocities]
    assert kin.angular_momentum(masses, positions, boosted) == pytest.approx(kin.angular_momentum(masses, positions, velocities), rel=1e-12)


# --- A generated planet and moon --------------------------------------------------------------------------------------------------------------------------------------

def a_system_with_a_planet_and_a_moon():
    for seed in range(1, 60):
        draw.set_run_seed(seed)
        config = SystemConfig()
        config.PLANETS = True
        system = StarSystem(system_config=config, galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            if planet.body_type != "a" and planet.moons and planet.eccentricity is not None and planet.star is system.star:
                return system, planet
    pytest.skip("no sampled system has a planet with a moon")


def test_a_planet_passed_to_a_second_star_keeps_its_place_in_space_and_what_it_is():
    system, planet = a_system_with_a_planet_and_a_moon()
    other = StarSystem(system_config=SystemConfig(), galactic_center_dist_ly=26000.0).star
    before = {"position": (planet.position_x, planet.position_y, planet.position_z), "velocity": (planet.velocity_x_kms, planet.velocity_y_kms, planet.velocity_z_kms),
              "mass": planet.mass, "radius": planet.radius, "birth": planet._birth_spin, "moons": len(planet.moons)}
    # The second star sits 0.3 of the planet's distance along x from the planet, and is moving so that the planet is on a near-circular orbit about it (0.9 of the circular speed).
    separation = 0.3 * planet.distance
    parent_position = tuple(c + d for c, d in zip(before["position"], (separation, 0.0, 0.0)))
    planet_velocity_au_yr = tuple(v / gen_reparent._AU_YR_TO_KM_S for v in before["velocity"])
    circular_speed = 0.9 * math.sqrt(kin.primary_mu_au3_per_year2(other.mass) / separation)
    parent_velocity = tuple(v - c for v, c in zip(planet_velocity_au_yr, (0.0, circular_speed, 0.0)))
    state = gen_reparent.reparent_planet(planet, parent_position, parent_velocity, other.mass, new_star=other)
    after = (planet.position_x, planet.position_y, planet.position_z)
    # Absolute position is kept to the roundoff of the elements the planet is rebuilt from (the old parent's frame is the sum of the new relative position and the new parent's).
    assert tuple(a + b for a, b in zip(after, parent_position)) == pytest.approx(before["position"], rel=1e-6, abs=1e-9)
    assert planet.star is other and planet.primary_mass_kg == other.mass
    assert planet.distance == pytest.approx(state.orbit["semi_major_axis"], rel=1e-9)
    assert planet.period == pytest.approx(planetPhysics.calculate_orbital_period_years(planet.distance, other.mass), rel=1e-12)
    assert planet.orbit_from_vector()["semi_major_axis"] == pytest.approx(planet.distance, rel=1e-6)
    assert planet.mass == before["mass"] and planet.radius == before["radius"] and planet._birth_spin == before["birth"] and len(planet.moons) == before["moons"]
    assert planet.zone == planetPhysics.zone_for(planet.habitable_zone, planet.distance, planet.star)


def test_a_planet_that_would_be_freed_is_raised_and_left_as_it_was():
    system, planet = a_system_with_a_planet_and_a_moon()
    before = (planet.distance, planet.eccentricity, planet.period, planet.star, planet.position_x)
    with pytest.raises(gen_reparent.Freed) as raised:
        gen_reparent.reparent_planet(planet, (planet.distance * 1.1, 0.0, 0.0), (0.0, 400.0, 0.0), 0.1 * MSUN)
    assert not raised.value.state.bound
    assert (planet.distance, planet.eccentricity, planet.period, planet.star, planet.position_x) == before


def test_a_moon_passed_to_another_planet_keeps_its_state_and_leaves_one_list_for_the_other():
    system, planet = a_system_with_a_planet_and_a_moon()
    moon = planet.moons[0]
    other = next((p for p in system.planets if p is not planet and p.body_type != "a"), None)
    if other is None:
        pytest.skip("the sampled system has one planet")
    moon_before = (moon.position_x, moon.position_y, moon.position_z)
    new_parent_position = tuple(0.5 * c for c in moon_before)       # a planet half way along the moon's line from the old planet, at rest relative to it
    other_mass, own_mass = other.mass, moon.mass
    state = gen_reparent.reparent_moon(moon, planet, other, new_parent_position, (0.0, 0.0, 0.0))
    assert moon in other.moons and moon not in planet.moons
    assert moon.primary_mass_kg == other_mass and moon.mass == own_mass
    assert tuple(a + b for a, b in zip((moon.position_x, moon.position_y, moon.position_z), new_parent_position)) == pytest.approx(moon_before, rel=1e-9)
    assert moon.distance == pytest.approx(state.orbit["semi_major_axis"], rel=1e-9) and moon._frame_axis is None


def _locked(body):
    """Whether the body turns once an orbit: its rotation period is its orbital period (hours)."""
    return body.rotation_period_hours and abs(body.rotation_period_hours - body.period * 8766.0) <= 1e-3 * body.rotation_period_hours


def test_a_locked_moons_day_follows_a_change_of_its_planet():
    for seed in range(1, 80):
        draw.set_run_seed(seed)
        config = SystemConfig()
        config.PLANETS = True
        system = StarSystem(system_config=config, galactic_center_dist_ly=26000.0)
        planets = [p for p in system.planets if p.body_type != "a"]
        for planet in planets:
            for moon in planet.moons:
                others = [p for p in planets if p is not planet]
                if others and _locked(moon):
                    other = others[0]
                    moon_position = (moon.position_x, moon.position_y, moon.position_z)
                    # the new planet sits at half the moon's distance on the far side of it, at rest relative to the old planet: the moon is now closer to it and its orbit is a different one
                    state = gen_reparent.reparent_moon(moon, planet, other, tuple(1.5 * c for c in moon_position), (0.0, 0.0, 0.0))
                    assert _locked(moon) or not state.bound
                    assert moon.period == pytest.approx(planetPhysics.calculate_orbital_period_years(moon.distance, other.mass), rel=1e-12)
                    return
    pytest.skip("no sampled moon is locked with a second planet in its system")
