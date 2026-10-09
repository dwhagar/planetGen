# tests/test_state_vectors.py

"""State vectors and orbital elements (planetgen/physics/state_vectors.py)."""

import math

import pytest

from planetgen.physics import orbits
from planetgen.physics import state_vectors as sv
from planetgen.physics.kepler import comet_orbital_state, gravitational_parameter_au3_yr2, true_anomaly_and_distance_elliptical

MU = 4.0 * math.pi ** 2  # AU^3/yr^2 for one solar mass


def _close(a, b, rel=1e-9, abs_=1e-9):
    return all(x == pytest.approx(y, rel=rel, abs=abs_) for x, y in zip(a, b))


def _angle_close(a, b, tol=1e-9):
    return abs((a - b + math.pi) % (2 * math.pi) - math.pi) < tol


@pytest.mark.parametrize("eccentricity", [0.0, 0.3, 0.97, 1.0, 1.8])
@pytest.mark.parametrize("inclination", [0.4, 1.3, 2.9])
def test_elements_state_elements_round_trip(eccentricity, inclination):
    node, argp, nu = 0.7, 2.1, 0.9
    position, velocity = sv.state_from_elements(MU, 1.7, eccentricity, inclination, node, argp, nu)
    back = sv.elements_from_state(position, velocity, MU)
    assert back["periapsis_distance"] == pytest.approx(1.7)
    assert back["eccentricity"] == pytest.approx(eccentricity, abs=1e-9)
    assert back["inclination"] == pytest.approx(inclination)
    assert _angle_close(back["ascending_node"], node)
    if eccentricity > 0.0:
        assert _angle_close(back["arg_periapsis"], argp)
        assert _angle_close(back["true_anomaly"], nu)
    else:  # a circle has no periapsis: the true anomaly counts from the node
        assert back["arg_periapsis"] == 0.0
        assert _angle_close(back["true_anomaly"], argp + nu)
    again = sv.state_from_elements(MU, back["periapsis_distance"], back["eccentricity"], back["inclination"],
                                   back["ascending_node"], back["arg_periapsis"], back["true_anomaly"])
    assert _close(again[0], position) and _close(again[1], velocity)


def test_the_kind_period_and_energy_follow_the_eccentricity():
    for e, kind in ((0.5, "elliptical"), (1.0, "parabolic"), (1.5, "hyperbolic")):
        position, velocity = sv.state_from_elements(MU, 2.0, e, 0.3, 0.2, 0.1, 0.4)
        elements = sv.elements_from_state(position, velocity, MU)
        assert elements["kind"] == kind
        assert (elements["period"] is not None) == (kind == "elliptical")
    position, velocity = sv.state_from_elements(MU, 2.0, 0.5, 0.3, 0.2, 0.1, 0.4)
    elements = sv.elements_from_state(position, velocity, MU)
    a = 2.0 / 0.5
    assert elements["semi_major_axis"] == pytest.approx(a)
    assert elements["period"] == pytest.approx(math.sqrt(a ** 3))  # years for one solar mass
    assert elements["specific_energy"] == pytest.approx(-MU / (2 * a))


def test_an_orbit_conserves_energy_and_angular_momentum_round_the_ellipse():
    energies, momenta = set(), set()
    for k in range(12):
        position, velocity = sv.state_from_elements(MU, 1.0, 0.8, 0.5, 0.3, 1.0, 2 * math.pi * k / 12)
        elements = sv.elements_from_state(position, velocity, MU)
        energies.add(round(elements["specific_energy"], 9))
        h = (position[1] * velocity[2] - position[2] * velocity[1], position[2] * velocity[0] - position[0] * velocity[2],
             position[0] * velocity[1] - position[1] * velocity[0])
        momenta.add(round(math.sqrt(sum(c * c for c in h)), 9))
    assert len(energies) == 1 and len(momenta) == 1


def test_a_circular_orbit_agrees_with_the_circular_position_and_velocity_helpers():
    distance, inclination, node, phase, period = 3.0, 17.0, 40.0, 123.0, 5.196
    mu = 4.0 * math.pi ** 2 * distance ** 3 / period ** 2
    position, velocity = sv.state_from_elements(mu, distance, 0.0, math.radians(inclination), math.radians(node),
                                                0.0, math.radians(phase))
    assert _close(position, orbits.orbital_position_au(distance, inclination, node, phase))
    assert _close(velocity, orbits.circular_orbital_velocity_au_per_year(distance, inclination, node, phase, period))
    assert math.sqrt(sum(c * c for c in velocity)) == pytest.approx(2 * math.pi * distance / period)


def test_the_circular_velocity_is_the_derivative_of_the_circular_position():
    distance, inclination, node, phase, period = 2.5, 33.0, 80.0, 200.0, 4.0
    step = 1.0e-6
    ahead = orbits.orbital_position_au(distance, inclination, node, phase + 360.0 * step / period)
    behind = orbits.orbital_position_au(distance, inclination, node, phase - 360.0 * step / period)
    numeric = tuple((a - b) / (2 * step) for a, b in zip(ahead, behind))
    assert _close(orbits.circular_orbital_velocity_au_per_year(distance, inclination, node, phase, period), numeric,
                  rel=1e-6, abs_=1e-6)
    with pytest.raises(ValueError):
        orbits.circular_orbital_velocity_au_per_year(distance, inclination, node, phase, 0.0)


def test_an_elliptical_comet_agrees_with_kepler_module():
    q, e, incl, argp, node, mass = 0.6, 0.9, 40.0, 70.0, 20.0, 1.0
    a = q / (1 - e)
    mean_anomaly = 1.1
    nu, distance = true_anomaly_and_distance_elliptical(mean_anomaly, e, a)
    state = comet_orbital_state("elliptical", q, e, incl, argp, node, mass, mean_anomaly_rad=mean_anomaly,
                                orbital_period_years=a ** 1.5)
    position, velocity = sv.state_from_elements(gravitational_parameter_au3_yr2(mass), q, e, math.radians(incl),
                                                math.radians(node), math.radians(argp), nu)
    assert _close(position, (state["position_x_au"], state["position_y_au"], state["position_z_au"]), rel=1e-7)
    assert math.sqrt(sum(c * c for c in velocity)) * 4.740470446 == pytest.approx(state["orbital_speed_kms"], rel=1e-5)
    assert sv.mean_anomaly_from_true(nu, e) == pytest.approx(mean_anomaly)


def test_an_untilted_or_retrograde_flat_orbit_keeps_its_direction():
    for incl in (0.0, math.pi):
        position, velocity = sv.state_from_elements(MU, 1.0, 0.4, incl, 0.0, 0.8, 0.5)
        back = sv.elements_from_state(position, velocity, MU)
        assert back["inclination"] == pytest.approx(incl)
        assert back["ascending_node"] == 0.0
        assert _angle_close(back["arg_periapsis"], 0.8)
        assert _angle_close(back["true_anomaly"], 0.5)


def test_closed_orbit_points_lie_on_the_ellipse():
    points = sv.closed_orbit_points(1.0, 0.5, 0.4, 0.3, 0.9, count=64)
    assert len(points) == 64
    distances = [math.sqrt(sum(c * c for c in p)) for p in points]
    assert min(distances) == pytest.approx(1.0) and max(distances) <= 3.0 + 1e-9
    assert distances[0] == pytest.approx(1.0)  # starts at periapsis
    with pytest.raises(ValueError):
        sv.closed_orbit_points(1.0, 1.0, 0.0, 0.0, 0.0)
    with pytest.raises(ValueError):
        sv.closed_orbit_points(1.0, 0.5, 0.0, 0.0, 0.0, count=2)


def test_bad_input_is_refused():
    with pytest.raises(ValueError):
        sv.state_from_elements(0.0, 1.0, 0.1, 0.0, 0.0, 0.0, 0.0)
    with pytest.raises(ValueError):
        sv.state_from_elements(MU, -1.0, 0.1, 0.0, 0.0, 0.0, 0.0)
    with pytest.raises(ValueError):
        sv.state_from_elements(MU, 1.0, -0.1, 0.0, 0.0, 0.0, 0.0)
    with pytest.raises(ValueError):  # beyond the asymptote of a hyperbola
        sv.state_from_elements(MU, 1.0, 2.0, 0.0, 0.0, 0.0, math.radians(170))
    with pytest.raises(ValueError):  # straight at the primary: no orbital plane
        sv.elements_from_state((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), MU)
    with pytest.raises(ValueError):
        sv.elements_from_state((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), MU)
    with pytest.raises(ValueError):
        sv.elements_from_state((1.0, float("nan"), 0.0), (0.0, 1.0, 0.0), MU)
    with pytest.raises(ValueError):
        sv.elements_from_state((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), -1.0)
    with pytest.raises(ValueError):
        sv.mean_anomaly_from_true(0.3, 1.2)
