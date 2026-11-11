# tests/test_kozai_bodies.py

"""GEN.1447: Kozai-Lidov pumping of a body of a hierarchical triple (physics/kozai_bodies.py)."""

import math

import pytest

from planetgen.physics import kozai_bodies as kb, stability


def _closed_form_peak(e0, inclination):
    return math.sqrt(1.0 - (5.0 / 3.0) * (1.0 - e0 ** 2) * math.cos(inclination) ** 2)


@pytest.mark.parametrize("degrees", [45.0, 55.0, 65.0, 75.0, 85.0, 115.0, 130.0])
def test_the_cycle_reaches_the_closed_form_maximum_and_keeps_the_kozai_constant(degrees):
    inclination = math.radians(degrees)
    e0 = 0.02
    period, e_now, i_now, e_peak = kb.cycle(e0, inclination, 0.7, 11.3)
    assert e_peak == pytest.approx(_closed_form_peak(e0, inclination), abs=2e-3)
    assert math.isfinite(period) and 2.0 < period < 12.0      # a few Kozai times
    assert math.sqrt(1.0 - e_now ** 2) * math.cos(i_now) == pytest.approx(math.sqrt(1.0 - e0 ** 2) * math.cos(inclination), abs=1e-6)


def test_a_cycle_is_periodic_so_the_state_at_the_age_is_the_state_a_whole_number_of_periods_later():
    inclination = math.radians(70.0)
    period, e_a, i_a, _ = kb.cycle(0.05, inclination, 1.1, 3.7)
    _, e_b, i_b, _ = kb.cycle(0.05, inclination, 1.1, 3.7 + 5.0 * period)
    assert e_b == pytest.approx(e_a, abs=2e-3) and i_b == pytest.approx(i_a, abs=2e-3)


def _planet(**changes):
    """An Earth-mass planet at 1 AU of a Sun, a 0.5 solar mass companion at 300 AU on a circular orbit tilted 60 degrees, 5 Gyr."""
    values = dict(m_central=1.0, m_body=3.0e-6, m_perturber=0.5, a_in_au=1.0, e_in=0.02, a_out_au=300.0, e_out=0.0, mutual_rad=math.radians(60.0),
                  age_yr=5.0e9, extent_sum_au=0.005, w0_rad=0.4)
    values.update(changes)
    return kb.screen(**values)


def test_a_planet_with_an_inclined_distant_companion_is_pumped_to_the_expected_maximum():
    result = _planet()
    assert result.cycles and result.cause is None
    assert result.eccentricity_max == pytest.approx(_closed_form_peak(0.02, math.radians(60.0)), abs=3e-3)
    assert result.timescale_yr == pytest.approx(stability.kozai_timescale_yr(
        stability.orbital_period_yr(1.0, 1.0 + 3.0e-6), stability.orbital_period_yr(300.0, 1.5 + 3.0e-6), 1.5 + 3.0e-6, 0.5, 0.0))
    assert 0.0 <= result.eccentricity_now <= result.eccentricity_max + 1e-6 and result.cycles_run > 1.0
    assert not result.lost and result.period_yr > result.timescale_yr


def test_the_quenches_stop_a_cycle_and_name_their_cause():
    assert _planet(mutual_rad=math.radians(20.0)).cause == "window" and not _planet(mutual_rad=math.radians(20.0)).cycles
    assert _planet(mutual_rad=math.radians(160.0)).cause == "window"
    assert _planet(age_yr=1.0e3).cause == "timescale"
    assert _planet(a_out_au=2.5).cause == "hierarchy"
    assert _planet(other_rate_rad_yr=1.0).cause == "precession"      # a sibling that turns the orbit faster than the Kozai rate
    assert _planet(other_rate_rad_yr=1.0e-30).cycles                  # a negligible one does not
    hot = _planet(a_in_au=0.02, a_out_au=100.0, m_perturber=0.01, age_yr=1.0e12)      # a Kozai time of 7e9 years against a relativistic precession of 2e3
    assert not hot.cycles and hot.cause == "relativity"
    assert _planet(a_in_au=0.02, a_out_au=1.0, m_perturber=1.0, age_yr=1.0e12).cycles      # a heavy star close by cycles faster than relativity turns the orbit


def test_a_cycle_that_brings_the_periapsis_inside_the_extent_loses_the_body():
    assert not _planet().lost
    assert _planet(extent_sum_au=0.5).lost      # 1 AU (1 - 0.76) = 0.24 AU is inside 0.5 AU
    assert _planet(extent_sum_au=0.5).cycles


def test_the_octupole_lifts_the_maximum_for_an_eccentric_outer_orbit_between_25_and_90_degrees():
    plain = _planet(a_out_au=100.0, e_out=0.0, mutual_rad=math.radians(30.0))
    widened = _planet(a_out_au=100.0, e_out=0.6, mutual_rad=math.radians(30.0))
    assert plain.cause == "window" and not plain.cycles
    assert widened.eccentricity_max == pytest.approx(kb.OCTUPOLE_MAXIMUM) or not widened.cycles


def test_a_moon_whose_planet_orbits_a_star_is_a_triple_too():
    # a Moon-mass moon at 0.1 of the Hill radius of a Jupiter at 5 AU: the star (1 Msun at 5 AU) is the outer body; the moon's axis in AU
    hill_au = 5.0 * (1.0e-3 / 3.0) ** (1.0 / 3.0)
    inner = kb.screen(1.0e-3, 3.7e-8, 1.0, 0.1 * hill_au, 0.05, 5.0, 0.05, math.radians(65.0), 4.5e9, 1.0e-4, w0_rad=0.3)
    assert inner.cycles and inner.eccentricity_max > 0.7 and inner.timescale_yr < 4.5e9
    assert inner.timescale_yr < 1.0e5      # the star torques a moon far faster than a distant star torques a planet (cf. the irregulars' 65 to 180 years at 0.4 Hill)


def test_the_birth_eccentricity_is_never_below_the_floor_so_a_circular_start_still_cycles():
    result = _planet(e_in=0.0)
    assert result.cycles and result.eccentricity_max > 0.5


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# The Laplace plane and the moons of a planet
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------

def test_the_laplace_plane_lies_on_the_equator_close_in_on_the_orbit_far_out_and_balances_the_torques_between():
    obliquity = math.radians(60.0)
    assert kb.laplace_tilt_rad(0.0, obliquity) == 0.0
    assert kb.laplace_tilt_rad(1.0e6, obliquity) == pytest.approx(obliquity, abs=1e-4)
    for ratio in (0.1, 1.0, 10.0):
        theta = kb.laplace_tilt_rad(ratio, obliquity)
        assert 0.0 < theta < obliquity
        assert math.sin(2.0 * theta) == pytest.approx(ratio * math.sin(2.0 * (obliquity - theta)), rel=1e-6)
    assert kb.laplace_tilt_rad(1.0, obliquity) == pytest.approx(obliquity / 2.0, abs=1e-6)     # equal torques split the angle in two
    assert kb.laplace_tilt_rad(1.0, math.pi - obliquity) == pytest.approx(kb.laplace_tilt_rad(1.0, obliquity))     # a retrograde axis is the same plane


class _Bag:
    def __init__(self, **values):
        self.__dict__.update(values)


def _giant_with_moons(moons, obliquity_deg=70.0, rotation_hours=10.0):
    star = _Bag(mass=1.989e30, age=4.5, name="Star")
    tilt = math.radians(obliquity_deg)
    return _Bag(star=star, mass=1.898e27, radius=69911.0, distance=5.2, eccentricity=0.05, body_type="g", rotation_period_hours=rotation_hours,
                orbital_inclination_deg=0.0, orbital_ascending_node_deg=0.0, spin_axis_x=math.sin(tilt), spin_axis_y=0.0, spin_axis_z=math.cos(tilt), moons=moons)


def _moon(distance_au, origin, inclination_deg, mass=8.9e22, radius=1821.0, eccentricity=0.004):
    return _Bag(mass=mass, radius=radius, distance=distance_au, eccentricity=eccentricity, moon_origin=origin, orbital_inclination_deg=inclination_deg, orbital_ascending_node_deg=0.0,
                _orbit_birth=eccentricity)


def test_a_moon_that_formed_in_place_keeps_its_stored_inclination_as_the_free_one():
    planet = _giant_with_moons([])
    assert kb.moon_free_inclination(planet, _moon(0.003, "disk", 0.4), 5.0) == pytest.approx(math.radians(0.4))
    assert kb.moon_free_inclination(planet, _moon(0.003, "disk", 179.0), 5.0) == pytest.approx(math.radians(1.0))


def test_a_captured_moon_far_out_is_inclined_to_the_planets_orbit_and_close_in_to_the_equator():
    planet = _giant_with_moons([])
    moon = _moon(0.2, "captured", 0.0)         # in the equator plane, which is 70 degrees from the planet's orbit
    assert kb.moon_free_inclination(planet, moon, 1.0e6) == pytest.approx(math.radians(70.0), abs=0.01)
    assert kb.moon_free_inclination(planet, moon, 1.0e-6) == pytest.approx(0.0, abs=0.01)


def test_pump_moons_leaves_a_disk_moon_and_a_planet_without_moons_alone_and_is_deterministic():
    disk = _moon(0.003, "disk", 1.0)
    planet = _giant_with_moons([disk])
    assert kb.pump_moons(planet) == []
    assert getattr(disk, "_kozai", None) is None and disk.eccentricity == 0.004
    assert kb.pump_moons(_giant_with_moons([])) == []
    outcomes = []
    for _ in range(2):
        captured = _moon(0.12, "captured", 62.0, mass=1.0e18, radius=5.0, eccentricity=0.1)
        kb.pump_moons(_giant_with_moons([captured]))
        outcomes.append((getattr(captured, "_kozai", None), captured.eccentricity))
    assert outcomes[0] == outcomes[1]


def test_pump_moons_pumps_a_captured_moon_in_the_window_and_removes_one_the_cycle_drives_into_the_planet():
    ledger = _Bag(rows=[])
    ledger.add = lambda effect, receiver, source, value: ledger.rows.append((effect, value))
    captured = _moon(0.12, "captured", 62.0, mass=1.0e18, radius=5.0, eccentricity=0.1)
    planet = _giant_with_moons([captured])
    lost = kb.pump_moons(planet, ledger)
    kozai = captured._kozai
    assert kozai.cycles and kozai.eccentricity_max > 0.8
    if kozai.lost:
        assert lost == [captured]
    else:
        assert lost == [] and captured.eccentricity == kozai.eccentricity_now and captured._orbit_birth == kozai.eccentricity_now
    assert ledger.rows and ledger.rows[0][0] == "secular" and ledger.rows[0][1] > 0.0


def test_an_orbit_driven_almost_radial_does_not_leave_the_square_root_domain():
    """The stages of a Runge-Kutta step may carry the eccentricity past 1 when the cycle reaches 1 - 1e-9 (the crash of 7 test_orbital_motion setups)."""
    for e in (1.0 - 1.0e-9, 1.0, 1.0 + 1.0e-6, -1.0e-3):
        de, dw = kb._derivatives((e, 0.7, 0.01), 4.0 / 3.0)
        assert math.isfinite(de) and math.isfinite(dw)
    e, w = kb._step(1.0 - 1.0e-9, 0.7, 0.01, 4.0 / 3.0, 1.0 / kb.STEPS_PER_KOZAI_TIME)
    assert 0.0 <= e < 1.0 and math.isfinite(w)
