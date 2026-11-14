# tests/test_starless_day.py

"""GEN.1496: a body no star lights (a rogue planet and the moon of one) stores no day: the set of suns is empty, so it is dark all year."""

from types import SimpleNamespace

import pytest

from planetgen.physics import day_night, daylight, habitability_world, latitude_bands


def _body(starless, body_type="t"):
    star = SimpleNamespace(radius=695_700.0, luminosity=0.0 if starless else 3.828e26, starless=starless)
    return SimpleNamespace(rotation_period_hours=20.0, period=1.0e6, axial_tilt_deg=23.0, distance=5.0, eccentricity=0.0, star=star, body_type=body_type,
                           envelope_mass_fraction=0.0, spin_axis_x=0.0, spin_axis_y=0.0, spin_axis_z=1.0, orbital_inclination_deg=0.0, orbital_ascending_node_deg=0.0)


def test_a_planet_under_no_star_is_dark_all_year():
    fields = daylight.compute(_body(True))
    assert fields["day_share_equator"] == 0.0 and fields["day_share_mean"] == 0.0 and fields["daylight_hours"] == 0.0
    assert fields["day_night_state"] == day_night.NO_SUN
    assert fields["solar_day_hours"] is None and fields["polar_circle_deg"] is None and fields["night_hours"] is None and fields["days_per_year"] is None


def test_a_planet_under_a_star_still_has_its_day():
    fields = daylight.compute(_body(False))
    assert fields["day_share_mean"] == pytest.approx(0.5, abs=0.01)
    assert fields["day_night_state"] != day_night.NO_SUN and fields["solar_day_hours"] is not None


def test_a_giant_under_no_star_still_has_no_surface():
    giant = _body(True, "g")
    giant.envelope_mass_fraction = 0.9
    fields = daylight.compute(giant)
    assert fields["day_night_state"] == day_night.NO_SURFACE and fields["day_share_mean"] == 0.0


def test_no_star_means_no_latitude_bands():
    values = {"surface_temperature_k": 60.0, "rotation_period_hours": 20.0, "period_years": 1.0e6, "axial_tilt_deg": 23.0, "atmospheric_pressure_pa": 0.0, "gravity_g": 1.0}
    assert len(latitude_bands.compute(values)) == latitude_bands.BAND_COUNT
    assert latitude_bands.compute({**values, "starless": True}) == []
    assert latitude_bands.compute({**values, "day_share_mean": 0.0}) == []


def test_the_values_of_a_body_say_whether_its_star_is_starless():
    assert habitability_world.values_from_body(_body(True))["starless"] is True
    assert habitability_world.values_from_body(_body(False))["starless"] is False


def test_a_built_rogue_planet_and_its_moon_store_no_day_and_no_bands():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.phenomena import starless_host
    from planetgen.generation.system import StarSystem
    from planetgen.physics import unit_registry
    from planetgen.util import draw

    draw.set_run_seed(1496)
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    donor = StarSystem(config).star
    rogue = starless_host.build_ejected_planet(config, donor, 2.0, unit_registry.converter("Mearth", "kg")(3.0), formation_au=5.0, moon_count=1)
    for body in [rogue, *getattr(rogue, "moons", [])]:
        assert body.day_share_mean == 0.0 and body.daylight_hours == 0.0
        assert body.day_night_state in (day_night.NO_SUN, day_night.NO_SURFACE)
        assert body.band1_mean_k is None and body.band4_state is None
