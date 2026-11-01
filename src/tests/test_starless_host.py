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
