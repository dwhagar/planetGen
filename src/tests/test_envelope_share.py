"""GEN.1226: the envelope share of every body is the true share of a hydrogen-rich air; GEN.1219 types a body from it; GEN.1220 stores the class digit."""

import copy

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import atmosphere, bulk_composition, constants, envelope_share, gases, planet_type
from tests import class_support


class Body:
    """A rocky body with a mass, radius, pressure and air, and the solid fractions of a drawn body."""

    def __init__(self, mass_earth=1.0, radius_km=6371.0, pressure_pa=101325.0, air=None):
        self.body_type = "t"
        self.mass = mass_earth * constants.EARTH_MASS_TO_KG
        self.radius = radius_km
        self.atmospheric_pressure = pressure_pa
        setattr(self, atmosphere.GAS_FIELD, air if air is not None else {"N2": 0.78 * pressure_pa, "O2": 0.21 * pressure_pa})
        self.metal_mass_fraction, self.rock_mass_fraction, self.water_mass_fraction, self.ice_mass_fraction = 0.3, 0.65, 0.05, 0.0
        self.envelope_mass_fraction = 0.0


def _hydrogen_air(pressure_pa):
    return {gases.stored("h2").formula: 0.9 * pressure_pa, gases.stored("ch4").formula: 0.1 * pressure_pa}


def test_a_nitrogen_air_has_no_envelope():
    body = Body()
    envelope_share.fill(body)
    assert body.envelope_mass_fraction == 0.0 and body.rock_mass_fraction == pytest.approx(0.65)


def test_a_hydrogen_air_is_an_envelope_of_the_air_mass_over_the_body_mass():
    body = Body(mass_earth=4.0, radius_km=11000.0, pressure_pa=5.0e9, air=None)
    setattr(body, atmosphere.GAS_FIELD, _hydrogen_air(5.0e9))
    share = envelope_share.air_share(body)
    assert 0.001 < share < 0.5
    envelope_share.fill(body)
    assert body.envelope_mass_fraction == pytest.approx(share)
    assert planet_type.planet_type(body) == ("transitional" if share >= planet_type.TERRESTRIAL_MAX_ENVELOPE else "terrestrial")


def test_the_five_fractions_still_add_up_to_one():
    body = Body(mass_earth=4.0, radius_km=11000.0, pressure_pa=5.0e9)
    setattr(body, atmosphere.GAS_FIELD, _hydrogen_air(5.0e9))
    envelope_share.fill(body)
    total = (body.metal_mass_fraction + body.rock_mass_fraction + body.water_mass_fraction + body.ice_mass_fraction + body.envelope_mass_fraction)
    assert total == pytest.approx(1.0, abs=bulk_composition.SUM_TOLERANCE)


def test_filling_twice_does_not_compound_and_restore_puts_the_drawn_fractions_back():
    body = Body(mass_earth=4.0, radius_km=11000.0, pressure_pa=5.0e9)
    setattr(body, atmosphere.GAS_FIELD, _hydrogen_air(5.0e9))
    envelope_share.fill(body)
    first = (body.envelope_mass_fraction, body.rock_mass_fraction)
    envelope_share.fill(body)
    assert (body.envelope_mass_fraction, body.rock_mass_fraction) == first
    envelope_share.restore(body)
    assert body.envelope_mass_fraction == 0.0 and body.rock_mass_fraction == pytest.approx(0.65)


def test_the_fractions_add_to_one_when_the_body_was_drawn_with_an_envelope_and_has_a_hydrogen_air_too():
    """The air's share comes out of the solids, which already hold only 1 less the envelope drawn: a body of 1.1 Earth masses with 0.1 percent drawn and a share of 7e-6
    summed to 1 + 8.6e-9, over the 1e-9 the bulk fractions are held to."""
    body = Body(mass_earth=1.1, radius_km=6600.0, pressure_pa=5.0e4, air=None)
    setattr(body, atmosphere.GAS_FIELD, _hydrogen_air(5.0e4))
    body.metal_mass_fraction, body.rock_mass_fraction, body.water_mass_fraction, body.ice_mass_fraction = 0.27 * 0.9988, 0.656 * 0.9988, 0.0617 * 0.9988, 0.011 * 0.9988
    body.envelope_mass_fraction = 1.0 - (body.metal_mass_fraction + body.rock_mass_fraction + body.water_mass_fraction + body.ice_mass_fraction)
    envelope_share.fill(body)
    assert envelope_share.air_share(body) > 0.0
    total = sum(getattr(body, name) for name in (*envelope_share.SOLIDS, "envelope_mass_fraction"))
    assert total == pytest.approx(1.0, abs=1e-12)


def test_an_air_with_little_hydrogen_is_secondary_and_has_no_envelope():
    body = Body(pressure_pa=1.0e7)
    setattr(body, atmosphere.GAS_FIELD, {gases.stored("h2").formula: 0.1 * 1.0e7, gases.stored("co2").formula: 0.9 * 1.0e7})
    assert envelope_share.air_share(body) == 0.0


def test_a_giant_keeps_the_envelope_its_interior_holds():
    body = Body()
    body.body_type, body.envelope_mass_fraction = "g", 0.9
    envelope_share.fill(body)
    assert body.envelope_mass_fraction == 0.9


def test_generated_planets_store_their_type_digit_and_whole_fractions():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    star = copy.copy(Star(config))
    star.age = 4.5
    star.luminosity = constants.SOLAR_LUMINOSITY
    kinds = set()
    for seed in range(12):
        planet = class_support.draw_planet(star, 0.6 + 0.4 * seed, seed, config, moon_count=1)
        for body in (planet, *planet.moons):
            assert body.planet_type in planet_type.TYPES
            assert body.class_variant == 0 or body.planet_class is None
            fractions = bulk_composition.bulk_fractions(body)
            if fractions is not None:
                assert sum(fractions.shares().values()) == pytest.approx(1.0, abs=bulk_composition.SUM_TOLERANCE)
            kinds.add(body.planet_type)
    assert "terrestrial" in kinds
