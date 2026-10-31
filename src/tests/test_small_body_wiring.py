"""GEN.442: the drawn composition of small bodies as the generator wires it (belts, comets and interstellar comets read
`physics/small_body_composition.py`; the belt reads its star's luminosity and aluminium-26)."""

import math

import pytest

from planetgen.admin import edits
from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.comet import Comet
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import InterstellarComet
from planetgen.generation.system import StarSystem
from planetgen.physics import condensate_shares, radiogenic, small_body_composition as sb
from planetgen.util import draw
from tests.test_systems import make_config


@pytest.fixture(scope="module")
def system():
    return StarSystem(system_config=make_config("G2V", PLANETS=False))


def test_a_belt_reads_its_stars_luminosity_and_aluminium_26(system):
    line = condensate_shares.ice_line_au(system.star.luminosity)
    with draw.bound(5):
        inner = AsteroidBelt(system.system_config, system.star, 0.5 * line, 0.5 * line, 0.6 * line)
        outer = AsteroidBelt(system.system_config, system.star, 2.0 * line, 2.0 * line, 2.4 * line)
    assert inner.makeup.detail["temperature_k"] == pytest.approx(condensate_shares.ICE_CONDENSATION_K["H2O"] * 2.0 ** 0.5, rel=1e-3)
    assert outer.makeup.detail["temperature_k"] == pytest.approx(condensate_shares.ICE_CONDENSATION_K["H2O"] * 0.5 ** 0.5, rel=1e-3)
    assert inner.makeup.detail["al26_ratio"] == outer.makeup.detail["al26_ratio"] == radiogenic.star_al26_ratio(system.star, system.system_config.GALAXY_SHAPE)
    assert inner.makeup.rock + inner.makeup.metal > outer.makeup.rock + outer.makeup.metal - 1.0  # both are valid mixtures
    assert all(level in ("high", "moderate", "small", "trace") for _label, level in inner.composition)
    assert 1 <= len(inner.composition) <= sb.TEXT_LINES


def test_a_belt_round_trips_with_its_makeup(system):
    with draw.bound(6):
        belt = AsteroidBelt(system.system_config, system.star, 2.0, 1.8, 2.2)
    again = AsteroidBelt.from_dict(belt.to_dict(), system.system_config)
    assert again.composition == belt.composition
    assert again.makeup.bulk_density_kg_m3 == pytest.approx(belt.makeup.bulk_density_kg_m3)
    assert again.makeup.components == pytest.approx(belt.makeup.components)
    assert again.get_composition_summary() == belt.get_composition_summary()


def test_a_belts_mass_and_density_word_come_from_its_discs_solids_ground_down_by_its_age(system):
    """GEN.477: the mass is `belt_mass.belt_state` at the star's mass, snow line, metals and age; it round-trips and the word is no longer a free draw."""
    from planetgen.physics import belt_mass, bulk_draw

    star = system.star
    with draw.bound(9):
        belt = AsteroidBelt(system.system_config, star, 2.5, 2.2, 3.4)
    largest = belt_mass.draw_largest_km(draw.Stream(f"belt-largest:{2.5!r}:{2.2!r}:{3.4!r}:{star.mass!r}"))
    mass_kg, word = belt_mass.belt_state(2.2, 3.4, star.mass, star.snow_line_au, bulk_draw.metal_ratio_to_solar(star), star.age, largest)
    assert belt.mass_kg == pytest.approx(mass_kg) and belt.density == word and belt.mass_kg > 0.0
    again = AsteroidBelt.from_dict(belt.to_dict(), system.system_config)
    assert again.mass_kg == belt.mass_kg and again.density == belt.density
    with draw.bound(10):
        twin = AsteroidBelt(system.system_config, star, 2.5, 2.2, 3.4)
    assert twin.mass_kg == belt.mass_kg, "the belt's own stream, not the run's, draws its largest body"


def test_a_belt_from_a_store_without_the_makeup_still_loads(system):
    with draw.bound(7):
        belt = AsteroidBelt(system.system_config, system.star, 2.0, 1.8, 2.2)
    data = belt.to_dict()
    data.pop("makeup")
    assert AsteroidBelt.from_dict(data, system.system_config).makeup is None


def test_regenerating_a_belt_rolls_a_new_makeup_around_the_same_star(system):
    with draw.bound(8):
        belt = AsteroidBelt(system.system_config, system.star, 2.0, 1.8, 2.2)
    belt.db_id = 1
    owner = [belt]
    original = system.planets
    system.planets = owner
    try:
        edits.regenerate_belt(system, belt, owner)
    finally:
        system.planets = original
    assert owner[0] is not belt and owner[0].makeup.detail["temperature_k"] == pytest.approx(
        condensate_shares.disk_temperature_k(2.0, system.star.luminosity))


def test_a_comet_mass_is_its_sphere_at_its_makeups_density():
    with draw.bound(9):
        comet = Comet(SystemConfig(), primary_mass_solar=1.0, orbit_type="elliptical")
    sphere = 4.0 / 3.0 * math.pi * (comet.nucleus_diameter_km * 500.0) ** 3
    assert comet.mass_kg == pytest.approx(sphere * comet.makeup.bulk_density_kg_m3)
    assert 100.0 < comet.makeup.bulk_density_kg_m3 < 1500.0
    again = Comet.from_dict(comet.to_dict(), SystemConfig())
    assert again.mass_kg == pytest.approx(comet.mass_kg) and again.composition == comet.composition
    assert again.makeup.porosity == pytest.approx(comet.makeup.porosity)


def test_an_interstellar_comet_has_a_makeup_and_a_mass():
    with draw.bound(10):
        comet = InterstellarComet(SystemConfig())
    assert comet.makeup.kind == sb.COMET and comet.mass_kg > 0.0
    assert InterstellarComet.from_dict(comet.to_dict(), SystemConfig()).mass_kg == pytest.approx(comet.mass_kg)


def test_a_standalone_field_has_a_family_and_a_makeup_from_its_taxa():
    from planetgen.generation.phenomena.asteroid_field import AsteroidField
    with draw.bound(11):
        fields = [AsteroidField(SystemConfig()) for _ in range(40)]
    assert {f.composition_family for f in fields} <= set(sb.FIELD_FAMILIES)
    for field in fields:
        shares = field.makeup.detail["taxon_shares"]
        assert sum(shares.values()) == pytest.approx(1.0)
        allowed = set(sb.FIELD_FAMILIES[field.composition_family])
        if allowed:
            assert {t for t, v in shares.items() if v > 0.0} <= allowed
        assert field.composition and all(level in ("high", "moderate", "small", "trace") for _label, level in field.composition)
        again = AsteroidField.from_dict(field.to_dict(), SystemConfig())
        assert again.composition == field.composition and again.makeup.bulk_density_kg_m3 == pytest.approx(field.makeup.bulk_density_kg_m3)


def test_the_field_family_table_is_gated_and_every_row_is_a_stand_in():
    table = sb.FIELD_FAMILY_TABLE
    assert [row.name for row in table.rows] == list(sb.FIELD_FAMILIES)
    assert all(row.state == "stand-in" and row.replacement for row in table.rows)
    assert all(sb.field_family_text(name) for name in sb.FIELD_FAMILIES)


def test_dust_has_no_pores_and_a_collisional_field_is_one_taxon():
    with draw.bound(12):
        dust = sb.draw_field("dust")
        broken = sb.draw_field("collisional")
    assert dust.porosity == 0.0
    assert sum(1 for v in broken.detail["taxon_shares"].values() if v > 0.0) == 1
    with pytest.raises(ValueError):
        sb.draw_field("volcanic")


def test_a_comet_of_a_star_takes_its_ice_from_the_zone_it_formed_in(system):
    """GEN.595: the median CO and CO2 follow the disc zone of the star's luminosity, snow line and edge, and the zone costs the run's stream nothing."""
    from planetgen import tuning
    from planetgen.physics import bulk_draw, comet_ice_budget

    star = system.star
    mass = star.mass / 1.989e30
    with draw.bound(21):
        comets = [Comet(system.system_config, mass, orbit_type="elliptical", star=star) for _ in range(60)]
    with draw.bound(21):
        plain = [Comet(system.system_config, mass, orbit_type="elliptical") for _ in range(60)]
    zone_co = sorted(c.makeup.detail["volatiles"]["CO"] for c in comets)[30]
    carbon, oxygen = bulk_draw.ice_forming_ratios(star)
    edge = tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER * star.snow_line_au
    zones = [comet_ice_budget.formation_ratios(star.luminosity, star.snow_line_au, edge, carbon, oxygen, (i + 0.5) / 200) for i in range(200)]
    assert min(z[1] for z in zones) * 100 / 3 < zone_co < max(z[1] for z in zones) * 100 * 3
    assert zone_co != sorted(c.makeup.detail["volatiles"]["CO"] for c in plain)[30]
    comet = comets[0]
    with draw.bound(22):
        first = draw.random()
    with draw.bound(22):
        assert comet._formation_ice_ratios(star) is not None
        assert draw.random() == first
    assert Comet(system.system_config, mass, orbit_type="elliptical", star=None)._formation_ice_ratios(None) is None


def test_a_star_with_no_abundances_keeps_the_population_ice():
    from types import SimpleNamespace
    from planetgen.physics import bulk_draw
    remnant = SimpleNamespace(metallicity_fe_h=None, age=1.0, mass=1e30)
    assert bulk_draw.ice_forming_ratios(remnant) is None
    assert Comet(SystemConfig(), 1.0, orbit_type="elliptical", star=remnant)._formation_ice_ratios(remnant) is None


def test_the_carbon_and_oxygen_of_a_star_are_ratios_to_the_suns_over_magnesium(system):
    from planetgen.physics import bulk_draw
    carbon, oxygen = bulk_draw.ice_forming_ratios(system.star)
    assert 0.1 < carbon < 10.0 and 0.1 < oxygen < 10.0
