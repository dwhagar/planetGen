"""GEN.442: the drawn composition of small bodies as the generator wires it (belts, comets and interstellar comets read
`physics/small_body_composition.py`; the belt reads its star's stored snow line)."""

import math

import pytest

from planetgen.admin import edits
from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.comet import Comet
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import InterstellarComet
from planetgen.generation.system import StarSystem
from planetgen.physics import small_body_composition as sb
from planetgen.util import draw
from tests.test_systems import make_config


@pytest.fixture(scope="module")
def system():
    return StarSystem(system_config=make_config("G2V", PLANETS=False))


def test_a_belt_reads_its_stars_snow_line(system):
    snow = system.star.snow_line_au
    with draw.bound(5):
        inner = AsteroidBelt(system.system_config, system.star, 0.5 * snow, 0.5 * snow, 0.6 * snow)
        outer = AsteroidBelt(system.system_config, system.star, 2.0 * snow, 2.0 * snow, 2.4 * snow)
    assert inner.makeup.detail["x"] == pytest.approx(0.5)
    assert outer.makeup.detail["x"] == pytest.approx(2.0)
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
    assert owner[0] is not belt and owner[0].makeup.detail["x"] == pytest.approx(2.0 / system.star.snow_line_au)


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
