"""GEN.466: the origin and composition of a moon (`planetgen/physics/moon_composition.py`; docs/design/small-bodies-composition.md
sections 5.4 and 5.5). Io, Europa, Ganymede and Callisto appear only as test targets for the order of the ice ramp."""

import math

import pytest

from planetgen.generation import possibility_table as pt
from planetgen.physics import compaction, moon_composition as mc
from planetgen.util import draw


def _names(kind, **kwargs):
    subject = pt.Subject(kind, {"inclination_deg": kwargs.get("inclination_deg", 1.0), "has_ring": int(kwargs.get("has_ring", False)),
                                "roche_ratio": kwargs.get("roche_ratio", math.inf)})
    return sorted(name for name, _p in [(row.name, p) for row, p in mc.ORIGIN_TABLE.probabilities(subject)])


def test_the_origin_table_is_registered_and_every_weight_is_open():
    assert pt.REGISTRY.get("moon-origin") is mc.ORIGIN_TABLE
    assert mc.ORIGIN_TABLE.state_counts()[pt.OPEN_STATE] == 5
    assert [row.name for row in mc.ORIGIN_TABLE.rows] == ["disk", "impact", "captured", "ring", "binary"]


def test_every_parent_kind_and_orbit_has_an_origin_and_each_origin_is_reachable():
    assert _names(mc.GIANT) == ["disk"]
    assert _names(mc.GIANT, has_ring=True, roche_ratio=2.0) == ["disk", "ring"]
    assert _names(mc.GIANT, has_ring=True, roche_ratio=20.0) == ["disk"]
    assert _names(mc.TERRESTRIAL) == ["impact"]
    assert _names(mc.DWARF_PLANET) == ["binary"]
    for kind in mc.PARENT_KINDS:
        assert _names(kind, inclination_deg=120.0) == ["captured"]


def test_open_rows_share_the_probability_equally():
    subject = pt.Subject(mc.GIANT, {"inclination_deg": 1.0, "has_ring": 1, "roche_ratio": 2.0})
    assert [p for _row, p in mc.ORIGIN_TABLE.probabilities(subject)] == pytest.approx([0.5, 0.5])
    with draw.bound(1):
        seen = {mc.draw_origin(mc.GIANT, 1.0, has_ring=True, roche_ratio=2.0).row.name for _ in range(60)}
    assert seen == {"disk", "ring"}


def test_a_bad_parent_kind_is_refused():
    with pytest.raises(ValueError):
        mc.draw_origin("moon", 1.0)


def test_the_parent_kind_from_envelope_and_mass():
    assert mc.parent_kind(0.8, 1.9e27) == mc.GIANT
    assert mc.parent_kind(0.0, 6e24) == mc.TERRESTRIAL
    assert mc.parent_kind(0.0, 1.3e22) == mc.DWARF_PLANET


def test_the_ice_line_is_the_sasaki_formula():
    assert mc.ice_line_km(mc.JUPITER_MASS_KG, 5.0) == pytest.approx(20.0 * mc.JUPITER_RADIUS_KM)
    assert mc.ice_line_km(mc.JUPITER_MASS_KG, 5.0) == pytest.approx(1.43e6, rel=0.01)
    assert mc.ice_line_km(0.3 * mc.JUPITER_MASS_KG, 5.0) == pytest.approx(6.4e5, rel=0.02)
    assert mc.ice_line_km(3.0 * mc.JUPITER_MASS_KG, 5.0) == pytest.approx(3.0e6, rel=0.02)
    assert mc.ice_line_km(mc.JUPITER_MASS_KG, 10.0) < mc.ice_line_km(mc.JUPITER_MASS_KG, 1.0)
    with pytest.raises(ValueError):
        mc.ice_line_km(0.0, 5.0)


def test_the_total_moon_mass_is_the_canup_and_ward_range():
    assert mc.total_moon_mass_ratio(3e-5) == pytest.approx(2.5e-4)
    assert mc.total_moon_mass_ratio(1e-6) == pytest.approx(0.8e-4, rel=0.02)
    assert mc.total_moon_mass_ratio(5e-4) == pytest.approx(6.4e-4, rel=0.02)


def test_a_disk_system_partitions_the_total_and_flags_an_extrapolated_ice_line():
    with draw.bound(2):
        systems = [mc.draw_disk_system(mc.JUPITER_MASS_KG) for _ in range(200)]
    for system in systems:
        assert 1 <= len(system.moon_masses_kg) <= 5
        total = sum(system.moon_masses_kg) / mc.JUPITER_MASS_KG
        assert mc.total_moon_mass_ratio(mc.ALPHA_OVER_F[0]) * 0.999 <= total <= mc.total_moon_mass_ratio(mc.ALPHA_OVER_F[1]) * 1.001
        assert 1.0 <= system.lifetime_myr <= 10.0 and 0.1 <= system.f_inf <= 0.6
        assert not system.extrapolated
    with draw.bound(3):
        assert mc.draw_disk_system(0.02 * mc.JUPITER_MASS_KG).extrapolated
    again = mc.DiskSystem.from_dict(systems[0].to_dict())
    assert again.ice_line_km == pytest.approx(systems[0].ice_line_km) and again.moon_masses_kg == systems[0].moon_masses_kg


def test_the_ice_fraction_rises_with_formation_distance_in_the_physical_direction():
    system = mc.DiskSystem(mc.JUPITER_MASS_KG, 5.0, 3e-5, 0.35, [1e22])
    line = system.ice_line_km
    ice = []
    with draw.bound(4):
        for ratio in (0.295, 0.47, 0.75, 1.3):  # the order of the four large moons of the Solar System's biggest planet
            ice.append(mc.draw_disk_moon(ratio * line, system).water)
    assert ice == sorted(ice) and ice[0] == 0.0 and ice[-1] == pytest.approx(0.35)
    with draw.bound(5):
        inner = mc.draw_disk_moon(0.2 * line, system)
        outer = mc.draw_disk_moon(2.0 * line, system)
    assert inner.detail["hydration"] == 0.0 and outer.detail["hydration"] == 1.0
    assert outer.components.get("serpentine", 0.0) > inner.components.get("serpentine", 0.0)


def test_every_origin_makes_fractions_that_add_up_to_one():
    system = mc.DiskSystem(mc.JUPITER_MASS_KG, 5.0, 3e-5, 0.35, [1e22])
    with draw.bound(6):
        moons = [mc.draw_disk_moon(f * system.ice_line_km, system) for f in (0.1, 0.5, 1.0, 3.0)]
        moons += [mc.draw_moon(name) for name in ("impact", "captured", "ring", "binary") for _ in range(25)]
    for moon in moons:
        assert moon.metal + moon.rock + moon.water + moon.ice == pytest.approx(1.0)
        assert math.fsum(moon.components.values()) == pytest.approx(1.0)
        assert moon.bulk_density_kg_m3 is None and 0.0 <= moon.porosity < 1.0


def test_an_impact_moon_has_little_metal_and_no_ice_and_a_captured_one_loses_ice():
    with draw.bound(7):
        impact = [mc.draw_impact_moon() for _ in range(60)]
        kuiper = [mc.draw_captured_moon() for _ in range(60)]
    assert all(0.01 <= m.metal <= 0.0501 and m.water <= 1.1e-3 and m.ice == 0.0 for m in impact)
    assert all(0.0 <= m.detail["ice_loss"] <= 0.5 for m in kuiper)
    assert all(m.water <= 1.0 - m.detail["rock_fraction"] + 1e-9 for m in kuiper)


def test_a_ring_moon_is_ice_rich_and_a_binary_partner_has_a_mass_ratio():
    with draw.bound(8):
        rings = [mc.draw_ring_moon() for _ in range(40)]
        partners = [mc.draw_binary_partner() for _ in range(40)]
    assert all(0.5 <= m.water <= 0.95 + 1e-9 for m in rings)
    assert all(5e-4 <= m.detail["mass_ratio"] <= 0.11 for m in partners)


def test_a_disk_moon_needs_its_distance_and_system_and_an_unknown_origin_is_refused():
    with pytest.raises(ValueError):
        mc.draw_moon("disk")
    with pytest.raises(ValueError):
        mc.draw_moon("volcanic")
    with pytest.raises(ValueError):
        mc.draw_disk_moon(0.0, mc.DiskSystem(mc.JUPITER_MASS_KG, 5.0, 3e-5, 0.35, [1e22]))


def test_the_density_comes_from_the_compaction_law_and_stops_at_its_edge():
    system = mc.DiskSystem(mc.JUPITER_MASS_KG, 5.0, 3e-5, 0.35, [1e22])
    with draw.bound(9):
        icy = mc.draw_disk_moon(1.3 * system.ice_line_km, system)
        dry = mc.draw_disk_moon(0.2 * system.ice_line_km, system)
    mass = 5e20
    assert mc.density_kg_m3(icy, mass) < mc.density_kg_m3(dry, mass) or icy.rock + icy.metal >= 0.95
    assert 800.0 < mc.density_kg_m3(icy, mass) < 3500.0
    with pytest.raises(compaction.OutOfRange):
        mc.density_kg_m3(icy, 1e23)


def test_the_same_seed_gives_the_same_moon():
    system = mc.DiskSystem(mc.JUPITER_MASS_KG, 5.0, 3e-5, 0.35, [1e22])
    with draw.bound(10):
        a = mc.draw_disk_moon(system.ice_line_km, system).to_dict()
    with draw.bound(10):
        b = mc.draw_disk_moon(system.ice_line_km, system).to_dict()
    assert a == b
