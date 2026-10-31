"""GEN.466: the origin and composition of a moon (`planetgen/physics/moon_composition.py`; docs/design/small-bodies-composition.md
sections 5.4 and 5.5). Io, Europa, Ganymede and Callisto appear only as test targets for the order of the ice ramp."""

import math

import pytest

from planetgen.generation import possibility_table as pt
from planetgen.physics import compaction, moon_composition as mc
from planetgen.util import draw


def _names(kind, **kwargs):
    subject = pt.Subject(kind, {"has_ring": int(kwargs.get("has_ring", False)), "roche_ratio": kwargs.get("roche_ratio", math.inf)})
    return sorted(name for name, _p in [(row.name, p) for row, p in mc.ORIGIN_TABLE.probabilities(subject)])


def test_the_origin_table_is_registered_and_every_weight_is_open():
    assert pt.REGISTRY.get("moon-origin") is mc.ORIGIN_TABLE
    assert mc.ORIGIN_TABLE.state_counts()[pt.OPEN_STATE] == 5
    assert [row.name for row in mc.ORIGIN_TABLE.rows] == ["disk", "impact", "captured", "ring", "binary"]


def test_every_parent_kind_and_orbit_has_an_origin_and_each_origin_is_reachable():
    assert _names(mc.GIANT) == ["captured", "disk"]
    assert _names(mc.GIANT, has_ring=True, roche_ratio=2.0) == ["captured", "disk", "ring"]
    assert _names(mc.GIANT, has_ring=True, roche_ratio=20.0) == ["captured", "disk"]
    assert _names(mc.TERRESTRIAL) == ["captured", "impact"]
    assert _names(mc.DWARF_PLANET) == ["binary", "captured"]


def test_open_rows_share_the_probability_equally():
    subject = pt.Subject(mc.GIANT, {"has_ring": 1, "roche_ratio": 2.0})
    assert [p for _row, p in mc.ORIGIN_TABLE.probabilities(subject)] == pytest.approx([1 / 3] * 3)
    with draw.bound(1):
        seen = {mc.draw_origin(mc.GIANT, has_ring=True, roche_ratio=2.0).row.name for _ in range(120)}
    assert seen == {"disk", "ring", "captured"}


def test_a_bad_parent_kind_is_refused():
    with pytest.raises(ValueError):
        mc.draw_origin("moon")


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


def test_a_regular_moon_lies_near_the_equator_and_a_captured_one_can_be_retrograde_outside_the_kozai_band():
    """GEN.714: the orbit follows from the origin."""
    from planetgen.physics import constants
    low, high = mc.KOZAI_FORBIDDEN_DEG
    with draw.bound(3):
        regular = [mc.draw_inclination_deg(origin) for origin in ("disk", "impact", "ring", "binary") for _ in range(100)]
        captured = [mc.draw_inclination_deg("captured") for _ in range(2000)]
    assert all(0.0 <= i <= constants.MOON_ORBITAL_INCLINATION_MAX_DEG for i in regular)
    assert all(0.0 <= i <= 180.0 for i in captured)
    assert not any(low < i < high for i in captured)
    assert any(i > 150.0 for i in captured) and any(i < 20.0 for i in captured)
    # prograde and retrograde in the ratio of the cosines each side of the band: 0.426 against 0.357
    retrograde = sum(i > 90.0 for i in captured) / len(captured)
    assert retrograde == pytest.approx(0.357 / (0.426 + 0.357), abs=0.04)
    with pytest.raises(ValueError):
        mc.draw_inclination_deg("moon")


def test_the_inclination_draw_is_one_number_and_repeats_with_the_seed():
    with draw.bound(5):
        first = [mc.draw_inclination_deg("captured") for _ in range(5)]
    with draw.bound(5):
        assert [mc.draw_inclination_deg("captured") for _ in range(5)] == first


def test_generated_moons_draw_their_origin_first_and_the_orbit_follows_from_it():
    """GEN.714: through `generate_moons`, every moon holds an origin, regular moons lie near the equator and captured ones, which exist, can be retrograde."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.physics import constants
    low, high = mc.KOZAI_FORBIDDEN_DEG
    origins = {}
    with draw.bound(21):
        for _ in range(60):
            config = SystemConfig()
            config.MOONS = True
            system = StarSystem(system_config=config)
            for planet in system.planets + system.secondary_planets:
                if planet.body_type == "a":
                    continue
                for moon in planet.moons:
                    assert moon.moon_origin in mc.ORIGINS
                    origins.setdefault(moon.moon_origin, []).append(moon.orbital_inclination_deg)
    assert "captured" in origins and ("disk" in origins or "impact" in origins)
    for origin, values in origins.items():
        if origin == "captured":
            assert not any(low < i < high for i in values)
        else:
            assert all(i <= constants.MOON_ORBITAL_INCLINATION_MAX_DEG for i in values)
    assert any(i > 130.0 for i in origins["captured"])


def _disk_ice_sides(moons):
    """`(near, far)`: the disk moons formed under a third of their disk's ice line and past one and a half times it."""
    disk = [moon for moon in moons if moon.moon_origin == "disk"]
    # a moon's makeup is set where it formed (`_birth_distance_au`); the tide has since moved it (`moon_history`, GEN.959)
    ratio = lambda moon: moon._birth_distance_au * 1.495978707e8 / moon.moon_disk_system.ice_line_km
    return [m for m in disk if ratio(m) < 0.3], [m for m in disk if ratio(m) > 1.5]


@pytest.fixture(scope="module")
def generated_moons():
    """The moons of seeded Sun-like single stars (GEN.466): at least 40 systems, and more (up to 120) until a few disk moons stand on each
    side of their disk's ice line, since only a giant makes a disk moon and few systems have one."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    moons = []
    with draw.bound(8):
        for count in range(120):
            if count >= 40 and all(len(side) >= 2 for side in _disk_ice_sides(moons)):
                break
            config = SystemConfig()
            config.MOONS = True
            config.STAR_TYPE = "G2V"
            config.BINARY_SYSTEM = False
            system = StarSystem(system_config=config)
            for planet in system.planets:
                if planet.body_type != "a":
                    moons.extend(getattr(planet, "_born_moons", None) or planet.moons)      # the moons as generated: the tide has since removed some (GEN.959)
    return moons


def test_every_generated_moon_has_a_composition_from_its_origin(generated_moons):
    origins = {moon.moon_origin for moon in generated_moons}
    assert {"disk", "impact", "captured"} <= origins
    for moon in generated_moons:
        shares = [moon.metal_mass_fraction, moon.rock_mass_fraction, moon.water_mass_fraction, moon.ice_mass_fraction,
                  moon.envelope_mass_fraction]
        assert sum(shares) == pytest.approx(1.0, abs=1e-9)
        assert moon.envelope_mass_fraction == 0.0
    impact = [moon for moon in generated_moons if moon.moon_origin == "impact"]
    assert all(moon.water_mass_fraction < 1e-3 and moon.ice_mass_fraction == 0.0 for moon in impact)
    assert all(0.0 < moon.metal_mass_fraction <= 0.06 for moon in impact)


def test_a_disk_moon_gains_ice_past_its_disks_ice_line(generated_moons):
    near, far = ([moon.water_mass_fraction for moon in side] for side in _disk_ice_sides(generated_moons))
    near, far = list(near), list(far)
    assert near and far
    assert max(near) < 0.05
    assert sum(far) / len(far) > sum(near) / len(near) + 0.1


def test_a_small_moons_density_follows_the_tapered_compaction_law(generated_moons):
    from planetgen.physics import compaction
    small = [moon for moon in generated_moons if moon.mass < compaction.TAPER_END_KG]
    assert small
    checked = 0
    for moon in small:
        rock = min(max(moon.metal_mass_fraction + moon.rock_mass_fraction, compaction.ROCK_FRACTION_RANGE[0]),
                   compaction.ROCK_FRACTION_RANGE[1])
        try:
            expected = compaction.compacted_density_kg_m3(moon.radius, rock) / 1000.0
        except compaction.OutOfRange:
            continue  # a moon just under 1e21 kg whose radius the table does not reach: the structure solver made it
        assert moon.density == pytest.approx(expected, rel=5e-3)  # the radius comes from the mass by the table's inverse (worst case 0.2 percent near the taper end)
        assert moon.porosity is not None and 0.0 <= moon.porosity < 1.0
        checked += 1
    assert checked
    big = [moon for moon in generated_moons if moon.mass >= compaction.TAPER_END_KG]
    assert big and all(moon.porosity is None for moon in big)


def test_the_same_seed_gives_the_same_moons():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem

    def masses():
        out = []
        with draw.bound(33):
            for _ in range(6):
                config = SystemConfig()
                config.MOONS = True
                system = StarSystem(system_config=config)
                out.extend((moon.moon_origin, moon.mass, moon.orbital_inclination_deg)
                           for planet in system.planets if planet.body_type != "a" for moon in planet.moons)
        return out
    assert masses() == masses()


# --- GEN.1081: the captured-origin weight depends on the moon's mass -----------------------------------------------------------

def test_a_small_moon_keeps_the_open_weights():
    assert mc.captured_weight(1e12) == 1.0 and mc.captured_weight(mc.CAPTURED_FULL_WEIGHT_MASS_KG) == 1.0
    assert [p for _name, p in mc.origin_probabilities(mc.GIANT, 1e15)] == pytest.approx([0.5, 0.5])
    assert [p for _name, p in mc.origin_probabilities(mc.GIANT, 1e15, has_ring=True, roche_ratio=2.0)] == pytest.approx([1 / 3] * 3)


@pytest.mark.parametrize("mass_kg,weight", [(1e20, 0.0464), (1e22, 0.01), (2.14e22, 0.0078)])
def test_a_heavy_moon_is_rarely_captured(mass_kg, weight):
    assert mc.captured_weight(mass_kg) == pytest.approx(weight, rel=0.02)
    shares = dict(mc.origin_probabilities(mc.GIANT, mass_kg))
    assert shares["captured"] < 0.05 and shares["disk"] > 0.95 and sum(shares.values()) == pytest.approx(1.0)


def test_the_weight_falls_with_mass_and_every_parent_kind_keeps_its_rows():
    masses = [1e16, 1e17, 1e18, 1e19, 1e20, 1e21]
    weights = [mc.captured_weight(m) for m in masses]
    assert weights == sorted(weights, reverse=True) and weights[-1] < 0.1
    for kind in mc.PARENT_KINDS:
        names = [name for name, _p in mc.origin_probabilities(kind, 1e20)]
        assert "captured" in names and sum(p for _name, p in mc.origin_probabilities(kind, 1e20)) == pytest.approx(1.0)
    assert dict(mc.origin_probabilities(mc.TERRESTRIAL, 7.3e22))["impact"] > 0.99


def test_the_weight_refuses_a_mass_that_is_not_positive():
    for mass in (0.0, -1.0):
        with pytest.raises(ValueError):
            mc.captured_weight(mass)
    with pytest.raises(ValueError):
        mc.origin_probabilities("moon", 1e20)
