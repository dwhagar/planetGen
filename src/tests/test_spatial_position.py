# tests/test_spatial_position.py

"""GEN.74: the one point-in-space object (planetgen/physics/position.py)."""

import math

import pytest

from planetgen.physics import constants
from planetgen.physics import position as pos
from planetgen.physics.position import SpatialPosition3D

SECTOR = (1.0e17, -2.0e17, 3.0e16)
STAR = (1.0e17 + 5.0e12, -2.0e17 + 1.0e12, 3.0e16 - 2.0e12)


def _planet():
    return SpatialPosition3D((STAR[0] + 1.5e11, STAR[1], STAR[2]), SECTOR, star_center_galactic=STAR, mass_kg=6.0e24)


def _close(a, b, rel=1e-9, abs_=1e-6):
    return all(x == pytest.approx(y, rel=rel, abs=abs_) for x, y in zip(a, b))


def test_a_new_position_has_every_frame_and_form():
    p = _planet()
    for frame in pos.FRAMES:
        for form in pos.FORMS:
            assert p.get_coordinates(frame, form) is not None
    assert _close(p.get_coordinates("system", "cartesian"), (1.5e11, 0.0, 0.0), abs_=1e5)
    assert p.get_coordinates("system", "cylindrical")[0] == pytest.approx(1.5e11, rel=1e-6)
    assert p.get_coordinates("system", "spherical")[2] == pytest.approx(math.pi / 2, abs=1e-6)


@pytest.mark.parametrize("frame", pos.FRAMES)
@pytest.mark.parametrize("form,values", [
    ("cartesian", (3.0e11, -4.0e11, 1.0e11)),
    ("cylindrical", (5.0e11, 0.7, -2.0e10)),
    ("spherical", (6.0e11, -2.1, 1.1)),
])
def test_setting_any_coordinate_in_any_frame_updates_every_other(frame, form, values):
    p = _planet()
    p.set_coordinates(frame, form, values)
    assert _close(p.get_coordinates(frame, form), values, rel=1e-9, abs_=1.0)
    cart = p.get_coordinates(frame, "cartesian")
    gal = p.get_coordinates("galactic", "cartesian")
    offsets = {"galactic": (0, 0, 0), "sector": SECTOR, "system": STAR}[frame]
    assert _close(gal, tuple(c + o for c, o in zip(cart, offsets)), rel=1e-12, abs_=1e5)
    # Every form agrees with its own Cartesian.
    for f in pos.FRAMES:
        c = p.get_coordinates(f, "cartesian")
        assert _close(pos.cylindrical_to_cartesian(*p.get_coordinates(f, "cylindrical")), c, abs_=1e5)
        assert _close(pos.spherical_to_cartesian(*p.get_coordinates(f, "spherical")), c, abs_=1e5)


def test_a_system_offset_keeps_its_own_metres_exactly():
    """Galactic coordinates (~1e20 m) can't hold a moon's metres; the frame it was set in does."""
    p = _planet()
    p.set_system_cartesian(3.844e8, 1.0, -2.0)
    assert p.get_coordinates("system", "cartesian") == (3.844e8, 1.0, -2.0)


def test_a_star_has_no_system_frame():
    star = SpatialPosition3D(STAR, SECTOR, is_star=True, star_center_galactic=STAR)
    assert star.get_coordinates("system", "cartesian") is None
    with pytest.raises(ValueError):
        star.set_system_cartesian(1.0, 2.0, 3.0)
    unbound = SpatialPosition3D(STAR, SECTOR)
    assert unbound.get_coordinates("system", "spherical") is None
    with pytest.raises(ValueError):
        unbound.set_system_cylindrical(1.0, 0.0, 0.0)


def test_moving_an_anchor_keeps_the_body_where_it_is():
    p = _planet()
    gal = p.get_coordinates("galactic", "cartesian")
    p.set_nearest_star_center((STAR[0] + 1.0e11, STAR[1], STAR[2]))
    assert _close(p.get_coordinates("galactic", "cartesian"), gal, abs_=1e5)
    assert _close(p.get_coordinates("system", "cartesian"), (5.0e10, 0.0, 0.0), abs_=1e5)
    p.set_sector_center((0.0, 0.0, 0.0))
    assert _close(p.get_coordinates("sector", "cartesian"), gal, abs_=1e5)
    p.set_nearest_star_center(None)
    assert p.get_coordinates("system", "cartesian") is None
    assert _close(p.get_coordinates("galactic", "cartesian"), gal, abs_=1e5)


def test_an_object_set_by_system_loses_its_frame_gracefully_when_the_star_is_dropped():
    p = _planet()
    p.set_system_cartesian(1.0e11, 0.0, 0.0)
    gal = p.get_coordinates("galactic", "cartesian")
    p.set_nearest_star_center(None)
    assert _close(p.get_coordinates("galactic", "cartesian"), gal, abs_=1e5)


def test_the_angles_at_the_origin_and_on_the_axis():
    p = SpatialPosition3D((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    assert p.get_coordinates("galactic", "spherical") == (0.0, 0.0, 0.0)
    assert p.get_coordinates("galactic", "cylindrical") == (0.0, 0.0, 0.0)
    p.set_galactic_cartesian(0.0, 0.0, -5.0)
    assert p.get_coordinates("galactic", "spherical") == (5.0, 0.0, math.pi)


@pytest.mark.parametrize("call", [
    lambda p: p.set_galactic_cartesian(float("nan"), 0, 0),
    lambda p: p.set_sector_cartesian(0, float("inf"), 0),
    lambda p: p.set_galactic_spherical(-1.0, 0.0, 1.0),
    lambda p: p.set_galactic_spherical(1.0, 0.0, 4.0),
    lambda p: p.set_galactic_cylindrical(-1.0, 0.0, 0.0),
    lambda p: p.set_velocity_cartesian(pos.SPEED_OF_LIGHT_MS, 0, 0),
    lambda p: p.set_mass(-1.0),
    lambda p: p.set_mass(float("nan")),
])
def test_values_the_maths_cannot_hold_are_refused(call):
    p = _planet()
    before = p.get_coordinates("galactic", "cartesian")
    with pytest.raises(ValueError):
        call(p)
    assert p.get_coordinates("galactic", "cartesian") == before


def test_unknown_frames_and_forms_raise_key_errors():
    p = _planet()
    with pytest.raises(KeyError):
        p.get_coordinates("planetary", "cartesian")
    with pytest.raises(KeyError):
        p.get_coordinates("galactic", "polar")
    with pytest.raises(KeyError):
        p.set_coordinates("sector", "polar", (1, 2, 3))


def test_velocity_speed_and_direction():
    p = _planet()
    assert p.get_speed() == 0.0 and p.get_velocity_direction() == (0.0, 0.0, 0.0)
    p.set_velocity_cartesian(3.0, 4.0, 0.0)
    assert p.get_speed() == 5.0
    assert p.get_velocity_vector() == (3.0, 4.0, 0.0)
    assert p.get_velocity_direction() == (0.6, 0.8, 0.0)


def test_time_to_observable_movement_follows_speed_and_the_scale_threshold():
    p = _planet()
    assert p.get_time_to_observable_movement("galactic") == pos.MAX_UPDATE_INTERVAL_S
    p.set_velocity_cartesian(1.0e4, 0.0, 0.0)
    assert p.get_time_to_observable_movement("system") == pytest.approx(0.01 * constants.AU_M / 1.0e4)
    assert p.get_time_to_observable_movement("planetary") == pytest.approx(1.0e8 / 1.0e4)
    assert p.get_time_to_observable_movement("galactic") == pytest.approx(pos.THRESHOLDS_M["galactic"] / 1.0e4)
    p.set_velocity_cartesian(1.0e-12, 0.0, 0.0)
    assert p.get_time_to_observable_movement("planetary") == pos.MAX_UPDATE_INTERVAL_S
    with pytest.raises(KeyError):
        p.get_time_to_observable_movement("sector")


def test_the_thresholds_are_the_design_documents():
    assert pos.THRESHOLDS_M["galactic"] == pytest.approx(0.01e-3 * constants.PARSEC_M)
    assert pos.THRESHOLDS_M["system"] == pytest.approx(0.01 * constants.AU_M)
    assert pos.THRESHOLDS_M["planetary"] == 1.0e8


def test_mass_and_mu_follow_each_other():
    p = SpatialPosition3D((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    assert p.mass_kg is None and p.mu is None
    p.set_mass(5.972e24)
    assert p.mu == pytest.approx(3.986e14, rel=1e-3)
    p.set_mu(1.32712440018e20)
    assert p.mass_kg == pytest.approx(1.98841e30, rel=1e-4)
    assert pos.mass_of(pos.mu_of(7.0)) == pytest.approx(7.0)


def test_coordinates_are_kept_in_the_unit_given_with_no_round_trip():
    ly = constants.LY_TO_M
    p = SpatialPosition3D((10.0, 20.0, 30.0), (9.0, 18.0, 27.0), length_unit_m=ly)
    assert p.get_coordinates("sector", "cartesian") == (1.0, 2.0, 3.0)
    p.set_sector_cartesian(0.1, 0.2, 0.3)
    assert p.get_coordinates("sector", "cartesian") == (0.1, 0.2, 0.3)
    assert p.get_coordinates("galactic", "cartesian") == pytest.approx((9.1, 18.2, 27.3))
    with pytest.raises(ValueError):
        SpatialPosition3D((0, 0, 0), (0, 0, 0), length_unit_m=0.0)


def test_carrying_an_anchor_moves_the_body_with_it():
    p = _planet()
    system = p.get_coordinates("system", "cartesian")
    sector = p.get_coordinates("sector", "cartesian")
    p.carry_star_center((0.0, 0.0, 0.0))
    assert p.get_coordinates("system", "cartesian") == system
    assert _close(p.get_coordinates("galactic", "cartesian"), system, abs_=1.0)
    p.carry_sector_center((5.0, 5.0, 5.0))
    assert p.get_coordinates("sector", "cartesian") == pytest.approx(
        tuple(g - 5.0 for g in p.get_coordinates("galactic", "cartesian")))
    star = SpatialPosition3D(STAR, SECTOR, is_star=True)
    keep = star.get_coordinates("sector", "cartesian")
    star.carry_sector_center((0.0, 0.0, 0.0))
    assert star.get_coordinates("sector", "cartesian") == keep
    with pytest.raises(ValueError):
        star.carry_star_center((0.0, 0.0, 0.0))


# --- Sector entries hold one (GEN.74 part 2) ------------------------------------------

def test_a_sector_entry_holds_one_with_its_mass_and_mu():
    from planetgen.galaxy.sector import SpaceSector
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.phenomena.compact_remnant import BlackHole
    from planetgen.generation.system import StarSystem

    sector = SpaceSector("Holders", edge_ly=11.5)
    system = StarSystem(system_config=SystemConfig())
    entry = sector.add_system(system, position=(1.0, -2.0, 0.5))
    assert entry.position == (1.0, -2.0, 0.5)
    assert entry.spatial.get_coordinates("sector", "cartesian") == (1.0, -2.0, 0.5)
    assert entry.spatial.mass_kg == pytest.approx(sum(star.mass for star in system.stars))
    assert entry.spatial.mu == pytest.approx(constants.G * entry.spatial.mass_kg)
    assert entry.spatial.get_coordinates("system", "cartesian") is None
    hole = sector.add_phenomenon(BlackHole(SystemConfig()), "black-hole", position=(-1.0, 1.0, 0.0))
    assert hole.spatial.mass_kg == hole.phenomenon.mass
    nebula_entry = sector.phenomena[0]
    assert nebula_entry is hole


def test_placing_a_sector_in_the_galaxy_carries_its_entries():
    from planetgen.galaxy.sector import SpaceSector
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem

    sector = SpaceSector("Placed", edge_ly=11.5)
    entry = sector.add_system(StarSystem(system_config=SystemConfig()), position=(1.0, 2.0, 3.0))
    sector.place_in_galaxy((100.0, 200.0, -50.0))
    assert entry.position == (1.0, 2.0, 3.0)
    assert entry.spatial.get_coordinates("galactic", "cartesian") == pytest.approx((101.0, 202.0, -47.0))
    later = sector.add_system(StarSystem(system_config=SystemConfig()), position=(-1.0, 0.0, 0.0))
    assert later.spatial.get_coordinates("galactic", "cartesian") == pytest.approx((99.0, 200.0, -50.0))


# --- Bodies in a system hold one (GEN.74 part 2b) --------------------------------------

def _system_with(predicate):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem

    for _ in range(200):
        system = StarSystem(system_config=SystemConfig())
        if predicate(system):
            return system
    pytest.fail("no system matching")


def _placed_system(predicate=lambda system: bool(system.planets)):
    from planetgen.galaxy.sector import SpaceSector

    system = _system_with(predicate)
    sector = SpaceSector("Bodies", edge_ly=11.5)
    sector.place_in_galaxy((26000.0, -40.0, 12.0))
    entry = sector.add_system(system, position=(1.5, -2.0, 0.25))
    return sector, entry, system


def _system_origin_au(entry, system):
    """Where the system frame's origin is, in AU in the galactic frame: the
    single or primary star, but for a close pair its barycenter (the system's
    center), which the primary star sits a fraction of an AU off. For a wide
    pair the primary star is thousands of AU from the system's center."""
    if getattr(system, "binary_type", None) == "close":
        return tuple(c * constants.LY_TO_AU for c in entry.spatial.get_coordinates("galactic", "cartesian"))
    return system.primary_star.spatial.get_coordinates("galactic", "cartesian")


def test_a_planet_holds_one_whose_system_frame_is_its_offset_from_the_star():
    _sector, entry, system = _placed_system(lambda s: any(getattr(p, "spatial", None) for p in s.planets))
    planet = next(p for p in system.planets if getattr(p, "spatial", None))
    system_frame = planet.spatial.get_coordinates("system", "cartesian")
    assert (planet.position_x, planet.position_y, planet.position_z) == pytest.approx(system_frame, rel=1e-9, abs=1e-9)
    assert planet.spatial.length_unit_m == constants.AU_M
    assert planet.spatial.mass_kg == planet.mass and planet.spatial.mu == pytest.approx(constants.G * planet.mass)
    # TEST.108: a random close binary used to fail this on the primary star.
    center_galactic = _system_origin_au(entry, system)
    galactic = planet.spatial.get_coordinates("galactic", "cartesian")
    assert galactic == pytest.approx(tuple(c + o for c, o in zip(center_galactic, system_frame)), rel=1e-9, abs=1e-3)


def test_placing_a_system_keeps_every_offset_and_follows_in_galactic_coordinates():
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.galaxy.sector import SpaceSector

    sector = SpaceSector("Move", edge_ly=11.5)
    system = _system_with(lambda s: any(getattr(p, "spatial", None) for p in s.planets))
    planet = next(p for p in system.planets if getattr(p, "spatial", None))
    offset = (planet.position_x, planet.position_y, planet.position_z)
    entry = sector.add_system(system, position=(1.0, 1.0, 1.0))
    sector.place_in_galaxy((100.0, 0.0, 0.0))
    assert (planet.position_x, planet.position_y, planet.position_z) == offset
    # The planet sits an orbit (AU) off its system's origin: a wide pair's is its
    # primary star, thousands of AU from the system's center (TEST.108).
    origin = _system_origin_au(entry, system)
    assert planet.spatial.get_coordinates("galactic", "cartesian")[0] == pytest.approx(origin[0] + offset[0], abs=1e-3)
    assert entry.spatial.get_coordinates("galactic", "cartesian")[0] == pytest.approx(101.0)
    assert entry.position == (1.0, 1.0, 1.0)


def test_a_moon_is_anchored_on_its_planet():
    _sector, _entry, system = _placed_system(lambda s: any(getattr(p, "moons", None) for p in s.planets))
    planet = next(p for p in system.planets if getattr(p, "moons", None))
    moon = planet.moons[0]
    planet_galactic = planet.spatial.get_coordinates("galactic", "cartesian")
    moon_galactic = moon.spatial.get_coordinates("galactic", "cartesian")
    moon_offset = (moon.position_x, moon.position_y, moon.position_z)
    assert moon_galactic == pytest.approx(tuple(p + o for p, o in zip(planet_galactic, moon_offset)), rel=1e-12)


def test_every_star_holds_one_and_a_binary_pair_straddles_the_system_center():
    _sector, entry, system = _placed_system(lambda s: len(s.stars) == 2)
    for star in system.stars:
        assert star.spatial.is_star() and star.spatial.mass_kg == star.mass
    first, second = (s.spatial.get_coordinates("galactic", "cartesian") for s in system.stars)
    assert first != second
    assert system.primary_star.spatial.get_coordinates("sector", "cartesian") != (0.0, 0.0, 0.0)


def test_a_comet_holds_one():
    _sector, _entry, system = _placed_system(lambda s: bool(s.comets))
    comet = system.comets[0]
    assert (comet.position_x_au, comet.position_y_au, comet.position_z_au) == comet.spatial.get_coordinates(
        "system", "cartesian")


def test_a_planet_rebuilt_from_its_dict_has_its_position_back():
    from planetgen.generation.planet import Planet

    _sector, _entry, system = _placed_system(lambda s: any(getattr(p, "spatial", None) for p in s.planets))
    planet = next(p for p in system.planets if getattr(p, "spatial", None))
    again = Planet.from_dict(planet.to_dict(), system.star, system.system_config)
    assert (again.position_x, again.position_y, again.position_z) == (
        planet.position_x, planet.position_y, planet.position_z)
    assert again.spatial.get_coordinates("system", "cartesian") == planet.spatial.get_coordinates("system", "cartesian")


def test_a_close_pairs_planets_are_offset_from_its_center_not_its_primary_star():
    """TEST.108: the system frame's origin is the barycenter, which a close
    binary's primary star sits off by a fraction of an AU (a wide pair's is
    its primary star, thousands of AU from the system center)."""
    _sector, entry, system = _placed_system(
        lambda s: getattr(s, "binary_type", None) == "close" and any(getattr(p, "spatial", None) for p in s.planets))
    planet = next(p for p in system.planets if getattr(p, "spatial", None))
    center = _system_origin_au(entry, system)
    star = system.primary_star.spatial.get_coordinates("galactic", "cartesian")
    assert max(abs(c - s) for c, s in zip(center, star)) > 1e-3
    offset = planet.spatial.get_coordinates("system", "cartesian")
    assert planet.spatial.get_coordinates("galactic", "cartesian") == pytest.approx(
        tuple(c + o for c, o in zip(center, offset)), rel=1e-9, abs=1e-3)


# --- Velocity in two frames, and the epoch -----------------------------------------

STAR_V = (2.0e5, 1.0e5, 0.0)


def _planet_moving():
    return SpatialPosition3D((STAR[0] + 1.5e11, STAR[1], STAR[2]), SECTOR, star_center_galactic=STAR,
                             star_velocity_galactic=STAR_V, mass_kg=6.0e24)


def test_a_system_velocity_is_the_stars_velocity_plus_the_relative_one():
    p = _planet_moving()
    p.set_velocity_cartesian(0.0, 3.0e4, 0.0, frame="system")
    assert p.get_velocity_vector("system") == (0.0, 3.0e4, 0.0)
    assert _close(p.get_velocity_vector("galactic"), (2.0e5, 1.3e5, 0.0))
    assert p.get_velocity_vector("sector") == p.get_velocity_vector("galactic")
    assert p.get_speed("system") == pytest.approx(3.0e4)
    assert p.get_velocity_direction("system") == (0.0, 1.0, 0.0)


def test_a_galactic_velocity_derives_the_system_one():
    p = _planet_moving()
    p.set_velocity_cartesian(2.0e5, 1.3e5, 4.0e3)
    assert _close(p.get_velocity_vector("system"), (0.0, 3.0e4, 4.0e3))


def test_the_stars_velocity_changing_keeps_the_galactic_velocity_or_the_relative_one():
    p = _planet_moving()
    p.set_velocity_cartesian(0.0, 3.0e4, 0.0, frame="system")
    p.set_star_velocity((1.0e5, 0.0, 0.0))
    assert _close(p.get_velocity_vector("galactic"), (2.0e5, 1.3e5, 0.0))
    assert _close(p.get_velocity_vector("system"), (1.0e5, 1.3e5, 0.0))
    p.set_velocity_cartesian(0.0, 3.0e4, 0.0, frame="system")
    p.carry_star_velocity((1.0e5, 0.0, 0.0))
    assert _close(p.get_velocity_vector("system"), (0.0, 3.0e4, 0.0))
    assert _close(p.get_velocity_vector("galactic"), (1.0e5, 3.0e4, 0.0))


def test_a_star_or_a_body_without_a_star_has_no_system_velocity():
    star = SpatialPosition3D((0.0, 0.0, 0.0), SECTOR, is_star=True)
    assert star.get_velocity_vector("system") is None and star.get_speed("system") is None
    assert star.get_velocity_direction("system") is None
    with pytest.raises(ValueError):
        star.set_velocity_cartesian(1.0, 0.0, 0.0, frame="system")
    with pytest.raises(ValueError):
        star.carry_star_velocity((1.0, 0.0, 0.0))
    with pytest.raises(KeyError):
        star.get_velocity_vector("orbital")


def test_dropping_the_star_keeps_the_galactic_velocity():
    p = _planet_moving()
    p.set_velocity_cartesian(0.0, 3.0e4, 0.0, frame="system")
    p.set_nearest_star_center(None)
    assert _close(p.get_velocity_vector("galactic"), (2.0e5, 1.3e5, 0.0))
    assert p.get_velocity_vector("system") is None


def test_the_speed_of_light_limit_covers_both_frames():
    p = _planet_moving()
    with pytest.raises(ValueError):
        p.set_velocity_cartesian(pos.SPEED_OF_LIGHT_MS, 0.0, 0.0, frame="system")
    with pytest.raises(ValueError):
        p.set_velocity_cartesian(pos.SPEED_OF_LIGHT_MS - 1.0e5, 0.0, 0.0, frame="system")


def test_the_time_to_observable_movement_uses_the_relative_speed_inside_a_system():
    p = _planet_moving()
    p.set_velocity_cartesian(1.0e4, 0.0, 0.0, frame="system")
    assert p.get_time_to_observable_movement("system") == pytest.approx(0.01 * constants.AU_M / 1.0e4)
    assert p.get_time_to_observable_movement("galactic") == pytest.approx(
        pos.THRESHOLDS_M["galactic"] / p.get_speed("galactic"))


def test_the_epoch_is_kept_and_checked():
    p = _planet_moving()
    assert p.epoch_unix is None
    p.set_epoch_unix(946_728_000)
    assert p.epoch_unix == 946_728_000.0
    p.set_epoch_unix(None)
    assert p.epoch_unix is None
    with pytest.raises(ValueError):
        p.set_epoch_unix(float("nan"))
    assert SpatialPosition3D((0.0, 0.0, 0.0), SECTOR, epoch_unix=5).epoch_unix == 5.0
