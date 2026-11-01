# tests/test_emission_members.py

"""
GEN.321: the ionizing members of an emission region (class D or E) are real stars (`nebula_contents.emission_members`):
the types the nebula's `hosts` text counts less the seed star, masses in each type's range and alive at the nebula's age,
points in a Plummer sphere of half-mass radius 1 pc (D) or 2 pc (E) clipped to the nebula's radius, all from the seed
star's point alone, so every sector that reaches a region builds the same members and each builds those in its cube.
"""

import math
import types

import numpy as np
import pytest

from planetgen import tuning
from planetgen.generation import bright_stars as brightStars, nebula_contents, run_sector, stored_trees
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import hii, stellar_evolution as evolution
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw
from tests.test_nebula_contents import _Sector, _args

SEED = bytes(range(16))
HOST_POINT = (8000.0, 20.0, 5.0)


def _star(number, type_="O5V"):
    return types.SimpleNamespace(type=f"{type_} Blue", age=0.002, mass=2.0e31 * (1.0 + number * 1e-6),
                                 temperature=40000.0)


def _regions(letter, count):
    """`count` hosted nebulae of class `letter`, with the host's type."""
    found = []
    number = 0
    while len(found) < count:
        number += 1
        nebula = nebula_contents.hosted_nebula_of(_star(number))
        if nebula is not None and nebula.nebula_class == letter:
            found.append(nebula)
    return found


D_REGIONS, E_REGIONS = _regions("D", 12), _regions("E", 6)


def test_a_types_mass_range_is_the_masses_that_name_it_and_the_ranges_tile_the_ionizing_stars():
    names = [anchor[0] for anchor in hii.TYPE_ANCHORS]
    ranges = [hii.mass_range_of_type(name, 150.0) for name in names]
    assert ranges[0][1] == 150.0 and ranges[-1][0] == hii.MIN_IONIZING_MASS_SOL
    for hotter, cooler in zip(ranges[:-1], ranges[1:]):
        assert cooler[1] == pytest.approx(hotter[0])        # each range ends where the next hotter one starts
    for name, (low, high) in zip(names, ranges):
        assert low < high
        for mass in np.geomspace(low * 1.001, high * 0.999, 9):
            assert hii.type_of_mass(float(mass)) == name
    with pytest.raises(ValueError):
        hii.mass_range_of_type("G2V", 150.0)


def test_the_members_are_the_hosts_text_less_the_seed_star():
    for nebula in D_REGIONS + E_REGIONS:
        members = nebula_contents.emission_members(SEED, nebula, "O5V", HOST_POINT)
        text = {name: int(count) for name, count in nebula_contents._HOST_COUNT.findall(nebula.hosts)}
        text["O5V"] = text.get("O5V", 0) - 1
        wanted = {name: count for name, count in text.items() if count}
        assert dict(nebula_contents.member_counts(nebula, "O5V")) == wanted
        assert dict(hii.host_counts([m.mass_sol for m in members])) == wanted     # the same types, star by star
        assert len(members) == sum(wanted.values())


def test_the_members_are_alive_at_the_nebulas_age_and_share_it():
    for nebula in D_REGIONS + E_REGIONS:
        for member in nebula_contents.emission_members(SEED, nebula, "O5V", HOST_POINT):
            assert member.mass_sol >= hii.MIN_IONIZING_MASS_SOL
            assert evolution.main_sequence_lifetime_gy(member.mass_sol) > nebula.age_myr / 1000.0
            assert member.age_yr == pytest.approx(nebula.age_myr * 1.0e6)


def test_the_members_lie_in_a_plummer_sphere_about_the_seed_inside_the_nebula():
    for letter, regions in (("D", D_REGIONS), ("E", E_REGIONS)):
        half = tuning.NEBULA_CLUSTER_HALF_MASS_PC[letter]
        radii = []
        for k, nebula in enumerate(regions):
            point = (HOST_POINT[0] + 11.0 * k, HOST_POINT[1], HOST_POINT[2])
            for member in nebula_contents.emission_members(SEED, nebula, "O5V", point):
                r = math.dist(member.position_pc, point)
                assert r <= ly_to_pc(nebula.radius_ly) * (1 + 1e-9)
                radii.append(r)
        assert len(radii) > 100
        assert np.median(radii) == pytest.approx(half, rel=0.2)        # half the members inside the half-mass radius
        assert max(radii) <= nebula_contents.emission_reach_pc() + 1e-9


def test_a_small_region_clips_its_members_to_its_radius():
    nebula = D_REGIONS[0]
    nebula.radius_ly = 0.6     # 0.2 pc, far inside the sphere's half-mass radius
    members = nebula_contents.emission_members(SEED, nebula, "O5V", HOST_POINT)
    assert members and all(math.dist(m.position_pc, HOST_POINT) <= ly_to_pc(0.6) * (1 + 1e-9) for m in members)


def test_the_members_come_from_the_galaxy_the_nebula_and_the_seed_stars_point_alone():
    nebula = E_REGIONS[0]
    first = nebula_contents.emission_members(SEED, nebula, "O5V", HOST_POINT)
    again = nebula_contents.emission_members(SEED, nebula, "O5V", HOST_POINT)
    key = lambda found: [(m.position_pc, m.mass_sol, m.age_yr, m.seed) for m in found]
    assert key(first) == key(again)
    draw.set_run_seed(99)        # other draws in between change nothing
    assert key(nebula_contents.emission_members(SEED, nebula, "O5V", HOST_POINT)) == key(first)
    assert [m.seed for m in nebula_contents.emission_members(bytes(16), nebula, "O5V", HOST_POINT)] != [m.seed for m in first]
    moved = (HOST_POINT[0] + 1.0, HOST_POINT[1], HOST_POINT[2])
    assert [m.seed for m in nebula_contents.emission_members(SEED, nebula, "O5V", moved)] != [m.seed for m in first]
    assert len({m.seed for m in first}) == len(first)


def test_a_region_of_another_class_has_no_members():
    other = types.SimpleNamespace(nebula_class="C", hosts="O5V x1", age_myr=0.5, radius_ly=10.0)
    assert nebula_contents.emission_members(SEED, other, "O5V", HOST_POINT) == []
    assert nebula_contents.member_counts(types.SimpleNamespace(nebula_class="D", hosts=None), "O5V") == []


class _Frame(_Sector):
    """A sector holding one side of the plane x = 0 (sector-local light-years), as `_Sector` does."""


def _cloud_field(centre_pc):
    return (SEED, None, centre_pc, 3.0)


def test_each_sector_builds_only_the_members_inside_it_and_together_they_make_every_one():
    """A region across two sectors: both reach it, the seed star's own and a neighbour's, and none is built twice."""
    nebula = D_REGIONS[1]
    centre_pc = (8000.0, 0.0, 0.0)
    field = _cloud_field(centre_pc)
    args = _args()
    built = []
    wanted = 0
    for k in range(3):
        from planetgen.galaxy.frame import local_to_galaxy_pc
        point = local_to_galaxy_pc(centre_pc, tuple(ly_to_pc(v) for v in (0.001 * k, 0.0, 0.0)))
        wanted += len(nebula_contents.emission_members(SEED, nebula, "O5V", point))
        for side in (True, False):
            sector = _Frame([], side)
            own = [(nebula, "O5V", (0.001 * k, 0.0, 0.0))] if side else []
            neighbour = [] if side else [(nebula, "O5V", point)]
            run_sector.add_emission_members(sector, args, field, own, neighbour, 26000.0)
            built += sector.systems
    assert wanted >= 30 and len(built) == wanted
    for entry in built:
        star = entry.star_system.primary_star
        assert star.mass / 1.98840987e30 >= hii.MIN_IONIZING_MASS_SOL


def test_the_sector_that_lights_a_region_and_the_one_beside_it_give_the_same_members():
    nebula = E_REGIONS[1]
    centre_pc = (8000.0, 0.0, 0.0)
    from planetgen.galaxy.frame import local_to_galaxy_pc
    local = (0.002, 0.0, 0.0)
    point = local_to_galaxy_pc(centre_pc, tuple(ly_to_pc(v) for v in local))
    field = _cloud_field(centre_pc)
    own = _Frame([], True)
    run_sector.add_emission_members(own, _args(), field, [(nebula, "O5V", local)], (), 26000.0)
    beside = _Frame([], False)
    run_sector.add_emission_members(beside, _args(), field, [], [(nebula, "O5V", point)], 26000.0)
    seeds = [e.system_seed for e in own.systems + beside.systems]
    expected = [m.seed for m in nebula_contents.emission_members(SEED, nebula, "O5V", point)]
    assert sorted(seeds) == sorted(expected) and len(own.systems) > 20 and len(beside.systems) > 20


def test_a_pre_placed_hot_star_rolls_the_same_region_as_the_star_its_sector_builds():
    """`row_host` reads only the row, and the sector that builds the row's star rolls the same region. The star is built as the
    fill builds a row with no stored members (`stored_trees.SINGLE`): a plain `StarSystem` draws a tree of its own (GEN.631),
    and a heavier companion drawn there took `stars[0]` (TEST.340)."""
    done = 0
    for k, mass in enumerate(np.linspace(14.0, 40.0, 400)):
        age_gy = 0.002
        if evolution.main_sequence_lifetime_gy(float(mass)) <= age_gy:
            continue
        params = evolution.star_params(float(mass), age_gy, evolution.evolve_star(float(mass), age_gy))
        row = {"star_type": params["type"], "yerkes_class": params["yerkes_class"], "mass_kg": params["mass_kg"],
               "radius_km": params["radius_km"], "temperature_k": params["temperature_k"],
               "luminosity_w": params["luminosity_w"], "age_gy": params["age_gy"], "lifespan_gy": params["lifespan_gy"],
               "initial_mass_sol": params["initial_mass_sol"], "phase_end_age_gy": params["phase_end_age_gy"],
               "metallicity_fe_h": -0.1, "position_x_mpc": 1.0e6 + k, "position_y_mpc": 2.0e3, "position_z_mpc": -50.0}
        config = SystemConfig()
        config.POPULATION = "young"
        with draw.bound(1000 + k):
            system = StarSystem(system_config=config, primary_star_params=brightStars.star_params(row),
                                member_tree=stored_trees.SINGLE)
        built = nebula_contents.hosted_nebula_of(system.stars[0])
        found = nebula_contents.row_host(row)
        if built is None or built.nebula_class not in "DE":
            assert found is None
            continue
        nebula, name, point = found
        assert (nebula.nebula_class, nebula.hosts, nebula.radius_ly, nebula.age_myr) == (
            built.nebula_class, built.hosts, built.radius_ly, built.age_myr)
        assert name == nebula_contents.host_name_of(system.stars[0])
        assert point == (1.0e3 + k / 1000.0, 2.0, -0.05)
        done += 1
    assert done >= 3


# --- The neighbours' hot stars, read from the stored rows -----------------------------------------------------------


def _hot_row(point_pc, address, k):
    """A `bright_stars` row (a dict of every column the scatter writes) of an O star at `point_pc` in `address`."""
    mass = 20.0 + 0.01 * k
    params = evolution.star_params(mass, 0.002, evolution.evolve_star(mass, 0.002))
    return {"shell_index": address[0], "band_index": address[1], "wedge_index": address[2],
            "position_x_mpc": point_pc[0] * 1000.0, "position_y_mpc": point_pc[1] * 1000.0,
            "position_z_mpc": point_pc[2] * 1000.0, "population": "young", "star_type": params["type"],
            "yerkes_class": params["yerkes_class"], "mass_kg": params["mass_kg"], "radius_km": params["radius_km"],
            "temperature_k": params["temperature_k"], "luminosity_w": params["luminosity_w"],
            "age_gy": params["age_gy"], "lifespan_gy": params["lifespan_gy"],
            "initial_mass_sol": params["initial_mass_sol"], "phase_end_age_gy": params["phase_end_age_gy"],
            "seed": 12345 + k, "metallicity_fe_h": -0.1,
            "rotation_period_hours": brightStars.rotation_period_hours(12345 + k, params),
            "nuclear_arrival_age_gy": params.get("nuclear_arrival_age_gy")}


@pytest.mark.db
def test_a_sector_reads_the_regions_of_its_neighbours_hot_stars_from_their_rows(mysql_config):
    from planetgen.db import store
    from planetgen.galaxy import spherical
    from planetgen.galaxy.profile import sector_profile
    from tests.test_nebula_field import SHAPE
    from tests.test_phenomenon_scatter import E_VALUE, EDGE_PC, EXTENTS
    from planetgen.generation import run_galaxy

    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8,
                            expected_system_count_at_density_1=E_VALUE, config=mysql_config, galaxy_seed=SEED)
    store.replace_galaxy_shells(EXTENTS, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        skeleton = store.get_galaxy_shape(conn)
        host_point = (SHAPE.disk_scale_length_pc * tuning.REFERENCE_RADIUS_SCALE_LENGTHS, 0.0, 0.0)
        host_address = spherical.sector_address_at(host_point, EDGE_PC)
        factor = sector_profile(host_address, skeleton).kind_factors["nebula"]
        assert factor > 0.0
        row = next(r for r in (_hot_row(host_point, host_address, k) for k in range(2000))
                   if nebula_contents.row_host(r, factor) is not None)
        store.insert_bright_stars(conn, [[row[column] for column in store.BRIGHT_STAR_COLUMNS]])
        conn.commit()
        nearby = (host_point[0], host_point[1] + 5.0, host_point[2])
        beside = spherical.sector_address_at(nearby, EDGE_PC)
        assert beside != host_address
        profile = sector_profile(beside, skeleton)
        found = run_galaxy._emission_hosts(conn, skeleton, beside, profile)
        want = nebula_contents.row_host(row, factor)
        assert len(found) == 1
        assert (found[0][0].hosts, found[0][1], found[0][2]) == (want[0].hosts, want[1], want[2])
        # The host's own sector makes its region itself, and a sector out of reach reads none.
        assert run_galaxy._emission_hosts(conn, skeleton, host_address, sector_profile(host_address, skeleton)) == []
        far = spherical.sector_address_at((host_point[0] + 200.0, 0.0, 0.0), EDGE_PC)
        assert run_galaxy._emission_hosts(conn, skeleton, far, sector_profile(far, skeleton)) == []
    finally:
        conn.close()


def test_a_hosted_nebula_reports_its_seed_star_for_the_members_to_be_placed_from():
    class _Host(_Frame):
        def __init__(self, stars):
            super().__init__([], True)
            self.entries = [types.SimpleNamespace(star_system=types.SimpleNamespace(stars=[star], placed_star=star, name=f"S{n}"),
                                                  position=(0.1 * n, 0.0, 0.0)) for n, star in enumerate(stars)]
            self.phenomena = []
            self.name = "Test"

        def add_phenomenon(self, nebula, kind, position=None):
            self.phenomena.append((nebula, kind, position))
            return types.SimpleNamespace(phenomenon=nebula)

    stars = [_star(number) for number in range(200)]
    sector = _Host(stars)
    reported = []
    entries = run_sector.add_star_hosted_nebulae(sector, _args(), hosts_out=reported)
    assert entries and len(reported) == len(entries) == len(sector.phenomena)
    for (nebula, name, position), (placed, _kind, where) in zip(reported, sector.phenomena):
        assert nebula is placed and name == "O5V" and position == where
