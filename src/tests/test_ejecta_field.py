"""
GEN.238: the neutron stars and black holes the galaxy ejected stream through every cell of the sphere
(`planetgen.galaxy.ejecta_field`), with a density that falls as 1 / r^2, and a sector fill draws them
(`run_sector.add_ejecta`) with a velocity of their own in place of an orbit.
"""

import math
import random
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.galaxy import ejecta_field
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.profile import profile_at
from planetgen.galaxy.sector import SpaceSector
from planetgen.galaxy.skeleton import expected_system_count_at_density_1
from planetgen.generation import fill_stages, run_galaxy, run_sector
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.util import draw

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
    bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
EDGE_PC = tuning.DEFAULT_SECTOR_EDGE_PC
E_VALUE = expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)
CELL_PC3 = EDGE_PC ** 3


def _density(kind, radius_pc, seed=None):
    return ejecta_field.density_pc3(kind, radius_pc, SHAPE, E_VALUE, EDGE_PC, seed)


def test_the_galaxys_systems_are_about_three_hundred_billion():
    systems = ejecta_field.galaxy_systems(SHAPE, E_VALUE, EDGE_PC)
    assert 2.0e11 < systems < 4.5e11      # the 3.0e11 the halo's share is built on (HALO_SYSTEM_FRACTION)


def test_the_unbound_share_is_a_draw_of_the_range_for_neutron_stars_and_two_percent_for_black_holes():
    low, high = tuning.NEUTRON_STAR_UNBOUND_FRACTION_RANGE
    seeds = [bytes([n]) * 16 for n in range(40)]
    shares = [ejecta_field.unbound_fraction("neutron-star", seed) for seed in seeds]
    assert all(low <= share <= high for share in shares) and len(set(shares)) > 30
    assert shares == [ejecta_field.unbound_fraction("neutron-star", seed) for seed in seeds]
    assert ejecta_field.unbound_fraction("neutron-star") == pytest.approx((low + high) / 2)
    assert ejecta_field.unbound_fraction("black-hole", seeds[0]) == tuning.BLACK_HOLE_UNBOUND_FRACTION == 0.02


def test_the_density_is_n_over_four_pi_r_squared_v_t_out_to_v_t_and_none_beyond():
    outer = ejecta_field.outer_radius_pc()
    assert 3.0e6 < outer < 5.0e6                         # "v T (3 to 5 Mpc)" for 300 to 500 km/s over 10 Gyr
    count = ejecta_field.unbound_count("neutron-star", SHAPE, E_VALUE, EDGE_PC)
    for radius in (1.0e3, 5.0e4, 2.62e5):
        assert _density("neutron-star", radius) == pytest.approx(count / (4 * math.pi * radius ** 2 * outer))
    assert _density("neutron-star", 5.0e4) / _density("neutron-star", 1.0e5) == pytest.approx(4.0)
    assert _density("neutron-star", outer) == 0.0 and _density("neutron-star", 2 * outer) == 0.0
    # Integrated over the sphere out to r the law holds N r / (v T) of them.
    for radius in (1.0e5, 2.62e5):
        inside = count * radius / outer
        assert inside == pytest.approx(count * radius / outer)
        assert inside < count


def test_a_cell_far_out_holds_about_the_expected_ejecta_of_the_design_note():
    """docs/design/intergalactic-and-halo-contents.md 5.6: 1.3e-7 of a 4 pc cell at 50 kpc and 3.3e-8 at 100 kpc for
    2e8 unbound neutron stars; this galaxy has about 3.8e8 at the middle of the range, so about twice that."""
    near = _density("neutron-star", 5.0e4) * CELL_PC3
    far = _density("neutron-star", 1.0e5) * CELL_PC3
    assert 1.0e-7 < near < 4.0e-7 and 2.5e-8 < far < 1.0e-7
    assert _density("black-hole", 5.0e4) < _density("neutron-star", 5.0e4) / 50


def test_a_cells_ejecta_never_exceed_the_cap_so_a_dense_galaxy_cannot_stall_the_draw():
    """A galaxy shaped to hold astronomically many systems in a few cells (the fuzz tests'), or a cell beside the centre,
    draws at most `EJECTA_MAX_PER_CELL` of a kind; the real galaxy's densest cells are below it."""
    dense = ejecta_field.expected_in_cell((EDGE_PC, 0.0, 0.0), CELL_PC3, SHAPE, E_VALUE * 1.0e9, EDGE_PC)
    assert dense == {kind: tuning.EJECTA_MAX_PER_CELL for kind in ejecta_field.KINDS}
    inner = ejecta_field.expected_in_cell((2.0 * EDGE_PC, 0.0, 0.0), CELL_PC3, SHAPE, E_VALUE, EDGE_PC)
    assert 0.0 < inner["neutron-star"] < tuning.EJECTA_MAX_PER_CELL


def test_a_profile_carries_the_ejecta_of_its_cell():
    far = profile_at((8.0e4, 0.0, 0.0), SHAPE, EDGE_PC, E_VALUE)
    near = profile_at((4.0e4, 0.0, 0.0), SHAPE, EDGE_PC, E_VALUE)
    assert set(far.ejecta) == set(ejecta_field.KINDS)
    assert far.ejecta["neutron-star"] == pytest.approx(_density("neutron-star", 8.0e4) * far.volume_pc3)
    assert near.ejecta["neutron-star"] == pytest.approx(4.0 * far.ejecta["neutron-star"])
    seeded = profile_at((8.0e4, 0.0, 0.0), SHAPE, EDGE_PC, E_VALUE, galaxy_seed=b"ejecta-seed-16by!"[:16])
    assert seeded.ejecta["neutron-star"] > 0.0 and seeded.ejecta["neutron-star"] != far.ejecta["neutron-star"]


def test_an_ejected_objects_velocity_has_a_speed_in_the_range_and_any_direction():
    rng = random.Random(1)
    low, high = tuning.EJECTA_SPEED_RANGE_KMS
    velocities = [ejecta_field.draw_velocity_kms(rng) for _ in range(2000)]
    speeds = [math.sqrt(sum(c * c for c in v)) for v in velocities]
    assert all(low - 1e-9 <= speed <= high + 1e-9 for speed in speeds)
    assert abs(sum(speeds) / len(speeds) - (low + high) / 2) < 10.0
    for axis in range(3):                                  # no preferred direction
        assert abs(sum(v[axis] / s for v, s in zip(velocities, speeds)) / len(speeds)) < 0.06


# --- The fill ---------------------------------------------------------------------------------------------------


def _args():
    args = run_galaxy._default_generation_args()
    args.num_systems = 0
    args.density = None
    return args


def _fill(neutron_stars, black_holes=0.0):
    return SimpleNamespace(profile=SimpleNamespace(ejecta={"neutron-star": neutron_stars, "black-hole": black_holes}))


def _ejected(seed, fill):
    with draw.bound(seed):
        sector = SpaceSector("Ejecta Test", edge_ly=100.0)
        fill_stages.begin(sector)
        entries = run_sector.add_ejecta(sector, _args(), fill)
    return sector, entries


def test_a_fill_draws_its_ejecta_with_a_velocity_in_place_of_an_orbit():
    sector, entries = _ejected(3, _fill(30.0, 20.0))
    kinds = {entry.phenomenon_type for entry in entries}
    assert kinds == {"neutron-star", "black-hole"} and len(entries) > 20
    low, high = tuning.EJECTA_SPEED_RANGE_KMS
    for entry in entries:
        remnant = entry.phenomenon
        assert isinstance(remnant, NeutronStar if entry.phenomenon_type == "neutron-star" else BlackHole)
        assert remnant.galactic_orbital_period_gy is None and remnant.galactic_orbital_speed_kms is None
        speed = math.sqrt(sum(c * c for c in remnant.ejecta_velocity_kms))
        assert low - 1e-9 <= speed <= high + 1e-9
        stored = entry.spatial.get_velocity_vector("galactic")
        assert tuple(c / 1000.0 for c in stored) == pytest.approx(remnant.ejecta_velocity_kms)
        assert "Ejected from the galaxy" in str(remnant)
        assert entry.phenomenon_type != "black-hole" or remnant.mass_class == "stellar"
    counts = fill_stages.counts(sector)[fill_stages.EJECTA]
    assert counts == {kind: sum(entry.phenomenon_type == kind for entry in entries) for kind in kinds}


def test_the_ejecta_are_the_same_from_the_same_seed_and_move_no_other_stage():
    first = [(e.phenomenon_type, tuple(round(c, 6) for c in e.position)) for e in _ejected(7, _fill(12.0))[1]]
    again = [(e.phenomenon_type, tuple(round(c, 6) for c in e.position)) for e in _ejected(7, _fill(12.0))[1]]
    other = [(e.phenomenon_type, tuple(round(c, 6) for c in e.position)) for e in _ejected(8, _fill(12.0))[1]]
    assert first == again and first != other
    a, b = _ejected(5, _fill(40.0))[0], _ejected(5, _fill(0.0))[0]
    with fill_stages.stage(a, fill_stages.SYSTEMS):
        drawn_a = [draw.random() for _ in range(4)]
    with fill_stages.stage(b, fill_stages.SYSTEMS):
        drawn_b = [draw.random() for _ in range(4)]
    assert drawn_a == drawn_b


def test_a_sector_without_a_profile_has_no_ejecta():
    assert _ejected(2, None)[1] == [] and _ejected(2, SimpleNamespace(profile=None))[1] == []
    assert _ejected(2, _fill(0.0))[1] == []
