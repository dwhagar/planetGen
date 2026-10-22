# tests/test_galaxy_scale.py

"""
The Scale dial (GEN.603): mass, size and density are three views of two numbers, so any two pinned fix the third, all
three are refused, and a pinned value replaces only its own draw.
"""

import math
import types

import pytest

from planetgen.galaxy import dials, shape_draw

SEED = bytes.fromhex("00112233445566778899AABBCCDDEEFF")
OTHER_SEED = bytes.fromhex("FFEEDDCCBBAA99887766554433221100")


def _args(**given):
    return types.SimpleNamespace(**{**{row.dest: None for row in dials.LEVERS}, **given})


def _plan(seed=SEED, **given):
    return shape_draw.plan_shape(_args(**given), seed)


def _density(shape, drawn):
    share = shape_draw.thin_disk_share(drawn.structure, drawn.bulge_to_total)
    return shape_draw.central_surface_density(drawn.chain.stellar_mass_solar, shape.disk_scale_length_pc, share)


def test_no_pin_draws_the_galaxy_the_seed_makes():
    shape, drawn = _plan()
    plain = shape_draw.draw_galaxy(SEED)
    assert shape.disk_scale_length_pc == plain.disk_scale_length_pc
    assert drawn.chain == plain.chain and shape.masses == (plain.chain.stellar_mass_solar, plain.chain.halo_mass_solar)
    assert shape_draw.pinned_stellar_mass(_args(), drawn) is None
    assert shape_draw.pins_of(_args()) == {}


def test_a_pinned_mass_is_the_galaxy_mass_and_moves_the_size_with_the_law():
    shape, drawn = _plan(mass=5.0e9)
    assert drawn.chain.stellar_mass_solar == 5.0e9 and shape.masses[0] == 5.0e9
    assert shape_draw.pinned_stellar_mass(_args(mass=5.0e9), drawn) == 5.0e9
    assert shape_draw.pins_of(_args(mass=5.0e9)) == {"mass": 5.0e9}


def test_the_same_seed_at_twice_the_mass_keeps_its_residuals():
    light, drawn_light = _plan(mass=4.0e9)
    heavy, drawn_heavy = _plan(mass=8.0e9)
    z_light = math.log(light.disk_scale_length_pc * shape_draw.EXPONENTIAL_HALF_LIGHT_OVER_SCALE / 1000.0
                       / shape_draw.size_median_kpc(4.0e9)) / shape_draw.size_sigma_ln(4.0e9)
    z_heavy = math.log(heavy.disk_scale_length_pc * shape_draw.EXPONENTIAL_HALF_LIGHT_OVER_SCALE / 1000.0
                       / shape_draw.size_median_kpc(8.0e9)) / shape_draw.size_sigma_ln(8.0e9)
    assert z_light == pytest.approx(z_heavy, abs=1e-9)
    assert heavy.disk_scale_length_pc > light.disk_scale_length_pc
    # the thick disc's residual in sigma is the seed's own too
    assert drawn_light.structure.thick_height_ratio == pytest.approx(drawn_heavy.structure.thick_height_ratio)


def test_mass_and_size_give_the_density():
    shape, drawn = _plan(mass=6.0e9, disk_scale_length=1500.0)
    assert shape.disk_scale_length_pc == 1500.0 and drawn.chain.stellar_mass_solar == 6.0e9
    scale = shape_draw.resolve_scale(SEED, mass=6.0e9, length_pc=1500.0)
    assert scale.derived == "density" and scale.density == pytest.approx(_density(shape, drawn), rel=1e-9)


def test_mass_and_density_give_the_size():
    reference = shape_draw.resolve_scale(SEED, mass=6.0e9, length_pc=1500.0)
    shape, drawn = _plan(mass=6.0e9, density=reference.density)
    assert shape.disk_scale_length_pc == pytest.approx(1500.0, rel=1e-6)
    assert _density(shape, drawn) == pytest.approx(reference.density, rel=1e-6)


def test_size_and_density_give_the_mass():
    reference = shape_draw.resolve_scale(SEED, mass=6.0e9, length_pc=1500.0)
    shape, drawn = _plan(disk_scale_length=1500.0, density=reference.density)
    assert drawn.chain.stellar_mass_solar == pytest.approx(6.0e9, rel=1e-6)
    assert shape_draw.pinned_stellar_mass(_args(disk_scale_length=1500.0, density=reference.density), drawn) == \
        drawn.chain.stellar_mass_solar


def test_each_pair_reproduces_the_third_for_a_different_seed():
    base = shape_draw.resolve_scale(OTHER_SEED, mass=2.0e10, length_pc=3000.0)
    by_size = shape_draw.resolve_scale(OTHER_SEED, mass=2.0e10, density=base.density)
    by_mass = shape_draw.resolve_scale(OTHER_SEED, length_pc=3000.0, density=base.density)
    assert by_size.disk_scale_length_pc == pytest.approx(3000.0, rel=1e-6)
    assert by_mass.stellar_mass_solar == pytest.approx(2.0e10, rel=1e-6)


def test_all_three_are_refused_naming_what_the_other_two_imply():
    with pytest.raises(shape_draw.ScaleError, match="give any two"):
        _plan(mass=6.0e9, disk_scale_length=1500.0, density=1.0)
    with pytest.raises(shape_draw.ScaleError, match=r"central surface density of [0-9.e+-]+ solar masses"):
        _plan(mass=6.0e9, compactness=0.5, density=1.0)


def test_the_size_set_twice_is_refused():
    with pytest.raises(shape_draw.ScaleError, match="set twice"):
        _plan(compactness=0.5, disk_scale_length=1500.0)


def test_a_solved_value_outside_its_range_is_refused():
    with pytest.raises(shape_draw.ScaleError, match="outside"):
        _plan(mass=6.0e9, density=1.0e-7)       # a disc far longer than the 8250 pc ceiling
    with pytest.raises(shape_draw.ScaleError):
        _plan(disk_scale_length=1500.0, density=1.0e9)   # a mass no galaxy of the table reaches


def test_compactness_is_the_size_in_sigma_and_survives_a_new_mass():
    short, _ = _plan(mass=6.0e9, compactness=-1.0)
    long, _ = _plan(mass=6.0e9, compactness=1.0)
    median, _ = _plan(mass=6.0e9, compactness=0.0)
    assert short.disk_scale_length_pc < median.disk_scale_length_pc < long.disk_scale_length_pc
    assert median.disk_scale_length_pc == pytest.approx(
        shape_draw.size_median_kpc(6.0e9) / shape_draw.EXPONENTIAL_HALF_LIGHT_OVER_SCALE * 1000.0)


def test_a_pinned_flatness_keeps_its_ratio_when_the_mass_moves():
    for mass in (3.0e9, 3.0e10):
        shape, _ = _plan(mass=mass, flatness=8.0)
        assert shape.disk_scale_length_pc / shape.disk_scale_height_pc == pytest.approx(8.0)


def test_the_scale_height_limit_gates_a_pinned_size():
    with pytest.raises(ValueError, match="scale height"):
        _plan(mass=6.0e9, disk_scale_length=8000.0, flatness=3.0)     # 2667 pc, over the 1155 pc ceiling


def test_the_stored_mass_is_drawn_again_by_the_loaded_shape():
    from planetgen.galaxy import mass_chain, structure_draw
    shape, drawn = _plan(mass=7.0e9)
    assert structure_draw.draw_structure(SEED, 7.0e9) == drawn.structure
    assert mass_chain.draw_mass_chain(SEED, 7.0e9) == drawn.chain


# --- The Halo dial (GEN.606) -------------------------------------------------------------------------------------

def test_a_pinned_halo_ratio_sets_the_halo_mass_and_keeps_the_stellar_mass():
    plain = _plan()[1].chain
    shape, drawn = _plan(halo_ratio=100.0)
    assert drawn.chain.stellar_mass_solar == plain.stellar_mass_solar
    assert drawn.chain.halo_mass_solar == shape_draw.halo_mass_of(plain.stellar_mass_solar, 100.0)
    assert shape.masses == (plain.stellar_mass_solar, drawn.chain.halo_mass_solar)
    assert drawn.chain.virial_radius_kpc > plain.virial_radius_kpc


def test_a_pinned_stellar_halo_share_is_the_structures_halo_ratio():
    _shape, drawn = _plan(stellar_halo_share=0.1)
    assert drawn.structure.halo_stellar_ratio == 0.1


def test_a_halo_pin_is_refused_for_a_galaxy_with_no_halo():
    with pytest.raises(ValueError, match="no halo to pin"):
        _plan(mass=1.0e6, halo_ratio=100.0)
    with pytest.raises(ValueError, match="no stellar halo to pin"):
        _plan(mass=1.0e8, stellar_halo_share=0.1)


def test_the_halo_report_names_the_halo_and_warns_of_a_disc_it_could_not_hold():
    _shape, drawn = _plan()
    line, warning = shape_draw.halo_report(drawn)
    assert line.startswith("Halo: ") and "virial radius" in line and warning is None
    _shape, pinned = _plan(mass=6.0e9, disk_scale_length=100.0)     # a disc a tenth of the usual size in its halo
    assert "below the usual" in shape_draw.halo_report(pinned)[1]
    _shape, light = _plan(mass=1.0e6, disk_scale_length=100.0)
    assert shape_draw.halo_report(light) == (
        "No halo: the galaxy is under the halo threshold, so no halo is drawn and no halo law runs", None)


# --- the hard limits of the density report -----------------------------------------------------------------------

def test_a_bulge_under_half_a_cell_is_refused():
    with pytest.raises(ValueError, match="bulge scale radius"):
        _plan(mass=6.0e9, disk_scale_length=60.0, bulge_size=0.02)


def test_the_warnings_name_a_shape_the_grid_holds_badly():
    notes = dials.shape_warnings(sigma0=0.5, height_pc=5.0, bulge_radius_pc=3.0)
    assert len(notes) == 3 and "percent of its stars sit in planned sectors" in notes[-1]
    assert dials.shape_warnings(sigma0=500.0, height_pc=180.0, bulge_radius_pc=400.0) == []
    assert dials.planned_share(0.5) == 0.0 and dials.planned_share(20.0) == pytest.approx(0.51)
    assert 0.2 < dials.planned_share(7.0) < 0.35 and dials.planned_share(1e6) == 0.94


def test_the_pitch_angle_is_refused_within_a_degree_of_flat():
    row = dials.lever("pitch_angle_deg")
    assert row.problem(0.5) and row.problem(-0.5) and row.problem(1.0) is None and row.problem(-90.0) is None
