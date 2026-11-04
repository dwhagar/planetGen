# tests/test_halo_threshold.py

"""A galaxy under the halo threshold (`dials.HALO_MIN_HALO_MASS_SOLAR`, Boss 7:33 PM EDT 10-21) runs no halo function."""

import types

import pytest

from planetgen.galaxy import density, dials, disk_history, enclosed_mass, mass_chain, shape_draw, structure_draw

SEED = bytes.fromhex("00112233445566778899AABBCCDDEEFF")
HALO_FUNCTIONS = (
    (mass_chain, "median_halo_mass"), (mass_chain, "conditional_halo_mass"), (mass_chain, "virial_radius_kpc"),
    (mass_chain, "median_concentration"), (disk_history, "halo_track"),
)
BELOW = mass_chain.HALO_MIN_STELLAR_MASS_SOLAR * 0.9
ABOVE = mass_chain.HALO_MIN_STELLAR_MASS_SOLAR * 1.1


def _watch(monkeypatch):
    """Replaces each halo function with one that records its call and runs the real one."""
    calls = []
    for module, name in HALO_FUNCTIONS:
        real = getattr(module, name)

        def spy(*args, _real=real, _name=name, **kwargs):
            calls.append(_name)
            return _real(*args, **kwargs)

        monkeypatch.setattr(module, name, spy)
    return calls


def test_the_threshold_is_the_stellar_mass_of_the_lowest_halo():
    assert dials.HALO_MIN_HALO_MASS_SOLAR == 1.0e10
    assert mass_chain.HALO_MIN_STELLAR_MASS_SOLAR == pytest.approx(
        mass_chain.median_stellar_mass(dials.HALO_MIN_HALO_MASS_SOLAR))
    assert 4.0e6 < mass_chain.HALO_MIN_STELLAR_MASS_SOLAR < 5.0e6
    assert mass_chain.has_halo(mass_chain.HALO_MIN_STELLAR_MASS_SOLAR)
    assert not mass_chain.has_halo(BELOW) and mass_chain.has_halo(ABOVE)


def test_below_the_threshold_no_halo_function_is_called(monkeypatch):
    calls = _watch(monkeypatch)
    chain = mass_chain.draw_mass_chain(SEED, BELOW)
    assert (chain.halo_mass_solar, chain.virial_radius_kpc, chain.concentration) == (None, None, None)
    structure = structure_draw.draw_structure(SEED, BELOW)
    assert structure.halo_stellar_ratio == 0.0
    inflow = disk_history.inflow_shape(None, 10.0, 20)
    assert list(inflow) == [1.0] * 21
    draw = shape_draw.draw_galaxy(SEED, stellar_mass_solar=BELOW)
    assert draw.chain.halo_mass_solar is None
    assert calls == []


def test_at_and_above_the_threshold_the_halo_functions_run(monkeypatch):
    for mass in (mass_chain.HALO_MIN_STELLAR_MASS_SOLAR, ABOVE):
        calls = _watch(monkeypatch)
        chain = mass_chain.draw_mass_chain(SEED, mass * (1.0 + 1e-9))
        assert chain.halo_mass_solar > 0 and chain.virial_radius_kpc > 0 and chain.concentration > 0
        assert {"conditional_halo_mass", "virial_radius_kpc", "median_concentration"} <= set(calls)
        monkeypatch.undo()


def test_a_galaxy_with_no_halo_has_only_its_stars_in_its_tides():
    def plan(mass):
        lever = {row.dest: None for row in dials.LEVERS}
        return shape_draw.plan_shape(types.SimpleNamespace(**{**lever, "mass": mass}), SEED)[0]

    below = enclosed_mass.mass_model(plan(BELOW), SEED)
    assert below.dark_solar == 0.0 and below.dark_halo_solar(1.0e5) == 0.0
    assert below.stars_solar(1.0e5) == pytest.approx(BELOW, rel=0.1)
    above = enclosed_mass.mass_model(plan(ABOVE), SEED)
    assert above.dark_solar > 0.0 and above.enclosed_solar(1.0e6) > ABOVE


def test_a_drawn_galaxy_is_unchanged_by_the_threshold():
    # the first family (3e9 and up) is far over it: the seed's own chain and halo ratio are drawn as before
    chain = mass_chain.draw_mass_chain(SEED)
    assert chain.stellar_mass_solar > 1.0e9 and chain.halo_mass_solar > 1.0e11
    assert structure_draw.draw_structure(SEED).halo_stellar_ratio > 0.0


def _args(**given):
    return type("Args", (), {**{row.dest: None for row in dials.LEVERS}, **given})()


def test_a_planned_galaxy_under_the_threshold_has_no_halo_anywhere():
    shape, drawn = shape_draw.plan_shape(_args(mass=BELOW, disk_scale_length=300.0), SEED)
    assert shape.masses == (BELOW, None)
    assert density.halo_mass(shape) == 0.0
    assert disk_history.stellar_and_halo_masses(shape) == (BELOW, None)
    assert disk_history.history_of(shape).rate.size > 0


def test_a_planned_galaxy_over_the_threshold_has_a_dark_halo_and_over_3e8_a_stellar_one():
    shape, drawn = shape_draw.plan_shape(_args(mass=ABOVE, disk_scale_length=300.0), SEED)
    assert shape.masses[1] > 0 and density.halo_mass(shape) == 0.0      # under 3e8: one disc plus a bulge
    shape, drawn = shape_draw.plan_shape(_args(mass=1.0e9, disk_scale_length=800.0), SEED)
    assert shape.masses[1] > 0 and density.halo_mass(shape) > 0.0


# --- one disc plus a bulge under 3e8 solar masses of stars (Boss, 7:58 PM EDT 10-21) --------------------------------

FULL = dials.FULL_STRUCTURE_MIN_STELLAR_MASS_SOLAR


def test_under_3e8_a_galaxy_is_one_disc_plus_a_bulge():
    for mass in (281.0, ABOVE, FULL * 0.99):
        structure = structure_draw.draw_structure(SEED, mass)
        assert structure.thick_mass_ratio == 0.0 and structure.halo_stellar_ratio == 0.0
        assert structure.nuclear_cluster_mass_solar == 0.0 and structure.nuclear_disc_mass_solar == 0.0
        assert not structure.bar_present
        assert structure.thick_length_ratio > 0 and structure.thick_height_ratio > 0   # nothing divides by zero


def test_from_3e8_up_the_full_structure_is_drawn_as_before():
    for mass in (FULL, 1.0e9, 3.0e9):
        structure = structure_draw.draw_structure(SEED, mass)
        assert structure.thick_mass_ratio > 0.0 and structure.halo_stellar_ratio > 0.0


def test_a_stellar_halo_pin_is_refused_under_3e8_and_a_dark_halo_pin_under_the_halo_threshold():
    with pytest.raises(ValueError, match="no stellar halo to pin"):
        structure_draw.draw_structure(SEED, 1.0e8, 0.1)
    assert structure_draw.draw_structure(SEED, FULL, 0.1).halo_stellar_ratio == 0.1
    # a galaxy with a dark halo (over 4.6e6) but under 3e8 may pin the dark halo only
    assert mass_chain.draw_mass_chain(SEED, 1.0e8, 5.0e10).halo_mass_solar == 5.0e10


def test_a_simple_galaxy_has_all_its_disc_stars_in_the_thin_disc_and_no_nuclear_terms():
    shape, drawn = shape_draw.plan_shape(_args(mass=1.0e8, disk_scale_length=300.0), SEED)
    terms = density.model_terms(shape)
    assert terms["thick_disk_amplitude"] == 0.0 and terms["nuclear_cluster_mass"] == 0.0
    assert terms["nuclear_disc_mass"] == 0.0 and density.halo_mass(shape) == 0.0
    assert shape_draw.thin_disk_share(drawn.structure, drawn.bulge_to_total) == pytest.approx(1.0 - drawn.bulge_to_total)


def test_a_galaxy_under_the_threshold_ejects_neutron_stars_against_its_stars_only(monkeypatch):
    calls = _watch(monkeypatch)
    share = structure_draw.draw_neutron_star_unbound_share(SEED, BELOW)
    assert 0.0 <= share <= 1.0
    assert calls == []
    # With no halo to hold them the stars' shallower potential lets out more than a halo of the same stars would.
    assert share >= structure_draw.draw_neutron_star_unbound_share(SEED, ABOVE)
