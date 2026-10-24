"""GEN.132: the regional factors for neutron stars and black holes."""

from planetgen.galaxy import remnant_distribution as rd, scales


def test_remnants_are_rarer_in_the_plane_and_commoner_high_up():
    for kind in ("neutron-star", "black-hole"):
        assert rd.vertical_factor(kind, 0.0, 300.0) < 0.6
        assert rd.vertical_factor(kind, 1500.0, 300.0) > 2.0


def test_vertical_factor_is_capped_and_unknown_kinds_are_flat():
    assert rd.vertical_factor("neutron-star", 1e6, 300.0) <= 4.0
    assert rd.vertical_factor("planetary-nebula", 0.0, 300.0) == 1.0


def test_black_holes_gain_toward_the_core_only():
    assert rd.radial_factor("black-hole", 0.0) > rd.radial_factor("black-hole", 8200.0) > 1.0
    assert rd.radial_factor("neutron-star", 0.0) == 1.0


def test_pulsar_profile_is_one_at_the_reference_radius_and_thin_in_the_inner_galaxy():
    reference = scales.reference_radius_pc()
    assert abs(rd.pulsar_radial_factor(reference) - 1.0) < 1e-9
    assert rd.pulsar_radial_factor(0.12 * reference) < 0.2
    assert rd.pulsar_radial_factor(0.0) == 0.0


def test_pulsar_profile_has_no_pole_far_out():
    """GEN.622: the disk's exponential underflows beyond about 580 kpc; the factor stays continuous and capped."""
    for radius_kpc in (579.0, 581.0, 10000.0):
        value = rd.pulsar_radial_factor(radius_kpc * 1000.0)
        assert 0.0 <= value <= 2.0
    assert abs(rd.pulsar_radial_factor(581000.0) - rd.pulsar_radial_factor(579000.0)) < 1e-9
