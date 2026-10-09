"""GEN.132: the regional factors for neutron stars and black holes."""

from planetgen.galaxy import remnant_distribution as rd


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


def test_pulsar_profile_is_one_at_the_sun_and_thin_in_the_inner_galaxy():
    assert abs(rd.pulsar_radial_factor(8200.0) - 1.0) < 1e-9
    assert rd.pulsar_radial_factor(1000.0) < 0.2
    assert rd.pulsar_radial_factor(0.0) == 0.0
