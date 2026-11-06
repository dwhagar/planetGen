# tests/test_cluster_profile.py

"""MAP.381: the projected density profile the Sector Map's cluster ball is drawn from."""

import pytest

from planetgen.physics import king
from planetgen.web.lib import cluster_profile

COMPONENT = {"mass_solar": 2.0e5, "half_mass_radius_pc": 4.0, "king_w0": 7}


def _component(w0=7):
    half_mass = COMPONENT["half_mass_radius_pc"]
    return dict(COMPONENT, king_w0=w0, tidal_radius_pc=king.tidal_radius_pc(w0, half_mass))


def test_the_profile_starts_at_one_falls_and_ends_at_zero_at_the_tidal_radius():
    profile = cluster_profile.surface_profile(_component())
    assert len(profile) == cluster_profile.PROFILE_STEPS + 1
    assert profile[0] == 1.0 and profile[-1] == 0.0
    assert all(b <= a for a, b in zip(profile, profile[1:]))
    assert profile[1] < 1.0


def test_a_more_concentrated_cluster_is_brighter_at_its_middle_than_at_its_edge_by_more():
    loose = cluster_profile.surface_profile(_component(3))
    tight = cluster_profile.surface_profile(_component(12))
    assert tight[2] < loose[2], "the denser core of a concentrated cluster falls off faster, in tidal radii"


def test_a_tidal_radius_cut_below_the_king_edge_shows_only_the_inner_part_of_the_profile():
    whole = cluster_profile.surface_profile(_component())
    cut = cluster_profile.surface_profile(dict(_component(), tidal_radius_pc=_component()["tidal_radius_pc"] / 2))
    assert cut[0] == 1.0 and cut[-1] == 0.0, "a ball cut short is still bright in the middle and has no limb light"
    assert cut[4] > whole[4], "the picture ends sooner, so the same fraction of its width is nearer the bright middle"


@pytest.mark.parametrize("broken", [{}, {"king_w0": None}, {"king_w0": 99}, {"mass_solar": 0}, {"half_mass_radius_pc": -1}])
def test_a_component_the_model_cannot_use_gives_no_profile(broken):
    component = dict(_component(), **broken) if broken else {}
    assert cluster_profile.surface_profile(component) is None
