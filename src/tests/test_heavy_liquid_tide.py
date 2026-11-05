# tests/test_heavy_liquid_tide.py

"""
GEN.1052: the fluid Love numbers of a body with a core under a deep global liquid layer (`planetgen.physics.heavy_liquid_tide`), held to their three limits, and the tide of a
non-water sea in `ocean_tides`.
"""

import pytest

from planetgen.physics import heavy_liquid_tide as hl
from planetgen.physics import ocean_tides


def test_a_homogeneous_fluid_body_has_the_love_numbers_three_halves_and_five_halves():
    k2, h2 = hl.two_layer_love(1000.0, 0.5, 1000.0)
    assert k2 == pytest.approx(1.5) and h2 == pytest.approx(2.5)


def test_a_fluid_body_has_its_surface_on_an_equipotential_so_one_plus_k2_minus_h2_is_zero():
    for core, fraction, liquid in ((11000.0, 0.55, 4500.0), (5500.0, 0.9, 1000.0), (3000.0, 0.3, 920.0)):
        k2, h2 = hl.fluid_love(core, fraction, liquid)
        assert 1.0 + k2 - h2 == pytest.approx(0.0, abs=1e-9)


def test_a_core_that_does_not_move_gives_the_thin_sea_loading_factor_at_any_depth():
    """A sea over a rigid floor: `h2 = 1 / (1 - 3 rho_l / (5 rho_mean))`, so the thin formula of `ocean_tides` is exact however deep the layer is."""
    for fraction in (0.99, 0.9, 0.5, 0.1):
        mean = hl.mean_density_kg_m3(5500.0, fraction, 1000.0)
        k2, h2 = hl.two_layer_love(5500.0, fraction, 1000.0, core_rigid=True)
        assert h2 == pytest.approx(ocean_tides.loading_factor(1000.0, mean))
        assert (k2, h2) == pytest.approx(hl.rigid_floor_love(mean, 1000.0))
        assert k2 == pytest.approx(3.0 * 1000.0 / (5.0 * mean) * h2)


def test_an_earth_like_body_with_a_molten_core_has_the_hydrostatic_love_numbers_to_about_a_tenth():
    """An iron core (55 percent of the radius, 11,000 kg/m3) under a 4,500 kg/m3 mantle: an incompressible fluid Earth, whose k2 is 0.93 (the hydrostatic Earth) to 1.04 here."""
    k2, h2 = hl.fluid_love(11000.0, 0.55, 4500.0)
    assert k2 == pytest.approx(0.93, rel=0.15) and h2 == pytest.approx(1.0 + k2)


def test_a_yielding_core_is_between_the_rigid_floor_and_the_homogeneous_fluid():
    rigid = hl.two_layer_love(5500.0, 0.9, 1000.0, core_rigid=True)
    yielding = hl.two_layer_love(5500.0, 0.9, 1000.0)
    assert rigid[1] < yielding[1] < 2.5 and rigid[0] < yielding[0] < 1.5


def test_a_layer_heavier_than_its_core_or_a_liquid_denser_than_the_body_is_refused():
    with pytest.raises(hl.HeavyLiquidError, match="would sink"):
        hl.two_layer_love(3000.0, 0.5, 13500.0)
    with pytest.raises(hl.HeavyLiquidError, match="denser than the body"):
        hl.rigid_floor_love(3000.0, 13500.0)
    with pytest.raises(hl.HeavyLiquidError):
        hl.two_layer_love(5500.0, 1.2, 1000.0)
    assert not hl.layer_is_stable(3000.0, 13500.0) and hl.layer_is_stable(11000.0, 4500.0)
