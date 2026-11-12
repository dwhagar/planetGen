# tests/test_hot_cloud_deck.py

"""GEN.994 (built ahead): where the clouds of a hot atmosphere sit along a profile, and the ceilings of its trace gases (`physics/hot_cloud_deck.py`)."""

import math

import pytest

from planetgen.physics import hot_cloud_deck as hcd
from planetgen.physics import hot_condensates as hc
from planetgen.physics import stellar_abundances as sa


def _host(fe_h=0.0, **offsets):
    return sa.composition(fe_h, offsets)


def _profile(t_at_1_bar, slope=0.1):
    """A rising profile: `T(P) = t_at_1_bar * P^slope`, K."""
    return lambda pressure: t_at_1_bar * pressure ** slope


def test_the_base_is_where_the_profile_meets_the_curve():
    host = _host()
    base = hcd.cloud_base("ZnS", host, _profile(700.0))
    assert base is not None and not base.reaches_bottom
    assert hc.cloud_temperature_k("ZnS", base.base_bar, host) == pytest.approx(base.base_k, abs=0.5)


def test_a_profile_too_hot_everywhere_makes_no_cloud_and_a_cold_one_reaches_the_bottom():
    host = _host()
    assert hcd.cloud_base("ZnS", host, lambda p: 5000.0) is None
    cold = hcd.cloud_base("ZnS", host, lambda p: 100.0)
    assert cold.reaches_bottom and cold.base_bar == hcd.GRID_BAR[-1]


def test_a_hotter_profile_puts_the_cloud_higher_or_removes_it():
    host = _host()
    bases = [hcd.cloud_base("MnS", host, _profile(t)) for t in (900.0, 1100.0, 1300.0)]
    pressures = [b.base_bar if b else 0.0 for b in bases]
    assert pressures == sorted(pressures, reverse=True) and pressures[0] > pressures[-1]


def test_the_deck_is_ordered_deepest_first_and_hot_species_sit_deeper():
    layers = hcd.deck(_host(), _profile(1000.0, 0.15))
    assert [layer.base_bar for layer in layers] == sorted((layer.base_bar for layer in layers), reverse=True)
    assert [layer.species for layer in layers] == ["Cr", "MnS", "Na2S", "ZnS", "KCl"]


def test_a_host_poor_in_zinc_moves_only_the_zinc_cloud():
    profile = _profile(800.0)
    base = {layer.species: layer.base_bar for layer in hcd.deck(_host(), profile)}
    poor = {layer.species: layer.base_bar for layer in hcd.deck(_host(Zn=-0.5), profile)}
    assert poor.get("ZnS", 0.0) != base.get("ZnS", 0.0)
    for species in ("Cr", "MnS", "Na2S", "KCl"):
        assert poor.get(species) == base.get(species)


def test_ceilings_cover_every_trace_gas_and_equal_gas_ceiling():
    host = _host()
    ceilings = hcd.trace_gas_ceilings(host)
    assert set(ceilings) == set(hc.HOT_TRACE_GASES)
    assert all(math.isfinite(value) and value > 0.0 for value in ceilings.values())
    assert ceilings["HCP"] == hc.gas_ceiling("HCP", host)


def test_the_column_plan_has_one_name_per_species_and_gas_and_the_stored_quantities_match():
    names = hcd.column_names()
    assert len(names) == len(set(names)) == len(hc.CLOUD_SPECIES) + len(hc.HOT_TRACE_GASES)
    stored = hcd.stored_quantities(_host(), lambda p: 5000.0)
    assert tuple(stored) == names
    assert all(stored[f"cloud_base_bar_{name.lower()}"] is None for name in hc.CLOUD_SPECIES)
