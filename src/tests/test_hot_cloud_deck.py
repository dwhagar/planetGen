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


def test_the_adiabat_floor_is_the_mixed_gas_exponent():
    assert hcd.ADIABAT_EXPONENT == pytest.approx(2.0 / 7.0)


def test_the_text_of_the_column_gives_back_the_deck_and_the_ceilings():
    host = _host()
    clouds = hcd.deck(host, _profile(1000.0, 0.15))
    ceilings = hcd.trace_gas_ceilings(host)
    record = hcd.unpack(hcd.pack(clouds, ceilings))
    assert [(c.species, c.reaches_bottom) for c in record.clouds] == [(c.species, c.reaches_bottom) for c in clouds]
    assert [c.base_bar for c in record.clouds] == pytest.approx([c.base_bar for c in clouds], rel=1e-3)
    assert record.ceilings == pytest.approx(ceilings, rel=1e-3)
    assert hcd.unpack(None) is None and hcd.unpack("") is None


def test_a_cloud_deeper_than_the_grid_is_kept_flagged():
    clouds = hcd.deck(_host(), lambda p: 100.0)
    record = hcd.unpack(hcd.pack(clouds, {}))
    assert record.clouds and all(c.reaches_bottom and c.base_bar == hcd.GRID_BAR[-1] for c in record.clouds)


def _giants(seeds):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from tests.fuzz_support import deterministic_entropy
    found = []
    for seed in seeds:
        config = SystemConfig()
        config.PLANETS = True
        with deterministic_entropy(seed):
            system = StarSystem(system_config=config)
        found += list(system.all_planets)
    return found


@pytest.fixture(scope="module")
def bodies():
    return _giants((111, 127))


def test_a_giant_of_a_star_stores_its_decks_and_nothing_else_does(bodies):
    giants = [body for body in bodies if body.body_type == "g" and not body.is_moon]
    assert giants
    for body in bodies:
        if body in giants:
            record = hcd.unpack(body.cloud_decks)
            assert record is not None and set(record.ceilings) == set(hc.HOT_TRACE_GASES)
            assert [c.base_bar for c in record.clouds] == sorted((c.base_bar for c in record.clouds), reverse=True)
        elif hasattr(body, "cloud_decks"):
            assert body.cloud_decks is None


def test_a_warmer_giant_has_its_clouds_higher_than_a_cooler_one():
    def deck_of(starlight):
        body = [b for b in _giants((111,)) if b.body_type == "g" and not b.is_moon][0]
        body.bond_albedo = 0.3
        hcd.generate(body, starlight)
        return {c.species: c.base_bar for c in hcd.unpack(body.cloud_decks).clouds}

    cool, warm = deck_of(2.0e3), deck_of(1.5e5)
    shared = set(cool) & set(warm)
    assert shared and all(warm[name] < cool[name] for name in shared)


def test_no_starlight_or_no_albedo_gives_no_decks(bodies):
    body = [b for b in bodies if b.body_type == "g" and not b.is_moon][0]
    hcd.generate(body, 0.0)
    assert body.cloud_decks is None
    albedo, body.bond_albedo = body.bond_albedo, None
    hcd.generate(body, 1.0e4)
    assert body.cloud_decks is None
    body.bond_albedo = albedo
