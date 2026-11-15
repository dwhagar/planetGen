# tests/test_comet_ice_formation.py

"""GEN.1022: one latent ice formation temperature per comet drives its CO / H2O, N2 and Ar together (`physics/comet_ice_formation.py`). The comets measured are tests, never inputs."""

import math
import random

import pytest

from planetgen.physics import comet_ice_formation as ice
from planetgen.physics import comet_species as cs

K5_MU, K5_SIGMA = math.log10(0.025), 0.65         # row k5: the CO / H2O of the population of comets, in decades


def latent_species_over_h2o(co_over_h2o, residuals=None, n_over_c=1.0, ar_over_c=1.0, s_over_c=1.0):
    """The chain row k7 draws, with the standard normals given: N2 and Ar about the temperature's value, Kr about the argon, Xe about the krypton."""
    z = {name: 0.0 for name in ice.LATENT_SPECIES}
    z.update(residuals or {})
    unknown = set(z) - set(ice.LATENT_SPECIES)
    if unknown:
        raise ValueError(sorted(unknown))
    n2 = ice.latent_log10("N2", co_over_h2o, n_over_c, ar_over_c, s_over_c) + ice.residual_dex("N2") * z["N2"]
    argon = ice.latent_log10("Ar", co_over_h2o, n_over_c, ar_over_c, s_over_c) + ice.residual_dex("Ar") * z["Ar"]
    krypton = ice.tie_log10("Kr", 10.0 ** argon) + ice.residual_dex("Kr") * z["Kr"]
    xenon = ice.tie_log10("Xe", 10.0 ** krypton) + ice.residual_dex("Xe") * z["Xe"]
    return {"N2": 10.0 ** n2, "Ar": 10.0 ** argon, "Kr": 10.0 ** krypton, "Xe": 10.0 ** xenon}


def _population(count=4000, seed=1022):
    rng = random.Random(seed)
    comets = []
    for _ in range(count):
        co = 10.0 ** rng.gauss(K5_MU, K5_SIGMA)
        z = {name: rng.gauss(0.0, 1.0) for name in ice.LATENT_SPECIES}
        comets.append((co, latent_species_over_h2o(co, z)))
    return comets


def _correlation(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    covariance = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return covariance / math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))


def test_the_laboratory_law_gives_halley_its_anchor_and_six_decades_in_75_kelvin():
    assert ice.formation_temperature_k(0.035) == pytest.approx(48.0)
    assert ice.formation_temperature_k(0.035 * 1e-6) - ice.formation_temperature_k(0.035) == pytest.approx(75.0)


def test_a_colder_ice_holds_more_co_and_the_median_and_hale_bopp_sit_where_the_notes_say():
    assert ice.formation_temperature_k(0.2) < ice.formation_temperature_k(0.025)
    assert ice.formation_temperature_k(0.025) == pytest.approx(49.8, abs=0.1)
    assert ice.formation_temperature_k(0.2) == pytest.approx(38.5, abs=0.1)


def test_the_position_is_0_at_the_median_comet_1_at_hale_bopp_and_held_at_the_ends():
    assert ice.position(cs.CO_OVER_H2O_MEDIAN) == pytest.approx(0.0, abs=1e-9)
    assert ice.position(0.2) == pytest.approx(1.0)
    assert ice.position(0.0005) == 0.0 and ice.position(0.9) == 1.0
    assert 0.0 < ice.position(0.07) < 1.0


def test_the_median_comet_has_the_warm_end_of_every_species_and_hale_bopp_the_cold_end():
    warm = latent_species_over_h2o(cs.CO_OVER_H2O_MEDIAN)
    cold = latent_species_over_h2o(0.2)
    for species in ice.LATENT_SPECIES:
        low, high = cs.bracket_over_h2o(species)
        if species == "N2":                                   # per CO at the ends, so per water with the comet's own CO
            assert warm[species] == pytest.approx(low) and cold[species] == pytest.approx(0.09 * 0.2)
        else:
            assert warm[species] == pytest.approx(low) and cold[species] == pytest.approx(high)
    assert cold["Ar"] == pytest.approx(8.7e-3)                # Hale-Bopp's own argon per water


def test_krypton_and_xenon_are_the_comets_argon_times_the_ratios_of_67p():
    value = latent_species_over_h2o(0.07)
    assert value["Kr"] == pytest.approx(value["Ar"] * 0.089)
    assert value["Xe"] == pytest.approx(value["Kr"] * 0.269)
    moved = latent_species_over_h2o(0.07, {"Ar": 1.0})
    assert moved["Kr"] / value["Kr"] == pytest.approx(10 ** 0.4) and moved["Xe"] / value["Xe"] == pytest.approx(10 ** 0.4)


def test_a_residual_moves_only_its_own_species_by_its_width():
    base = latent_species_over_h2o(0.07)
    moved = latent_species_over_h2o(0.07, {"N2": 2.0})
    assert moved["N2"] / base["N2"] == pytest.approx(10 ** (2.0 * ice.RESIDUAL_DEX["N2"]))
    assert moved["Ar"] == base["Ar"]
    kr = latent_species_over_h2o(0.07, {"Kr": -1.0})
    assert kr["Kr"] / base["Kr"] == pytest.approx(10 ** -ice.NOBLE_TIE_DEX) and kr["Ar"] == base["Ar"]


def test_the_host_ratios_scale_the_species_they_belong_to():
    base = latent_species_over_h2o(0.07)
    assert latent_species_over_h2o(0.07, n_over_c=2.0)["N2"] == pytest.approx(2.0 * base["N2"])
    rich = latent_species_over_h2o(0.07, ar_over_c=3.0)
    assert rich["Ar"] == pytest.approx(3.0 * base["Ar"]) and rich["Kr"] == pytest.approx(3.0 * base["Kr"]) and rich["N2"] == base["N2"]


def test_the_species_rise_together_with_the_co_the_way_the_data_say():
    comets = _population()
    co = [math.log10(c) for c, _v in comets]
    for species in ("N2", "Ar", "Kr", "Xe"):
        assert _correlation(co, [math.log10(v[species]) for _c, v in comets]) > 0.7, species
    argon_poor_with_co = sum(1 for c, v in comets if c > 0.1 and v["Ar"] < 1e-5)
    argon_rich_without_co = sum(1 for c, v in comets if c < 0.025 and v["Ar"] > 1e-4)
    assert argon_rich_without_co / len(comets) < 0.01 and argon_poor_with_co / len(comets) < 0.01


def test_the_spread_across_comets_still_covers_67p_to_hale_bopp_for_argon_and_the_two_comets_for_n2():
    comets = _population()
    argon = sorted(v["Ar"] for _c, v in comets)
    assert argon[int(0.02 * len(argon))] < 1.5e-5 and argon[int(0.98 * len(argon))] > 5e-3        # 67P (5.8e-6) to Hale-Bopp (8.7e-3)
    n2_per_co = sorted(v["N2"] / c for c, v in comets)
    assert n2_per_co[int(0.05 * len(n2_per_co))] < 0.04 and n2_per_co[int(0.95 * len(n2_per_co))] > 0.07        # 67P 0.0287 to C/2016 R2 0.09


def test_the_co_of_the_comets_is_the_k5_draw_unchanged():
    """The latent is read off the CO / H2O, never drawn over it: the population's spread is the one k5 gives."""
    logs = [math.log10(c) for c, _v in _population()]
    mean = sum(logs) / len(logs)
    assert mean == pytest.approx(K5_MU, abs=0.03)
    assert math.sqrt(sum((x - mean) ** 2 for x in logs) / len(logs)) == pytest.approx(K5_SIGMA, abs=0.03)


def test_refusals():
    with pytest.raises(ValueError):
        ice.formation_temperature_k(0.0)
    with pytest.raises(ValueError):
        ice.latent_log10("OCS", 0.05)
    with pytest.raises(ValueError):
        ice.tie_log10("Ar", 1e-5)
    with pytest.raises(ValueError):
        ice.residual_dex("O2")
    with pytest.raises(ValueError):
        ice.latent_log10("Ar", 0.05, ar_over_c=0.0)
