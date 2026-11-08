"""
GEN.1362: the seven minor oxides of a rocky body's rock (`planetgen.physics.mantle_oxides`, called by `bulk_draw.dry_partition`).
"""

import dataclasses
from types import SimpleNamespace

import pytest

from planetgen.physics import bulk_composition, bulk_draw, mantle_oxides, rock_types, surface_composition, volatile_depletion

SUN_RATIOS = bulk_draw.HostRatios(fe=0.81, si=0.91, al=0.07, ca=0.06, ni=0.045, na=0.0468, k=0.0033, ti=0.0026, cr=0.0117, mn=0.0074, p=0.0072,
                                  c_over_o=0.55, y_envelope=0.27, z=0.014)

BSE = dict(Na2O=0.35, K2O=0.0285, TiO2=0.196, Cr2O3=0.377, MnO=0.134, NiO=0.242, P2O5=0.0189)
"""dict: The bulk silicate Earth, wt percent (the 12-oxide concordance of [R4608])."""

MARS = dict(Na2O=0.59, K2O=0.043, TiO2=0.17, Cr2O3=0.88, MnO=0.37, NiO=0.046, P2O5=0.17)
"""dict: The bulk silicate Mars, wt percent ([R4833: Table 4])."""


def _earth():
    return bulk_draw.dry_partition(SUN_RATIOS, 0.0805, 0.07, False)


def test_the_earth_is_within_a_factor_of_one_and_a_half_of_its_bulk_silicate_in_every_oxide():
    part = _earth()
    for oxide, wt in BSE.items():
        assert 100.0 * part.oxides[oxide] == pytest.approx(wt, rel=0.5), oxide
        assert 1 / 1.5 < 100.0 * part.oxides[oxide] / wt < 1.5, oxide


def test_mars_is_within_a_factor_of_two_and_a_half_in_all_but_its_nickel():
    cut = volatile_depletion.cut_temperature_k(1.524, 3.828e26)
    part = bulk_draw.dry_partition(SUN_RATIOS, 0.147, 0.07, False, 0.0, cut)
    for oxide, wt in MARS.items():
        ratio = 100.0 * part.oxides[oxide] / wt
        if oxide == "NiO":
            assert ratio > 2.5, "the open item: the leading-term partition leaves out the pressure and sulfur that take Mars's nickel (GEN.725 MO3)"
        else:
            assert 1 / 2.5 < ratio < 2.5, oxide


def test_the_core_is_the_dry_mass_less_the_rock_and_both_hold_what_the_star_gave():
    part = _earth()
    assert part.metal + part.rock == pytest.approx(1.0)
    assert 0.27 < part.metal < 0.34
    assert part.rock_feo == pytest.approx(0.0805)


def test_a_hot_zone_loses_the_volatile_oxides_and_a_cool_one_keeps_them():
    hot = bulk_draw.dry_partition(SUN_RATIOS, 0.0805, 0.07, False, 0.0, 2100.0).oxides
    cool = bulk_draw.dry_partition(SUN_RATIOS, 0.0805, 0.07, False, 0.0, 900.0).oxides
    assert hot["Na2O"] < 0.1 * _earth().oxides["Na2O"] < 0.1 * cool["Na2O"]
    assert hot["K2O"] < _earth().oxides["K2O"] < cool["K2O"]


def test_a_more_oxidised_mantle_keeps_more_of_its_phosphorus_and_a_reduced_one_less():
    reduced = bulk_draw.dry_partition(SUN_RATIOS, 0.01, 0.07, False).oxides["P2O5"]
    oxidised = bulk_draw.dry_partition(SUN_RATIOS, 0.18, 0.07, False).oxides["P2O5"]
    assert reduced < _earth().oxides["P2O5"] < oxidised
    assert oxidised / reduced > 20.0


def test_titanium_is_lithophile_and_sodium_and_potassium_stay_in_the_mantle():
    for symbol in ("Ti", "Na", "K"):
        assert mantle_oxides.core_share(symbol, 0.325) <= 0.05 + 1e-12
    assert mantle_oxides.core_share("Ti", 1.0) == 0.0


def test_the_earth_anchor_gives_back_the_earths_core_shares():
    for symbol, share in mantle_oxides.EARTH_CORE_SHARE.items():
        assert mantle_oxides.core_share(symbol, mantle_oxides.EARTH_CORE_FRACTION, mantle_oxides.EARTH_ROCK_FEO) == pytest.approx(share)


def test_a_bigger_core_takes_a_bigger_share():
    assert mantle_oxides.core_share("Cr", 0.5) > mantle_oxides.core_share("Cr", 0.325) > mantle_oxides.core_share("Cr", 0.1)
    assert mantle_oxides.core_share("Cr", 0.0) == 0.0 and mantle_oxides.core_share("Cr", 1.0) == 1.0
    with pytest.raises(ValueError):
        mantle_oxides.core_share("Cr", 1.5)


def test_a_carbon_rock_holds_none_of_the_seven():
    part = bulk_draw.dry_partition(SUN_RATIOS, 0.0, 0.07, True, 0.2)
    assert part.oxides is None


def test_the_rock_oxides_hold_the_star_atoms_the_zone_kept_less_what_the_core_took():
    """For titanium (no core share) the rock holds all the atoms the volatility law kept, per Mg atom: the TiO2 over the MgO is the star's Ti over Mg."""
    part = _earth()
    kept = volatile_depletion.retained_fraction(volatile_depletion.t50_k("Ti"), volatile_depletion.CUT_K)
    grams_tio2 = SUN_RATIOS.ti * kept * 79.866
    grams_mgo = 40.304
    # MgO is a fixed share of the rock: its mass fraction times the rock is the grams per Mg atom
    rock_grams = grams_mgo / (grams_mgo / (grams_mgo + SUN_RATIOS.si * 60.084 + SUN_RATIOS.ca * 56.077 + 0.5 * SUN_RATIOS.al * 101.961))
    assert part.oxides["TiO2"] > 0.0
    assert part.oxides["TiO2"] < grams_tio2 / grams_mgo * 1.0 and part.oxides["TiO2"] > 0.5 * grams_tio2 / rock_grams


def test_the_fractions_are_stored_on_the_bulk_composition_and_the_rock_of_a_body_reads_them_back():
    part = _earth()
    fields = bulk_draw._minor_oxides(part)
    assert set(fields) == {f"rock_{name.lower()}" for name in mantle_oxides.OXIDE_NAMES}
    bulk = bulk_composition.BulkFractions(metal=0.3, rock=0.7, water=0.0, ice=0.0, envelope=0.0, rock_feo=part.rock_feo, rock_al2o3=part.rock_al2o3,
                                          rock_mg_si=part.rock_mg_si, **fields)
    attributes = bulk.to_attributes()
    assert attributes["rock_na2o_fraction"] == part.oxides["Na2O"] and attributes["rock_p2o5_fraction"] == part.oxides["P2O5"]
    star = SimpleNamespace(metallicity_fe_h=0.0)
    body = SimpleNamespace(star=star, **attributes)
    others = {oxide: getattr(body, f"rock_{oxide.lower()}_fraction") for oxide in mantle_oxides.OXIDE_NAMES}
    oxides = rock_types.rock_oxides(part.rock_feo, part.rock_al2o3, part.rock_mg_si, 0.08, others)
    assert oxides["Na2O"] == pytest.approx(100.0 * part.oxides["Na2O"])
    assert sum(oxides.values()) == pytest.approx(100.0)
    norm = rock_types.cipw_norm(oxides)
    assert norm["albite"] > 0.0 and norm["orthoclase"] > 0.0 and norm["ilmenite"] > 0.0 and norm["chromite"] > 0.0 and norm["apatite"] > 0.0


def test_the_earth_ground_is_not_felsic_with_the_alkalis_and_stays_ultramafic():
    part = _earth()
    others = {oxide: part.oxides[oxide] for oxide in mantle_oxides.OXIDE_NAMES}
    oxides = rock_types.rock_oxides(part.rock_feo, part.rock_al2o3, part.rock_mg_si, 0.08, others)
    assert rock_types.ground_type_of(oxides) == "ultramafic"


def test_the_ratios_of_a_star_carry_the_seven_elements():
    fields = {f.name for f in dataclasses.fields(bulk_draw.HostRatios)}
    assert {"na", "k", "ti", "cr", "mn", "ni", "p"} <= fields
